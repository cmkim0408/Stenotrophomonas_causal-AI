"""Shared utilities for the iScience revision workstreams.

Reuses the existing acetate_xai pipeline outputs (regime_dataset.parquet) and
applies a consistent feature-panel / SHAP / cross-validation interface.

XGBoost-native SHAP via Booster.predict(pred_contribs=True) is used to bypass
shap 0.49.1 ↔ xgboost 3.1.x incompatibility.
"""
from __future__ import annotations

import json
import logging
import platform
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import (
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    mean_absolute_error,
    mean_squared_error,
    precision_recall_fscore_support,
    r2_score,
)
from sklearn.model_selection import StratifiedKFold, KFold

ROOT = Path(__file__).resolve().parents[2]
LEGACY_PARQUET = ROOT / "results" / "regime_dataset.parquet"
EXTENDED_PARQUET = ROOT / "revision_runs" / "iscience_rev1" / "regime_dataset_extended.parquet"
# Auto-promote to extended once the FVA campaign extension exists.
DEFAULT_PARQUET = EXTENDED_PARQUET if EXTENDED_PARQUET.exists() else LEGACY_PARQUET
DEFAULT_ANCHORS_YAML = ROOT / "acetate_xai" / "configs" / "anchors.yaml"
REVISION_OUT = ROOT / "revision_runs" / "iscience_rev1"

logger = logging.getLogger("revision.utils")

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_regime_dataset(path: Path | str = DEFAULT_PARQUET) -> pd.DataFrame:
    df = pd.read_parquet(path)
    if "label" not in df.columns:
        raise ValueError(f"'label' column missing in {path}")
    return df


def width_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith("width__")]


def all_feature_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith(("width__", "mid__", "signchange__"))]


def severity_target(df: pd.DataFrame) -> pd.Series:
    """Normalized growth-potential index G_i = obj / obj_max (paper definition).

    Paper computes obj_max per simulation run; here we normalize against the
    dataset-level max, which is monotonic and adequate for benchmarking.
    """
    obj = df["objective_value"].astype(float)
    return obj / obj.max()


# ---------------------------------------------------------------------------
# Feature panel resolution
# ---------------------------------------------------------------------------

def _rxn_id(col: str) -> str:
    return col.split("__", 1)[1] if "__" in col else col


def load_anchors_yaml(path: Path | str = DEFAULT_ANCHORS_YAML) -> list[dict[str, Any]]:
    with open(path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    return cfg.get("anchors", [])


# Curated paper-aligned reaction IDs that exist in the iSO1_933 model.
# After the extended FVA campaign (revision_runs/.../extended_fva/), the
# deployed dataset's width__ universe expands from 120 to ~300 columns and
# now covers the full paper-narrative TCA / glyoxylate / glycolysis /
# respiration anchors (MDH/ICDH/ICL/MALS/PYK/PPC/NADH16pp/FUM, etc.).
# curated_paper_panel(df) intersects this list with whatever width__ columns
# are present, so the same constant works on both the legacy 120-width and
# the extended ~300-width parquets.
PAPER_CURATED_RXNS = (
    # Limiting / uptake exchanges (paper anchors)
    "EX_o2_e", "EX_nh4_e", "EX_pi_e", "EX_co2_e", "EX_h_e", "EX_h2o_e",
    # Energy / respiration (paper SHAP top)
    "ATPS4rpp", "ADK1", "ATPM", "AKGDH", "CYO1_KT", "NADH16pp",
    # TCA cycle (paper Fig 4 narrative — added via extended FVA)
    "CS", "ACONT", "ACONTa", "ACONTb",
    "ICDHyr", "ICDHx", "ICL", "MALS",
    "MDH", "MDH2", "MDH3", "FUM",
    # Acetate uptake / activation
    "ACS", "ACSERL",
    # Glycolytic / anaplerotic
    "ENO", "ENOPH", "PYK", "PYK3", "PPC",
    # Acetolactate synthase variants (paper severity SHAP)
    "ACLS", "ACLSa", "ACLSb",
    # ADC synthase / APS reductase (paper severity TOP-2)
    "ADCS", "APSR", "APSR2",
    # N/aa biosynthesis touched by paper text
    "GLNS", "GLUDy", "ACGS", "ARGSL", "ARGSS", "ASPTA",
)


def curated_paper_panel(df: pd.DataFrame) -> list[str]:
    """Hand-curated panel: paper-named reactions that actually exist in the
    120-width parquet. Skips IDs that are absent so the panel stays clean.
    """
    have = set(width_cols(df))
    chosen = [f"width__{rid}" for rid in PAPER_CURATED_RXNS
              if f"width__{rid}" in have]
    return chosen


def curated_panel_from_anchors(df: pd.DataFrame, anchors: list[dict[str, Any]]) -> list[str]:
    """Legacy: match anchors.yaml keywords against width__ column ids.

    Kept for diagnostic comparison; due to substring collisions and the
    truncated width__ universe this returns a noisy ~19-column set, so the
    canonical curated panel for the revision is curated_paper_panel(df).
    """
    cols = width_cols(df)
    rxns = {c: _rxn_id(c).lower() for c in cols}
    chosen: list[str] = []
    seen: set[str] = set()
    for anchor in anchors:
        keywords = [k.lower() for k in anchor.get("keywords", [])]
        for col, rid in rxns.items():
            if col in seen:
                continue
            if any(_kw_match(rid, kw) for kw in keywords):
                chosen.append(col)
                seen.add(col)
                break
    return chosen


def _kw_match(rid: str, kw: str) -> bool:
    """Match a reaction id against an anchors.yaml keyword (case-insensitive).

    Uses substring matching: anchors.yaml keywords are specific enough
    (e.g. 'EX_o2_e', 'ATPM', 'MDH') that substring is reliable; full-word
    boundaries fail on COBRA-style underscored ids like 'ex_o2_e'.
    """
    if len(kw.strip()) < 2:
        return False
    return kw.strip().lower() in rid.lower()


@dataclass
class Panel:
    name: str
    columns: list[str]
    note: str = ""

    @property
    def size(self) -> int:
        return len(self.columns)


def build_panels(
    df: pd.DataFrame,
    *,
    seeds_for_random: Sequence[int] = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    random_size: int = 30,
) -> list[Panel]:
    """Construct the static feature panels (curated, all, random controls).

    SHAP-ranked panels are produced separately by shap_ranked_panels().
    """
    all_widths = width_cols(df)
    curated = curated_paper_panel(df)
    panels: list[Panel] = [
        Panel("curated_paper", curated,
              note=f"hand-curated paper-aligned reactions present in 120-width set "
                   f"(n={len(curated)})"),
        Panel("all_widths", all_widths, note="every width__ column (superset)"),
    ]
    for seed in seeds_for_random:
        sub_rng = np.random.default_rng(seed)
        idx = sub_rng.choice(len(all_widths), size=min(random_size, len(all_widths)),
                             replace=False)
        cols = sorted(all_widths[i] for i in idx)
        panels.append(Panel(f"random_30_seed{seed}", cols, note=f"seed={seed}"))
    return panels


def shap_ranked_panels(df: pd.DataFrame, base_columns: list[str], y: np.ndarray,
                       *, sizes: Sequence[int] = (10, 20, 50)) -> list[Panel]:
    """Rank features by mean(|SHAP|) using an XGBoost classifier on `base_columns`,
    then build top-K panels.
    """
    import xgboost as xgb
    X = df[base_columns].to_numpy()
    clf = xgb.XGBClassifier(
        n_estimators=300, max_depth=4, learning_rate=0.1,
        eval_metric="mlogloss", random_state=42, n_jobs=-1, tree_method="hist",
    )
    clf.fit(X, y)
    booster = clf.get_booster()
    contribs = booster.predict(xgb.DMatrix(X), pred_contribs=True)
    # Multi-class: shape (n_samples, n_classes, n_features+1). Last = bias.
    if contribs.ndim == 3:
        feat_contrib = contribs[:, :, :-1]                  # (n, K, p)
        importance = np.abs(feat_contrib).mean(axis=(0, 1)) # (p,)
    else:
        importance = np.abs(contribs[:, :-1]).mean(axis=0)  # (p,)
    order = np.argsort(-importance)
    ranked = [base_columns[i] for i in order]
    out = []
    for k in sizes:
        cols = ranked[:k]
        out.append(Panel(
            name=f"top_{k}",
            columns=cols,
            note=f"top-{k} by mean(|SHAP|) on the {len(base_columns)}-width superset",
        ))
    return out


# ---------------------------------------------------------------------------
# Cross-validation evaluators
# ---------------------------------------------------------------------------

def cv_classify(
    X: np.ndarray, y: np.ndarray, model_factory: Callable[[], Any],
    *, n_splits: int = 5, seed: int = 42,
) -> dict[str, Any]:
    """Stratified K-fold classification. Returns metrics + OOF predictions."""
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    n_classes = int(y.max() + 1)
    oof = np.zeros_like(y)
    fold_macro_f1: list[float] = []
    fold_bal_acc: list[float] = []
    for tr, te in cv.split(X, y):
        model = model_factory()
        model.fit(X[tr], y[tr])
        pred = model.predict(X[te])
        oof[te] = pred
        fold_macro_f1.append(f1_score(y[te], pred, average="macro", zero_division=0))
        fold_bal_acc.append(balanced_accuracy_score(y[te], pred))
    macro_f1 = f1_score(y, oof, average="macro", zero_division=0)
    bal_acc = balanced_accuracy_score(y, oof)
    p, r, f, _ = precision_recall_fscore_support(y, oof, average=None,
                                                 labels=list(range(n_classes)),
                                                 zero_division=0)
    cm = confusion_matrix(y, oof, labels=list(range(n_classes)))
    return {
        "macro_f1": float(macro_f1),
        "balanced_accuracy": float(bal_acc),
        "macro_f1_fold_mean": float(np.mean(fold_macro_f1)),
        "macro_f1_fold_std": float(np.std(fold_macro_f1)),
        "balanced_accuracy_fold_mean": float(np.mean(fold_bal_acc)),
        "balanced_accuracy_fold_std": float(np.std(fold_bal_acc)),
        "per_class_precision": p.tolist(),
        "per_class_recall": r.tolist(),
        "per_class_f1": f.tolist(),
        "confusion_matrix": cm.tolist(),
        "oof": oof.tolist(),
    }


def cv_regress(
    X: np.ndarray, y: np.ndarray, model_factory: Callable[[], Any],
    *, n_splits: int = 5, seed: int = 42,
) -> dict[str, Any]:
    """K-fold regression. Returns RMSE / MAE / R² + OOF predictions."""
    cv = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    oof = np.full(len(y), np.nan, dtype=float)
    for tr, te in cv.split(X):
        model = model_factory()
        model.fit(X[tr], y[tr])
        oof[te] = model.predict(X[te])
    mse = mean_squared_error(y, oof)
    rmse = float(np.sqrt(mse))
    mae = float(mean_absolute_error(y, oof))
    r2 = float(r2_score(y, oof))
    return {
        "rmse": rmse,
        "mae": mae,
        "r2": r2,
        "oof": oof.tolist(),
    }


# ---------------------------------------------------------------------------
# Default model factories
# ---------------------------------------------------------------------------

def xgb_clf_factory(*, n_estimators: int = 300, max_depth: int = 4,
                    learning_rate: float = 0.1, seed: int = 42) -> Callable[[], Any]:
    import xgboost as xgb
    def _make():
        return xgb.XGBClassifier(
            n_estimators=n_estimators, max_depth=max_depth,
            learning_rate=learning_rate, eval_metric="mlogloss",
            random_state=seed, n_jobs=-1, tree_method="hist", verbosity=0,
        )
    return _make


def xgb_reg_factory(*, n_estimators: int = 300, max_depth: int = 4,
                    learning_rate: float = 0.1, seed: int = 42) -> Callable[[], Any]:
    import xgboost as xgb
    def _make():
        return xgb.XGBRegressor(
            n_estimators=n_estimators, max_depth=max_depth,
            learning_rate=learning_rate, random_state=seed,
            n_jobs=-1, tree_method="hist", verbosity=0,
        )
    return _make


# ---------------------------------------------------------------------------
# Misc
# ---------------------------------------------------------------------------

def env_summary() -> dict[str, str]:
    import importlib.metadata as md
    pkgs = ("numpy", "pandas", "scipy", "scikit-learn", "xgboost", "shap",
            "cobra", "causal-learn", "matplotlib", "pyarrow", "numba")
    versions: dict[str, str] = {}
    for p in pkgs:
        try:
            versions[p] = md.version(p)
        except md.PackageNotFoundError:
            versions[p] = "missing"
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": str(__import__("os").cpu_count() or 0),
        **{f"pkg.{k}": v for k, v in versions.items()},
    }


def ensure_outdirs(workstream: str) -> dict[str, Path]:
    base = REVISION_OUT
    paths = {
        "base": base,
        "ws": base / workstream,
        "figures": base / "figures",
    }
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)
    return paths


def write_md(path: Path, lines: Iterable[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@dataclass
class Timer:
    """Context-manager / decorator timer that records seconds."""
    label: str
    seconds: float = field(default=0.0, init=False)
    _t0: float = field(default=0.0, init=False)

    def __enter__(self) -> "Timer":
        self._t0 = time.perf_counter()
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.seconds = time.perf_counter() - self._t0


__all__ = [
    "ROOT", "DEFAULT_PARQUET", "DEFAULT_ANCHORS_YAML", "REVISION_OUT",
    "load_regime_dataset", "width_cols", "all_feature_cols", "severity_target",
    "Panel", "build_panels", "shap_ranked_panels",
    "load_anchors_yaml", "curated_panel_from_anchors",
    "cv_classify", "cv_regress",
    "xgb_clf_factory", "xgb_reg_factory",
    "env_summary", "ensure_outdirs", "write_md", "Timer",
]
