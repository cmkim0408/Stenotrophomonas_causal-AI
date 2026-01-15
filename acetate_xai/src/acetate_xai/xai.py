from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from acetate_xai.config import load_config
from acetate_xai.io import load_conditions_csv

logger = logging.getLogger(__name__)


def _feature_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith(("width__", "mid__", "signchange__"))]


def _width_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith("width__")]


def _reaction_id_from_col(col: str) -> str:
    if "__" not in col:
        return col
    return col.split("__", 1)[1]


def _top_k_from_widths(row: pd.Series, width_cols: list[str], k: int, smallest: bool) -> str:
    vals = row[width_cols]
    vals = vals.dropna()
    if vals.empty:
        return ""
    vals = vals.sort_values(ascending=smallest)
    top = vals.head(k)
    rxns = [_reaction_id_from_col(c) for c in top.index.tolist()]
    return ";".join(rxns)


def _get_bound_from_medium_cfg(
    *,
    medium_cfg: dict[str, Any],
    conditions_row: pd.Series,
    key: str,
) -> float | None:
    """
    Return the condition-specific lower bound (uptake cap) for a named exchange.

    For acetate and ammonium, uses scaling with condition concentrations.
    For oxygen/phosphate, uses base_bounds if present.
    """
    exchanges = medium_cfg.get("exchanges", {})
    base_bounds = medium_cfg.get("base_bounds", {})
    scaling = medium_cfg.get("scaling", {})

    rxn_id = exchanges.get(key, None) if isinstance(exchanges, dict) else None
    if not rxn_id:
        return None
    rxn_id = str(rxn_id)

    # Condition-specific: acetate/ammonium
    if key == "acetate":
        k_ac = float(scaling.get("k_ac", 0.0))
        acetate_mM = conditions_row.get("acetate_mM", np.nan)
        if pd.isna(acetate_mM):
            return None
        uptake_max = max(0.0, k_ac * float(acetate_mM))
        return -uptake_max
    if key == "ammonium":
        k_nh4 = float(scaling.get("k_nh4", 0.0))
        nh4_gL = conditions_row.get("nh4cl_gL", np.nan)
        if pd.isna(nh4_gL):
            return None
        uptake_max = max(0.0, k_nh4 * float(nh4_gL))
        return -uptake_max

    # Base-only: oxygen/phosphate/etc.
    if isinstance(base_bounds, dict) and rxn_id in base_bounds and isinstance(base_bounds[rxn_id], dict):
        lb = base_bounds[rxn_id].get("lb", None)
        if lb is None:
            return None
        return float(lb)

    return None


def _is_saturated(
    *,
    features_row: pd.Series,
    rxn_id: str,
    lb: float | None,
    sat_tol: float = 1e-3,
    width_tol: float = 1e-3,
) -> pd._libs.missing.NAType | bool:
    """
    Saturation heuristic (v0):
    - needs mid__RXN and width__RXN
    - checks (mid ~ lb) AND (width small)
    """
    if lb is None:
        return pd.NA
    mid_col = f"mid__{rxn_id}"
    width_col = f"width__{rxn_id}"
    if mid_col not in features_row.index or width_col not in features_row.index:
        return pd.NA
    mid = features_row.get(mid_col, np.nan)
    width = features_row.get(width_col, np.nan)
    if pd.isna(mid) or pd.isna(width):
        return pd.NA
    return (abs(float(mid) - float(lb)) <= sat_tol) and (abs(float(width)) <= width_tol)


def build_regime_table(
    *,
    features_df: pd.DataFrame,
    conditions_df: pd.DataFrame,
    medium_cfg: dict[str, Any],
    objective_value_series: pd.Series | None = None,
) -> pd.DataFrame:
    """
    Build regime_table.csv with:
    - saturation flags for acetate/o2/nh4/pi
    - top 10 narrow/wide reactions by width
    """
    width_cols = _width_cols(features_df)

    # Ensure we have the condition metadata columns (from Step 8 join)
    base_cols = ["condition_id", "set_name", "measured_OD"]
    for c in base_cols:
        if c not in features_df.columns:
            raise ValueError(f"features.parquet missing required column: {c}")

    # Join with conditions table for bound calculations (safe even if same columns exist)
    cond_meta = conditions_df.set_index("condition_id")
    feat = features_df.set_index("condition_id")

    out_rows: list[dict[str, Any]] = []
    for cid, row in feat.iterrows():
        if cid not in cond_meta.index:
            continue
        c_row = cond_meta.loc[cid]

        ex = medium_cfg.get("exchanges", {}) if isinstance(medium_cfg.get("exchanges", {}), dict) else {}
        rxn_ac = str(ex.get("acetate", "")) if ex else ""
        rxn_o2 = str(ex.get("oxygen", "")) if ex else ""
        rxn_nh4 = str(ex.get("ammonium", "")) if ex else ""
        rxn_pi = str(ex.get("phosphate", "")) if ex else ""

        lb_ac = _get_bound_from_medium_cfg(medium_cfg=medium_cfg, conditions_row=c_row, key="acetate")
        lb_o2 = _get_bound_from_medium_cfg(medium_cfg=medium_cfg, conditions_row=c_row, key="oxygen")
        lb_nh4 = _get_bound_from_medium_cfg(medium_cfg=medium_cfg, conditions_row=c_row, key="ammonium")
        lb_pi = _get_bound_from_medium_cfg(medium_cfg=medium_cfg, conditions_row=c_row, key="phosphate")

        out_rows.append(
            {
                "condition_id": cid,
                "set_name": row.get("set_name"),
                "measured_OD": row.get("measured_OD"),
                "objective_value": (objective_value_series.loc[cid] if objective_value_series is not None and cid in objective_value_series.index else np.nan),
                "acetate_sat": _is_saturated(features_row=row, rxn_id=rxn_ac, lb=lb_ac) if rxn_ac else pd.NA,
                "o2_sat": _is_saturated(features_row=row, rxn_id=rxn_o2, lb=lb_o2) if rxn_o2 else pd.NA,
                "nh4_sat": _is_saturated(features_row=row, rxn_id=rxn_nh4, lb=lb_nh4) if rxn_nh4 else pd.NA,
                "pi_sat": _is_saturated(features_row=row, rxn_id=rxn_pi, lb=lb_pi) if rxn_pi else pd.NA,
                "top_10_narrow_reactions": _top_k_from_widths(row, width_cols, k=10, smallest=True),
                "top_10_wide_reactions": _top_k_from_widths(row, width_cols, k=10, smallest=False),
            }
        )

    return pd.DataFrame(out_rows)


def fit_elasticnet_coefficients(
    *,
    features_df: pd.DataFrame,
    target_col: str = "measured_OD",
    top_n: int = 30,
    random_state: int = 0,
) -> pd.DataFrame:
    """
    Train an ElasticNetCV (interpretable linear model) and return top |coef| features.
    NaN targets are excluded. If not enough rows, returns an empty table with headers.
    """
    from sklearn.linear_model import ElasticNetCV
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    df = features_df.copy()
    df = df.loc[~df[target_col].isna()].copy()
    if len(df) < 3:
        return pd.DataFrame(columns=["feature", "coef"])

    X = df[_feature_cols(df)].copy()
    # booleans -> 0/1
    for c in X.columns:
        if X[c].dtype == bool:
            X[c] = X[c].astype(int)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    y = pd.to_numeric(df[target_col], errors="coerce").astype(float)

    # CV folds cannot exceed n_samples
    cv = min(5, len(df))
    model = Pipeline(
        steps=[
            ("scaler", StandardScaler(with_mean=True, with_std=True)),
            ("enet", ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9, 0.95, 0.99], cv=cv, random_state=random_state)),
        ]
    )
    model.fit(X, y)
    coefs = model.named_steps["enet"].coef_
    coef_df = pd.DataFrame({"feature": X.columns, "coef": coefs})
    coef_df["abs"] = coef_df["coef"].abs()
    coef_df = coef_df.sort_values("abs", ascending=False).head(top_n).drop(columns=["abs"])
    return coef_df.reset_index(drop=True)


def fit_tree_rules(
    *,
    features_df: pd.DataFrame,
    target_col: str = "measured_OD",
    max_depth: int = 3,
    random_state: int = 0,
) -> str:
    """
    Train a shallow DecisionTreeRegressor and export rules as text.
    If not enough rows, returns a short note.
    """
    from sklearn.tree import DecisionTreeRegressor, export_text

    df = features_df.copy()
    df = df.loc[~df[target_col].isna()].copy()
    if len(df) < 3:
        return "Not enough non-NaN measured_OD rows to fit a tree (need >= 3)."

    X = df[_feature_cols(df)].copy()
    for c in X.columns:
        if X[c].dtype == bool:
            X[c] = X[c].astype(int)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    y = pd.to_numeric(df[target_col], errors="coerce").astype(float)

    tree = DecisionTreeRegressor(max_depth=max_depth, random_state=random_state)
    tree.fit(X, y)
    return export_text(tree, feature_names=list(X.columns))


def default_conditions_and_medium(
    *,
    conditions_path: str | Path = "data/conditions_experiment.csv",
    medium_path: str | Path = "configs/medium_base.yaml",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    return load_conditions_csv(conditions_path), load_config(medium_path)

