from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Train XGBoost + SHAP on Random LHS campaign outputs.")
    p.add_argument("--labels", default="results/campaigns/C_random_LHS/regime_labels.parquet", help="Labels parquet")
    p.add_argument("--features", default="results/campaigns/C_random_LHS/features.parquet", help="Features parquet")
    p.add_argument("--outdir", default="results/campaigns/C_random_LHS/xai_xgb", help="Output directory")
    p.add_argument("--task", choices=["regime", "severity"], default="regime", help="Task type (default: regime)")
    p.add_argument("--topk", type=int, default=20, help="Top-K features to display in SHAP plots (default: 20)")
    p.add_argument("--test-size", type=float, default=0.2, help="Test size fraction (default: 0.2)")
    p.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    p.add_argument("--n-jobs", type=int, default=-1, help="XGBoost n_jobs (default: -1)")
    p.add_argument(
        "--class-weight",
        default="balanced",
        choices=["balanced", "none"],
        help="Class imbalance handling (default: balanced).",
    )
    p.add_argument("--skip-dependence", action="store_true", help="Skip dependence plots")
    p.add_argument(
        "--drop-features",
        default="",
        help=(
            "Comma-separated feature names to drop from X (e.g., 'objective_value,pi_sat'). "
            "Applied after feature construction; missing names are ignored."
        ),
    )
    p.add_argument(
        "--severity-drop-objective",
        action="store_true",
        help="For --task severity, additionally drop objective_value from X (convenience flag for bias check).",
    )
    return p


def _safe_filename(s: str, max_len: int = 120) -> str:
    s2 = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_")
    return (s2 or "feature")[:max_len]


def _feature_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith(("width__", "mid__", "signchange__"))]


def _read_parquet(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"File not found: {p}")
    return pd.read_parquet(p)


def _merge_labels_features(labels: pd.DataFrame, features: pd.DataFrame, *, features_path: Path) -> pd.DataFrame:
    """
    Merge on sample_id; if features lacks sample_id, try to attach via sibling design.parquet.
    """
    if "sample_id" not in labels.columns:
        raise ValueError("labels parquet must contain sample_id")

    lab = labels.copy()
    lab["sample_id"] = lab["sample_id"].astype(str)

    feat = features.copy()
    if "sample_id" not in feat.columns:
        design = features_path.parent / "design.parquet"
        if not design.exists():
            raise ValueError("features parquet missing sample_id and design.parquet not found to recover it")
        d = pd.read_parquet(design)
        if "sample_id" not in d.columns:
            raise ValueError("design.parquet missing sample_id")
        # try to merge by design knobs
        keys = [c for c in ["acetate_mM", "o2_uptake_max", "nh4_uptake_max", "atpm"] if c in d.columns and c in feat.columns]
        if len(keys) < 2:
            raise ValueError("Cannot recover sample_id for features: insufficient shared columns with design.parquet")
        feat = feat.merge(d[["sample_id"] + keys], on=keys, how="left")
        if feat["sample_id"].isna().any():
            raise ValueError("Failed to recover sample_id for some feature rows via design.parquet merge")

    feat["sample_id"] = feat["sample_id"].astype(str)
    merged = lab.merge(feat, on="sample_id", how="left", suffixes=("", "_feat"))
    return merged


def _prepare_X_y_regime(df: pd.DataFrame, *, drop_features: set[str]) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    # 1) status optimal only
    if "status" in df.columns:
        df = df.loc[df["status"].astype(str) == "optimal"].copy()

    # 2) primary_regime filter
    if "primary_regime" not in df.columns:
        raise ValueError("labels must contain primary_regime")
    df["primary_regime"] = df["primary_regime"].astype(str).str.strip()
    keep = {"O2_limited", "N_limited", "Ac_limited"}
    df = df.loc[df["primary_regime"].isin(keep)].copy()
    if len(df) == 0:
        raise ValueError("No rows left after filtering primary_regime to {O2_limited,N_limited,Ac_limited}")

    y_raw = df["primary_regime"].astype(str).to_numpy()

    # 3) build X
    base_cols = ["acetate_mM", "o2_uptake_max", "nh4_uptake_max", "atpm", "objective_value"]
    sat_cols = [c for c in ["acetate_sat", "o2_sat", "nh4_sat", "pi_sat"] if c in df.columns]
    fva_cols = _feature_cols(df)

    X = pd.DataFrame(index=df.index)
    for c in base_cols:
        if c in df.columns:
            X[c] = pd.to_numeric(df[c], errors="coerce")
        else:
            X[c] = np.nan
    for c in sat_cols:
        # bool -> int
        if df[c].dtype == bool:
            X[c] = df[c].astype(int)
        else:
            X[c] = pd.to_numeric(df[c], errors="coerce")

    # has_fva indicator: any non-null in FVA cols
    has_fva = np.zeros(len(df), dtype=int)
    if fva_cols:
        tmp = df[fva_cols].copy()
        # bool -> int, then to numeric
        for c in tmp.columns:
            if tmp[c].dtype == bool:
                tmp[c] = tmp[c].astype(int)
        tmp = tmp.apply(pd.to_numeric, errors="coerce")
        has_fva = (tmp.notna().any(axis=1)).astype(int).to_numpy()
        tmp = tmp.fillna(0.0)
        for c in tmp.columns:
            X[c] = tmp[c]
    X["has_fva"] = has_fva

    # 4) numeric-only
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # Drop requested features
    drop = [c for c in X.columns if c in drop_features]
    if drop:
        X = X.drop(columns=drop)
    return X, y_raw, df["sample_id"].astype(str).tolist()


def _prepare_X_y_severity(df: pd.DataFrame, *, drop_features: set[str]) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    # 1) status optimal only
    if "status" in df.columns:
        df = df.loc[df["status"].astype(str) == "optimal"].copy()
    if len(df) == 0:
        raise ValueError("No optimal rows available for severity training")

    # severity y
    if "maintenance_severity" in df.columns:
        y = pd.to_numeric(df["maintenance_severity"], errors="coerce")
    else:
        if "objective_value" not in df.columns:
            raise ValueError("Need objective_value to compute maintenance_severity")
        obj = pd.to_numeric(df["objective_value"], errors="coerce")
        # no run_id in LHS campaign; normalize by global max
        denom = float(np.nanmax(obj.to_numpy())) if np.isfinite(obj.to_numpy()).any() else np.nan
        y = obj / denom
    y = y.replace([np.inf, -np.inf], np.nan).dropna()
    df = df.loc[y.index].copy()
    yv = y.astype(float).to_numpy()

    base_cols = ["acetate_mM", "o2_uptake_max", "nh4_uptake_max", "atpm", "objective_value"]
    sat_cols = [c for c in ["acetate_sat", "o2_sat", "nh4_sat", "pi_sat"] if c in df.columns]
    fva_cols = _feature_cols(df)

    X = pd.DataFrame(index=df.index)
    for c in base_cols:
        if c in df.columns:
            X[c] = pd.to_numeric(df[c], errors="coerce")
        else:
            X[c] = np.nan
    for c in sat_cols:
        if df[c].dtype == bool:
            X[c] = df[c].astype(int)
        else:
            X[c] = pd.to_numeric(df[c], errors="coerce")

    has_fva = np.zeros(len(df), dtype=int)
    if fva_cols:
        tmp = df[fva_cols].copy()
        for c in tmp.columns:
            if tmp[c].dtype == bool:
                tmp[c] = tmp[c].astype(int)
        tmp = tmp.apply(pd.to_numeric, errors="coerce")
        has_fva = (tmp.notna().any(axis=1)).astype(int).to_numpy()
        tmp = tmp.fillna(0.0)
        for c in tmp.columns:
            X[c] = tmp[c]
    X["has_fva"] = has_fva

    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    drop = [c for c in X.columns if c in drop_features]
    if drop:
        X = X.drop(columns=drop)
    return X, yv, df["sample_id"].astype(str).tolist()


def _shap_global_importance_and_dep_values(shap_values) -> tuple[np.ndarray, np.ndarray | None]:
    """
    Return:
      - global_importance: (F,) mean(|SHAP|) aggregated over samples + classes (if multi-class)
      - dep_values: (N,F) SHAP values for dependence plots (averaged over classes if needed), or None
    Handles:
      - list[C] of (N,F)
      - ndarray (N,F)
      - ndarray (C,N,F)
      - ndarray (N,F,C)
    """
    if isinstance(shap_values, list) and shap_values:
        arr = np.stack([np.asarray(v) for v in shap_values], axis=0)  # (C,N,F)
        global_importance = np.mean(np.abs(arr), axis=(0, 1))  # (F,)
        dep_vals = np.mean(arr, axis=0)  # (N,F)
        return global_importance, dep_vals

    arr = np.asarray(shap_values)
    if arr.ndim == 2:
        return np.mean(np.abs(arr), axis=0), arr
    if arr.ndim == 3:
        # (N,F,C) or (C,N,F)
        if arr.shape[0] < 20 and arr.shape[0] != arr.shape[1]:
            # likely (C,N,F)
            global_importance = np.mean(np.abs(arr), axis=(0, 1))
            dep_vals = np.mean(arr, axis=0)
            return global_importance, dep_vals
        # assume (N,F,C)
        global_importance = np.mean(np.abs(arr), axis=(0, 2))
        dep_vals = np.mean(arr, axis=2)
        return global_importance, dep_vals
    raise TypeError(f"Unsupported shap_values shape: {getattr(arr, 'shape', None)}")


def _shap_to_n_f_k(shap_values, *, n_features: int) -> np.ndarray:
    """
    Coerce multiclass SHAP outputs into a unified ndarray (N, F, K).

    Accepts common shapes:
      - list[K] of (N, F)
      - ndarray (N, F) -> treated as (N, F, 1)
      - ndarray (K, N, F)
      - ndarray (N, F, K)
    """
    if isinstance(shap_values, list) and shap_values:
        arr = np.stack([np.asarray(v) for v in shap_values], axis=0)  # (K,N,F)
        if arr.ndim != 3:
            raise TypeError(f"Unexpected list->stack shap_values ndim: {arr.ndim}")
        if arr.shape[2] != n_features:
            raise ValueError(f"shap_values feature dim mismatch: got {arr.shape[2]}, expected {n_features}")
        return np.transpose(arr, (1, 2, 0))  # (N,F,K)

    arr = np.asarray(shap_values)
    if arr.ndim == 2:
        if arr.shape[1] != n_features:
            raise ValueError(f"shap_values feature dim mismatch: got {arr.shape[1]}, expected {n_features}")
        return arr[:, :, None]

    if arr.ndim == 3:
        # (K,N,F)
        if arr.shape[2] == n_features and arr.shape[0] < 50:
            return np.transpose(arr, (1, 2, 0))
        # (N,F,K)
        if arr.shape[1] == n_features:
            return arr
        raise ValueError(f"Unsupported 3D shap_values shape (n_features={n_features}): {arr.shape}")

    raise TypeError(f"Unsupported shap_values shape: {getattr(arr, 'shape', None)}")


def _write_shap_richness_reports(
    *,
    outdir: Path,
    task: str,
    seed: int,
    topk: int,
    feature_names: list[str],
    shap_n_f_k: np.ndarray,
    class_names: list[str] | None = None,
    class_distribution: dict[str, int] | None = None,
) -> float:
    """
    Write:
      - shap_richness_report.json
      - shap_richness_report.csv
    and optionally:
      - shap_by_class_importance.csv (regime only)

    Also returns domination_index and writes warning into run_metadata.json if >=0.7.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    n_used = int(shap_n_f_k.shape[0])
    n_features = int(shap_n_f_k.shape[1])
    k = int(shap_n_f_k.shape[2])
    n_topk = int(min(max(1, int(topk)), n_features))

    # global importance: mean(|shap|) over samples+classes -> (F,)
    mean_abs = np.mean(np.abs(shap_n_f_k), axis=(0, 2))
    order = np.argsort(mean_abs)[::-1]
    top_idx = order[:n_topk]
    top_vals = mean_abs[top_idx]
    sum_top = float(np.sum(top_vals)) if float(np.sum(top_vals)) > 0 else float("nan")

    domination_index = float("nan")
    if np.isfinite(sum_top) and sum_top > 0 and len(top_vals) > 0:
        domination_index = float(top_vals[0] / sum_top)

    # nonzero ratio per feature among topk (flatten N and K)
    nz = np.abs(shap_n_f_k[:, top_idx, :]) > 1e-8
    nonzero_ratio = np.mean(nz, axis=(0, 2))  # (topk,)

    # summary stats for topk mean_abs_shap
    top_stats = {
        "mean": float(np.mean(top_vals)) if len(top_vals) else float("nan"),
        "std": float(np.std(top_vals, ddof=1)) if len(top_vals) >= 2 else float("nan"),
        "min": float(np.min(top_vals)) if len(top_vals) else float("nan"),
        "max": float(np.max(top_vals)) if len(top_vals) else float("nan"),
    }

    report: dict[str, Any] = {
        "n_used": n_used,
        "n_features": n_features,
        "n_topk": n_topk,
        "seed": int(seed),
        "task": str(task),
        "topk_mean_abs_shap_stats": top_stats,
        "domination_index": domination_index,
    }
    if task == "regime" and class_distribution is not None:
        report["class_distribution"] = class_distribution

    # write CSV per-feature for topk
    rows = []
    for i, f_idx in enumerate(top_idx):
        f = str(feature_names[int(f_idx)])
        rows.append(
            {
                "feature": f,
                "mean_abs_shap": float(mean_abs[int(f_idx)]),
                "nonzero_ratio": float(nonzero_ratio[i]),
                "domination_contrib": float(mean_abs[int(f_idx)] / sum_top) if np.isfinite(sum_top) and sum_top > 0 else float("nan"),
            }
        )
    pd.DataFrame(rows).to_csv(outdir / "shap_richness_report.csv", index=False)

    # Nonzero ratios full vector for JSON (topk only, compact)
    report["nonzero_shap_ratio_topk"] = {str(feature_names[int(i)]): float(r) for i, r in zip(top_idx, nonzero_ratio)}

    (outdir / "shap_richness_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    # Per-class importance (regime)
    if task == "regime" and class_names is not None and k == len(class_names):
        recs = []
        # mean(|shap|) over samples per class -> (F,) for each class
        for ci, cname in enumerate(class_names):
            imp_c = np.mean(np.abs(shap_n_f_k[:, :, ci]), axis=0)
            ord_c = np.argsort(imp_c)[::-1]
            for rank, j in enumerate(ord_c[:10], start=1):
                recs.append(
                    {
                        "class_name": str(cname),
                        "feature": str(feature_names[int(j)]),
                        "mean_abs_shap_class": float(imp_c[int(j)]),
                        "rank_in_class": int(rank),
                    }
                )
        pd.DataFrame(recs).to_csv(outdir / "shap_by_class_importance.csv", index=False)

    # warning propagation into run_metadata.json
    if np.isfinite(domination_index) and domination_index >= 0.7:
        meta_path = outdir / "run_metadata.json"
        if meta_path.exists():
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
            except Exception:  # noqa: BLE001
                meta = {}
            meta["warning"] = f"SHAP domination_index >= 0.7 (domination_index={domination_index:.3g})"
            meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return domination_index


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Merge with existing JSON if present (so warnings from other steps are not overwritten).
    merged: dict[str, Any] = {}
    if path.exists():
        try:
            merged = json.loads(path.read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001
            merged = {}
    merged.update(payload)
    path.write_text(json.dumps(merged, indent=2), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # deps
    try:
        import shap
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Missing dependency: shap ({e})", file=sys.stderr)
        return 2
    try:
        from xgboost import XGBClassifier, XGBRegressor
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Missing dependency: xgboost ({e})", file=sys.stderr)
        return 2

    labels_path = Path(args.labels)
    features_path = Path(args.features)
    labels = _read_parquet(labels_path)
    features = _read_parquet(features_path)
    df = _merge_labels_features(labels, features, features_path=features_path)

    # Parse drop-features
    drop_features: set[str] = set()
    if str(args.drop_features).strip():
        for tok in re.split(r"[,\s]+", str(args.drop_features).strip()):
            if tok:
                drop_features.add(tok.strip())
    if args.task == "severity" and bool(args.severity_drop_objective):
        drop_features.add("objective_value")

    # Train task
    if args.task == "regime":
        from sklearn.metrics import classification_report
        from sklearn.model_selection import train_test_split
        from sklearn.preprocessing import LabelEncoder

        X, y_raw, sample_ids = _prepare_X_y_regime(df, drop_features=drop_features)
        le = LabelEncoder()
        y = le.fit_transform(y_raw)
        classes = le.classes_.tolist()

        X_train, X_test, y_train, y_test = train_test_split(
            X,
            y,
            test_size=float(args.test_size),
            random_state=int(args.seed),
            stratify=y,
        )

        # class weights (balanced)
        sample_weight = None
        if args.class_weight == "balanced":
            counts = np.bincount(y_train)
            weights = {i: (len(y_train) / (len(counts) * max(1, int(c)))) for i, c in enumerate(counts)}
            sample_weight = np.asarray([weights[int(c)] for c in y_train], dtype=float)

        model = XGBClassifier(
            objective="multi:softprob",
            eval_metric="mlogloss",
            num_class=int(len(classes)),
            max_depth=4,
            n_estimators=800,
            learning_rate=0.03,
            subsample=0.8,
            colsample_bytree=0.8,
            n_jobs=int(args.n_jobs),
            random_state=int(args.seed),
            verbosity=0,
        )
        model.fit(X_train, y_train, sample_weight=sample_weight)

        y_pred = model.predict(X_test)
        rep = classification_report(y_test, y_pred, target_names=classes)
        (outdir / "classification_report.txt").write_text(rep + "\n", encoding="utf-8")

        cm = pd.crosstab(
            pd.Series(le.inverse_transform(y_test), name="true"),
            pd.Series(le.inverse_transform(y_pred), name="pred"),
            dropna=False,
        )
        for lab in classes:
            if lab not in cm.index:
                cm.loc[lab] = 0
            if lab not in cm.columns:
                cm[lab] = 0
        cm = cm.loc[classes, classes]
        cm.to_csv(outdir / "confusion_matrix.csv")

        pd.DataFrame({"label_int": list(range(len(classes))), "label": classes}).to_csv(
            outdir / "label_mapping.csv", index=False
        )

        model.save_model(str(outdir / "model.json"))

        # SHAP
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X_test)
        global_importance, dep_vals = _shap_global_importance_and_dep_values(shap_values)
        shap_n_f_k = _shap_to_n_f_k(shap_values, n_features=len(X_test.columns))

        shap_imp = (
            pd.DataFrame({"feature": list(X_test.columns), "mean_abs_shap": global_importance})
            .sort_values("mean_abs_shap", ascending=False)
            .reset_index(drop=True)
        )
        shap_imp.to_csv(outdir / "shap_importance.csv", index=False)

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # noqa: E402

        plt.figure(figsize=(10, 6))
        shap.summary_plot(shap_values, X_test, show=False, max_display=int(args.topk))
        plt.tight_layout()
        plt.savefig(outdir / "shap_beeswarm.png", dpi=200, bbox_inches="tight")
        plt.close()

        plt.figure(figsize=(10, 6))
        shap.summary_plot(shap_values, X_test, show=False, plot_type="bar", max_display=int(args.topk))
        plt.tight_layout()
        plt.savefig(outdir / "shap_bar.png", dpi=200, bbox_inches="tight")
        plt.close()

        if not args.skip_dependence:
            top3_idx = np.argsort(global_importance)[-3:][::-1]
            top3_features = [str(X_test.columns[int(i)]) for i in top3_idx]
            if dep_vals is None:
                logger.warning("Skipping dependence plots (no dep_vals)")
            else:
                for feat in top3_features:
                    safe = _safe_filename(feat)
                    plt.figure(figsize=(7, 4.5))
                    shap.dependence_plot(feat, dep_vals, X_test, show=False)
                    plt.tight_layout()
                    plt.savefig(outdir / f"shap_dependence_{safe}.png", dpi=200, bbox_inches="tight")
                    plt.close()

        # Per-class beeswarm (subset by true class label) + richness reports
        class_dist = pd.Series(y_raw).value_counts().to_dict()
        domination_index = _write_shap_richness_reports(
            outdir=outdir,
            task="regime",
            seed=int(args.seed),
            topk=int(args.topk),
            feature_names=list(X_test.columns),
            shap_n_f_k=shap_n_f_k,
            class_names=classes,
            class_distribution={str(k): int(v) for k, v in class_dist.items()},
        )

        # class-specific beeswarm: use SHAP for that class output, but subset rows where true label == class
        for ci, cname in enumerate(classes):
            mask = (le.inverse_transform(y_test) == cname)
            if not np.any(mask):
                continue
            sv_c = shap_n_f_k[mask, :, ci]
            X_c = X_test.loc[mask].copy()
            plt.figure(figsize=(10, 6))
            shap.summary_plot(sv_c, X_c, show=False, max_display=int(args.topk))
            plt.tight_layout()
            safe_c = _safe_filename(cname, max_len=80)
            plt.savefig(outdir / f"shap_beeswarm_class_{safe_c}.png", dpi=200, bbox_inches="tight")
            plt.close()

        _write_json(
            outdir / "run_metadata.json",
            {
                "task": "regime",
                "labels": str(labels_path),
                "features": str(features_path),
                "n_samples_used": int(len(X)),
                "class_distribution": pd.Series(y_raw).value_counts().to_dict(),
                "seed": int(args.seed),
                "test_size": float(args.test_size),
                "topk": int(args.topk),
                "n_jobs": int(args.n_jobs),
                "class_weight": str(args.class_weight),
                "skip_dependence": bool(args.skip_dependence),
                "drop_features": sorted(drop_features),
                "domination_index": domination_index,
            },
        )
        print(f"[OK] Wrote outputs under: {outdir}")
        return 0

    # severity regression
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    from sklearn.model_selection import train_test_split

    X, y, _sample_ids = _prepare_X_y_severity(df, drop_features=drop_features)
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=float(args.test_size),
        random_state=int(args.seed),
    )

    model = XGBRegressor(
        objective="reg:squarederror",
        eval_metric="rmse",
        max_depth=4,
        n_estimators=800,
        learning_rate=0.03,
        subsample=0.8,
        colsample_bytree=0.8,
        n_jobs=int(args.n_jobs),
        random_state=int(args.seed),
        verbosity=0,
    )
    model.fit(X_train, y_train)

    pred = model.predict(X_test)
    r2 = float(r2_score(y_test, pred))
    mae = float(mean_absolute_error(y_test, pred))
    rmse = float(np.sqrt(mean_squared_error(y_test, pred)))
    pd.DataFrame([{"r2": r2, "mae": mae, "rmse": rmse, "n_train": len(y_train), "n_test": len(y_test)}]).to_csv(
        outdir / "regression_metrics.csv", index=False
    )

    model.save_model(str(outdir / "model.json"))

    import shap
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402

    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_test)
    global_importance, dep_vals = _shap_global_importance_and_dep_values(shap_values)
    shap_n_f_k = _shap_to_n_f_k(shap_values, n_features=len(X_test.columns))

    shap_imp = (
        pd.DataFrame({"feature": list(X_test.columns), "mean_abs_shap": global_importance})
        .sort_values("mean_abs_shap", ascending=False)
        .reset_index(drop=True)
    )
    shap_imp.to_csv(outdir / "shap_importance.csv", index=False)

    plt.figure(figsize=(10, 6))
    shap.summary_plot(shap_values, X_test, show=False, max_display=int(args.topk))
    plt.tight_layout()
    plt.savefig(outdir / "shap_beeswarm.png", dpi=200, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(10, 6))
    shap.summary_plot(shap_values, X_test, show=False, plot_type="bar", max_display=int(args.topk))
    plt.tight_layout()
    plt.savefig(outdir / "shap_bar.png", dpi=200, bbox_inches="tight")
    plt.close()

    if not args.skip_dependence:
        top3_idx = np.argsort(global_importance)[-3:][::-1]
        top3_features = [str(X_test.columns[int(i)]) for i in top3_idx]
        if dep_vals is None:
            logger.warning("Skipping dependence plots (no dep_vals)")
        else:
            for feat in top3_features:
                safe = _safe_filename(feat)
                plt.figure(figsize=(7, 4.5))
                shap.dependence_plot(feat, dep_vals, X_test, show=False)
                plt.tight_layout()
                plt.savefig(outdir / f"shap_dependence_{safe}.png", dpi=200, bbox_inches="tight")
                plt.close()

    domination_index = _write_shap_richness_reports(
        outdir=outdir,
        task="severity",
        seed=int(args.seed),
        topk=int(args.topk),
        feature_names=list(X_test.columns),
        shap_n_f_k=shap_n_f_k,
    )

    if np.isfinite(domination_index) and domination_index >= 0.7:
        msg = f"[WARN] severity SHAP domination_index>=0.7 (domination_index={domination_index:.3g})"
        if bool(args.severity_drop_objective):
            msg += " while --severity-drop-objective is ON (compare with OFF to assess bias)."
        else:
            msg += " (try --severity-drop-objective to see if domination decreases)."
        print(msg)

    _write_json(
        outdir / "run_metadata.json",
        {
            "task": "severity",
            "labels": str(labels_path),
            "features": str(features_path),
            "n_samples_used": int(len(X)),
            "seed": int(args.seed),
            "test_size": float(args.test_size),
            "topk": int(args.topk),
            "n_jobs": int(args.n_jobs),
            "skip_dependence": bool(args.skip_dependence),
            "drop_features": sorted(drop_features),
            "severity_drop_objective": bool(args.severity_drop_objective),
            "domination_index": domination_index,
        },
    )
    print(f"[OK] Wrote outputs under: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

