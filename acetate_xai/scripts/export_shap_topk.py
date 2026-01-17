from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Export SHAP top-k features for causal discovery (regime + severity).")
    p.add_argument("--data", default="results/causal_dataset.parquet", help="Input causal dataset parquet")
    p.add_argument("--outdir", default="results/causal_topk", help="Output directory")
    p.add_argument("--topk", type=int, default=15, help="Top-K features to export")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def _feature_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith(("width__", "mid__", "signchange__"))]


def _safe_filename(s: str, max_len: int = 120) -> str:
    s2 = re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_")
    if not s2:
        s2 = "feature"
    return s2[:max_len]


def _mean_abs_shap_multiclass(shap_values) -> np.ndarray:
    # shap_values for multi-class is often list[C] of (N,F)
    if isinstance(shap_values, list) and shap_values:
        arr = np.stack([np.asarray(v) for v in shap_values], axis=0)  # (C,N,F)
        return np.mean(np.abs(arr), axis=(0, 1))  # (F,)
    arr = np.asarray(shap_values)
    return np.mean(np.abs(arr), axis=0)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    data_path = Path(args.data)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        from xgboost import XGBClassifier, XGBRegressor
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Missing dependency: xgboost ({e})")
        print("Install on server: pip install xgboost")
        return 2

    try:
        import shap
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Missing dependency: shap ({e})")
        print("Install on server: pip install shap")
        return 2

    if not data_path.exists():
        print(f"[ERROR] data not found: {data_path}")
        return 2

    df = pd.read_parquet(data_path)
    for c in ("primary_regime", "maintenance_severity"):
        if c not in df.columns:
            print(f"[ERROR] missing required column: {c}")
            return 2

    feat_cols = _feature_cols(df)
    if not feat_cols:
        print("[ERROR] No feature columns found (expected width__/mid__/signchange__).")
        return 2

    X = df[feat_cols].copy()
    for c in X.columns:
        if X[c].dtype == bool:
            X[c] = X[c].astype(int)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # 1) Classification: primary_regime
    y_class_raw = df["primary_regime"].astype(str).str.strip()
    mask_c = y_class_raw.notna() & (y_class_raw != "") & (y_class_raw.str.lower() != "nan")
    Xc = X.loc[mask_c].copy()
    yc_raw = y_class_raw.loc[mask_c].copy()

    from sklearn.preprocessing import LabelEncoder

    le = LabelEncoder()
    yc = le.fit_transform(yc_raw.values)
    classes = le.classes_.tolist()

    clf = XGBClassifier(
        objective="multi:softprob",
        eval_metric="mlogloss",
        num_class=int(len(classes)),
        max_depth=4,
        n_estimators=500,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        n_jobs=-1,
        random_state=int(args.seed),
        verbosity=0,
    )
    clf.fit(Xc, yc)
    expl_c = shap.TreeExplainer(clf)
    shap_c = expl_c.shap_values(Xc)
    imp_c = _mean_abs_shap_multiclass(shap_c)

    imp_c_df = pd.DataFrame({"feature": Xc.columns, "mean_abs_shap": imp_c}).sort_values(
        "mean_abs_shap", ascending=False
    )
    imp_c_df.to_csv(outdir / "shap_importance_regime.csv", index=False)
    imp_c_df.head(int(args.topk)).to_csv(outdir / "top_features_regime.csv", index=False)

    # 2) Regression: maintenance_severity
    y_reg = pd.to_numeric(df["maintenance_severity"], errors="coerce")
    mask_r = y_reg.notna() & np.isfinite(y_reg.values)
    Xr = X.loc[mask_r].copy()
    yr = y_reg.loc[mask_r].astype(float).values

    reg = XGBRegressor(
        objective="reg:squarederror",
        eval_metric="rmse",
        max_depth=4,
        n_estimators=800,
        learning_rate=0.03,
        subsample=0.8,
        colsample_bytree=0.8,
        n_jobs=-1,
        random_state=int(args.seed),
        verbosity=0,
    )
    reg.fit(Xr, yr)
    expl_r = shap.TreeExplainer(reg)
    shap_r = expl_r.shap_values(Xr)
    imp_r = np.mean(np.abs(np.asarray(shap_r)), axis=0)

    imp_r_df = pd.DataFrame({"feature": Xr.columns, "mean_abs_shap": imp_r}).sort_values(
        "mean_abs_shap", ascending=False
    )
    imp_r_df.to_csv(outdir / "shap_importance_severity.csv", index=False)
    imp_r_df.head(int(args.topk)).to_csv(outdir / "top_features_severity.csv", index=False)

    print(f"[INFO] regime classes={len(classes)} rows={len(Xc)} features={len(feat_cols)}")
    print(f"[INFO] severity rows={len(Xr)} features={len(feat_cols)}")
    print(f"[OK] Wrote: {outdir / 'top_features_regime.csv'}")
    print(f"[OK] Wrote: {outdir / 'top_features_severity.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

