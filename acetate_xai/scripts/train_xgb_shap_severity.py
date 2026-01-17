from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Train XGBoost regressor for maintenance severity + SHAP explanations.")
    p.add_argument("--data", default="results/regime_dataset.parquet", help="Input parquet path")
    p.add_argument("--outdir", default="results/xai_xgb_maintenance", help="Output directory")
    p.add_argument("--objective-col", default="objective_value", help="Objective column name")
    p.add_argument("--run-col", default="run_id", help="Run identifier column name")
    p.add_argument("--topk", type=int, default=20, help="Top-K features for SHAP bar plot focus")
    p.add_argument("--test-size", type=float, default=0.2, help="Test size (default 0.2)")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    p.add_argument("--max-depth", type=int, default=4, help="XGBoost tree depth")
    p.add_argument("--n-estimators", type=int, default=800, help="Number of boosting rounds")
    p.add_argument("--learning-rate", type=float, default=0.03, help="Learning rate")
    p.add_argument("--subsample", type=float, default=0.8, help="Subsample")
    p.add_argument("--colsample-bytree", type=float, default=0.8, help="Colsample bytree")
    p.add_argument("--n-jobs", type=int, default=-1, help="Parallel jobs for XGBoost (-1 = all)")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def _feature_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith(("width__", "mid__", "signchange__"))]


def _safe_filename(s: str, max_len: int = 120) -> str:
    s2 = re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_")
    if not s2:
        s2 = "feature"
    return s2[:max_len]


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    data_path = Path(args.data)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        from xgboost import XGBRegressor
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

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402

    if not data_path.exists():
        print(f"[ERROR] data not found: {data_path}")
        return 2

    df = pd.read_parquet(data_path)
    for c in (args.objective_col, args.run_col):
        if c not in df.columns:
            print(f"[ERROR] missing required column: {c}")
            return 2

    # severity construction
    df[args.objective_col] = pd.to_numeric(df[args.objective_col], errors="coerce")
    run_max = df.groupby(args.run_col)[args.objective_col].max().rename("max_objective_in_same_run")
    df = df.merge(run_max.reset_index(), on=args.run_col, how="left")
    den = df["max_objective_in_same_run"].replace({0.0: np.nan})
    df["severity"] = df[args.objective_col] / den

    # Clean severity: keep finite and within [0, 1] (with small tolerance)
    sev = pd.to_numeric(df["severity"], errors="coerce")
    ok = sev.notna() & np.isfinite(sev.values) & (sev >= -1e-12) & (sev <= 1.0 + 1e-12)
    df = df.loc[ok].copy()
    df["severity"] = df["severity"].clip(lower=0.0, upper=1.0)

    feat_cols = _feature_cols(df)
    if not feat_cols:
        print("[ERROR] No feature columns found (expected width__/mid__/signchange__).")
        return 2

    X = df[feat_cols].copy()
    for c in X.columns:
        if X[c].dtype == bool:
            X[c] = X[c].astype(int)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    y = pd.to_numeric(df["severity"], errors="coerce").astype(float).values

    from sklearn.model_selection import train_test_split
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=float(args.test_size),
        random_state=int(args.seed),
    )

    print(f"[INFO] data rows={len(df)}, features={len(feat_cols)}")
    print(
        f"[INFO] severity distribution: min={float(np.min(y)):.6g} mean={float(np.mean(y)):.6g} max={float(np.max(y)):.6g}"
    )

    model = XGBRegressor(
        objective="reg:squarederror",
        eval_metric="rmse",
        max_depth=int(args.max_depth),
        n_estimators=int(args.n_estimators),
        learning_rate=float(args.learning_rate),
        subsample=float(args.subsample),
        colsample_bytree=float(args.colsample_bytree),
        n_jobs=int(args.n_jobs),
        random_state=int(args.seed),
        verbosity=0,
    )
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    rmse = float(np.sqrt(mean_squared_error(y_test, y_pred)))
    metrics = pd.DataFrame(
        [
            {
                "n_train": int(len(X_train)),
                "n_test": int(len(X_test)),
                "r2": float(r2_score(y_test, y_pred)),
                "mae": float(mean_absolute_error(y_test, y_pred)),
                "rmse": rmse,
            }
        ]
    )
    metrics.to_csv(outdir / "regression_metrics.csv", index=False)

    model.save_model(str(outdir / "model.json"))

    # SHAP
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_test)

    # Beeswarm
    plt.figure()
    shap.summary_plot(shap_values, X_test, show=False, max_display=int(args.topk))
    plt.tight_layout()
    plt.savefig(outdir / "shap_beeswarm.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Bar
    plt.figure()
    shap.summary_plot(shap_values, X_test, show=False, plot_type="bar", max_display=int(args.topk))
    plt.tight_layout()
    plt.savefig(outdir / "shap_bar.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Dependence: top-3 by mean(|SHAP|)
    mean_abs = np.mean(np.abs(np.asarray(shap_values)), axis=0)
    top3_idx = np.argsort(-mean_abs)[:3]
    top3_features = [X_test.columns[int(i)] for i in top3_idx]
    for feat in top3_features:
        safe = _safe_filename(feat)
        plt.figure()
        shap.dependence_plot(feat, shap_values, X_test, show=False)
        plt.tight_layout()
        plt.savefig(outdir / f"shap_dependence_{safe}.png", dpi=200, bbox_inches="tight")
        plt.close()

    # Console logs
    print("[INFO] regression metrics:")
    print(metrics.to_string(index=False))
    print(f"[OK] Wrote: {outdir / 'regression_metrics.csv'}")
    print(f"[OK] Wrote: {outdir / 'model.json'}")
    print(f"[OK] Wrote: {outdir / 'shap_beeswarm.png'}")
    print(f"[OK] Wrote: {outdir / 'shap_bar.png'}")
    for feat in top3_features:
        print(f"[OK] Wrote: {outdir / ('shap_dependence_' + _safe_filename(feat) + '.png')}")

    meta = {
        "data": str(data_path),
        "rows": int(len(df)),
        "n_features": int(len(feat_cols)),
        "severity": {"objective_col": args.objective_col, "run_col": args.run_col},
        "params": {
            "max_depth": int(args.max_depth),
            "n_estimators": int(args.n_estimators),
            "learning_rate": float(args.learning_rate),
            "subsample": float(args.subsample),
            "colsample_bytree": float(args.colsample_bytree),
            "seed": int(args.seed),
        },
    }
    (outdir / "run_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

