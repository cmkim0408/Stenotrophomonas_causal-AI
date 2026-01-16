from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd


def _feature_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith(("width__", "mid__", "signchange__"))]


def _safe_div(num: pd.Series, den: pd.Series) -> pd.Series:
    den2 = den.replace({0.0: np.nan})
    return num / den2


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Train maintenance/severity regression XAI models (tree + elasticnet).")
    p.add_argument("--data", required=True, help="Input parquet (e.g., results/regime_dataset.parquet)")
    p.add_argument("--outdir", required=True, help="Output directory (e.g., results/xai_platform_maintenance)")
    p.add_argument("--tree-depth", type=int, default=4, help="Decision tree max_depth")
    p.add_argument("--test-size", type=float, default=0.2, help="Test split size (default 0.2)")
    p.add_argument("--random-state", type=int, default=0, help="Random state")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    # Avoid oversubscription on servers
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    data_path = Path(args.data)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(data_path)
    if "run_id" not in df.columns:
        raise ValueError("Dataset missing required column: run_id")
    if "objective_value" not in df.columns:
        raise ValueError("Dataset missing required column: objective_value")

    # severity = objective_value / max(objective_value) within run_id
    df["objective_value"] = pd.to_numeric(df["objective_value"], errors="coerce")
    run_max = df.groupby("run_id")["objective_value"].max().rename("objective_value_max")
    df = df.merge(run_max.reset_index(), on="run_id", how="left")
    df["severity"] = _safe_div(df["objective_value"], df["objective_value_max"])

    # Basic distribution output
    sev = df["severity"].dropna()
    dist = pd.DataFrame(
        {
            "count": [int(sev.shape[0])],
            "min": [float(sev.min()) if not sev.empty else np.nan],
            "p25": [float(sev.quantile(0.25)) if not sev.empty else np.nan],
            "median": [float(sev.quantile(0.5)) if not sev.empty else np.nan],
            "p75": [float(sev.quantile(0.75)) if not sev.empty else np.nan],
            "max": [float(sev.max()) if not sev.empty else np.nan],
            "mean": [float(sev.mean()) if not sev.empty else np.nan],
            "std": [float(sev.std()) if not sev.empty else np.nan],
        }
    )
    dist.to_csv(outdir / "severity_distribution.csv", index=False)

    feat_cols = _feature_columns(df)
    if not feat_cols:
        raise ValueError("No feature columns found (expected width__/mid__/signchange__).")

    train_df = df.loc[~df["severity"].isna()].copy()
    if len(train_df) < 5:
        raise ValueError("Not enough rows with computable severity to train (need >= 5).")

    X = train_df[feat_cols].copy()
    for c in X.columns:
        if X[c].dtype == bool:
            X[c] = X[c].astype(int)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    y = pd.to_numeric(train_df["severity"], errors="coerce").astype(float).values

    from sklearn.model_selection import train_test_split

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=float(args.test_size),
        random_state=int(args.random_state),
    )

    # 1) DecisionTreeRegressor
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    from sklearn.tree import DecisionTreeRegressor, export_text

    tree = DecisionTreeRegressor(max_depth=int(args.tree_depth), random_state=int(args.random_state))
    tree.fit(X_train, y_train)
    pred_tree = tree.predict(X_test)
    (outdir / "tree_rules.txt").write_text(export_text(tree, feature_names=list(X.columns)) + "\n", encoding="utf-8")

    # 2) ElasticNetCV (with scaling)
    from sklearn.linear_model import ElasticNetCV
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    cv = min(5, len(X_train))
    enet = Pipeline(
        steps=[
            ("scaler", StandardScaler(with_mean=True, with_std=True)),
            ("enet", ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9, 0.95, 0.99], cv=cv, random_state=int(args.random_state))),
        ]
    )
    enet.fit(X_train, y_train)
    pred_enet = enet.predict(X_test)
    coef = enet.named_steps["enet"].coef_
    coef_df = pd.DataFrame({"feature": X.columns, "coef": coef})
    coef_df["abs"] = coef_df["coef"].abs()
    coef_df = coef_df.sort_values("abs", ascending=False).head(30).drop(columns=["abs"]).reset_index(drop=True)
    coef_df.to_csv(outdir / "linear_top_features.csv", index=False)

    # Metrics
    def _metrics(name: str, yhat: np.ndarray) -> dict:
        rmse = float(np.sqrt(mean_squared_error(y_test, yhat)))
        return {
            "model": name,
            "n_train": int(len(X_train)),
            "n_test": int(len(X_test)),
            "r2": float(r2_score(y_test, yhat)),
            "mae": float(mean_absolute_error(y_test, yhat)),
            "rmse": rmse,
        }

    metrics = pd.DataFrame([_metrics("DecisionTreeRegressor", pred_tree), _metrics("ElasticNetCV", pred_enet)])
    metrics.to_csv(outdir / "regression_metrics.csv", index=False)

    print(f"[OK] Wrote: {outdir / 'severity_distribution.csv'}")
    print(f"[OK] Wrote: {outdir / 'regression_metrics.csv'}")
    print(f"[OK] Wrote: {outdir / 'tree_rules.txt'}")
    print(f"[OK] Wrote: {outdir / 'linear_top_features.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

