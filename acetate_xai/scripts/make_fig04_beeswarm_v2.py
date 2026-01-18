from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def _find_results_root(override: str | None) -> Path:
    if override:
        p = Path(override)
        if not p.exists():
            raise FileNotFoundError(f"--results-root not found: {p}")
        return p
    for p in [Path("results_all_v2_extracted") / "results", Path("results")]:
        if p.exists() and p.is_dir():
            return p
    raise FileNotFoundError("No results root found (tried results_all_v2_extracted/results and results/).")


def _pretty_feature_name(col: str) -> str:
    """
    Keep mapping consistent with Fig03 v2 script.
    """
    rxn_map = {
        "EX_h2o_e": "Water transport",
        "EX_h_e": "Proton exchange",
        "EX_o2_e": "Oxygen uptake",
        "EX_nh4_e": "Ammonium uptake",
        "EX_pi_e": "Phosphate uptake",
        "EX_co2_e": "CO2 exchange",
    }

    prefix = None
    rid = col
    for pfx in ["width__", "mid__", "signchange__"]:
        if col.startswith(pfx):
            prefix = pfx[:-2]  # width/mid/signchange
            rid = col[len(pfx) :]
            break

    base = rxn_map.get(rid, rid)
    if prefix == "width":
        return f"{base} flexibility"
    if prefix == "mid":
        return f"{base} flux midpoint"
    if prefix == "signchange":
        return f"{base} sign change"
    return col


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Rebuild Fig04 beeswarm (v2) using shap.plots.beeswarm for severity regression.")
    p.add_argument("--results-root", default=None, help="Results root containing xai_xgb_maintenance/ and regime_dataset.parquet")
    p.add_argument("--xgb-dir", default="xai_xgb_maintenance", help="Folder under results-root (default: xai_xgb_maintenance)")
    p.add_argument("--data", default="regime_dataset.parquet", help="Parquet under results-root (default: regime_dataset.parquet)")
    p.add_argument("--run-col", default="run_id", help="Run id column (default: run_id)")
    p.add_argument("--objective-col", default="objective_value", help="Objective column (default: objective_value)")
    p.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    p.add_argument("--test-size", type=float, default=0.2, help="Test size (default: 0.2)")
    p.add_argument("--max-display", type=int, default=15, help="max_display for beeswarm (default: 15)")
    p.add_argument(
        "--out",
        default=str(Path("results") / "figures_draft" / "Fig04_shap_severity_beeswarm_v2.png"),
        help="Output PNG path",
    )
    args = p.parse_args(argv)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import shap
    import xgboost as xgb
    from xgboost import XGBRegressor
    from sklearn.model_selection import train_test_split

    results_root = _find_results_root(args.results_root)
    xgb_dir = results_root / args.xgb_dir
    model_json = xgb_dir / "model.json"
    if not model_json.exists():
        raise FileNotFoundError(f"Missing model.json: {model_json}")

    data_path = results_root / args.data
    if not data_path.exists():
        raise FileNotFoundError(f"Missing data parquet: {data_path}")

    df = pd.read_parquet(data_path)

    feat_cols = [c for c in df.columns if c.startswith(("width__", "mid__", "signchange__"))]
    if not feat_cols:
        raise ValueError("No feature columns found (expected width__/mid__/signchange__).")

    if args.run_col not in df.columns:
        raise ValueError(f"Missing run column: {args.run_col}")
    if args.objective_col not in df.columns:
        raise ValueError(f"Missing objective column: {args.objective_col}")

    # Build severity target
    obj = pd.to_numeric(df[args.objective_col], errors="coerce")
    run = df[args.run_col].astype(str)
    max_by_run = obj.groupby(run).transform("max")
    y = obj / max_by_run
    y = y.replace([np.inf, -np.inf], np.nan)

    X = df[feat_cols].copy()
    for c in X.columns:
        if X[c].dtype == bool:
            X[c] = X[c].astype(int)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    mask = y.notna()
    X = X.loc[mask].copy()
    y = y.loc[mask].astype(float)

    # For plotting, we use the test split to mirror training scripts.
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y.to_numpy(),
        test_size=float(args.test_size),
        random_state=int(args.seed),
    )

    # Load trained model
    model = XGBRegressor()
    model.load_model(str(model_json))

    # Compute SHAP values via xgboost pred_contribs (robust to SHAP TreeExplainer JSON parsing issues)
    booster = model.get_booster()
    dmat = xgb.DMatrix(X_test, feature_names=list(X_test.columns))
    try:
        contrib = booster.predict(dmat, pred_contribs=True, strict_shape=True)
    except TypeError:
        contrib = booster.predict(dmat, pred_contribs=True)
    arr = np.asarray(contrib)
    # pred_contribs shapes:
    #  - regression: (N, F+1)
    #  - strict_shape regression: (N, 1, F+1)
    if arr.ndim == 3 and arr.shape[1] == 1:
        arr = arr[:, 0, :]
    if arr.ndim != 2 or arr.shape[1] != (len(X_test.columns) + 1):
        raise TypeError(f"Unexpected pred_contribs shape: {arr.shape}")
    values = arr[:, :-1]  # drop bias

    pretty_names = [_pretty_feature_name(c) for c in X_test.columns]
    expl = shap.Explanation(values=values, data=X_test.to_numpy(), feature_names=pretty_names)

    # Paper style: landscape + larger fonts
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    )
    plt.figure(figsize=(10, 6))
    shap.plots.beeswarm(expl, max_display=int(args.max_display), show=False, plot_size=None)

    out_png = Path(args.out)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=260, bbox_inches="tight")
    plt.close("all")

    print(f"[OK] Wrote: {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

