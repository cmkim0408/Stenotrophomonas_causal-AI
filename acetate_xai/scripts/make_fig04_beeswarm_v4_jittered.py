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
    p = argparse.ArgumentParser(description="Fig04 v4: jittered swarm-style SHAP plot for visual density (no data edits).")
    p.add_argument("--results-root", default=None, help="Results root containing xai_xgb_maintenance/ and regime_dataset.parquet")
    p.add_argument("--xgb-dir", default="xai_xgb_maintenance", help="Folder under results-root (default: xai_xgb_maintenance)")
    p.add_argument("--data", default="regime_dataset.parquet", help="Parquet under results-root (default: regime_dataset.parquet)")
    p.add_argument("--run-col", default="run_id", help="Run id column (default: run_id)")
    p.add_argument("--objective-col", default="objective_value", help="Objective column (default: objective_value)")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    p.add_argument("--test-size", type=float, default=0.2, help="Test size (default 0.2)")
    p.add_argument("--max-display", type=int, default=12, help="Top features to display (default 12)")
    p.add_argument("--alpha", type=float, default=0.5, help="Point alpha (default 0.5)")
    p.add_argument("--s", type=float, default=40, help="Point size (default 40)")
    p.add_argument("--jitter", type=float, default=0.42, help="Base y-jitter strength (default 0.42)")
    p.add_argument("--xlim-lo", type=float, default=1.0, help="Lower percentile for xlim (default 1)")
    p.add_argument("--xlim-hi", type=float, default=99.0, help="Upper percentile for xlim (default 99)")
    p.add_argument(
        "--out",
        default=str(Path("results") / "figures_draft" / "Fig04_shap_severity_beeswarm_v4_jittered.png"),
        help="Output PNG path",
    )
    args = p.parse_args(argv)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
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

    # Severity target (not used for plotting, but keeps consistent split filtering)
    obj = pd.to_numeric(df[args.objective_col], errors="coerce")
    run = df[args.run_col].astype(str)
    max_by_run = obj.groupby(run).transform("max")
    y = (obj / max_by_run).replace([np.inf, -np.inf], np.nan)

    X = df[feat_cols].copy()
    for c in X.columns:
        if X[c].dtype == bool:
            X[c] = X[c].astype(int)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    mask = y.notna()
    X = X.loc[mask].copy()
    y = y.loc[mask].astype(float)

    _X_train, X_test, _y_train, _y_test = train_test_split(
        X,
        y.to_numpy(),
        test_size=float(args.test_size),
        random_state=int(args.seed),
    )

    # Load trained regressor + SHAP via pred_contribs
    model = XGBRegressor()
    model.load_model(str(model_json))
    booster = model.get_booster()
    dmat = xgb.DMatrix(X_test, feature_names=list(X_test.columns))
    try:
        contrib = booster.predict(dmat, pred_contribs=True, strict_shape=True)
    except TypeError:
        contrib = booster.predict(dmat, pred_contribs=True)
    arr = np.asarray(contrib)
    if arr.ndim == 3 and arr.shape[1] == 1:
        arr = arr[:, 0, :]
    if arr.ndim != 2 or arr.shape[1] != (len(X_test.columns) + 1):
        raise TypeError(f"Unexpected pred_contribs shape: {arr.shape}")
    shap_vals = arr[:, :-1]  # (N,F)

    # Choose top features by mean(|SHAP|)
    mean_abs = np.mean(np.abs(shap_vals), axis=0)
    order = np.argsort(mean_abs)[::-1]
    topk = int(min(max(1, args.max_display), len(order)))
    idx = order[:topk]

    # Aesthetics: larger fonts + landscape
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    )
    fig, ax = plt.subplots(figsize=(10, 6))

    rng = np.random.default_rng(int(args.seed))
    cmap = plt.get_cmap("coolwarm")

    # Swarm-like: stronger jitter near center, weaker at extremes to create oval density.
    for row, j in enumerate(idx):
        x = shap_vals[:, j].astype(float)
        fv = X_test.iloc[:, int(j)].to_numpy().astype(float)

        # Normalize feature values for color (robust percentiles)
        lo_f = np.nanpercentile(fv, 5.0) if np.isfinite(fv).any() else 0.0
        hi_f = np.nanpercentile(fv, 95.0) if np.isfinite(fv).any() else 1.0
        denom = (hi_f - lo_f) if (hi_f > lo_f) else 1.0
        cval = np.clip((fv - lo_f) / denom, 0.0, 1.0)

        # Center-weighted jitter magnitude: larger near 0, smaller at tails (oval look)
        # scale by robust x-range to be stable against outliers
        x_lo = np.percentile(x, 1.0)
        x_hi = np.percentile(x, 99.0)
        x_scale = max(1e-9, (x_hi - x_lo) / 2.0)
        w = np.exp(-0.5 * (x / x_scale) ** 2)  # ~1 at center, ->0 at tails
        sigma = float(args.jitter) * (0.15 + 0.85 * w)
        yj = row + rng.normal(0.0, sigma, size=len(x))

        ax.scatter(
            x,
            yj,
            s=float(args.s),
            alpha=float(args.alpha),
            c=cval,
            cmap=cmap,
            edgecolors="none",
            rasterized=True,
        )

    # Y labels
    names = [_pretty_feature_name(X_test.columns[int(j)]) for j in idx]
    ax.set_yticks(range(topk))
    ax.set_yticklabels(names)
    ax.invert_yaxis()

    ax.set_xlabel("SHAP value (impact on maintenance severity)")
    ax.set_title("Maintenance severity: SHAP swarm (jittered for visual density)")
    ax.grid(True, axis="x", alpha=0.18)

    # Percentile zoom-in on x (style-only)
    flat = shap_vals[:, idx].reshape(-1)
    flat = flat[np.isfinite(flat)]
    if flat.size > 0:
        lo = float(np.percentile(flat, float(args.xlim_lo)))
        hi = float(np.percentile(flat, float(args.xlim_hi)))
        if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
            pad = 0.03 * (hi - lo)
            ax.set_xlim(lo - pad, hi + pad)

    # Colorbar (minimal)
    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Feature value (low→high)", fontsize=11)

    fig.tight_layout()
    out_png = Path(args.out)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] Wrote: {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

