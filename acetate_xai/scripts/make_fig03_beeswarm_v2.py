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
    Convert raw feature column names (e.g., width__EX_h2o_e) into paper-friendly labels.
    Keeps this conservative to avoid wrong biology; uses simple reaction-ID mapping for common exchanges.
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
    p = argparse.ArgumentParser(description="Rebuild Fig03 beeswarm (v2) using shap.plots.beeswarm for multiclass.")
    p.add_argument("--results-root", default=None, help="Results root containing xai_xgb/ and regime_dataset.parquet")
    p.add_argument("--xgb-dir", default="xai_xgb", help="Folder under results-root (default: xai_xgb)")
    p.add_argument("--data", default="regime_dataset.parquet", help="Parquet under results-root (default: regime_dataset.parquet)")
    p.add_argument("--label-col", default="label", help="Label column for stratified split (default: label)")
    p.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    p.add_argument("--test-size", type=float, default=0.2, help="Test size (default: 0.2)")
    p.add_argument("--max-display", type=int, default=15, help="max_display for beeswarm (default: 15)")
    p.add_argument(
        "--multiclass-layout",
        choices=["single", "panels"],
        default="panels",
        help=(
            "How to render multiclass beeswarm. "
            "'panels' draws one beeswarm per class in a single figure (recommended; shap.plots.beeswarm expects 2D). "
            "'single' averages SHAP values over classes into a single 2D explanation."
        ),
    )
    p.add_argument("--debug", action="store_true", help="Print shapes for debugging.")
    p.add_argument(
        "--out",
        default=str(Path("results") / "figures_draft" / "Fig03_shap_regime_beeswarm_v2.png"),
        help="Output PNG path",
    )
    args = p.parse_args(argv)

    # heavy deps: import only inside main (for venv use)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import shap
    from sklearn.model_selection import train_test_split
    from xgboost import XGBClassifier

    results_root = _find_results_root(args.results_root)
    xgb_dir = results_root / args.xgb_dir
    model_json = xgb_dir / "model.json"
    if not model_json.exists():
        raise FileNotFoundError(f"Missing model.json: {model_json}")

    data_path = results_root / args.data
    if not data_path.exists():
        raise FileNotFoundError(f"Missing data parquet: {data_path}")

    df = pd.read_parquet(data_path)
    if args.label_col not in df.columns:
        raise ValueError(f"Label column not found: {args.label_col}")

    feat_cols = [c for c in df.columns if c.startswith(("width__", "mid__", "signchange__"))]
    if not feat_cols:
        raise ValueError("No feature columns found (expected width__/mid__/signchange__).")

    X = df[feat_cols].copy()
    for c in X.columns:
        if X[c].dtype == bool:
            X[c] = X[c].astype(int)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    y_raw = df[args.label_col].astype(str)
    mask = y_raw.notna() & (y_raw.str.strip() != "") & (y_raw.str.lower() != "nan")
    X = X.loc[mask].copy()
    y_raw = y_raw.loc[mask].astype(str).str.strip()

    # We only need X_test for beeswarm; use stratified split to mirror training script.
    X_train, X_test, _y_train, _y_test = train_test_split(
        X,
        y_raw.values,
        test_size=float(args.test_size),
        random_state=int(args.seed),
        stratify=y_raw.values,
    )

    # Load the trained model from model.json
    model = XGBClassifier()
    model.load_model(str(model_json))

    # Compute shap_values.
    # Prefer shap.TreeExplainer, but multiclass XGBoost JSON can store base_score as a list string
    # (e.g., "[5E-1,5E-1,5E-1]") which breaks some SHAP versions. Fall back to xgboost-native
    # pred_contribs (SHAP values) in that case.
    try:
        explainer = shap.TreeExplainer(model)
        sv = explainer.shap_values(X_test)
    except Exception as e:  # noqa: BLE001
        import xgboost as xgb

        print(f"[WARN] shap.TreeExplainer failed ({e}). Falling back to xgboost pred_contribs.")
        booster = model.get_booster()
        dmat = xgb.DMatrix(X_test, feature_names=list(X_test.columns))
        try:
            contrib = booster.predict(dmat, pred_contribs=True, strict_shape=True)
        except TypeError:
            contrib = booster.predict(dmat, pred_contribs=True)
        sv = np.asarray(contrib)

    # Normalize multiclass SHAP to (N, F, C) as requested.
    if isinstance(sv, list) and sv:
        # list[C] of (N,F) -> (N,F,C)
        arr = np.stack([np.asarray(v) for v in sv], axis=0)  # (C,N,F)
        values = np.transpose(arr, (1, 2, 0))  # (N,F,C)
    else:
        arr = np.asarray(sv)
        F = int(len(X_test.columns))
        # Regression/single output: (N, F) or (N, F+1)
        if arr.ndim == 2:
            # drop bias if present
            if arr.shape[1] == (F + 1):
                arr = arr[:, :-1]
            values = arr
        # Multiclass: could be (N, C, F+1)/(N, C, F) or (N, F, C)
        elif arr.ndim == 3:
            # drop bias column if present (xgboost pred_contribs yields ...,(F+1))
            if arr.shape[-1] == (F + 1):
                arr = arr[:, :, :-1]
            if arr.shape[1] == F:
                # (N, F, C)
                values = arr
            elif arr.shape[2] == F:
                # (N, C, F) -> (N, F, C)
                values = np.transpose(arr, (0, 2, 1))
            else:
                raise TypeError(f"Unexpected 3D shap_values shape (F={F}): {arr.shape}")
        else:
            raise TypeError(f"Unexpected shap_values shape: {arr.shape}")

    # Make an Explanation and pass it to shap.plots.beeswarm.
    pretty_names = [_pretty_feature_name(c) for c in X_test.columns]
    shap_values = shap.Explanation(values=values, data=X_test.to_numpy(), feature_names=pretty_names)

    # NOTE: shap.plots.beeswarm currently expects a 2D Explanation (N,F).
    # For multiclass (N,F,C), we provide two options:
    #  - panels: one beeswarm per class using slicing shap_values[:,:,k]
    #  - single: average over classes into a single (N,F)
    if values.ndim == 3 and args.multiclass_layout == "panels":
        # Try to load class names for nicer titles
        class_names = None
        mapping_csv = xgb_dir / "label_mapping.csv"
        if mapping_csv.exists():
            m = pd.read_csv(mapping_csv)
            if {"label_int", "label"}.issubset(m.columns):
                m = m.sort_values("label_int")
                class_names = [str(x) for x in m["label"].tolist()]

        n_classes = int(values.shape[2])
        if class_names is None or len(class_names) != n_classes:
            class_names = [f"class_{i}" for i in range(n_classes)]

        fig, axes = plt.subplots(1, n_classes, figsize=(5.2 * n_classes, 5.6), sharey=True)
        if n_classes == 1:
            axes = [axes]
        for i in range(n_classes):
            ax = axes[i]
            sv_i = shap.Explanation(
                values=values[:, :, i],
                data=X_test.to_numpy(),
                feature_names=pretty_names,
            )
            if args.debug:
                import numpy as _np

                print(
                    f"[DEBUG] class={i} values={_np.asarray(sv_i.values).shape} data={_np.asarray(sv_i.data).shape} "
                    f"n_features={len(sv_i.feature_names)}"
                )
            shap.plots.beeswarm(
                sv_i,
                max_display=int(args.max_display),
                show=False,
                ax=ax,
                color_bar=(i == n_classes - 1),
                plot_size=None,
            )
            ax.set_title(class_names[i], fontsize=12)
        fig.tight_layout()
    else:
        # single: average over classes if needed
        if values.ndim == 3:
            v2 = np.mean(values, axis=2)
            shap_values = shap.Explanation(values=v2, data=X_test.to_numpy(), feature_names=pretty_names)
        plt.figure(figsize=(10, 6))
        shap.plots.beeswarm(shap_values, max_display=int(args.max_display), show=False, plot_size=None)
        plt.tight_layout()

    out_png = Path(args.out)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=260, bbox_inches="tight")
    plt.close("all")

    print(f"[OK] Wrote: {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

