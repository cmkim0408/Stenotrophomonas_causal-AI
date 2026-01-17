from __future__ import annotations

import argparse
import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Train XGBoost regime classifier + SHAP explanations.")
    p.add_argument("--data", default="results/regime_dataset.parquet", help="Input parquet path")
    p.add_argument("--outdir", default="results/xai_xgb", help="Output directory")
    p.add_argument("--label-col", default="label", help="Label column name")
    p.add_argument("--topk", type=int, default=20, help="Top-K features for SHAP bar plot focus (informational)")
    p.add_argument("--test-size", type=float, default=0.2, help="Test size (default 0.2)")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    p.add_argument("--max-depth", type=int, default=4, help="XGBoost tree depth")
    p.add_argument("--n-estimators", type=int, default=500, help="Number of boosting rounds")
    p.add_argument("--learning-rate", type=float, default=0.05, help="Learning rate")
    p.add_argument("--subsample", type=float, default=0.8, help="Subsample")
    p.add_argument("--colsample-bytree", type=float, default=0.8, help="Colsample bytree")
    p.add_argument("--n-jobs", type=int, default=-1, help="Parallel jobs for XGBoost (-1 = all)")
    p.add_argument("--no-dependence", action="store_true", help="Skip SHAP dependence plots.")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def _feature_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith(("width__", "mid__", "signchange__"))]


def _safe_filename(s: str, max_len: int = 120) -> str:
    s2 = re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_")
    if not s2:
        s2 = "feature"
    return s2[:max_len]


def _shap_global_importance_and_dep_values(shap_values) -> tuple[np.ndarray, np.ndarray | None]:
    """
    Return:
      - global_importance: (F,) mean(|SHAP|) aggregated over samples + classes (if multi-class)
      - dep_values: (N,F) SHAP values for dependence plots (averaged over classes if needed), or None

    Handles common SHAP return shapes:
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
        # (N,F)
        global_importance = np.mean(np.abs(arr), axis=0)  # (F,)
        return global_importance, arr

    if arr.ndim == 3:
        # Could be (C,N,F) or (N,F,C)
        if arr.shape[0] < 20 and arr.shape[0] != arr.shape[1]:
            # likely (C,N,F)
            global_importance = np.mean(np.abs(arr), axis=(0, 1))  # (F,)
            dep_vals = np.mean(arr, axis=0)  # (N,F)
            return global_importance, dep_vals
        # assume (N,F,C)
        global_importance = np.mean(np.abs(arr), axis=(0, 2))  # (F,)
        dep_vals = np.mean(arr, axis=2)  # (N,F)
        return global_importance, dep_vals

    raise TypeError(f"Unsupported shap_values shape: {getattr(arr, 'shape', None)}")


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    data_path = Path(args.data)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Hard requirement: we rely on xgboost + shap, but we don't modify project dependencies here.
    try:
        from xgboost import XGBClassifier
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

    # Matplotlib (Agg backend) for headless servers
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402

    if not data_path.exists():
        print(f"[ERROR] data not found: {data_path}")
        return 2

    df = pd.read_parquet(data_path)
    if args.label_col not in df.columns:
        print(f"[ERROR] label column not found: {args.label_col}")
        return 2

    # y cleaning: drop NaN/empty labels
    y_raw = df[args.label_col].astype(str)
    mask = y_raw.notna() & (y_raw.str.strip() != "") & (y_raw.str.lower() != "nan")
    df = df.loc[mask].copy()
    y_raw = df[args.label_col].astype(str).str.strip()

    feat_cols = _feature_cols(df)
    if not feat_cols:
        print("[ERROR] No feature columns found (expected width__/mid__/signchange__).")
        return 2

    X = df[feat_cols].copy()
    # bool -> int, fill NaN
    for c in X.columns:
        if X[c].dtype == bool:
            X[c] = X[c].astype(int)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # label encoding
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import LabelEncoder
    from sklearn.metrics import classification_report

    le = LabelEncoder()
    y = le.fit_transform(y_raw.values)
    classes = le.classes_.tolist()

    # class distribution
    dist = pd.Series(y_raw).value_counts()
    print(f"[INFO] data rows={len(df)}, features={len(feat_cols)}, classes={len(classes)}")
    print("[INFO] class distribution:")
    for k, v in dist.items():
        print(f"  - {k}: {int(v)}")

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=float(args.test_size),
        random_state=int(args.seed),
        stratify=y,
    )

    model = XGBClassifier(
        objective="multi:softprob",
        eval_metric="mlogloss",
        num_class=int(len(classes)),
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

    # Predictions + evaluation
    y_pred = model.predict(X_test)
    rep = classification_report(y_test, y_pred, target_names=classes)
    (outdir / "classification_report.txt").write_text(rep + "\n", encoding="utf-8")

    # Confusion matrix (true x pred)
    cm = pd.crosstab(
        pd.Series(le.inverse_transform(y_test), name="true"),
        pd.Series(le.inverse_transform(y_pred), name="pred"),
        dropna=False,
    )
    # ensure all labels appear
    for lab in classes:
        if lab not in cm.index:
            cm.loc[lab] = 0
        if lab not in cm.columns:
            cm[lab] = 0
    cm = cm.loc[classes, classes]
    cm.to_csv(outdir / "confusion_matrix.csv")

    # Label mapping
    pd.DataFrame({"label_int": list(range(len(classes))), "label": classes}).to_csv(
        outdir / "label_mapping.csv", index=False
    )

    # Save model
    model.save_model(str(outdir / "model.json"))

    # SHAP
    explainer = shap.TreeExplainer(model)
    # shap values for multi-class are typically a list: [n_classes] of (n_samples, n_features)
    shap_values = explainer.shap_values(X_test)
    global_importance, dep_vals = _shap_global_importance_and_dep_values(shap_values)

    # Save global SHAP importance
    shap_imp = (
        pd.DataFrame({"feature": list(X_test.columns), "mean_abs_shap": global_importance})
        .sort_values("mean_abs_shap", ascending=False)
        .reset_index(drop=True)
    )
    shap_imp.to_csv(outdir / "shap_importance.csv", index=False)

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

    # Top-3 dependence plots using global importance (mean over samples + classes)
    top3_idx = np.argsort(global_importance)[-3:][::-1]
    top3_features = [str(X_test.columns[int(i)]) for i in top3_idx]

    if args.no_dependence:
        print("[INFO] --no-dependence set: skipping dependence plots")
    else:
        if dep_vals is None:
            print("[WARN] dependence plots skipped (no suitable dep_vals)")
        else:
            for feat in top3_features:
                safe = _safe_filename(feat)
                plt.figure()
                shap.dependence_plot(feat, dep_vals, X_test, show=False)
                plt.tight_layout()
                plt.savefig(outdir / f"shap_dependence_{safe}.png", dpi=200, bbox_inches="tight")
                plt.close()

    # Log a short metrics snippet
    acc = float(np.mean(y_pred == y_test))
    print(f"[INFO] test accuracy={acc:.4f} (quick check)")
    print(f"[OK] Wrote: {outdir / 'confusion_matrix.csv'}")
    print(f"[OK] Wrote: {outdir / 'classification_report.txt'}")
    print(f"[OK] Wrote: {outdir / 'label_mapping.csv'}")
    print(f"[OK] Wrote: {outdir / 'model.json'}")
    print(f"[OK] Wrote: {outdir / 'shap_beeswarm.png'}")
    print(f"[OK] Wrote: {outdir / 'shap_bar.png'}")
    print(f"[OK] Wrote: {outdir / 'shap_importance.csv'}")
    if not args.no_dependence:
        for feat in top3_features:
            print(f"[OK] Wrote: {outdir / ('shap_dependence_' + _safe_filename(feat) + '.png')}")

    # Keep a small metadata JSON (optional, but useful)
    meta = {
        "data": str(data_path),
        "rows": int(len(df)),
        "n_features": int(len(feat_cols)),
        "classes": classes,
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

