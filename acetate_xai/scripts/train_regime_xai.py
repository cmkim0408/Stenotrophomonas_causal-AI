from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd


def _feature_columns(df: pd.DataFrame) -> list[str]:
    # Prefer the FVA-derived features; keep this conservative for interpretability.
    return [c for c in df.columns if c.startswith(("width__", "mid__", "signchange__"))]


def _write_confusion_matrix(y_true: np.ndarray, y_pred: np.ndarray, labels: list[str], out_csv: Path) -> None:
    cm = pd.crosstab(pd.Series(y_true, name="true"), pd.Series(y_pred, name="pred"), dropna=False)
    # ensure all labels appear
    for lab in labels:
        if lab not in cm.index:
            cm.loc[lab] = 0
        if lab not in cm.columns:
            cm[lab] = 0
    cm = cm.loc[labels, labels]
    cm.to_csv(out_csv)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Train regime classification models (tree + linear) from regime_dataset.parquet.")
    p.add_argument("--data", required=True, help="Input parquet (e.g., results/regime_dataset.parquet)")
    p.add_argument("--outdir", required=True, help="Output directory (e.g., results/xai_platform)")
    p.add_argument("--tree-depth", type=int, default=4, help="Decision tree max_depth")
    p.add_argument("--test-size", type=float, default=0.2, help="Test split size (default 0.2)")
    p.add_argument("--random-state", type=int, default=0, help="Random state")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    # Avoid oversubscription
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    data_path = Path(args.data)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(data_path)
    if "label" not in df.columns:
        raise ValueError("Dataset missing required column: label")

    # label distribution (all data)
    label_dist = df["label"].value_counts().rename_axis("label").reset_index(name="count")
    label_dist.to_csv(outdir / "label_distribution.csv", index=False)

    feat_cols = _feature_columns(df)
    if not feat_cols:
        raise ValueError("No feature columns found (expected width__/mid__/signchange__).")

    X = df[feat_cols].copy()
    # booleans -> 0/1, fill missing
    for c in X.columns:
        if X[c].dtype == bool:
            X[c] = X[c].astype(int)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    y = df["label"].astype(str).values

    # Train/test split (stratified)
    from sklearn.model_selection import train_test_split

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=float(args.test_size),
        random_state=int(args.random_state),
        stratify=y,
    )

    labels_sorted = sorted(pd.unique(y).tolist())

    # 1) DecisionTreeClassifier
    from sklearn.tree import DecisionTreeClassifier, export_text

    tree = DecisionTreeClassifier(max_depth=int(args.tree_depth), random_state=int(args.random_state))
    tree.fit(X_train, y_train)
    y_pred_tree = tree.predict(X_test)
    (outdir / "tree_rules.txt").write_text(export_text(tree, feature_names=list(X.columns)) + "\n", encoding="utf-8")
    _write_confusion_matrix(y_test, y_pred_tree, labels_sorted, outdir / "confusion_matrix.csv")

    # 2) LogisticRegression (saga, elasticnet)
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    clf = Pipeline(
        steps=[
            ("scaler", StandardScaler(with_mean=True, with_std=True)),
            (
                "logreg",
                LogisticRegression(
                    solver="saga",
                    penalty="elasticnet",
                    l1_ratio=0.5,
                    max_iter=5000,
                    random_state=int(args.random_state),
                    n_jobs=1,
                    multi_class="auto",
                ),
            ),
        ]
    )
    clf.fit(X_train, y_train)
    logreg: LogisticRegression = clf.named_steps["logreg"]
    coef = logreg.coef_
    # multi-class: aggregate |coef| across classes
    if coef.ndim == 2:
        score = np.mean(np.abs(coef), axis=0)
    else:
        score = np.abs(coef)
    top_idx = np.argsort(-score)[:30]
    top = pd.DataFrame(
        {
            "feature": [X.columns[i] for i in top_idx],
            "mean_abs_coef": [float(score[i]) for i in top_idx],
        }
    )
    top.to_csv(outdir / "linear_top_features.csv", index=False)

    print(f"[OK] Wrote: {outdir / 'label_distribution.csv'}")
    print(f"[OK] Wrote: {outdir / 'confusion_matrix.csv'}")
    print(f"[OK] Wrote: {outdir / 'tree_rules.txt'}")
    print(f"[OK] Wrote: {outdir / 'linear_top_features.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

