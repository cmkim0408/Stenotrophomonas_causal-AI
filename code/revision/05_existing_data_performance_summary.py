"""Workstream #3: existing-data performance summary.

Computes from already-collected files (no new simulations):
  - Confusion matrix (5-fold CV oof on regime_dataset.parquet)
  - Per-class precision / recall / F1
  - Severity regression residual summary
  - Top mismatch conditions (holdout: predicted vs measured)

Outputs:
    revision_runs/iscience_rev1/05_existing_data/performance_metrics.csv
    revision_runs/iscience_rev1/05_existing_data/residual_summary.csv
    revision_runs/iscience_rev1/05_existing_data/top_mismatch_conditions.csv
    revision_runs/iscience_rev1/figures/confusion_matrix.{png,pdf}
    revision_runs/iscience_rev1/figures/regression_residuals.{png,pdf}
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xgboost as xgb

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from revision import utils as u

OUT = u.ensure_outdirs("05_existing_data")


def _confusion_matrix_plot(cm: np.ndarray, classes: list[str], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.5, 4.0))
    im = ax.imshow(cm, cmap="Blues")
    ax.set_xticks(range(len(classes)))
    ax.set_yticks(range(len(classes)))
    ax.set_xticklabels(classes, rotation=20, ha="right")
    ax.set_yticklabels(classes)
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    ax.set_title("Regime classification — 5-fold CV confusion matrix")
    for i in range(len(classes)):
        for j in range(len(classes)):
            ax.text(j, i, int(cm[i, j]), ha="center", va="center",
                    color="white" if cm[i, j] > cm.max() / 2 else "black")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=180, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def _residual_plot(y_true: np.ndarray, y_pred: np.ndarray, path: Path) -> None:
    res = y_pred - y_true
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    axes[0].scatter(y_true, y_pred, s=12, alpha=0.7, color="#0072B2")
    lo, hi = float(min(y_true.min(), y_pred.min())), float(max(y_true.max(), y_pred.max()))
    axes[0].plot([lo, hi], [lo, hi], "k--", lw=0.7)
    axes[0].set_xlabel("Measured severity (G = obj/obj_max)")
    axes[0].set_ylabel("Predicted (5-fold OOF)")
    axes[0].set_title("Predicted vs measured")
    axes[0].grid(alpha=0.3)

    axes[1].hist(res, bins=30, color="#E69F00", edgecolor="black", lw=0.4)
    axes[1].axvline(0, color="black", lw=0.7)
    axes[1].set_xlabel("Residual (predicted − measured)")
    axes[1].set_ylabel("Count")
    axes[1].set_title(f"Residuals (mean={res.mean():.3f}, sd={res.std():.3f})")
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=180, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    print("[05 existing-data] starting...")
    df = u.load_regime_dataset()
    y, classes = pd.factorize(df["label"])
    classes = list(classes)
    sev = u.severity_target(df).to_numpy()

    # Use ALL widths (matches paper's intracellular_only filter scope)
    widths = u.width_cols(df)
    X = df[widths].to_numpy()

    clf_res = u.cv_classify(X, y, u.xgb_clf_factory(), n_splits=5)
    reg_res = u.cv_regress(X, sev, u.xgb_reg_factory(), n_splits=5)

    # ----- Per-class metrics CSV -----
    pcm = pd.DataFrame({
        "class": classes,
        "precision": clf_res["per_class_precision"],
        "recall": clf_res["per_class_recall"],
        "f1": clf_res["per_class_f1"],
    })
    pcm.loc[len(pcm)] = ["macro_avg",
                          float(np.mean(clf_res["per_class_precision"])),
                          float(np.mean(clf_res["per_class_recall"])),
                          clf_res["macro_f1"]]
    pcm.to_csv(OUT["ws"] / "performance_metrics.csv", index=False)
    print("  wrote performance_metrics.csv")

    # ----- Confusion matrix figure -----
    cm = np.array(clf_res["confusion_matrix"])
    _confusion_matrix_plot(cm, classes, OUT["figures"] / "confusion_matrix")
    print("  wrote confusion_matrix.{png,pdf}")

    # ----- Severity regression residual summary -----
    pred = np.array(reg_res["oof"])
    res = pred - sev
    res_summary = pd.DataFrame({
        "metric": ["rmse", "mae", "r2", "residual_mean", "residual_std",
                   "residual_p10", "residual_p90"],
        "value": [reg_res["rmse"], reg_res["mae"], reg_res["r2"],
                  float(res.mean()), float(res.std()),
                  float(np.percentile(res, 10)),
                  float(np.percentile(res, 90))],
    })
    res_summary.to_csv(OUT["ws"] / "residual_summary.csv", index=False)
    _residual_plot(sev, pred, OUT["figures"] / "regression_residuals")
    print("  wrote residual_summary.csv + regression_residuals.{png,pdf}")

    # ----- Top mismatch conditions from C1..C10 holdout -----
    pred_csv = u.ROOT / "results" / "holdout_predictions.csv"
    od_csv = u.ROOT / "data" / "holdout_od_results.csv"
    mismatch_path = OUT["ws"] / "top_mismatch_conditions.csv"
    if pred_csv.exists() and od_csv.exists():
        preds = pd.read_csv(pred_csv)
        od_raw = pd.read_csv(od_csv)
        od = (od_raw.groupby("condition_id")["od600_32h"]
                    .agg(od_mean="mean", od_sd="std", n="count")
                    .reset_index())
        merged = preds.merge(od, on="condition_id", how="left")
        # Map predicted_severity to OD-equivalent scale: rank both, compute
        # absolute Spearman-style residual (0..1 normalized).
        merged["sev_rank"] = merged["predicted_severity"].rank(pct=True)
        merged["od_rank"] = merged["od_mean"].rank(pct=True)
        merged["rank_residual"] = (merged["sev_rank"] - merged["od_rank"]).abs()
        merged = merged.sort_values("rank_residual", ascending=False)
        merged.to_csv(mismatch_path, index=False)
        print(f"  wrote top_mismatch_conditions.csv (n={len(merged)})")
    else:
        print(f"  SKIP top_mismatch_conditions: missing {pred_csv} or {od_csv}")

    # ----- Summary md -----
    md = ["# Existing-data performance summary",
          "",
          f"- N conditions: {len(df)}, classes: {classes}",
          f"- Universe: {len(widths)} `width__` columns (all widths)",
          "",
          "## 5-fold CV classification",
          "",
          f"- macro-F1: **{clf_res['macro_f1']:.4f}**",
          f"- balanced accuracy: **{clf_res['balanced_accuracy']:.4f}**",
          f"- per-class F1: {dict(zip(classes, [round(x, 3) for x in clf_res['per_class_f1']]))}",
          "",
          "## 5-fold CV severity regression",
          "",
          f"- RMSE = {reg_res['rmse']:.5f}",
          f"- MAE = {reg_res['mae']:.5f}",
          f"- R² = {reg_res['r2']:.4f}",
          "",
          "## Top mismatch (holdout C1..C10)",
          "",
          "See `top_mismatch_conditions.csv` (sorted by rank residual).",
          "",
          ]
    (OUT["ws"] / "existing_data_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("[05 existing-data] done")


if __name__ == "__main__":
    main()
