#!/usr/bin/env python3
"""
Local validation script for holdout predictions vs experimental OD.

Run on local PC after experimental OD data is available.
Reads: results/holdout_predictions.csv + data/holdout_od_results.csv
Outputs: local_validation/ (summary_metrics.csv, severity_vs_od.png, rank_plot.png, confusion_matrix.png)

Lightweight - no server computation needed.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Validate holdout predictions against experimental OD (local run)."
    )
    ap.add_argument(
        "--pred",
        default="results/holdout_predictions.csv",
        help="Holdout predictions CSV",
    )
    ap.add_argument(
        "--od",
        default="data/holdout_od_results.csv",
        help="Experimental OD results CSV (condition_id, replicate, od600_32h)",
    )
    ap.add_argument(
        "--out",
        default="local_validation",
        help="Output directory",
    )
    args = ap.parse_args()

    pred_path = Path(args.pred)
    od_path = Path(args.od)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not pred_path.exists():
        print(f"[ERROR] Predictions not found: {pred_path}", file=sys.stderr)
        return 2

    pred_df = pd.read_csv(pred_path)
    if "condition_id" not in pred_df.columns or "predicted_severity" not in pred_df.columns:
        print("[ERROR] Predictions CSV must have condition_id, predicted_severity", file=sys.stderr)
        return 2

    if not od_path.exists():
        print(f"[WARN] OD results not found: {od_path}", file=sys.stderr)
        print("[INFO] Create from template: cp data/holdout_od_results_template.csv data/holdout_od_results.csv", file=sys.stderr)
        print("[INFO] Filling with metrics only (no OD-based plots)", file=sys.stderr)
        od_df = None
    else:
        od_df = pd.read_csv(od_path)
        if "condition_id" not in od_df.columns or "od600_32h" not in od_df.columns:
            print("[ERROR] OD CSV must have condition_id, od600_32h", file=sys.stderr)
            return 2

    # Summary metrics (always)
    metrics = []
    metrics.append({"metric": "n_conditions", "value": len(pred_df)})

    if od_df is not None:
        od_df["od600_32h"] = pd.to_numeric(od_df["od600_32h"], errors="coerce")
        od_valid = od_df.dropna(subset=["od600_32h"])
        if len(od_valid) > 0:
            od_mean = od_valid.groupby("condition_id")["od600_32h"].mean()
            merged = pred_df.merge(
                od_mean.rename("mean_od600_32h").reset_index(),
                on="condition_id",
                how="inner",
            )
            if len(merged) >= 2:
                try:
                    from scipy.stats import spearmanr
                except ImportError:
                    rho, pval = float("nan"), float("nan")
                else:
                    rho, pval = spearmanr(merged["predicted_severity"], merged["mean_od600_32h"])
                metrics.append({"metric": "spearman_rho", "value": float(rho)})
                metrics.append({"metric": "spearman_pval", "value": float(pval)})
                metrics.append({"metric": "n_matched", "value": len(merged)})

                # Plots
                try:
                    import matplotlib
                    matplotlib.use("Agg")
                    import matplotlib.pyplot as plt

                    # severity vs OD
                    fig, ax = plt.subplots(figsize=(6, 5))
                    ax.scatter(merged["mean_od600_32h"], merged["predicted_severity"], s=60)
                    for _, r in merged.iterrows():
                        ax.annotate(r["condition_id"], (r["mean_od600_32h"], r["predicted_severity"]),
                                   xytext=(5, 5), textcoords="offset points", fontsize=9)
                    ax.set_xlabel("Mean OD600 @ 32h")
                    ax.set_ylabel("Predicted severity")
                    ax.set_title("Holdout: Predicted severity vs experimental OD")
                    fig.tight_layout()
                    fig.savefig(out_dir / "severity_vs_od.png", dpi=150, bbox_inches="tight")
                    plt.close()

                    # rank plot (OD rank vs severity rank)
                    merged = merged.copy()
                    merged["od_rank"] = merged["mean_od600_32h"].rank()
                    merged["sev_rank"] = merged["predicted_severity"].rank()
                    fig, ax = plt.subplots(figsize=(6, 5))
                    ax.scatter(merged["od_rank"], merged["sev_rank"], s=60)
                    for _, r in merged.iterrows():
                        ax.annotate(r["condition_id"], (r["od_rank"], r["sev_rank"]),
                                   xytext=(5, 5), textcoords="offset points", fontsize=9)
                    ax.plot([0.5, merged["od_rank"].max() + 0.5], [0.5, merged["od_rank"].max() + 0.5], "k--", alpha=0.5)
                    ax.set_xlabel("OD600 rank (experimental)")
                    ax.set_ylabel("Predicted severity rank")
                    ax.set_title("Rank correlation (Spearman rho=%.3f)" % rho)
                    fig.tight_layout()
                    fig.savefig(out_dir / "rank_plot.png", dpi=150, bbox_inches="tight")
                    plt.close()

                    # confusion-style: regime vs OD-based regime (if we map OD to regime)
                    # Simple: OD low/high by median split -> binary
                    median_od = merged["mean_od600_32h"].median()
                    merged["od_low"] = merged["mean_od600_32h"] < median_od
                    merged["sev_low"] = merged["predicted_severity"] < merged["predicted_severity"].median()
                    try:
                        from sklearn.metrics import confusion_matrix as sk_confusion_matrix
                        cm = sk_confusion_matrix(merged["od_low"].astype(int), merged["sev_low"].astype(int))
                    except ImportError:
                        cm = np.array([[0, 0], [0, 0]])
                    fig, ax = plt.subplots(figsize=(5, 4))
                    im = ax.imshow(cm, cmap="Blues")
                    ax.set_xticks([0, 1])
                    ax.set_yticks([0, 1])
                    ax.set_xticklabels(["High OD", "Low OD"])
                    ax.set_yticklabels(["High sev", "Low sev"])
                    for i in range(2):
                        for j in range(2):
                            ax.text(j, i, str(cm[i, j]), ha="center", va="center")
                    ax.set_title("OD (median-split) vs Predicted severity (median-split)")
                    fig.colorbar(im, ax=ax, label="count")
                    fig.tight_layout()
                    fig.savefig(out_dir / "confusion_matrix.png", dpi=150, bbox_inches="tight")
                    plt.close()

                except ImportError:
                    print("[WARN] matplotlib/scipy not available, skipping plots", file=sys.stderr)

    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(out_dir / "summary_metrics.csv", index=False)
    print(f"[OK] Wrote {out_dir / 'summary_metrics.csv'}")
    if od_df is not None and len(od_valid) > 0:
        print(f"[OK] Wrote {out_dir / 'severity_vs_od.png'}")
        print(f"[OK] Wrote {out_dir / 'rank_plot.png'}")
        print(f"[OK] Wrote {out_dir / 'confusion_matrix.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
