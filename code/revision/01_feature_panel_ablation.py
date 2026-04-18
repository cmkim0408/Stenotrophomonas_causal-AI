"""Workstream #1: feature panel ablation / expansion.

Compares classification + regression performance across:
  - curated_paper (paper-aligned reactions present in 120-width set)
  - top_10, top_20, top_50 (SHAP-ranked from all_widths)
  - all_widths (120)
  - random_30 controls × 10 seeds

Outputs:
    revision_runs/iscience_rev1/01_feature_panel_ablation/ablation_metrics.csv
    revision_runs/iscience_rev1/01_feature_panel_ablation/ablation_summary.md
    revision_runs/iscience_rev1/figures/feature_panel_ablation.{png,pdf}
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from revision import utils as u

OUT = u.ensure_outdirs("01_feature_panel_ablation")


def main() -> None:
    print("[01 ablation] starting...")
    df = u.load_regime_dataset()
    y, classes = pd.factorize(df["label"])
    classes = list(classes)
    sev = u.severity_target(df).to_numpy()

    # Static + SHAP-ranked panels
    static = u.build_panels(df)
    ranked = u.shap_ranked_panels(df, u.width_cols(df), y, sizes=(10, 20, 50))

    # Order: curated, top_10, top_20, top_50, all_widths, random_30 controls
    panels = [static[0]] + list(ranked) + [static[1]] + static[2:]

    rows: list[dict] = []
    for p in panels:
        if p.size == 0:
            print(f"  - {p.name}: SKIP (empty)")
            continue
        X = df[p.columns].to_numpy()
        clf_res = u.cv_classify(X, y, u.xgb_clf_factory(), n_splits=5)
        reg_res = u.cv_regress(X, sev, u.xgb_reg_factory(), n_splits=5)
        rows.append({
            "panel": p.name,
            "panel_kind": ("curated" if "curated" in p.name
                           else "top_k" if p.name.startswith("top_")
                           else "random" if p.name.startswith("random_")
                           else "full"),
            "n_features": p.size,
            "macro_f1": round(clf_res["macro_f1"], 4),
            "balanced_accuracy": round(clf_res["balanced_accuracy"], 4),
            "macro_f1_fold_mean": round(clf_res["macro_f1_fold_mean"], 4),
            "macro_f1_fold_std": round(clf_res["macro_f1_fold_std"], 4),
            **{f"per_class_f1__{cls}": round(v, 4)
               for cls, v in zip(classes, clf_res["per_class_f1"])},
            "rmse": round(reg_res["rmse"], 5),
            "mae": round(reg_res["mae"], 5),
            "r2": round(reg_res["r2"], 4),
            "note": p.note,
        })
        print(f"  - {p.name:25s} n={p.size:3d}  macro_F1={clf_res['macro_f1']:.4f}  "
              f"R²={reg_res['r2']:.4f}")

    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT["ws"] / "ablation_metrics.csv", index=False)

    # ---------- Figure: panel size vs metric, with random envelope ----------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    rand = metrics[metrics["panel_kind"] == "random"]
    others = metrics[metrics["panel_kind"] != "random"]

    for ax, metric, ylabel, ylim in (
        (axes[0], "macro_f1", "Macro-F1 (regime classification)", (0.0, 1.02)),
        (axes[1], "r2",       "R² (severity regression)",        (-0.1, 1.02)),
    ):
        # Random envelope (mean ± std at n=30)
        if not rand.empty:
            mean, std = rand[metric].mean(), rand[metric].std()
            ax.axhspan(mean - std, mean + std, color="lightgray", alpha=0.7,
                       label=f"random_30 ±1σ (mean={mean:.3f})")
            ax.axhline(mean, color="gray", lw=0.8, ls="--")

        for kind, marker, color in (("curated", "D", "#0072B2"),
                                    ("top_k",   "o", "#E69F00"),
                                    ("full",    "s", "#009E73")):
            sub = others[others["panel_kind"] == kind]
            if sub.empty: continue
            sub = sub.sort_values("n_features")
            ax.plot(sub["n_features"], sub[metric], marker=marker,
                    color=color, lw=1.4, ms=8, label=kind)
            for _, r in sub.iterrows():
                ax.annotate(r["panel"], (r["n_features"], r[metric]),
                            xytext=(4, 4), textcoords="offset points",
                            fontsize=7, color=color)
        ax.set_xlabel("# features in panel")
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="lower right")

    fig.suptitle("Feature panel ablation: classification & severity-regression performance",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT["figures"] / "feature_panel_ablation.png", dpi=180,
                bbox_inches="tight")
    fig.savefig(OUT["figures"] / "feature_panel_ablation.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[01 ablation] wrote figure → {OUT['figures'] / 'feature_panel_ablation.png'}")

    # ---------- Markdown summary ----------
    md = ["# Feature panel ablation — summary",
          "",
          f"- N conditions: {len(df)}  ·  classes: {classes}",
          f"- Universe: {len(u.width_cols(df))} `width__` columns "
          f"(parquet alphabetically truncated at 'FACOAL161'; M/I/N-prefixed paper anchors absent)",
          f"- Random-30 controls: {(metrics['panel_kind']=='random').sum()} seeds",
          "",
          "## Headline numbers",
          ""]
    headline = metrics[metrics["panel_kind"] != "random"][
        ["panel", "n_features", "macro_f1", "balanced_accuracy", "r2", "rmse"]
    ].to_string(index=False)
    md.append("```")
    md.append(headline)
    md.append("```")
    md.append("")
    if not rand.empty:
        md.append("## Random-30 controls")
        md.append("")
        md.append(f"- macro_F1 mean ± std: {rand['macro_f1'].mean():.3f} ± {rand['macro_f1'].std():.3f}")
        md.append(f"- R²        mean ± std: {rand['r2'].mean():.3f} ± {rand['r2'].std():.3f}")
        md.append("")
    md.append("## Curated paper-aligned panel composition")
    md.append("")
    curated = next(p for p in static if p.name == "curated_paper")
    md.append("```")
    md.extend(curated.columns)
    md.append("```")
    (OUT["ws"] / "ablation_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"[01 ablation] wrote summary → {OUT['ws'] / 'ablation_summary.md'}")


if __name__ == "__main__":
    main()
