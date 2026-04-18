"""Workstream consolidator: build metrics_summary.csv + REPORT.md from per-workstream outputs."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from revision import utils as u

OUT = u.REVISION_OUT


def _safe_read(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path)
    except Exception:
        return None


def main() -> None:
    print("[06 consolidate] starting...")

    # Collect headline numbers from each workstream into a unified summary table
    rows: list[dict] = []

    abl = _safe_read(OUT / "01_feature_panel_ablation" / "ablation_metrics.csv")
    if abl is not None:
        non_random = abl[abl["panel_kind"] != "random"]
        random = abl[abl["panel_kind"] == "random"]
        for _, r in non_random.iterrows():
            rows.append({
                "workstream": "01_ablation",
                "key": r["panel"],
                "n_features": int(r["n_features"]),
                "macro_f1": r["macro_f1"],
                "r2": r["r2"],
                "rmse": r.get("rmse"),
            })
        if not random.empty:
            rows.append({
                "workstream": "01_ablation",
                "key": "random_30 (mean ± std)",
                "n_features": 30,
                "macro_f1": f"{random['macro_f1'].mean():.4f} ± {random['macro_f1'].std():.4f}",
                "r2": f"{random['r2'].mean():.4f} ± {random['r2'].std():.4f}",
                "rmse": f"{random['rmse'].mean():.5f}",
            })

    bench = _safe_read(OUT / "02_benchmarking" / "benchmark_metrics.csv")
    if bench is not None:
        for _, r in bench.iterrows():
            rows.append({
                "workstream": "02_benchmark",
                "key": f"{r['feature_set']} + {r['clf_model']}",
                "n_features": int(r["n_features"]),
                "macro_f1": r["macro_f1"],
                "r2": r["r2"],
                "rmse": r["rmse"],
            })

    transfer = _safe_read(OUT / "03_transfer" / "transfer_metrics.csv")
    if transfer is not None:
        for _, r in transfer.iterrows():
            rows.append({
                "workstream": "03_transfer",
                "key": r["system"],
                "n_features": r.get("n_features"),
                "macro_f1": r.get("macro_f1"),
                "r2": r.get("r2_severity"),
                "rmse": r.get("rmse_severity"),
            })

    runtime = _safe_read(OUT / "04_runtime" / "runtime_summary.csv")
    if runtime is not None:
        for _, r in runtime.iterrows():
            rows.append({
                "workstream": "04_runtime",
                "key": r.get("stage"),
                "n_features": r.get("n", ""),
                "macro_f1": "",
                "r2": "",
                "rmse": f"{r.get('wall_seconds')} s, {r.get('peak_rss_mb')} MB",
            })

    perf = _safe_read(OUT / "05_existing_data" / "performance_metrics.csv")
    if perf is not None:
        macro = perf[perf["class"] == "macro_avg"]
        if not macro.empty:
            r = macro.iloc[0]
            rows.append({
                "workstream": "05_existing_data",
                "key": "5-fold CV (all_widths, XGBoost)",
                "n_features": 120,
                "macro_f1": round(r["f1"], 4),
                "r2": "(see residual_summary.csv)",
                "rmse": "",
            })
    res = _safe_read(OUT / "05_existing_data" / "residual_summary.csv")
    if res is not None:
        d = dict(zip(res["metric"], res["value"]))
        rows.append({
            "workstream": "05_existing_data",
            "key": "severity regression",
            "n_features": 120,
            "macro_f1": "",
            "r2": round(d.get("r2", float("nan")), 4),
            "rmse": round(d.get("rmse", float("nan")), 5),
        })

    metrics_summary = pd.DataFrame(rows,
        columns=["workstream", "key", "n_features", "macro_f1", "r2", "rmse"])
    metrics_summary.to_csv(OUT / "metrics_summary.csv", index=False)
    print(f"  wrote metrics_summary.csv ({len(metrics_summary)} rows)")

    # ----- REPORT.md -----
    md: list[str] = []
    md.append("# iScience revision — REPORT")
    md.append("")
    md.append("**Branch:** `revision/iscience-rev1`")
    md.append("**Date:** 2026-04-18")
    md.append("")
    md.append("## At-a-glance")
    md.append("")
    md.append("| Workstream | Key result |")
    md.append("|---|---|")
    if abl is not None:
        cur = abl[abl["panel"] == "curated_paper"].iloc[0]
        all_w = abl[abl["panel"] == "all_widths"].iloc[0]
        md.append(f"| #1 Feature ablation | curated_paper (n={int(cur['n_features'])}): "
                  f"F1={cur['macro_f1']:.4f}, R²={cur['r2']:.4f}; all_widths (120): "
                  f"F1={all_w['macro_f1']:.4f}, R²={all_w['r2']:.4f}; random_30 envelope ≈ same |")
    if bench is not None:
        ours = bench[(bench["feature_set"] == "C_curated_widths") & (bench["clf_model"] == "XGBoost")].iloc[0]
        a = bench[bench["feature_set"] == "A_inputs_only"]["macro_f1"].max()
        md.append(f"| #2 Benchmark | curated+XGBoost F1={ours['macro_f1']:.4f}, R²={ours['r2']:.4f}; "
                  f"best inputs-only F1={a:.4f} |")
    if perf is not None and res is not None:
        macro = perf[perf["class"] == "macro_avg"].iloc[0]
        d = dict(zip(res["metric"], res["value"]))
        md.append(f"| #3 Existing-data | 5-fold CV macro-F1={macro['f1']:.4f}, "
                  f"severity R²={d.get('r2', 0):.4f} |")
    if transfer is not None:
        e = transfer[transfer["system"] == "iML1515"].iloc[0]
        md.append(f"| #4 iML1515 transfer | macro-F1={e['macro_f1']:.4f}, "
                  f"R²={e['r2_severity']:.4f} on {e['n_conditions']} LHS conditions |")
    if runtime is not None:
        total = runtime["wall_seconds"].sum()
        md.append(f"| #5 Runtime | total stages = {total:.1f} s on a single core |")
    md.append("")

    # Per-workstream sections — pull each summary.md as-is
    for sub, label in (("01_feature_panel_ablation", "1. Feature panel ablation"),
                       ("02_benchmarking", "2. Baseline benchmarking"),
                       ("05_existing_data", "3. Existing-data performance summary"),
                       ("03_transfer", "4. External transfer (iML1515)"),
                       ("04_runtime", "5. Runtime profiling")):
        md.append(f"## {label}")
        md.append("")
        candidates = list((OUT / sub).glob("*summary.md"))
        if candidates:
            text = candidates[0].read_text(encoding="utf-8")
            md.append(text.strip())
        else:
            md.append("_(no summary.md found)_")
        md.append("")
        md.append(f"Outputs: `revision_runs/iscience_rev1/{sub}/` + `figures/`")
        md.append("")

    md.append("## Known limitations")
    md.append("")
    md.append("- `regime_dataset.parquet` width__ universe is alphabetically truncated at "
              "`FACOAL161`, so paper-named TCA enzymes beyond 'F' (MDH, ICDH, ICL, MALS, "
              "PFK, PYK, NDH, PPC, PCK) are absent. The ablation thus operates within "
              "the 120-width superset that the deployed model actually uses; results "
              "still answer Reviewer 2 (panel-size robustness) but qualitative cross-"
              "feature comparisons against paper Fig 4 should be qualified.")
    md.append("- The B_gem_summary baseline includes `objective_value`, which is the "
              "numerator of the severity target G = obj/obj_max — its R²≈0.998 reflects "
              "construction overlap, not new predictive power.")
    md.append("- iML1515 transfer is in-silico only; experimental matching deferred.")
    md.append("")

    (OUT / "REPORT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"  wrote REPORT.md ({sum(len(l) for l in md)} chars)")

    print("[06 consolidate] done")


if __name__ == "__main__":
    main()
