"""Workstream #5 (numbered 04 in the plan order): runtime / scalability profiling.

Wall-clock and environment capture for every major pipeline stage:
    LHS (synthetic) → FBA batch → targeted FVA → classifier+SHAP → regressor+SHAP → PC bootstrap.

Produces:
    revision_runs/iscience_rev1/04_runtime/runtime_summary.csv
    revision_runs/iscience_rev1/04_runtime/environment_summary.txt
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from revision import utils as u

OUT = u.ensure_outdirs("04_runtime")
ROOT = u.ROOT


def _peak_rss_mb() -> float:
    try:
        import psutil  # type: ignore
        return psutil.Process().memory_info().rss / 1024 / 1024
    except Exception:
        return float("nan")


def stage_lhs_sample(n: int = 500) -> dict:
    from scipy.stats.qmc import LatinHypercube
    with u.Timer("lhs") as t:
        sampler = LatinHypercube(d=4, seed=42)
        sample = sampler.random(n=n)
    return {"stage": "lhs_sample", "n": n, "wall_seconds": round(t.seconds, 4),
            "peak_rss_mb": round(_peak_rss_mb(), 1)}


def stage_fba_batch(n: int = 50) -> dict:
    """Run n FBA solves on the iSO1_933 model with random uptake bounds.

    Sets bounds via reaction.bounds = (lb, ub) tuple to avoid asymmetric-update
    constraint errors when sweeping aggressively.
    """
    import cobra
    model = cobra.io.read_sbml_model(str(ROOT / "BaseModel.xml"))
    rng = np.random.default_rng(42)
    rxn_ids = {r.id for r in model.reactions}
    targets = [(rid, lo, 0.0) for rid, lo in (("EX_o2_e", 100), ("EX_ac_e", 50),
                                              ("EX_nh4_e", 20))
               if rid in rxn_ids]
    with u.Timer("fba") as t:
        for _ in range(n):
            for rid, lo_max, ub in targets:
                lb = -float(rng.uniform(0.5 * lo_max, lo_max))
                model.reactions.get_by_id(rid).bounds = (lb, ub)
            model.optimize()
    return {"stage": "fba_batch", "n": n, "wall_seconds": round(t.seconds, 4),
            "peak_rss_mb": round(_peak_rss_mb(), 1)}


def stage_targeted_fva(n_conditions: int = 10, n_targets: int = 30) -> dict:
    import cobra
    from cobra.flux_analysis import flux_variability_analysis as fva
    model = cobra.io.read_sbml_model(str(ROOT / "BaseModel.xml"))
    rxn_ids = [r.id for r in model.reactions[:n_targets]]
    rng = np.random.default_rng(42)
    o2_rxn = model.reactions.get_by_id("EX_o2_e") if "EX_o2_e" in [r.id for r in model.reactions] else None
    with u.Timer("fva") as t:
        for _ in range(n_conditions):
            if o2_rxn: o2_rxn.lower_bound = -float(rng.uniform(40, 100))
            try:
                fva(model, reaction_list=rxn_ids, fraction_of_optimum=0.95, processes=1)
            except Exception:
                pass
    return {"stage": "targeted_fva", "n": n_conditions, "n_targets": n_targets,
            "wall_seconds": round(t.seconds, 4), "peak_rss_mb": round(_peak_rss_mb(), 1)}


def stage_classifier_shap() -> dict:
    import xgboost as xgb
    df = u.load_regime_dataset()
    cols = u.width_cols(df)
    X = df[cols].to_numpy()
    y, _ = pd.factorize(df["label"])
    with u.Timer("clf") as t:
        clf = u.xgb_clf_factory()()
        clf.fit(X, y)
    with u.Timer("shap_clf") as t2:
        contribs = clf.get_booster().predict(xgb.DMatrix(X), pred_contribs=True)
    return {"stage": "classifier+shap", "n": len(X), "n_features": len(cols),
            "wall_seconds": round(t.seconds + t2.seconds, 4),
            "train_seconds": round(t.seconds, 4), "shap_seconds": round(t2.seconds, 4),
            "peak_rss_mb": round(_peak_rss_mb(), 1)}


def stage_regressor_shap() -> dict:
    import xgboost as xgb
    df = u.load_regime_dataset()
    cols = u.width_cols(df)
    X = df[cols].to_numpy()
    y = u.severity_target(df).to_numpy()
    with u.Timer("reg") as t:
        reg = u.xgb_reg_factory()()
        reg.fit(X, y)
    with u.Timer("shap_reg") as t2:
        reg.get_booster().predict(xgb.DMatrix(X), pred_contribs=True)
    return {"stage": "regressor+shap", "n": len(X), "n_features": len(cols),
            "wall_seconds": round(t.seconds + t2.seconds, 4),
            "train_seconds": round(t.seconds, 4), "shap_seconds": round(t2.seconds, 4),
            "peak_rss_mb": round(_peak_rss_mb(), 1)}


def stage_pc_bootstrap(n_bootstrap: int = 25) -> dict:
    """PC algorithm bootstrap stability (small n_bootstrap to keep runtime modest)."""
    try:
        from causallearn.search.ConstraintBased.PC import pc
    except Exception as e:
        return {"stage": "pc_bootstrap", "n": 0, "wall_seconds": float("nan"),
                "peak_rss_mb": float("nan"), "note": f"causal-learn unavailable: {e}"}
    df = u.load_regime_dataset()
    cols = u.curated_paper_panel(df)
    sev = u.severity_target(df).to_numpy()
    X = df[cols].to_numpy()
    # Stack severity as last column for joint structure
    data = np.column_stack([X, sev]).astype(float)
    rng = np.random.default_rng(42)
    with u.Timer("pc") as t:
        for _ in range(n_bootstrap):
            idx = rng.choice(len(data), size=len(data), replace=True)
            try:
                pc(data[idx], alpha=0.05, indep_test="fisherz", show_progress=False)
            except Exception:
                continue
    return {"stage": "pc_bootstrap", "n": n_bootstrap,
            "wall_seconds": round(t.seconds, 4),
            "peak_rss_mb": round(_peak_rss_mb(), 1)}


def main() -> None:
    print("[04 runtime] starting...")
    rows = []
    for fn in (stage_lhs_sample, stage_fba_batch, stage_targeted_fva,
               stage_classifier_shap, stage_regressor_shap, stage_pc_bootstrap):
        print(f"  - {fn.__name__}...")
        try:
            row = fn()
        except Exception as e:
            row = {"stage": fn.__name__, "wall_seconds": float("nan"),
                   "peak_rss_mb": float("nan"), "error": str(e)}
        print(f"    {row}")
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(OUT["ws"] / "runtime_summary.csv", index=False)
    print(f"[04 runtime] wrote {OUT['ws'] / 'runtime_summary.csv'}")

    env = u.env_summary()
    lines = ["# Environment summary", ""]
    for k, v in env.items():
        lines.append(f"{k}: {v}")
    (OUT["ws"] / "environment_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[04 runtime] wrote {OUT['ws'] / 'environment_summary.txt'}")

    # short summary md
    md = ["# Runtime profiling summary", "",
          f"- environment: {env['python']} on {env['platform']}",
          f"- CPUs: {env['cpu_count']}", ""]
    md.append("| Stage | n | wall (s) | peak RSS (MB) |")
    md.append("|---|---|---|---|")
    for r in rows:
        md.append(f"| {r.get('stage','')} | {r.get('n','')} | "
                  f"{r.get('wall_seconds','')} | {r.get('peak_rss_mb','')} |")
    (OUT["ws"] / "runtime_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
