"""Workstream #7: Point-flux vs flexibility-interval baseline (novelty defense).

Directly tests the manuscript's central novelty claim — that the **width**
of the FVA-derived feasible interval (flexibility representation) carries
more diagnostic signal than the **point flux** (whether parsimonious FBA
or FVA midpoint).

For each organism we evaluate XGBoost classifier + regressor on the SAME
curated panel of reactions, varying only the feature representation:

  - W   = width__rxn         (vmax - vmin)              [our method]
  - M   = mid__rxn           ((vmax + vmin) / 2)        [FVA midpoint = point center proxy]
  - PF  = pFBA flux on rxn   (parsimonious FBA point)   [true point-flux baseline]
  - PFA = |pFBA flux| on rxn (magnitude only)           [point-magnitude baseline]
  - OBJ = FBA objective_value only                       [scalar floor]

Both iSO1_933 (this study, n=242) and iML1515 (E. coli transfer, n=244)
are evaluated on their own curated panels (42 and 45 reactions).

Outputs:
    revision_runs/iscience_rev1/07_pointflux/pointflux_metrics.csv
    revision_runs/iscience_rev1/07_pointflux/pointflux_summary.md
    revision_runs/iscience_rev1/07_pointflux/pfba_fluxes_iso1.parquet
    revision_runs/iscience_rev1/07_pointflux/pfba_fluxes_iml1515.parquet
    revision_runs/iscience_rev1/figures/pointflux_vs_width.{png,pdf}
"""
from __future__ import annotations

import re
import sys
import time
import warnings
from pathlib import Path

import cobra
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from cobra.flux_analysis import pfba

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from revision import utils as u

OUT = u.ensure_outdirs("07_pointflux")
ROOT = u.ROOT
MODEL_PATH = ROOT / "acetate_xai" / "models" / "model.xml"
MEDIUM_YAML = ROOT / "acetate_xai" / "configs" / "medium_base.yaml"


# ---------------------------------------------------------------------------
# iSO1_933 — apply medium + LHS condition + ATPM/PI overrides (verified replay)
# ---------------------------------------------------------------------------

def _load_medium() -> dict:
    with open(MEDIUM_YAML) as f:
        return yaml.safe_load(f)


def replay_iso1_row(model: cobra.Model, row: pd.Series, med: dict) -> None:
    mids = {r.id for r in model.reactions}
    base = med.get("base_bounds", {})
    ye_cfg = med.get("yeast_extract", {})
    ye_open = ye_cfg.get("open_exchanges_when_enabled", [])
    ye_lb = float(ye_cfg.get("open_uptake_lb", -1.0))
    ye_th = float(ye_cfg.get("enabled_if_gL_gt", 0.0))

    for rid, b in base.items():
        if rid in mids:
            model.reactions.get_by_id(rid).bounds = (float(b["lb"]), float(b["ub"]))
    if float(row.get("yeast_extract_gL", 0) or 0) > ye_th:
        for rid in ye_open:
            if rid in mids:
                r_obj = model.reactions.get_by_id(rid)
                r_obj.bounds = (ye_lb, max(r_obj.upper_bound, 1000.0))
    for pre in ("acetate", "oxygen", "ammonium", "phosphate"):
        rid = row[f"{pre}_rid"]
        if rid in mids:
            model.reactions.get_by_id(rid).bounds = (
                float(row[f"{pre}_lb"]), float(row[f"{pre}_ub"]))
    rf = str(row["run_folder"])
    m_atpm = re.search(r"ATPM(\d+)", rf)
    if m_atpm and "ATPM" in mids:
        v = float(m_atpm.group(1))
        atpm = model.reactions.ATPM
        atpm.bounds = (v, max(atpm.upper_bound, 1000.0))


def run_pfba_iso1(df: pd.DataFrame, panel_rids: list[str]) -> pd.DataFrame:
    """Run pFBA on each of the 242 stored conditions and extract flux on the
    curated panel. Returns wide-format frame indexed by row index.
    """
    print(f"[pfba/iso1] solving {len(df)} conditions, panel={len(panel_rids)} rxns ...")
    model = cobra.io.read_sbml_model(str(MODEL_PATH))
    med = _load_medium()
    mids = {r.id for r in model.reactions}
    panel_have = [r for r in panel_rids if r in mids]
    rows = []
    failures = 0
    t0 = time.perf_counter()
    for idx, row in df.iterrows():
        with model as mm:
            try:
                replay_iso1_row(mm, row, med)
                sol = mm.optimize()
                if sol.status != "optimal" or sol.objective_value is None:
                    failures += 1
                    continue
                psol = pfba(mm)
                fluxes = psol.fluxes
            except Exception as e:
                failures += 1
                print(f"  row {idx}: pFBA failed ({e})")
                continue
        rec = {"row_idx": idx, "pfba_obj": float(psol.objective_value)}
        for rid in panel_have:
            rec[f"pfba__{rid}"] = float(fluxes.get(rid, np.nan))
        rows.append(rec)
        if (idx + 1) % 25 == 0:
            print(f"    {idx+1}/{len(df)}  elapsed={time.perf_counter()-t0:.0f}s")
    print(f"  done: {len(rows)} succeeded, {failures} failed, "
          f"wall={time.perf_counter()-t0:.0f}s")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# iML1515 — reuse 03_transfer condition matrix and re-solve with pFBA
# ---------------------------------------------------------------------------

IML_UPTAKE = ("EX_glc__D_e", "EX_ac_e", "EX_o2_e", "EX_nh4_e", "EX_pi_e")


def run_pfba_iml1515(transfer_df: pd.DataFrame, panel_rids: list[str]) -> pd.DataFrame:
    print(f"[pfba/iml1515] solving {len(transfer_df)} conditions, "
          f"panel={len(panel_rids)} rxns ...")
    model = cobra.io.load_model("iML1515")
    mids = {r.id for r in model.reactions}
    saved_bounds = {rid: model.reactions.get_by_id(rid).bounds for rid in IML_UPTAKE}
    panel_have = [r for r in panel_rids if r in mids]

    rows = []
    failures = 0
    t0 = time.perf_counter()
    for i, row in transfer_df.iterrows():
        for rid in IML_UPTAKE:
            lb = float(row[f"{rid}_lb"])
            ub = saved_bounds[rid][1]
            model.reactions.get_by_id(rid).bounds = (lb, ub)
        try:
            psol = pfba(model)
            fluxes = psol.fluxes
        except Exception as e:
            failures += 1
            continue
        rec = {"row_idx": i,
               "condition_id": int(row["condition_id"]),
               "pfba_obj": float(psol.objective_value)}
        for rid in panel_have:
            rec[f"pfba__{rid}"] = float(fluxes.get(rid, np.nan))
        rows.append(rec)
        if (i + 1) % 25 == 0:
            print(f"    {i+1}/{len(transfer_df)}  elapsed={time.perf_counter()-t0:.0f}s")
    # restore
    for rid, b in saved_bounds.items():
        model.reactions.get_by_id(rid).bounds = b
    print(f"  done: {len(rows)} succeeded, {failures} failed, "
          f"wall={time.perf_counter()-t0:.0f}s")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Evaluation harness
# ---------------------------------------------------------------------------

def eval_panel(X: np.ndarray, y: np.ndarray, sev: np.ndarray, n_splits: int = 5,
               clf_min: int = 2) -> dict:
    """Run both classification and regression CV; pad nan with zero."""
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    # Need at least clf_min folds based on smallest class
    cmin = pd.Series(y).value_counts().min()
    splits = max(2, min(n_splits, int(cmin)))
    clf_res = u.cv_classify(X, y, u.xgb_clf_factory(), n_splits=splits)
    reg_res = u.cv_regress(X, sev, u.xgb_reg_factory(), n_splits=n_splits)
    return {
        "macro_f1": clf_res["macro_f1"],
        "balanced_accuracy": clf_res["balanced_accuracy"],
        "rmse": reg_res["rmse"],
        "r2": reg_res["r2"],
        "per_class_f1": clf_res["per_class_f1"],
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("[07 pointflux] starting ...")

    # ============ iSO1_933 ============
    iso1 = u.load_regime_dataset()
    print(f"  iSO1 dataset: {iso1.shape}")
    curated_w = u.curated_paper_panel(iso1)         # ['width__EX_o2_e', ...]
    curated_rids = [c.replace("width__", "") for c in curated_w]
    print(f"  iSO1 curated panel: {len(curated_rids)} reactions")

    # Run / load pFBA fluxes (cache to parquet to avoid re-solving on rerun)
    pfba_iso1_path = OUT["ws"] / "pfba_fluxes_iso1.parquet"
    if pfba_iso1_path.exists():
        pfba_iso1 = pd.read_parquet(pfba_iso1_path)
        print(f"  iSO1 pFBA fluxes cached: {pfba_iso1.shape} ← {pfba_iso1_path}")
    else:
        pfba_iso1 = run_pfba_iso1(iso1, curated_rids)
        pfba_iso1.to_parquet(pfba_iso1_path, index=False)
        print(f"  iSO1 pFBA fluxes saved → {pfba_iso1_path}")

    # Build joined matrix (parquet row order preserved)
    pfba_iso1 = pfba_iso1.set_index("row_idx").reindex(iso1.index)
    pfba_cols_iso1 = [f"pfba__{r}" for r in curated_rids if f"pfba__{r}" in pfba_iso1.columns]

    # Feature matrices (same row count = 242)
    y_iso1, classes_iso1 = pd.factorize(iso1["label"])
    sev_iso1 = u.severity_target(iso1).to_numpy()

    feat_iso1: dict[str, tuple[np.ndarray, int]] = {}
    # W: width__
    w_present = [c for c in curated_w if c in iso1.columns]
    feat_iso1["W_widths"] = (iso1[w_present].astype(float).fillna(0).to_numpy(),
                              len(w_present))
    # M: mid__   (overlap subset)
    mid_cols_iso1 = [f"mid__{r}" for r in curated_rids if f"mid__{r}" in iso1.columns]
    feat_iso1["M_midpoints"] = (iso1[mid_cols_iso1].astype(float).fillna(0).to_numpy(),
                                 len(mid_cols_iso1))
    # PF: pFBA point flux
    feat_iso1["PF_pfba"] = (pfba_iso1[pfba_cols_iso1].astype(float).fillna(0).to_numpy(),
                             len(pfba_cols_iso1))
    # PFA: |pFBA|
    feat_iso1["PFA_pfba_abs"] = (np.abs(pfba_iso1[pfba_cols_iso1].astype(float)
                                          .fillna(0).to_numpy()),
                                  len(pfba_cols_iso1))
    # OBJ: scalar floor
    obj_only = iso1[["objective_value"]].astype(float).to_numpy()
    feat_iso1["OBJ_only"] = (obj_only, 1)

    rows_out: list[dict] = []
    for name, (X, n_feat) in feat_iso1.items():
        if X.shape[1] == 0:
            print(f"  iSO1 SKIP {name} — no features available")
            continue
        res = eval_panel(X, y_iso1, sev_iso1)
        rows_out.append({
            "system": "iSO1_933",
            "representation": name,
            "n_conditions": len(iso1),
            "n_features": n_feat,
            "macro_f1": round(res["macro_f1"], 4),
            "balanced_accuracy": round(res["balanced_accuracy"], 4),
            "rmse": round(res["rmse"], 5),
            "r2": round(res["r2"], 4),
            **{f"per_class_f1__{c}": round(v, 4)
               for c, v in zip(classes_iso1, res["per_class_f1"])},
        })
        print(f"  iSO1 {name:18s} (n_feat={n_feat:3d}) → "
              f"F1={res['macro_f1']:.4f}  R²={res['r2']:.4f}")

    # ============ iML1515 ============
    transfer_path = ROOT / "revision_runs" / "iscience_rev1" / "03_transfer" / "transfer_dataset.parquet"
    transfer = pd.read_parquet(transfer_path)
    print(f"\n  iML1515 transfer dataset: {transfer.shape}")
    iml_width_cols = [c for c in transfer.columns if c.startswith("width__")]
    iml_panel = [c.replace("width__", "") for c in iml_width_cols]

    pfba_iml_path = OUT["ws"] / "pfba_fluxes_iml1515.parquet"
    pfba_iml = None
    if pfba_iml_path.exists():
        pfba_iml = pd.read_parquet(pfba_iml_path)
        print(f"  iML1515 pFBA fluxes cached: {pfba_iml.shape} ← {pfba_iml_path}")
    else:
        # Try to compute; gracefully skip if iML1515 is unreachable (sandbox/offline)
        try:
            pfba_iml = run_pfba_iml1515(transfer, iml_panel)
            pfba_iml.to_parquet(pfba_iml_path, index=False)
            print(f"  iML1515 pFBA fluxes saved → {pfba_iml_path}")
        except Exception as e:
            print(f"  iML1515 pFBA UNAVAILABLE ({e.__class__.__name__}: {e})")
            print(f"    → run on machine with iML1515 cached:")
            print(f"      python3 code/revision/_pfba_runner.py iml 0 {len(transfer)}")
            pfba_iml = None

    y_iml, classes_iml = pd.factorize(transfer["regime"])
    classes_iml = list(classes_iml)
    sev_iml = (transfer["objective_value"] / transfer["objective_value"].max()).to_numpy()

    # Always evaluate W_widths and OBJ_only (no pFBA needed)
    feat_iml = {}
    feat_iml["W_widths"] = (transfer[iml_width_cols].astype(float).fillna(0).to_numpy(),
                              len(iml_width_cols))
    feat_iml["OBJ_only"] = (transfer[["objective_value"]].astype(float).to_numpy(), 1)
    if pfba_iml is not None:
        pfba_iml_aligned = pfba_iml.set_index("row_idx").reindex(transfer.index)
        pfba_cols_iml = [f"pfba__{r}" for r in iml_panel
                         if f"pfba__{r}" in pfba_iml_aligned.columns]
        feat_iml["PF_pfba"] = (pfba_iml_aligned[pfba_cols_iml].astype(float)
                                 .fillna(0).to_numpy(), len(pfba_cols_iml))
        feat_iml["PFA_pfba_abs"] = (np.abs(pfba_iml_aligned[pfba_cols_iml]
                                              .astype(float).fillna(0).to_numpy()),
                                     len(pfba_cols_iml))

    for name, (X, n_feat) in feat_iml.items():
        if X.shape[1] == 0:
            print(f"  iML1515 SKIP {name} — no features")
            continue
        res = eval_panel(X, y_iml, sev_iml)
        rows_out.append({
            "system": "iML1515",
            "representation": name,
            "n_conditions": len(transfer),
            "n_features": n_feat,
            "macro_f1": round(res["macro_f1"], 4),
            "balanced_accuracy": round(res["balanced_accuracy"], 4),
            "rmse": round(res["rmse"], 5),
            "r2": round(res["r2"], 4),
            **{f"per_class_f1__{c}": round(v, 4)
               for c, v in zip(classes_iml, res["per_class_f1"])},
        })
        print(f"  iML1515 {name:18s} (n_feat={n_feat:3d}) → "
              f"F1={res['macro_f1']:.4f}  R²={res['r2']:.4f}")

    metrics = pd.DataFrame(rows_out)
    metrics.to_csv(OUT["ws"] / "pointflux_metrics.csv", index=False)
    print(f"\n[07 pointflux] wrote metrics → {OUT['ws'] / 'pointflux_metrics.csv'}")

    # ----- Figure: grouped bar chart -----
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    REP_ORDER = ["OBJ_only", "PF_pfba", "PFA_pfba_abs", "M_midpoints", "W_widths"]
    REP_PALETTE = {
        "OBJ_only":     "#999999",
        "PF_pfba":      "#0072B2",   # point flux
        "PFA_pfba_abs": "#56B4E9",
        "M_midpoints":  "#E69F00",   # midpoint (center proxy)
        "W_widths":     "#009E73",   # interval width (OUR method)
    }
    REP_LABEL = {
        "OBJ_only":     "FBA obj only",
        "PF_pfba":      "pFBA point flux",
        "PFA_pfba_abs": "|pFBA flux|",
        "M_midpoints":  "FVA midpoint",
        "W_widths":     "FVA width (ours)",
    }

    for ax_idx, (metric, ylabel, ylim) in enumerate([
        ("macro_f1", "Macro-F1 (regime classification)", (0.0, 1.05)),
        ("r2",       "R² (severity regression)",         (-0.1, 1.05)),
    ]):
        ax = axes[ax_idx]
        systems = ["iSO1_933", "iML1515"]
        n_reps = len(REP_ORDER)
        bar_w = 0.16
        x_base = np.arange(len(systems))
        for j, rep in enumerate(REP_ORDER):
            sub = metrics[metrics["representation"] == rep].set_index("system")
            vals = [sub.loc[s, metric] if s in sub.index else np.nan for s in systems]
            offsets = (j - (n_reps - 1) / 2) * bar_w
            xs = x_base + offsets
            ax.bar(xs, vals, width=bar_w, color=REP_PALETTE[rep], edgecolor="black",
                   lw=0.4, label=REP_LABEL[rep])
            for xi, v in zip(xs, vals):
                if v is not None and not np.isnan(v):
                    ax.text(xi, v + 0.015, f"{v:.3f}", ha="center", fontsize=6.5,
                            rotation=90)
        ax.set_xticks(x_base)
        ax.set_xticklabels(systems)
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        ax.grid(axis="y", alpha=0.3)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=n_reps,
               bbox_to_anchor=(0.5, 1.04), fontsize=8.5, frameon=False)
    fig.suptitle("Novelty defense: point-flux vs flexibility-interval representation "
                 "(5-fold CV, same curated panel, same XGBoost)",
                 fontsize=10.5, y=1.10)
    fig.tight_layout()
    fig.savefig(OUT["figures"] / "pointflux_vs_width.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT["figures"] / "pointflux_vs_width.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[07 pointflux] wrote figure → {OUT['figures'] / 'pointflux_vs_width.png'}")

    # ----- Markdown summary -----
    md = [
        "# Point-flux vs flexibility-interval baseline — summary",
        "",
        "**Question:** does the *width* of the FVA-derived feasible interval "
        "(our flexibility representation) carry more diagnostic signal than the "
        "*point flux* (pFBA solution, FVA midpoint, FBA objective)?",
        "",
        "Each representation is evaluated on the SAME curated reaction panel, "
        "with the SAME XGBoost classifier and regressor under 5-fold CV. The "
        "only thing that changes is how the panel's reactions are encoded as "
        "features.",
        "",
        "## Feature representations",
        "",
        "- **W (widths)** — `width__rxn = vmax − vmin` from targeted FVA "
        "(`fraction_of_optimum=0.95`). This is the manuscript's flexibility "
        "feature.",
        "- **M (midpoints)** — `mid__rxn = (vmax + vmin) / 2` from the same FVA "
        "solve. A point-center proxy that comes from the same data as W.",
        "- **PF (pFBA flux)** — parsimonious FBA flux for each reaction. The "
        "gold-standard \"point flux\" baseline.",
        "- **PFA (|pFBA flux|)** — magnitude of pFBA flux (sign-agnostic, "
        "directly comparable in scale to widths).",
        "- **OBJ** — biomass objective only (1 feature). Scalar floor.",
        "",
        "## Headline results (5-fold CV, same XGBoost hyperparameters)",
        "",
        "```",
        metrics.to_string(index=False),
        "```",
        "",
        "## Interpretation (conservative framing)",
        "",
        ("- **W vs PF on the same panel.** When the FVA-derived width "
         "representation is replaced by the pFBA point-flux (or its magnitude), "
         "the macro-F1 and severity R² values change — directly testing whether "
         "the diagnostic signal lives in the *interval* or in the *point*. The "
         "comparison is fair because the underlying reactions and conditions "
         "are identical."),
        ("- **W vs M (within the FVA solve).** Width and midpoint use the same "
         "FVA solve but encode complementary information: the midpoint is a "
         "single representative flux value, whereas the width is the size of "
         "the feasible interval. Comparing them isolates the contribution of "
         "*the interval itself*, controlling for the FVA solve and the "
         "underlying network."),
        ("- **OBJ floor.** Reporting the FBA objective alone establishes how "
         "much of the signal is trivially encoded in biomass rate."),
        "",
        "## Caveats (mandatory transparency)",
        "",
        ("- On the iSO1 dataset, the `mid__` column universe is the original "
         "120-reaction FVA campaign (alphabetically truncated at FACOAL161), so "
         "the M_midpoints panel is necessarily smaller than the W_widths panel; "
         "the dominant comparison is therefore W vs PF/PFA (full panel match), "
         "with M reported as a control."),
        ("- The pFBA fluxes are deterministic — pFBA returns a unique flux "
         "vector that minimizes total absolute flux subject to optimal biomass. "
         "This eliminates non-uniqueness as a confound and gives the *strongest "
         "possible* point-flux baseline."),
        ("- The same hyperparameters and CV split are used across all "
         "representations, so any performance differences are attributable to "
         "feature representation, not learner tuning."),
        "",
        "## Files",
        "",
        "- `pointflux_metrics.csv` -- full headline table",
        "- `pfba_fluxes_iso1.parquet` -- cached iSO1 pFBA solutions (242 conditions x panel)",
        "- `pfba_fluxes_iml1515.parquet` -- cached iML1515 pFBA solutions (244 conditions x panel)",
        "- `../figures/pointflux_vs_width.{png,pdf}` -- grouped bar comparison",
        "",
    ]
    (OUT["ws"] / "pointflux_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"[07 pointflux] wrote summary -> {OUT['ws'] / 'pointflux_summary.md'}")


if __name__ == "__main__":
    main()
