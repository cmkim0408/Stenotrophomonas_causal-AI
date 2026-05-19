"""Per-row pFBA runner with resume support.

Usage:
  python3 _pfba_runner.py iso1 START END
  python3 _pfba_runner.py iml  START END

Writes per-row parquets to revision_runs/iscience_rev1/07_pointflux/parts_<sys>/.
"""
from __future__ import annotations

import re
import sys
import time
import warnings
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
import yaml
from cobra.flux_analysis import pfba

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from revision import utils as u

ROOT = u.ROOT
OUT_BASE = ROOT / "revision_runs" / "iscience_rev1" / "07_pointflux"
MODEL_PATH = ROOT / "acetate_xai" / "models" / "model.xml"
MEDIUM_YAML = ROOT / "acetate_xai" / "configs" / "medium_base.yaml"


def replay_iso1_row(model, row, med):
    mids = {r.id for r in model.reactions}
    for rid, b in med.get("base_bounds", {}).items():
        if rid in mids:
            model.reactions.get_by_id(rid).bounds = (float(b["lb"]), float(b["ub"]))
    ye = med.get("yeast_extract", {})
    if float(row.get("yeast_extract_gL", 0) or 0) > float(ye.get("enabled_if_gL_gt", 0)):
        for rid in ye.get("open_exchanges_when_enabled", []):
            if rid in mids:
                r_obj = model.reactions.get_by_id(rid)
                r_obj.bounds = (float(ye.get("open_uptake_lb", -1.0)),
                                max(r_obj.upper_bound, 1000.0))
    for pre in ("acetate", "oxygen", "ammonium", "phosphate"):
        rid = row[f"{pre}_rid"]
        if rid in mids:
            model.reactions.get_by_id(rid).bounds = (
                float(row[f"{pre}_lb"]), float(row[f"{pre}_ub"]))
    m_atpm = re.search(r"ATPM(\d+)", str(row["run_folder"]))
    if m_atpm and "ATPM" in mids:
        v = float(m_atpm.group(1))
        atpm = model.reactions.ATPM
        atpm.bounds = (v, max(atpm.upper_bound, 1000.0))


def run_iso1(start, end):
    out_dir = OUT_BASE / "parts_iso1"
    out_dir.mkdir(parents=True, exist_ok=True)
    iso1 = u.load_regime_dataset()
    panel_rids = [c.replace("width__", "") for c in u.curated_paper_panel(iso1)]
    print(f"[iso1] panel={len(panel_rids)} rxns; rows {start}..{end-1}")

    model = cobra.io.read_sbml_model(str(MODEL_PATH))
    with open(MEDIUM_YAML) as f:
        med = yaml.safe_load(f)
    mids = {r.id for r in model.reactions}
    panel_have = [r for r in panel_rids if r in mids]

    t0 = time.perf_counter()
    n_new = n_skip = n_fail = 0
    for idx in range(start, min(end, len(iso1))):
        part = out_dir / f"row_{idx:04d}.parquet"
        if part.exists():
            n_skip += 1
            continue
        row = iso1.iloc[idx]
        with model as mm:
            try:
                replay_iso1_row(mm, row, med)
                sol = mm.optimize()
                if sol.status != "optimal" or sol.objective_value is None:
                    n_fail += 1
                    continue
                psol = pfba(mm)
                fluxes = psol.fluxes
            except Exception as e:
                n_fail += 1
                print(f"  row {idx}: {e}")
                continue
        rec = {"row_idx": idx, "pfba_obj": float(psol.objective_value)}
        for rid in panel_have:
            rec[f"pfba__{rid}"] = float(fluxes.get(rid, np.nan))
        pd.DataFrame([rec]).to_parquet(part, index=False)
        n_new += 1
    print(f"[iso1] {start}..{end-1}: {n_new} new, {n_skip} skip, {n_fail} fail, "
          f"wall={time.perf_counter()-t0:.1f}s")


IML_UPTAKE = ("EX_glc__D_e", "EX_ac_e", "EX_o2_e", "EX_nh4_e", "EX_pi_e")


def run_iml(start, end):
    out_dir = OUT_BASE / "parts_iml"
    out_dir.mkdir(parents=True, exist_ok=True)
    transfer = pd.read_parquet(ROOT / "revision_runs" / "iscience_rev1" /
                                "03_transfer" / "transfer_dataset.parquet")
    iml_panel = [c.replace("width__", "") for c in transfer.columns
                 if c.startswith("width__")]
    print(f"[iml] panel={len(iml_panel)} rxns; rows {start}..{end-1}")

    model = cobra.io.load_model("iML1515")
    mids = {r.id for r in model.reactions}
    saved = {rid: model.reactions.get_by_id(rid).bounds for rid in IML_UPTAKE}
    panel_have = [r for r in iml_panel if r in mids]

    t0 = time.perf_counter()
    n_new = n_skip = n_fail = 0
    for idx in range(start, min(end, len(transfer))):
        part = out_dir / f"row_{idx:04d}.parquet"
        if part.exists():
            n_skip += 1
            continue
        row = transfer.iloc[idx]
        for rid in IML_UPTAKE:
            lb = float(row[f"{rid}_lb"])
            ub = saved[rid][1]
            model.reactions.get_by_id(rid).bounds = (lb, ub)
        try:
            psol = pfba(model)
            fluxes = psol.fluxes
        except Exception:
            n_fail += 1
            continue
        rec = {"row_idx": idx, "condition_id": int(row["condition_id"]),
               "pfba_obj": float(psol.objective_value)}
        for rid in panel_have:
            rec[f"pfba__{rid}"] = float(fluxes.get(rid, np.nan))
        pd.DataFrame([rec]).to_parquet(part, index=False)
        n_new += 1
    # restore
    for rid, b in saved.items():
        model.reactions.get_by_id(rid).bounds = b
    print(f"[iml] {start}..{end-1}: {n_new} new, {n_skip} skip, {n_fail} fail, "
          f"wall={time.perf_counter()-t0:.1f}s")


if __name__ == "__main__":
    which = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    if which == "iso1":
        run_iso1(start, end)
    elif which == "iml":
        run_iml(start, end)
    else:
        raise SystemExit(f"unknown system: {which}")
