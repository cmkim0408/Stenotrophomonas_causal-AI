"""Extended FVA campaign — close the width__ universe truncation.

The deployed `regime_dataset.parquet` width__ columns are alphabetically
truncated at FACOAL161, so paper-named TCA / glyoxylate / glycolysis /
respiration anchors beyond 'F' (MDH, ICDH, ICL, MALS, PYK, PPC, NADH16pp, …)
are absent. This script:

  1. Re-applies each of the 242 (campaign × run_folder × condition_id) rows
     using stored uptake bounds + ATPM override (parsed from run_folder).
     Replay is 100% verified against stored objective_value (tested 2026-04-18).
  2. Runs targeted FVA (fraction_of_optimum=0.95, GLPK) on the 174 reactions
     that are in `targets_300.json` but absent from the deployed parquet.
  3. Saves per-row partial parquets (resumable) under
     revision_runs/iscience_rev1/extended_fva/parts/.
  4. Final merge → revision_runs/iscience_rev1/extended_fva/widths_extended.parquet.

After this script, run `06_build_extended_dataset.py` to assemble the merged
`regime_dataset_extended.parquet`.
"""
from __future__ import annotations

import json
import re
import sys
import time
import warnings
from pathlib import Path

import cobra
import pandas as pd
import yaml
from cobra.flux_analysis import flux_variability_analysis as fva

warnings.filterwarnings("ignore", category=UserWarning)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from revision import utils as u

ROOT = u.ROOT
MODEL_PATH = ROOT / "acetate_xai" / "models" / "model.xml"
MEDIUM_YAML = ROOT / "acetate_xai" / "configs" / "medium_base.yaml"
TARGETS_300 = ROOT / "acetate_xai" / "configs" / "targets_300.json"

OUT_BASE = ROOT / "revision_runs" / "iscience_rev1" / "extended_fva"
PARTS_DIR = OUT_BASE / "parts"
PARTS_DIR.mkdir(parents=True, exist_ok=True)


def replay_row(model: cobra.Model, row: pd.Series, base_bounds: dict,
               ye_open_list: list[str], ye_open_lb: float, ye_threshold: float) -> None:
    """Apply medium + per-row uptake bounds + ATPM override (in place)."""
    mids = {r.id for r in model.reactions}
    # 1. Base medium bounds
    for rid, b in base_bounds.items():
        if rid in mids:
            model.reactions.get_by_id(rid).bounds = (float(b["lb"]), float(b["ub"]))
    # 2. Yeast extract opens (vitamin/cofactor exchanges)
    if float(row.get("yeast_extract_gL", 0) or 0) > ye_threshold:
        for rid in ye_open_list:
            if rid in mids:
                r_obj = model.reactions.get_by_id(rid)
                r_obj.bounds = (ye_open_lb, max(r_obj.upper_bound, 1000.0))
    # 3. Per-row stored uptake bounds (FINAL values; bypass scaling logic)
    for prefix in ("acetate", "oxygen", "ammonium", "phosphate"):
        rid = row[f"{prefix}_rid"]
        if rid in mids:
            model.reactions.get_by_id(rid).bounds = (
                float(row[f"{prefix}_lb"]), float(row[f"{prefix}_ub"]))
    # 4. ATPM lower-bound override (parsed from run_folder name)
    rf = str(row["run_folder"])
    m_atpm = re.search(r"ATPM(\d+)", rf)
    if m_atpm and "ATPM" in mids:
        v = float(m_atpm.group(1))
        atpm = model.reactions.ATPM
        atpm.bounds = (v, max(atpm.upper_bound, 1000.0))


def main() -> None:
    print("[extend_fva] starting...")
    df = u.load_regime_dataset()
    print(f"  parent dataset: {df.shape}")

    # Load model + medium + 300-target list
    model = cobra.io.read_sbml_model(str(MODEL_PATH))
    mids = {r.id for r in model.reactions}
    print(f"  loaded model: {len(mids)} reactions")

    with open(MEDIUM_YAML) as f:
        med = yaml.safe_load(f)
    base_bounds = med.get("base_bounds", {})
    ye_cfg = med.get("yeast_extract", {})
    ye_open = ye_cfg.get("open_exchanges_when_enabled", [])
    ye_lb = float(ye_cfg.get("open_uptake_lb", -1.0))
    ye_th = float(ye_cfg.get("enabled_if_gL_gt", 0.0))

    with open(TARGETS_300) as f:
        targets = json.load(f)

    # Existing widths in parquet → identify NEW reactions to compute
    present = {c.split("__", 1)[1] for c in u.width_cols(df)}
    new_rxns = [r for r in targets if r in mids and r not in present]
    print(f"  NEW reactions to add: {len(new_rxns)}")
    print(f"  examples: {new_rxns[:8]} ...")

    # Per-row FVA loop (resumable: skip rows whose part-parquet already exists)
    t0 = time.perf_counter()
    n_done = 0
    n_skipped = 0
    n_failed = 0
    for idx, row in df.iterrows():
        part_path = PARTS_DIR / f"row_{idx:04d}.parquet"
        if part_path.exists():
            n_skipped += 1
            continue
        with model as m:
            try:
                replay_row(m, row, base_bounds, ye_open, ye_lb, ye_th)
                # Verify FBA still feasible
                sol = m.optimize()
                if sol.status != "optimal" or sol.objective_value is None:
                    n_failed += 1
                    print(f"  row {idx}: FBA failed (status={sol.status})")
                    continue
                obj_replay = float(sol.objective_value)
                obj_stored = float(row["objective_value"])
                if abs(obj_replay - obj_stored) > 0.05:
                    print(f"  row {idx}: WARN obj mismatch replay={obj_replay:.4f} "
                          f"stored={obj_stored:.4f}")
                # Targeted FVA on new reactions
                fva_df = fva(m, reaction_list=new_rxns, fraction_of_optimum=0.95,
                             processes=1, loopless=False)
            except Exception as e:
                n_failed += 1
                print(f"  row {idx}: ERR {e}")
                continue
        # Save partial parquet: long format (row_idx, reaction_id, fva_min, fva_max)
        part = pd.DataFrame({
            "row_idx": idx,
            "reaction_id": fva_df.index.tolist(),
            "fva_min": fva_df["minimum"].astype(float).tolist(),
            "fva_max": fva_df["maximum"].astype(float).tolist(),
        })
        part.to_parquet(part_path, index=False)
        n_done += 1

        if (n_done + n_skipped) % 10 == 0:
            elapsed = time.perf_counter() - t0
            eta = elapsed / max(n_done, 1) * (len(df) - n_done - n_skipped)
            print(f"  progress: {n_done + n_skipped}/{len(df)}  "
                  f"({n_done} new, {n_skipped} skipped, {n_failed} failed)  "
                  f"elapsed={elapsed:.0f}s  ETA={eta:.0f}s")

    elapsed = time.perf_counter() - t0
    print(f"\n[extend_fva] FVA loop done: {n_done} computed, {n_skipped} skipped, "
          f"{n_failed} failed.  Total wall: {elapsed:.0f}s")

    # Merge all parts into a single wide-format parquet
    parts = sorted(PARTS_DIR.glob("row_*.parquet"))
    long = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
    long["width"] = long["fva_max"] - long["fva_min"]

    # Pivot to wide: row_idx × reaction_id → width
    wide = (long.pivot(index="row_idx", columns="reaction_id", values="width")
                 .add_prefix("width__")
                 .reset_index())
    out_wide = OUT_BASE / "widths_extended.parquet"
    wide.to_parquet(out_wide, index=False)
    print(f"[extend_fva] wrote wide widths: {wide.shape} → {out_wide}")
    print(f"  new width__ columns: "
          f"{[c for c in wide.columns if c.startswith('width__')][:8]} ...")


if __name__ == "__main__":
    main()
