from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd

from acetate_xai.config import load_config
from acetate_xai.fva import FVAError, run_targeted_fva
from acetate_xai.io import load_conditions_csv, load_sbml_model, save_table
from acetate_xai.medium import MediumConfigError, apply_condition_to_model


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run FBA + targeted FVA for one condition.")
    p.add_argument("--model", required=True, help="SBML model path (e.g., models/model.xml)")
    p.add_argument("--conditions", required=True, help="Conditions CSV (long format)")
    p.add_argument("--condition-id", required=True, help="condition_id to run (must exist in CSV)")
    p.add_argument("--targets", required=True, help="targets JSON (array of reaction ids)")
    p.add_argument("--medium", required=True, help="medium YAML config (e.g., configs/medium_base.yaml)")
    p.add_argument("--out", required=True, help="Output parquet path (e.g., results/fva_one.parquet)")
    p.add_argument("--fraction", type=float, default=0.95, help="fraction_of_optimum for FVA")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def _load_targets_json(path: str | Path) -> list[str]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Targets JSON not found: {p}")
    with p.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, list) or not data:
        raise ValueError("Targets JSON must be a non-empty JSON array of reaction ids.")
    return [str(x) for x in data]


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(message)s",
    )

    try:
        model = load_sbml_model(args.model)
        cfg_medium = load_config(args.medium)
        targets = _load_targets_json(args.targets)
        cond_df = load_conditions_csv(args.conditions)
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Failed to load inputs: {e}", file=sys.stderr)
        return 2

    rows = cond_df.loc[cond_df["condition_id"] == args.condition_id]
    if rows.empty:
        available = ", ".join(cond_df["condition_id"].astype(str).head(10).tolist())
        print(
            f"[ERROR] condition_id not found: {args.condition_id}. Example available ids: {available}",
            file=sys.stderr,
        )
        return 2
    row = rows.iloc[0].to_dict()

    try:
        apply_condition_to_model(model, row, cfg_medium)
    except MediumConfigError as e:
        print(f"[ERROR] medium config/model mismatch: {e}", file=sys.stderr)
        return 2
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] failed to apply medium: {e}", file=sys.stderr)
        return 2

    sol = model.optimize()
    if sol.status != "optimal":
        print(f"[ERROR] FBA failed: status={sol.status}", file=sys.stderr)
        return 2

    objective_value = float(sol.objective_value)
    logging.info("FBA optimal objective_value=%.6g", objective_value)

    try:
        fva_df = run_targeted_fva(model, targets=targets, fraction_of_optimum=args.fraction)
    except FVAError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 2

    out = fva_df.copy()
    out.insert(0, "condition_id", str(args.condition_id))
    out.insert(1, "objective_value", objective_value)

    # Ensure schema order
    out = out[["condition_id", "objective_value", "reaction_id", "fva_min", "fva_max"]]

    save_table(out, args.out, fmt="parquet")
    print(f"[OK] Wrote {len(out)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

