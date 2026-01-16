from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from acetate_xai.config import load_config
from acetate_xai.io import load_conditions_csv, load_sbml_model, save_table
from acetate_xai.medium import apply_condition_to_model
from acetate_xai.regime import compute_saturation_for_reaction, pick_first_existing_reaction_id


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Run FBA-based regime table (flux+bound saturation) per condition with fixed-flux constraints."
    )
    # Same args as scripts/run_fba_regime.py
    p.add_argument("--model", required=True, help="SBML model path (e.g., acetate_xai/models/model.xml)")
    p.add_argument("--conditions", required=True, help="Conditions CSV (long format)")
    p.add_argument("--medium", required=True, help="Medium YAML config")
    p.add_argument("--regime-config", required=True, help="Regime exchanges YAML config")
    p.add_argument("--out", required=True, help="Output parquet path")
    p.add_argument("--limit", type=int, default=None, help="Run only first N conditions")
    p.add_argument("--condition-ids", nargs="+", default=None, help="Run only specified condition_id values")
    p.add_argument("--eps", type=float, default=1e-6, help="Saturation tolerance")
    p.add_argument("--infty-bound", type=float, default=999.0, help="|bound| >= this considered open/infinite")

    # New
    p.add_argument(
        "--rxn-fix",
        action="append",
        default=[],
        help='Fix reaction flux: "RXNID=value". Repeatable. Example: --rxn-fix ATPM=20',
    )
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def _parse_rxn_fix(items: list[str]) -> list[tuple[str, float]]:
    out: list[tuple[str, float]] = []
    for raw in items:
        if raw is None:
            continue
        s = str(raw).strip()
        if not s:
            continue
        if "=" not in s:
            raise ValueError(f"--rxn-fix must be in form RXNID=value, got: {s!r}")
        rid, val = s.split("=", 1)
        rid = rid.strip()
        if not rid:
            raise ValueError(f"--rxn-fix has empty reaction id: {s!r}")
        try:
            v = float(val.strip())
        except Exception as e:  # noqa: BLE001
            raise ValueError(f"--rxn-fix has non-numeric value: {s!r}") from e
        out.append((rid, v))
    return out


def _apply_rxn_fix(model, rxn_fix: list[tuple[str, float]]) -> None:
    present = set(r.id for r in model.reactions)
    for rid, value in rxn_fix:
        if rid not in present:
            logging.warning("Reaction not found for --rxn-fix (skipped): %s", rid)
            continue
        rxn = model.reactions.get_by_id(rid)
        # Make sure we don't hit transient bound invalidity
        rxn.upper_bound = 1000.0
        rxn.lower_bound = float(value)
        rxn.upper_bound = float(value)
        logging.info("Fixed flux %s = %s", rid, value)


def _candidates(regime_cfg: dict, key: str) -> list[str]:
    c = regime_cfg.get(key, [])
    if isinstance(c, list):
        return [str(x) for x in c]
    return []


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    if not Path(args.model).exists():
        print(f"[ERROR] model not found: {args.model}", file=sys.stderr)
        return 2
    if not Path(args.medium).exists():
        print(f"[ERROR] medium yaml not found: {args.medium}", file=sys.stderr)
        return 2
    if not Path(args.regime_config).exists():
        print(f"[ERROR] regime-config yaml not found: {args.regime_config}", file=sys.stderr)
        return 2

    try:
        conditions_df = load_conditions_csv(args.conditions)
        medium_cfg = load_config(args.medium)
        regime_cfg = load_config(args.regime_config)
        rxn_fix = _parse_rxn_fix(list(args.rxn_fix or []))
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Failed to load inputs: {e}", file=sys.stderr)
        return 2

    df = conditions_df.copy()
    if args.condition_ids:
        wanted = set(str(x) for x in args.condition_ids)
        df = df.loc[df["condition_id"].astype(str).isin(wanted)].copy()
        if df.empty:
            example = ", ".join(conditions_df["condition_id"].astype(str).head(10).tolist())
            print(f"[ERROR] None of --condition-ids found. Example ids: {example}", file=sys.stderr)
            return 2
    elif args.limit is not None:
        if args.limit <= 0:
            print("[ERROR] --limit must be a positive integer", file=sys.stderr)
            return 2
        df = df.head(int(args.limit)).copy()

    nutrients = {
        "acetate": _candidates(regime_cfg, "acetate"),
        "oxygen": _candidates(regime_cfg, "oxygen"),
        "ammonium": _candidates(regime_cfg, "ammonium"),
        "phosphate": _candidates(regime_cfg, "phosphate"),
    }

    out_rows: list[dict] = []
    for _, row in df.iterrows():
        cid = str(row["condition_id"])

        model = load_sbml_model(args.model)
        apply_condition_to_model(model, row.to_dict(), medium_cfg)
        _apply_rxn_fix(model, rxn_fix)

        sol = model.optimize()
        if sol.status != "optimal":
            rec: dict = {"condition_id": cid, "objective_value": np.nan}
            # Still keep nutrient columns stable
            for nutrient in nutrients.keys():
                rec.update(
                    {
                        f"{nutrient}_rid": "",
                        f"{nutrient}_flux": np.nan,
                        f"{nutrient}_lb": np.nan,
                        f"{nutrient}_ub": np.nan,
                        f"{nutrient}_is_constrained": False,
                        f"{nutrient}_sat": False,
                        f"{nutrient}_sat_side": "missing",
                    }
                )
            out_rows.append(rec)
            continue

        rec = {"condition_id": cid, "objective_value": float(sol.objective_value)}
        for nutrient, cand in nutrients.items():
            rid_used = pick_first_existing_reaction_id(model, cand) if cand else None
            if rid_used is None:
                rec.update(
                    {
                        f"{nutrient}_rid": "",
                        f"{nutrient}_flux": np.nan,
                        f"{nutrient}_lb": np.nan,
                        f"{nutrient}_ub": np.nan,
                        f"{nutrient}_is_constrained": False,
                        f"{nutrient}_sat": False,
                        f"{nutrient}_sat_side": "missing",
                    }
                )
                continue

            sat = compute_saturation_for_reaction(
                rid=rid_used,
                model=model,
                solution=sol,
                eps=float(args.eps),
                infty_bound=float(args.infty_bound),
            )
            rec.update(
                {
                    f"{nutrient}_rid": sat.rid,
                    f"{nutrient}_flux": sat.flux,
                    f"{nutrient}_lb": sat.lb,
                    f"{nutrient}_ub": sat.ub,
                    f"{nutrient}_is_constrained": sat.is_constrained,
                    f"{nutrient}_sat": sat.saturated,
                    f"{nutrient}_sat_side": sat.sat_side,
                }
            )

        out_rows.append(rec)

    out_df = pd.DataFrame(out_rows)
    save_table(out_df, args.out, fmt="parquet")
    print(f"[OK] Wrote {len(out_df)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

