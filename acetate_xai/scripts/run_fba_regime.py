from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from acetate_xai.config import load_config
from acetate_xai.io import load_conditions_csv, save_table
from acetate_xai.regime import run_fba_regime_table


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run FBA-based regime table (flux+bound saturation) per condition.")
    p.add_argument("--model", required=True, help="SBML model path (e.g., acetate_xai/models/model.xml)")
    p.add_argument("--conditions", required=True, help="Conditions CSV (long format)")
    p.add_argument("--medium", required=True, help="Medium YAML config")
    p.add_argument("--regime-config", required=True, help="Regime exchanges YAML config")
    p.add_argument("--out", required=True, help="Output parquet path (e.g., results/xai/regime_fba.parquet)")
    p.add_argument("--limit", type=int, default=None, help="Run only first N conditions")
    p.add_argument("--condition-ids", nargs="+", default=None, help="Run only specified condition_id values")
    p.add_argument("--eps", type=float, default=1e-6, help="Saturation tolerance")
    p.add_argument("--infty-bound", type=float, default=999.0, help="|bound| >= this considered open/infinite")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    if not Path(args.model).exists():
        print(f"[ERROR] model not found: {args.model}", file=sys.stderr)
        return 2

    try:
        conditions_df = load_conditions_csv(args.conditions)
        medium_cfg = load_config(args.medium)
        regime_cfg = load_config(args.regime_config)
        out = run_fba_regime_table(
            model_path=str(args.model),
            conditions_df=conditions_df,
            medium_cfg=medium_cfg,
            regime_cfg=regime_cfg,
            limit=args.limit,
            condition_ids=args.condition_ids,
            eps=float(args.eps),
            infty_bound=float(args.infty_bound),
        )
        save_table(out, args.out, fmt="parquet")
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] {e}", file=sys.stderr)
        return 2

    print(f"[OK] Wrote {len(out)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

