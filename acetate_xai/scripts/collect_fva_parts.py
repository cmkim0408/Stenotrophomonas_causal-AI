from __future__ import annotations

import argparse
import logging
import sys

from acetate_xai.collect import collect_and_build_features
from acetate_xai.io import load_conditions_csv, save_table


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Collect distributed FVA parts and build XAI feature table.")
    p.add_argument("--parts-dir", required=True, help="Directory containing per-condition parquet parts")
    p.add_argument("--conditions", required=True, help="Conditions CSV (with measured_OD)")
    p.add_argument("--out", required=True, help="Output parquet for concatenated long table (fva_all)")
    p.add_argument("--features-out", required=True, help="Output parquet for wide feature matrix joined with conditions")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(message)s",
    )

    try:
        conditions_df = load_conditions_csv(args.conditions)
        fva_all, features = collect_and_build_features(parts_dir=args.parts_dir, conditions_df=conditions_df)
        save_table(fva_all, args.out, fmt="parquet")
        save_table(features, args.features_out, fmt="parquet")
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] {e}", file=sys.stderr)
        return 2

    n_conditions_parts = fva_all["condition_id"].nunique()
    print(f"[OK] Wrote: {args.out} (rows={len(fva_all)}, conditions={n_conditions_parts})")
    print(f"[OK] Wrote: {args.features_out} (rows={len(features)}, cols={len(features.columns)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

