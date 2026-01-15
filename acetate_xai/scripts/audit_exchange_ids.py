from __future__ import annotations

import argparse
import logging
import sys

from acetate_xai.config import load_config
from acetate_xai.io import load_sbml_model
from acetate_xai.audit import audit_exchange_ids, write_audit_csv


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Audit exchange IDs referenced in medium config vs model.")
    p.add_argument("--model", required=True, help="SBML model path (e.g., models/model.xml)")
    p.add_argument("--medium", required=True, help="Medium YAML config (e.g., configs/medium_base.yaml)")
    p.add_argument(
        "--out",
        default="results/audit_missing_exchanges.csv",
        help="Output CSV path (default: results/audit_missing_exchanges.csv)",
    )
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(message)s",
    )

    try:
        model = load_sbml_model(args.model)
        medium_cfg = load_config(args.medium)
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Failed to load inputs: {e}", file=sys.stderr)
        # Diagnostic tool: still exit 0 per requirement? This is "can't run" though.
        # We'll follow requirement strictly: always exit 0.
        return 0

    rows = audit_exchange_ids(model, medium_cfg)
    present = sum(1 for r in rows if r.status == "present")
    missing = [r for r in rows if r.status == "missing"]

    print(f"[REPORT] Requested exchange ids: {len(rows)}")
    print(f"[REPORT] Present: {present}")
    print(f"[REPORT] Missing: {len(missing)}")
    if missing:
        print("[REPORT] Missing IDs:")
        for r in missing:
            sug = ", ".join([x for x in [r.suggestion_1, r.suggestion_2, r.suggestion_3] if x])
            print(f"  - {r.requested_id}  | suggestions: {sug if sug else '(none)'}")

    out_path = write_audit_csv(rows, args.out)
    print(f"[OK] Wrote audit CSV: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

