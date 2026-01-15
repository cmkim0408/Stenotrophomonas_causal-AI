from __future__ import annotations

import argparse
import sys

from acetate_xai.io import load_sbml_model
from acetate_xai.targets import (
    TargetSelectionError,
    load_anchors_yaml,
    save_targets_json,
    select_targets_anchored,
)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Select targeted FVA reactions (anchors + auto-fill).")
    p.add_argument("--model", required=True, help="SBML model path (e.g., models/model.xml)")
    p.add_argument("--anchors", required=True, help="anchors.yaml path (e.g., configs/anchors.yaml)")
    p.add_argument("--out", required=True, help="Output JSON (e.g., configs/targets_120.json)")
    p.add_argument("--target-count", type=int, default=120, help="Final number of targets.")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        model = load_sbml_model(args.model)
    except FileNotFoundError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 2
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Failed to load model: {e}", file=sys.stderr)
        return 2

    try:
        anchors = load_anchors_yaml(args.anchors)
        targets = select_targets_anchored(model=model, anchors=anchors, target_count=args.target_count)
        save_targets_json(targets, args.out)
    except (FileNotFoundError, ValueError, TargetSelectionError) as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 2
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Unexpected failure: {e}", file=sys.stderr)
        return 2

    print(f"[OK] Wrote {len(targets)} targets to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

