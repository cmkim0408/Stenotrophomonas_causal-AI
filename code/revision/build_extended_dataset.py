"""Merge widths_extended.parquet into regime_dataset_extended.parquet.

Joins by row_idx (= original parquet's positional index). The result has the
same 242 rows but ~300 width__ columns instead of 120.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from revision import utils as u

ROOT = u.ROOT
PARENT = ROOT / "results" / "regime_dataset.parquet"
NEW = ROOT / "revision_runs" / "iscience_rev1" / "extended_fva" / "widths_extended.parquet"
OUT = ROOT / "revision_runs" / "iscience_rev1" / "regime_dataset_extended.parquet"


def main() -> None:
    print("[merge] starting...")
    parent = pd.read_parquet(PARENT)
    new = pd.read_parquet(NEW)
    print(f"  parent: {parent.shape}  ({len([c for c in parent.columns if c.startswith('width__')])} width__)")
    print(f"  new:    {new.shape}     ({len([c for c in new.columns if c.startswith('width__')])} new width__)")

    # parent's positional index aligns with new.row_idx
    parent = parent.reset_index(drop=False).rename(columns={"index": "row_idx"})
    merged = parent.merge(new, on="row_idx", how="left", validate="one_to_one")
    merged = merged.drop(columns=["row_idx"])
    n_widths = len([c for c in merged.columns if c.startswith("width__")])
    print(f"  merged: {merged.shape}  ({n_widths} total width__)")
    merged.to_parquet(OUT, index=False)
    print(f"[merge] wrote {OUT}")

    # Sanity check: paper-named anchors now present
    paper_named = ["MDH", "MDH2", "MDH3", "ICDHyr", "ICDHx", "ICL", "MALS",
                   "PYK", "PYK3", "PPC", "NADH16pp", "FUM"]
    present = {c.split("__", 1)[1] for c in merged.columns if c.startswith("width__")}
    print(f"\n  Paper-named anchors now present:")
    for p in paper_named:
        mark = "[ok]" if p in present else "[--]"
        print(f"    {mark} {p}")


if __name__ == "__main__":
    main()
