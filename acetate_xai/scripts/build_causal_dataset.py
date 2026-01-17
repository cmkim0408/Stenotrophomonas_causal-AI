from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Build causal_dataset.parquet from regime_dataset.parquet + parsed interventions.")
    p.add_argument("--in", dest="in_path", default="results/regime_dataset.parquet", help="Input parquet path")
    p.add_argument("--out", dest="out_path", default="results/causal_dataset.parquet", help="Output parquet path")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def _parse_atpm_fixed(run_id: str) -> float | None:
    # Examples: ...__ATPM20..., ATPM5, run__baseC3__ATPM10
    m = re.search(r"ATPM(\d+)", run_id, flags=re.IGNORECASE)
    if not m:
        return None
    return float(m.group(1))


def _parse_o2_lb(run_id: str) -> float | None:
    # Examples: o2lb-100, o2lb-50, o2lb-20 (meaning lower bound = -100)
    m = re.search(r"o2lb(-\d+)", run_id, flags=re.IGNORECASE)
    if not m:
        return None
    return float(m.group(1))


def _parse_frac_opt(run_id: str) -> float | None:
    # Example: f099 -> 0.99, f095 -> 0.95
    m = re.search(r"\bf(\d{3})\b", run_id, flags=re.IGNORECASE)
    if not m:
        return None
    v = int(m.group(1))
    return v / 100.0


def _parse_targets_n(run_id: str) -> float | None:
    # Example: 300targets
    m = re.search(r"(\d+)\s*targets", run_id, flags=re.IGNORECASE)
    if not m:
        return None
    return float(int(m.group(1)))


def _parse_acetate_mode(run_id: str) -> str | None:
    # Examples: acFree, acScale1.0, acetate cap tokens (heuristic)
    if re.search(r"\bacfree\b", run_id, flags=re.IGNORECASE):
        return "acFree"
    m = re.search(r"\bacscale([0-9.]+)\b", run_id, flags=re.IGNORECASE)
    if m:
        return f"acScale{m.group(1)}"
    if re.search(r"accap|acetatecap|acetate_cap", run_id, flags=re.IGNORECASE):
        return "acCap"
    return None


def _parse_nh4_mode(run_id: str) -> str | None:
    # Examples: nh4cap-100, k_nh4 tokens
    m = re.search(r"\bnh4cap(-?\d+)\b", run_id, flags=re.IGNORECASE)
    if m:
        return f"nh4cap{m.group(1)}"
    if re.search(r"\bk_nh4\b", run_id, flags=re.IGNORECASE):
        return "k_nh4"
    return None


def _compute_severity(df: pd.DataFrame, *, run_col: str = "run_id", obj_col: str = "objective_value") -> pd.Series:
    obj = pd.to_numeric(df[obj_col], errors="coerce")
    run_max = obj.groupby(df[run_col].astype(str)).max()
    max_obj = df[run_col].astype(str).map(run_max)
    den = pd.to_numeric(max_obj, errors="coerce").replace({0.0: np.nan})
    sev = obj / den
    # keep in [0,1] when possible
    sev = sev.where(np.isfinite(sev), np.nan)
    sev = sev.clip(lower=0.0, upper=1.0)
    return sev


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    in_path = Path(args.in_path)
    out_path = Path(args.out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not in_path.exists():
        print(f"[ERROR] input not found: {in_path}")
        return 2

    df = pd.read_parquet(in_path)
    for col in ("run_id", "objective_value", "label"):
        if col not in df.columns:
            print(f"[ERROR] required column missing in input: {col}")
            return 2

    run_id = df["run_id"].astype(str)
    df_out = df.copy()
    df_out["atpm_fixed"] = run_id.map(_parse_atpm_fixed)
    df_out["o2_lb"] = run_id.map(_parse_o2_lb)
    df_out["acetate_mode"] = run_id.map(_parse_acetate_mode)
    df_out["nh4_mode"] = run_id.map(_parse_nh4_mode)
    df_out["frac_opt"] = run_id.map(_parse_frac_opt)
    df_out["targets_n"] = run_id.map(_parse_targets_n)

    df_out["maintenance_severity"] = _compute_severity(df_out, run_col="run_id", obj_col="objective_value")

    # primary_regime is the existing label
    df_out["primary_regime"] = df_out["label"]

    df_out.to_parquet(out_path, index=False)
    print(f"[OK] Wrote {out_path} (rows={len(df_out)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

