from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from acetate_xai.xai import (
    build_regime_table,
    default_conditions_and_medium,
    fit_elasticnet_coefficients,
    fit_tree_rules,
)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Generate XAI reports from features.parquet.")
    p.add_argument("--features", required=True, help="Path to results/features.parquet")
    p.add_argument("--outdir", required=True, help="Output directory (e.g., results/xai)")
    p.add_argument(
        "--conditions",
        default="data/conditions_experiment.csv",
        help="Conditions CSV path (default: data/conditions_experiment.csv)",
    )
    p.add_argument(
        "--medium",
        default="configs/medium_base.yaml",
        help="Medium YAML path (default: configs/medium_base.yaml)",
    )
    p.add_argument(
        "--regime-fba",
        default=None,
        help="Optional parquet with FBA-based regime saturation to merge/override (e.g., results/xai/regime_fba.parquet)",
    )
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def _try_load_objective_value(features_path: Path, features_df: pd.DataFrame) -> pd.Series | None:
    """
    objective_value is not in Step 8 features by default. If we can find a nearby
    fva_all.parquet, extract objective_value per condition_id.
    """
    if "objective_value" in features_df.columns:
        s = pd.to_numeric(features_df["objective_value"], errors="coerce")
        s.index = features_df["condition_id"].astype(str)
        return s

    candidate = features_path.parent / "fva_all.parquet"
    if candidate.exists():
        fva = pd.read_parquet(candidate)
        if {"condition_id", "objective_value"}.issubset(fva.columns):
            s = fva.groupby("condition_id")["objective_value"].first()
            s.index = s.index.astype(str)
            return s
    return None


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    features_path = Path(args.features)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    features_df = pd.read_parquet(features_path)
    conditions_df, medium_cfg = default_conditions_and_medium(
        conditions_path=args.conditions,
        medium_path=args.medium,
    )
    objective_series = _try_load_objective_value(features_path, features_df)

    # 1) Regime table
    regime = build_regime_table(
        features_df=features_df,
        conditions_df=conditions_df,
        medium_cfg=medium_cfg,
        objective_value_series=objective_series,
    )

    if args.regime_fba:
        rf = pd.read_parquet(args.regime_fba)
        if "condition_id" not in rf.columns:
            raise ValueError("--regime-fba parquet must contain condition_id")
        # Merge & override sat columns; also keep diagnostic columns (flux/lb/ub/rid/etc.)
        regime = regime.merge(rf, on="condition_id", how="left", suffixes=("", "_fba"))

        # objective_value: prefer FBA-regime value if provided
        if "objective_value_fba" in regime.columns:
            regime["objective_value"] = regime["objective_value_fba"]
            regime = regime.drop(columns=["objective_value_fba"])

        # Saturation overrides
        if "acetate_sat_fba" in regime.columns:
            regime["acetate_sat"] = regime["acetate_sat_fba"]
            regime = regime.drop(columns=["acetate_sat_fba"])

        if "oxygen_sat" in regime.columns:
            regime["o2_sat"] = regime["oxygen_sat"]

        if "ammonium_sat" in regime.columns:
            regime["nh4_sat"] = regime["ammonium_sat"]

        if "phosphate_sat" in regime.columns:
            regime["pi_sat"] = regime["phosphate_sat"]

    regime_out = outdir / "regime_table.csv"
    regime.to_csv(regime_out, index=False)

    # 2) ElasticNet coefficients (top 30 by |coef|)
    coef_df = fit_elasticnet_coefficients(features_df=features_df, target_col="measured_OD", top_n=30)
    coef_out = outdir / "lasso_coefficients.csv"
    coef_df.to_csv(coef_out, index=False)

    # 3) Tree rules (depth 3)
    rules = fit_tree_rules(features_df=features_df, target_col="measured_OD", max_depth=3)
    rules_out = outdir / "tree_rules.txt"
    rules_out.write_text(rules + "\n", encoding="utf-8")

    # Small console report
    n = len(features_df)
    n_train = int((~features_df["measured_OD"].isna()).sum()) if "measured_OD" in features_df.columns else 0
    print(f"[OK] Wrote: {regime_out}")
    print(f"[OK] Wrote: {coef_out} (n_nonzero_shown={len(coef_df)})")
    print(f"[OK] Wrote: {rules_out}")
    print(f"[INFO] features rows={n}, trainable(non-NaN measured_OD)={n_train}")
    if objective_series is None:
        print("[INFO] objective_value not found (kept as NaN). If you have results/fva_all.parquet next to features, it will be picked up.")
    else:
        print(f"[INFO] objective_value available for {int(objective_series.notna().sum())} conditions")

    # always exit 0 for report generation unless a hard exception occurred earlier
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

