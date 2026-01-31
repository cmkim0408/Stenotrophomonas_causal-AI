#!/usr/bin/env python3
"""
Sensitivity analysis script.

Runs the existing pipeline with:
  - near-optimal threshold (fraction_of_optimum) = 0.95
  - co-limitation epsilon = 0.5x and 2x

Saves outputs under results/sensitivity/. NO retraining. NO model modification.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

EPS_DEFAULT = 1e-6
FRAC_DEFAULT = 0.95


def _load_targets(path: str | Path) -> list[str]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Targets not found: {p}")
    with p.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, list) or not data:
        raise ValueError("Targets must be a non-empty JSON array")
    return [str(x) for x in data]


def _load_atpm_fit(path: str | Path | None) -> tuple[float, float]:
    if not path or not Path(path).exists():
        return 10.0, 0.1
    with Path(path).open("r", encoding="utf-8") as f:
        fit = json.load(f)
    return float(fit.get("a", 10.0)), float(fit.get("b", 0.1))


def _map_holdout_to_condition(row: pd.Series) -> dict:
    o2lb = row.get("o2lb")
    if o2lb is None or (isinstance(o2lb, float) and np.isnan(o2lb)):
        o2lb = -100.0
    else:
        o2lb = float(o2lb)
    return {
        "condition_id": str(row.get("condition_id", "")),
        "acetate_mM": float(row.get("acetate_mM", row.get("acetate", 0.0)) or 0.0),
        "nh4cl_gL": float(row.get("NH4_gL", row.get("nh4cl_gL", row.get("NH4Cl", 1.0))) or 1.0),
        "yeast_extract_gL": float(row.get("YE_gL", row.get("yeast_extract_gL", row.get("YE", 0.0))) or 0.0),
        "pH0": float(row.get("pH", row.get("pH0", 7.0)) or 7.0),
        "o2lb": o2lb,
    }


def _apply_o2lb(model, condition: dict, medium_cfg: dict) -> None:
    o2lb = condition.get("o2lb")
    if o2lb is None or not np.isfinite(o2lb):
        return
    exchanges = medium_cfg.get("exchanges", {})
    o2_rid = exchanges.get("oxygen") or "EX_o2_e"
    try:
        rxn = model.reactions.get_by_id(o2_rid)
    except KeyError:
        return
    rxn.upper_bound = 1000.0
    rxn.lower_bound = float(o2lb)




def _ensure_xgbclassifier_metadata(model, label_mapping_path="results/xai_xgb/label_mapping.csv"):
    """
    Restore sklearn metadata (classes_/n_classes_) that may be missing after
    XGBoost JSON model loading (xgboost>=2.0).

    Note: In some xgboost versions, `classes_` is a read-only property.
    We therefore restore via `model._le.classes_` (LabelEncoder).
    """
    import os
    import numpy as np

    # If already present and consistent, nothing to do
    if hasattr(model, "n_classes_") and getattr(model, "n_classes_", None):
        # also ensure _le exists
        if hasattr(model, "_le") and getattr(model._le, "classes_", None) is not None:
            return

    # Build class list from label_mapping.csv if possible
    cls = None
    try:
        import pandas as pd
        if os.path.exists(label_mapping_path):
            lm = pd.read_csv(label_mapping_path)
            # prefer integer-like column as class id
            cand = [c for c in lm.columns if getattr(lm[c].dtype, "kind", "") in ("i", "u")]
            if cand:
                cls = sorted(lm[cand[0]].dropna().astype(int).unique().tolist())
            else:
                cls = list(range(len(lm)))
    except Exception:
        cls = None

    if not cls:
        cls = [0, 1]  # fallback binary

    # Restore via LabelEncoder
    try:
        from sklearn.preprocessing import LabelEncoder
        le = LabelEncoder()
        le.classes_ = np.array(cls, dtype=int)
        model._le = le
    except Exception:
        # last resort: attach minimal _le-like object
        class _LE: pass
        le = _LE()
        le.classes_ = np.array(cls, dtype=int)
        model._le = le

    # Restore n_classes_
    try:
        model.n_classes_ = len(cls)
    except Exception:
        pass


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run sensitivity: frac=0.95, eps=0.5x and 2x."
    )
    ap.add_argument(
        "--conditions",
        default="data/holdout_10_conditions.csv",
        help="Conditions CSV",
    )
    ap.add_argument(
        "--model",
        default="acetate_xai/models/model.xml",
        help="SBML model path",
    )
    ap.add_argument(
        "--targets",
        default="acetate_xai/configs/targets_120.json",
        help="Targets JSON for FVA",
    )
    ap.add_argument(
        "--medium",
        default="acetate_xai/configs/medium_base.yaml",
        help="Medium YAML config",
    )
    ap.add_argument(
        "--regime-config",
        default="acetate_xai/configs/regime_exchanges.yaml",
        help="Regime exchanges YAML",
    )
    ap.add_argument(
        "--atpm-fit",
        default=None,
        help="atpm_fit.json path",
    )
    ap.add_argument(
        "--outdir",
        default="results/sensitivity",
        help="Output base directory",
    )
    ap.add_argument(
        "--fraction",
        type=float,
        default=FRAC_DEFAULT,
        help="FVA fraction_of_optimum (fixed at 0.95)",
    )
    args = ap.parse_args()

    atpm_fit = args.atpm_fit
    if not atpm_fit:
        for c in [
            "results/campaigns/C6_atpm_calibrated/calibration/atpm_fit.json",
            "acetate_xai/configs/atpm_fit.json",
        ]:
            if Path(c).exists():
                atpm_fit = c
                break
    a, b = _load_atpm_fit(atpm_fit)

    from acetate_xai.config import load_config
    from acetate_xai.fva import FVAError, run_targeted_fva
    from acetate_xai.io import load_sbml_model
    from acetate_xai.medium import apply_condition_to_model
    from acetate_xai.regime import (
        compute_saturation_for_reaction,
        pick_first_existing_reaction_id,
    )
    from acetate_xai.collect import build_fva_long_features, build_wide_feature_matrix

    conditions_path = Path(args.conditions)
    if not conditions_path.exists():
        logger.error("Conditions not found: %s", conditions_path)
        return 2

    df_cond = pd.read_csv(conditions_path)
    cond_list = [_map_holdout_to_condition(r) for _, r in df_cond.iterrows()]

    medium_cfg = load_config(args.medium)
    regime_cfg = load_config(args.regime_config)
    targets = _load_targets(args.targets)

    nutrients = {
        "acetate": regime_cfg.get("acetate") or ["EX_ac_e"],
        "oxygen": regime_cfg.get("oxygen") or ["EX_o2_e"],
        "ammonium": regime_cfg.get("ammonium") or ["EX_nh4_e"],
        "phosphate": regime_cfg.get("phosphate") or ["EX_pi_e"],
    }

    eps_multipliers = [0.5, 2.0]
    frac = float(args.fraction)

    for eps_mult in eps_multipliers:
        eps = EPS_DEFAULT * eps_mult
        subdir = Path(args.outdir) / f"frac_{frac:.2f}_eps_{eps_mult}x"
        subdir.mkdir(parents=True, exist_ok=True)

        regime_rows = []
        fva_parts = []

        for cond in cond_list:
            cid = cond["condition_id"]
            try:
                model = load_sbml_model(args.model)
                apply_condition_to_model(model, cond, medium_cfg)
                _apply_o2lb(model, cond, medium_cfg)

                atpm_eff = max(0.0, min(200.0, a + b * cond["acetate_mM"]))
                try:
                    rxn = model.reactions.get_by_id("ATPM")
                    rxn.upper_bound = atpm_eff
                    rxn.lower_bound = atpm_eff
                except KeyError:
                    pass

                sol = model.optimize()
                if sol.status != "optimal":
                    regime_rows.append({"condition_id": cid, "objective_value": np.nan})
                    continue

                rec = {
                    "condition_id": cid,
                    "objective_value": float(sol.objective_value),
                    "eps_used": eps,
                }

                for nutrient, cand in nutrients.items():
                    rid = pick_first_existing_reaction_id(
                        model, cand if isinstance(cand, list) else [cand]
                    )
                    if rid is None:
                        rec.update({
                            f"{nutrient}_rid": "",
                            f"{nutrient}_sat": False,
                        })
                        continue
                    sat = compute_saturation_for_reaction(
                        rid=rid, model=model, solution=sol, eps=eps
                    )
                    rec[f"{nutrient}_rid"] = rid
                    rec[f"{nutrient}_sat"] = sat.saturated

                regime_rows.append(rec)

                fva_df = run_targeted_fva(
                    model, targets=targets, fraction_of_optimum=frac
                )
                fva_df.insert(0, "condition_id", cid)
                fva_df.insert(1, "objective_value", rec["objective_value"])
                fva_df = fva_df[["condition_id", "objective_value", "reaction_id", "fva_min", "fva_max"]]
                fva_parts.append(fva_df)

            except (FVAError, Exception) as e:
                logger.warning("Condition %s failed: %s", cid, e)
                regime_rows.append({"condition_id": cid, "objective_value": np.nan})

        regime_df = pd.DataFrame(regime_rows)
        regime_df.to_csv(subdir / "regime_table_sensitivity.csv", index=False)

        if fva_parts:
            fva_all = pd.concat(fva_parts, ignore_index=True)
            long_feat = build_fva_long_features(fva_all)
            wide_feat = build_wide_feature_matrix(long_feat)
            wide_feat.to_csv(subdir / "features_sensitivity.csv", index=False)

        meta = {
            "fraction_of_optimum": frac,
            "eps_default": EPS_DEFAULT,
            "eps_multiplier": eps_mult,
            "eps_used": eps,
            "n_conditions": len(cond_list),
        }
        (subdir / "sensitivity_metadata.json").write_text(
            json.dumps(meta, indent=2), encoding="utf-8"
        )
        logger.info("Wrote %s", subdir)

    # Summary CSV: compare baseline vs sensitivity (delta in predicted_severity)
    _write_sensitivity_summary(
        outdir=Path(args.outdir),
        holdout_pred_path=Path("results/holdout_predictions.csv"),
        severity_model_path=Path("results/xai_xgb_maintenance/model.json"),
    )

    return 0


def _write_sensitivity_summary(
    outdir: Path,
    holdout_pred_path: Path,
    severity_model_path: Path,
) -> None:
    """Write sensitivity_summary.csv comparing baseline vs eps variants."""
    if not holdout_pred_path.exists():
        logger.warning("Holdout predictions not found, skipping sensitivity summary")
        return
    if not severity_model_path.exists():
        logger.warning("Severity model not found, skipping sensitivity summary")
        return

    base_pred = pd.read_csv(holdout_pred_path)
    if "predicted_severity" not in base_pred.columns:
        return

    base_sev = base_pred.set_index("condition_id")["predicted_severity"]

    rows = []
    for cid in base_sev.index:
        row = {"condition_id": cid, "baseline_severity": base_sev.get(cid, np.nan)}
        for eps_name in ["frac_0.95_eps_0.5x", "frac_0.95_eps_2.0x"]:
            feat_path = outdir / eps_name / "features_sensitivity.csv"
            if feat_path.exists():
                try:
                    from xgboost import XGBRegressor

                    df = pd.read_csv(feat_path)
                    sub = df.loc[df["condition_id"] == cid]
                    if not sub.empty:
                        model = XGBRegressor()
                        model.load_model(str(severity_model_path))
                        feat_names = list(model.get_booster().feature_names)
                        X = pd.DataFrame(index=[0])
                        r = sub.iloc[0]
                        for col in feat_names:
                            X[col] = r[col] if col in r.index else 0.0
                        sev = float(np.clip(model.predict(X[feat_names])[0], 0, 1))
                        row[f"{eps_name}_severity"] = sev
                        row[f"delta_{eps_name}"] = sev - row["baseline_severity"]
                except Exception as e:
                    logger.warning("Could not predict for %s %s: %s", cid, eps_name, e)
        rows.append(row)

    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "sensitivity_summary.csv", index=False)
    logger.info("Wrote %s", outdir / "sensitivity_summary.csv")


if __name__ == "__main__":
    sys.exit(main())
