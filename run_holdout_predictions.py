#!/usr/bin/env python3
"""
Holdout prediction script.

Loads the existing trained digital-twin pipeline (GEM + pre-trained ML models),
applies each of the 10 holdout conditions, runs prediction ONLY (no retraining),
and outputs results/holdout_predictions.csv.

NO model modification. NO retraining.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def _get_git_commit() -> str:
    try:
        r = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=5, check=False,
        )
        return (r.stdout or "").strip() or "unknown"
    except Exception:
        return "unknown"


def _timestamp_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


def _load_targets(path: str | Path) -> list[str]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Targets not found: {p}")
    with p.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, list) or not data:
        raise ValueError("Targets must be a non-empty JSON array")
    return [str(x) for x in data]


def _load_atpm_fit(path: str | Path) -> tuple[float, float]:
    p = Path(path)
    if not p.exists():
        logger.warning("atpm_fit.json not found, using fallback a=10, b=0.1")
        return 10.0, 0.1
    with p.open("r", encoding="utf-8") as f:
        fit = json.load(f)
    a = float(fit.get("a", 10.0))
    b = float(fit.get("b", 0.1))
    return a, b


def _map_holdout_to_condition(row: pd.Series) -> dict:
    """Map holdout CSV columns to pipeline condition dict.
    CSV columns: condition_id, acetate_mM, NH4_gL, YE_gL, pH, o2_level, o2lb
    - o2lb: direct model input (EX_o2_e.lb); -100=high, -50=mid, -10=low
    - NH4_gL -> nh4cl_gL, YE_gL -> yeast_extract_gL for medium
    """
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
        "o2_level": str(row.get("o2_level", "")),
        "o2lb": o2lb,
    }


def _apply_o2lb(model, condition: dict, medium_cfg: dict) -> None:
    """Apply o2lb (EX_o2_e lower_bound) from condition. o2lb is direct model input."""
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


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run holdout predictions (no retraining)."
    )
    ap.add_argument(
        "--conditions",
        default="data/holdout_10_conditions.csv",
        help="Holdout conditions CSV",
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
        "--atpm-fit",
        default=None,
        help="atpm_fit.json path (optional; uses fallback if missing)",
    )
    ap.add_argument(
        "--regime-model",
        default="results/xai_xgb/model.json",
        help="Pre-trained regime classifier path",
    )
    ap.add_argument(
        "--severity-model",
        default="results/xai_xgb_maintenance/model.json",
        help="Pre-trained severity regressor path",
    )
    ap.add_argument(
        "--regime-labels",
        default="results/xai_xgb/label_mapping.csv",
        help="Label mapping for regime (int -> label)",
    )
    ap.add_argument(
        "--out",
        default="results/holdout_predictions.csv",
        help="Output CSV path",
    )
    ap.add_argument(
        "--fraction",
        type=float,
        default=0.95,
        help="FVA fraction_of_optimum",
    )
    args = ap.parse_args()

    # ----- PREREQUISITES CHECK (fail immediately) -----
    conditions_path = Path(args.conditions)
    regime_model_path = Path(args.regime_model)
    severity_model_path = Path(args.severity_model)
    regime_labels_path = Path(args.regime_labels)

    missing = []
    if not conditions_path.exists():
        missing.append(str(conditions_path))
    if not regime_model_path.exists():
        missing.append(str(regime_model_path))
    if not severity_model_path.exists():
        missing.append(str(severity_model_path))
    if not regime_labels_path.exists():
        missing.append(str(regime_labels_path))

    if missing:
        logger.error("PREREQUISITES MISSING - cannot proceed. Train models first.")
        for m in missing:
            logger.error("  - %s", m)
        return 2

    # Resolve atpm_fit (try common locations if not given)
    atpm_fit = args.atpm_fit
    if not atpm_fit:
        candidates = [
            "results/campaigns/C6_atpm_calibrated/calibration/atpm_fit.json",
            "acetate_xai/configs/atpm_fit.json",
        ]
        for c in candidates:
            if Path(c).exists():
                atpm_fit = c
                break
    a, b = _load_atpm_fit(atpm_fit) if atpm_fit else (10.0, 0.1)

    out_path = Path(args.out)

    # Load pipeline components
    from acetate_xai.config import load_config
    from acetate_xai.fva import FVAError, run_targeted_fva
    from acetate_xai.io import load_sbml_model
    from acetate_xai.medium import apply_condition_to_model

    from acetate_xai.collect import (
        build_fva_long_features,
        build_wide_feature_matrix,
    )

    targets = _load_targets(args.targets)
    medium_cfg = load_config(args.medium)
    df_cond = pd.read_csv(conditions_path)

    if len(df_cond) != 10:
        logger.warning("Holdout CSV has %d rows (expected 10)", len(df_cond))

    required_cols = {"condition_id", "acetate_mM", "o2lb"}
    missing_cols = required_cols - set(df_cond.columns)
    if missing_cols:
        logger.error("Holdout CSV missing required columns: %s", sorted(missing_cols))
        return 2

    # Map CSV columns to pipeline condition dict
    rows = []
    for _, r in df_cond.iterrows():
        rows.append(_map_holdout_to_condition(r))
    cond_list = rows

    # Regime model: load and get feature names
    from xgboost import XGBClassifier, XGBRegressor

    regime_model = XGBClassifier()
    regime_model.load_model(str(regime_model_path))
    regime_feat_names = list(regime_model.get_booster().feature_names)

    severity_model = XGBRegressor()
    severity_model.load_model(str(severity_model_path))
    severity_feat_names = list(severity_model.get_booster().feature_names)

    # Label mapping for regime (already checked in prerequisites)
    label_map: dict[int, str] = {}
    if regime_labels_path.exists():
        lm = pd.read_csv(regime_labels_path)
        for _, r in lm.iterrows():
            label_map[int(r["label_int"])] = str(r["label"])

    # FVA parts collector
    fva_parts: list[pd.DataFrame] = []
    obj_vals: dict[str, float] = {}

    exchanges = medium_cfg.get("exchanges", {})
    ac_rid = exchanges.get("acetate") or "EX_ac_e"
    nh4_rid = exchanges.get("ammonium") or "EX_nh4_e"
    o2_rid = exchanges.get("oxygen") or "EX_o2_e"
    k_ac = float(medium_cfg.get("scaling", {}).get("k_ac", 0.1))
    k_nh4 = float(medium_cfg.get("scaling", {}).get("k_nh4", 10.0))

    for cond in cond_list:
        cid = cond["condition_id"]
        ac = cond["acetate_mM"]
        nh4 = cond["nh4cl_gL"]
        ye = cond["yeast_extract_gL"]
        ph = cond["pH0"]
        o2lb = cond["o2lb"]
        # Mapping log (pre-run)
        ex_ac = -k_ac * ac if ac else 0
        ex_nh4 = -k_nh4 * nh4 if nh4 else 0
        logger.info(
            "[%s] acetate=%.1fmM NH4=%.2fg/L YE=%.2fg/L pH=%.1f o2lb=%.0f -> constraints: %s=%.1f %s=%.1f %s=%.0f",
            cid, ac, nh4, ye, ph, o2lb, ac_rid, ex_ac, nh4_rid, ex_nh4, o2_rid, o2lb,
        )
        try:
            model = load_sbml_model(args.model)
            apply_condition_to_model(model, cond, medium_cfg)
            _apply_o2lb(model, cond, medium_cfg)

            acetate_mM = cond["acetate_mM"]
            atpm_eff = max(0.0, min(200.0, a + b * acetate_mM))

            atpm_rid = "ATPM"
            try:
                rxn = model.reactions.get_by_id(atpm_rid)
                rxn.upper_bound = atpm_eff
                rxn.lower_bound = atpm_eff
            except KeyError:
                pass

            sol = model.optimize()
            if sol.status != "optimal":
                logger.warning("Condition %s: FBA status=%s, skipping", cid, sol.status)
                obj_vals[cid] = float("nan")
                continue

            obj_vals[cid] = float(sol.objective_value)
            fva_df = run_targeted_fva(
                model,
                targets=targets,
                fraction_of_optimum=float(args.fraction),
            )
            fva_df.insert(0, "condition_id", cid)
            fva_df.insert(1, "objective_value", obj_vals[cid])
            fva_df = fva_df[["condition_id", "objective_value", "reaction_id", "fva_min", "fva_max"]]
            fva_parts.append(fva_df)

        except (FVAError, Exception) as e:
            logger.warning("Condition %s failed: %s", cid, e)
            obj_vals[cid] = float("nan")

    if not fva_parts:
        logger.error("No successful FVA runs; cannot predict")
        return 2

    fva_all = pd.concat(fva_parts, ignore_index=True)
    long_feat = build_fva_long_features(fva_all)
    wide_feat = build_wide_feature_matrix(long_feat)

    # Predict per condition
    out_rows = []
    run_max = max(obj_vals.values()) if obj_vals else 1.0
    if not np.isfinite(run_max) or run_max <= 0:
        run_max = 1.0

    for cid in df_cond["condition_id"].astype(str):
        sub = wide_feat.loc[wide_feat["condition_id"] == cid]
        if sub.empty:
            out_rows.append({
                "condition_id": cid,
                "predicted_regime": "unknown",
                "predicted_severity": float("nan"),
            })
            continue

        row = sub.iloc[0]
        # Build X for regime
        X_reg = pd.DataFrame(index=[0])
        for col in regime_feat_names:
            X_reg[col] = row[col] if col in row.index else 0.0

        # Build X for severity (may differ)
        X_sev = pd.DataFrame(index=[0])
        for col in severity_feat_names:
            X_sev[col] = row[col] if col in row.index else 0.0

        pred_reg_int = regime_model.predict(X_reg[regime_feat_names])[0]
        pred_reg = label_map.get(int(pred_reg_int), f"class_{pred_reg_int}")

        pred_sev_raw = float(severity_model.predict(X_sev[severity_feat_names])[0])
        pred_sev = float(np.clip(pred_sev_raw, 0.0, 1.0))

        out_rows.append({
            "condition_id": cid,
            "predicted_regime": pred_reg,
            "predicted_severity": pred_sev,
        })

    git_commit = _get_git_commit()
    timestamp_utc = _timestamp_utc()
    model_path_regime = str(regime_model_path.resolve())
    model_path_severity = str(severity_model_path.resolve())

    out_df = pd.DataFrame(out_rows)
    out_df["git_commit"] = git_commit
    out_df["timestamp_utc"] = timestamp_utc
    out_df["model_path_regime"] = model_path_regime
    out_df["model_path_severity"] = model_path_severity
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)
    logger.info("Wrote %s", out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
