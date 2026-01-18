from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from acetate_xai.config import load_config
from acetate_xai.fva import FVAError, run_targeted_fva
from acetate_xai.io import load_sbml_model, save_table
from acetate_xai.medium import MediumConfigError, apply_condition_to_model
from acetate_xai.regime import compute_saturation_for_reaction, pick_first_existing_reaction_id

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Failure:
    condition_id: str
    error_type: str
    error_message: str


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run a large Latin Hypercube (LHS) random campaign (FBA + targeted FVA).")
    p.add_argument("--model", required=True, help="SBML model path (e.g., acetate_xai/models/model.xml)")
    p.add_argument("--targets", default="acetate_xai/configs/targets_120.json", help="targets JSON (array of reaction ids)")
    p.add_argument("--medium", default="acetate_xai/configs/medium_base.yaml", help="medium YAML config")
    p.add_argument("--regime-config", default="acetate_xai/configs/regime_exchanges.yaml", help="regime exchanges YAML config")
    p.add_argument("--outdir", default="results/campaigns/C_random_LHS", help="Campaign output directory")
    p.add_argument("--n-samples", type=int, default=2000, help="Number of LHS samples (default: 2000)")
    p.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    p.add_argument("--fraction", type=float, default=0.95, help="FVA fraction_of_optimum (default: 0.95)")
    p.add_argument("--n-jobs", type=int, default=-1, help="Parallel workers (joblib). -1 = all cores.")
    p.add_argument("--backend", default="loky", choices=["loky", "threading"], help="joblib backend.")
    p.add_argument("--keep-parts", action="store_true", help="Keep per-condition fva_parts parquet files.")
    p.add_argument("--atpm-rxn", default="ATPM", help="Reaction id for ATP maintenance (default: ATPM)")
    return p


def _load_targets_json(path: str | Path) -> list[str]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Targets JSON not found: {p}")
    with p.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, list) or not data:
        raise ValueError("Targets JSON must be a non-empty JSON array of reaction ids.")
    return [str(x) for x in data]


def _safe_filename_component(s: str) -> str:
    return s.strip().replace("/", "_").replace("\\", "_")


def _lhs_samples(n: int, *, seed: int) -> pd.DataFrame:
    """
    LatinHypercube sampling for 4 variables:
      - acetate_mM: 0..200
      - o2_lb: 0..20 (uptake max, mmol/gDW/h)
      - nh4_lb: 0..10 (uptake max)
      - atpm_fixed: 0..25
    """
    try:
        from scipy.stats import qmc
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(f"Missing dependency: scipy ({e}). Install: pip install scipy") from e

    n = int(n)
    if n <= 0:
        raise ValueError("--n-samples must be positive")

    sampler = qmc.LatinHypercube(d=4, seed=int(seed))
    u = sampler.random(n=n)
    lo = np.array([0.0, 0.0, 0.0, 0.0])
    hi = np.array([200.0, 20.0, 10.0, 25.0])
    x = qmc.scale(u, lo, hi)
    df = pd.DataFrame(x, columns=["acetate_mM", "o2_lb", "nh4_lb", "atpm_fixed"])
    return df


def _set_exchange_lb_uptake(model, rid: str, uptake_max: float) -> None:
    """
    Set exchange uptake cap with standard COBRA sign convention: uptake is negative.
      lower_bound = -uptake_max, upper_bound >= 0
    """
    rxn = model.reactions.get_by_id(rid)
    u = max(0.0, float(uptake_max))
    if rxn.upper_bound < 0.0:
        rxn.upper_bound = 0.0
    rxn.lower_bound = -u


def _fix_flux(model, rid: str, value: float) -> None:
    rxn = model.reactions.get_by_id(rid)
    v = float(value)
    rxn.upper_bound = v
    rxn.lower_bound = v


def _primary_regime_label(rec: dict[str, Any]) -> str:
    # Align with earlier label logic
    if bool(rec.get("acetate_sat", False)):
        return "Ac_limited"
    if bool(rec.get("ammonium_sat", False)):
        return "N_limited"
    if bool(rec.get("phosphate_sat", False)):
        return "Pi_limited"
    if bool(rec.get("oxygen_sat", False)):
        return "O2_limited"
    return "Unconstrained"


def _run_one(
    *,
    model_path: str,
    medium_cfg: dict[str, Any],
    regime_cfg: dict[str, Any],
    targets: list[str],
    atpm_rxn: str,
    parts_dir: Path,
    cond: dict[str, Any],
    fraction: float,
) -> tuple[dict[str, Any] | None, Failure | None]:
    """
    One worker:
      - load model
      - apply base medium + acetate scaling via apply_condition_to_model
      - override oxygen/ammonium exchange LBs (uptake caps)
      - fix ATPM flux
      - run FBA + targeted FVA
      - write fva part parquet
      - compute saturation table (acetate/oxygen/ammonium/phosphate) from FBA solution
    """
    cid = str(cond["condition_id"])
    out_path = parts_dir / f"condition_id={_safe_filename_component(cid)}.parquet"

    try:
        model = load_sbml_model(model_path)
        apply_condition_to_model(model, cond, medium_cfg)

        # pick exchange rids from regime config (fallback to medium config where available)
        cand_o2 = regime_cfg.get("oxygen", [])
        cand_nh4 = regime_cfg.get("ammonium", [])
        rid_o2 = pick_first_existing_reaction_id(model, [str(x) for x in cand_o2]) if isinstance(cand_o2, list) else None
        rid_nh4 = pick_first_existing_reaction_id(model, [str(x) for x in cand_nh4]) if isinstance(cand_nh4, list) else None

        if rid_o2:
            _set_exchange_lb_uptake(model, rid_o2, uptake_max=float(cond["o2_lb"]))
        if rid_nh4:
            _set_exchange_lb_uptake(model, rid_nh4, uptake_max=float(cond["nh4_lb"]))

        # Fix ATPM
        try:
            _fix_flux(model, str(atpm_rxn), float(cond["atpm_fixed"]))
        except KeyError:
            logger.warning("ATPM rxn not found in model (skipped): %s", atpm_rxn)

        sol = model.optimize()
        if sol.status != "optimal":
            raise FVAError(f"FBA failed: status={sol.status}")

        objective_value = float(sol.objective_value)

        fva_df = run_targeted_fva(model, targets=targets, fraction_of_optimum=float(fraction))
        out = fva_df.copy()
        out.insert(0, "condition_id", cid)
        out.insert(1, "objective_value", objective_value)
        out = out[["condition_id", "objective_value", "reaction_id", "fva_min", "fva_max"]]
        save_table(out, out_path, fmt="parquet")

        # saturation rows (wide)
        nutrients = {
            "acetate": regime_cfg.get("acetate", []),
            "oxygen": regime_cfg.get("oxygen", []),
            "ammonium": regime_cfg.get("ammonium", []),
            "phosphate": regime_cfg.get("phosphate", []),
        }
        rec: dict[str, Any] = {
            "condition_id": cid,
            "objective_value": objective_value,
            "acetate_mM": float(cond["acetate_mM"]),
            "o2_lb": float(cond["o2_lb"]),
            "nh4_lb": float(cond["nh4_lb"]),
            "atpm_fixed": float(cond["atpm_fixed"]),
        }
        for nutrient, cand in nutrients.items():
            cand_list = [str(x) for x in cand] if isinstance(cand, list) else []
            rid = pick_first_existing_reaction_id(model, cand_list) if cand_list else None
            if rid is None:
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
            sat = compute_saturation_for_reaction(rid=rid, model=model, solution=sol)
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

        rec["primary_regime"] = _primary_regime_label(rec)
        return rec, None
    except (FileNotFoundError, MediumConfigError, FVAError, ValueError) as e:
        return None, Failure(condition_id=cid, error_type=type(e).__name__, error_message=str(e))
    except Exception as e:  # noqa: BLE001
        return None, Failure(condition_id=cid, error_type=type(e).__name__, error_message=str(e))


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    model_path = str(Path(args.model))
    if not Path(model_path).exists():
        print(f"[ERROR] Model file not found: {model_path}", file=sys.stderr)
        return 2

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    parts_dir = outdir / "fva_parts"
    parts_dir.mkdir(parents=True, exist_ok=True)

    try:
        medium_cfg = load_config(args.medium)
        regime_cfg = load_config(args.regime_config)
        targets = _load_targets_json(args.targets)
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Failed to load inputs: {e}", file=sys.stderr)
        return 2

    # LHS design table (one row per condition)
    design = _lhs_samples(int(args.n_samples), seed=int(args.seed))
    design.insert(0, "condition_id", [f"LHS_{i:05d}" for i in range(len(design))])
    # Fixed metadata expected by apply_condition_to_model
    design["pH0"] = 7.0
    design["yeast_extract_gL"] = 0.0
    design["nh4cl_gL"] = 0.0  # we override ammonium by nh4_lb directly
    design["set_name"] = "random_lhs"
    design.to_parquet(outdir / "lhs_conditions.parquet", index=False)

    # Avoid oversubscription on servers
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    rows = design.to_dict(orient="records")
    logger.info(
        "Running LHS campaign: n=%d, n_targets=%d, fraction=%.3f, n_jobs=%d",
        len(rows),
        len(targets),
        float(args.fraction),
        int(args.n_jobs),
    )

    results = Parallel(n_jobs=int(args.n_jobs), backend=args.backend)(
        delayed(_run_one)(
            model_path=model_path,
            medium_cfg=medium_cfg,
            regime_cfg=regime_cfg,
            targets=targets,
            atpm_rxn=str(args.atpm_rxn),
            parts_dir=parts_dir,
            cond=row,
            fraction=float(args.fraction),
        )
        for row in rows
    )

    recs: list[dict[str, Any]] = []
    failures: list[Failure] = []
    for rec, fail in results:
        if rec is not None:
            recs.append(rec)
        if fail is not None:
            failures.append(fail)

    # Save failures
    pd.DataFrame([asdict(f) for f in failures]).to_csv(outdir / "failed_samples.csv", index=False)

    # Aggregate FVA parts -> wide features
    from acetate_xai.collect import build_fva_long_features, build_wide_feature_matrix, load_fva_parts

    fva_all = load_fva_parts(parts_dir)
    long_feat = build_fva_long_features(fva_all)
    wide = build_wide_feature_matrix(long_feat)

    # Join design metadata (no measured_OD)
    meta_cols = ["condition_id", "acetate_mM", "o2_lb", "nh4_lb", "atpm_fixed", "pH0", "yeast_extract_gL", "set_name"]
    meta = design[meta_cols].copy()
    # Keep ALL sampled conditions (even if FVA failed for some); failed ones will have NaNs in feature columns.
    if wide["condition_id"].duplicated().any():
        raise ValueError("Duplicate condition_id rows in wide feature matrix.")
    features = meta.merge(wide, on="condition_id", how="left")
    save_table(features, outdir / "features.parquet", fmt="parquet")

    # Regime labels
    regime_labels = pd.DataFrame(recs)
    if not regime_labels.empty:
        # ensure stable order
        regime_labels = regime_labels.merge(meta[["condition_id"]], on="condition_id", how="right")
    save_table(regime_labels, outdir / "regime_labels.parquet", fmt="parquet")

    # Optionally cleanup parts
    if not args.keep_parts:
        try:
            for f in parts_dir.glob("*.parquet"):
                f.unlink()
            parts_dir.rmdir()
        except Exception:
            pass

    print(f"[OK] features: {outdir / 'features.parquet'}")
    print(f"[OK] regime_labels: {outdir / 'regime_labels.parquet'}")
    print(f"[OK] failed_samples: {outdir / 'failed_samples.csv'} (n_failed={len(failures)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

