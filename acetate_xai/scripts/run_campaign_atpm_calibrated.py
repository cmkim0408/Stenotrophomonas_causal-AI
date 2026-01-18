from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from acetate_xai.collect import collect_and_build_features
from acetate_xai.config import load_config
from acetate_xai.fva import FVAError, run_targeted_fva
from acetate_xai.io import load_conditions_csv, load_sbml_model, save_table
from acetate_xai.medium import MediumConfigError, apply_condition_to_model
from acetate_xai.regime import compute_saturation_for_reaction, pick_first_existing_reaction_id


@dataclass(frozen=True)
class Failure:
    condition_id: str
    error_type: str
    error_message: str


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Run full calibrated campaign: per-condition ATPM_eff from f(Ac) + FBA + targeted FVA + features + regime_fba + xai."
    )
    p.add_argument("--model", required=True, help="SBML model path (e.g., acetate_xai/models/model.xml)")
    p.add_argument("--conditions", required=True, help="Conditions CSV (long format)")
    p.add_argument("--targets", required=True, help="targets JSON (array of reaction ids)")
    p.add_argument("--medium", required=True, help="Medium YAML config")
    p.add_argument("--regime-config", default="acetate_xai/configs/regime_exchanges.yaml", help="Regime exchanges YAML config")
    p.add_argument("--atpm-fit", required=True, help="Path to atpm_fit.json from calibration step")
    p.add_argument("--outdir", required=True, help="Output run directory (e.g., results/campaigns/C6_atpm_calibrated/run__atpmCalib_linear/)")
    p.add_argument("--n-jobs", type=int, default=1, help="Parallel workers (joblib). Use 1 to disable.")
    p.add_argument("--backend", default="loky", choices=["loky", "threading"], help="joblib backend.")
    p.add_argument("--fraction", type=float, default=0.95, help="FVA fraction_of_optimum (default: 0.95)")
    p.add_argument("--atpm-rid", default="ATPM", help="Reaction id for ATP maintenance (default: ATPM)")
    p.add_argument("--clip-min", type=float, default=0.0, help="ATPM_eff clip min (default: 0)")
    p.add_argument("--clip-max", type=float, default=200.0, help="ATPM_eff clip max (default: 200)")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def _load_targets_json(path: str | Path) -> list[str]:
    import json as _json

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Targets JSON not found: {p}")
    with p.open("r", encoding="utf-8") as f:
        data = _json.load(f)
    if not isinstance(data, list) or not data:
        raise ValueError("Targets JSON must be a non-empty JSON array of reaction ids.")
    return [str(x) for x in data]


def _safe_filename_component(s: str) -> str:
    return s.strip().replace("/", "_").replace("\\", "_")


def _apply_fixed_flux(model, rid: str, value: float) -> None:
    present = set(r.id for r in model.reactions)
    if rid not in present:
        raise ValueError(f"Reaction not found for fixed flux: {rid}")
    rxn = model.reactions.get_by_id(rid)
    rxn.upper_bound = 1000.0
    rxn.lower_bound = float(value)
    rxn.upper_bound = float(value)


def _compute_atpm_eff(*, a: float, b: float, acetate_mM: float, clip_min: float, clip_max: float) -> float:
    v = float(a + b * float(acetate_mM))
    if not np.isfinite(v):
        v = float("nan")
    v = max(float(clip_min), v)
    v = min(float(clip_max), v)
    return float(v)


def _candidates(regime_cfg: dict, key: str) -> list[str]:
    c = regime_cfg.get(key, [])
    if isinstance(c, list):
        return [str(x) for x in c]
    return []


def _run_one_condition_fva(
    *,
    model_path: str,
    medium_cfg: dict,
    targets: list[str],
    atpm_rid: str,
    a: float,
    b: float,
    clip_min: float,
    clip_max: float,
    row: dict,
    out_parts_dir: Path,
    fraction: float,
) -> Failure | None:
    cid = str(row.get("condition_id", ""))
    if not cid:
        return Failure(condition_id="", error_type="ValueError", error_message="Missing condition_id in row")

    out_parts_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_parts_dir / f"condition_id={_safe_filename_component(cid)}.parquet"

    try:
        model = load_sbml_model(model_path)
        apply_condition_to_model(model, row, medium_cfg)

        acetate_mM = float(row.get("acetate_mM", 0.0) or 0.0)
        atpm_eff = _compute_atpm_eff(a=a, b=b, acetate_mM=acetate_mM, clip_min=clip_min, clip_max=clip_max)
        _apply_fixed_flux(model, atpm_rid, atpm_eff)
        logging.info("ATPM_eff(%s): acetate_mM=%.6g => %.6g", cid, acetate_mM, atpm_eff)

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
        return None
    except (FileNotFoundError, MediumConfigError, FVAError, ValueError) as e:
        return Failure(condition_id=cid, error_type=type(e).__name__, error_message=str(e))
    except Exception as e:  # noqa: BLE001
        return Failure(condition_id=cid, error_type=type(e).__name__, error_message=str(e))


def _run_regime_table_per_condition(
    *,
    model_path: str,
    conditions_df: pd.DataFrame,
    medium_cfg: dict,
    regime_cfg: dict,
    atpm_rid: str,
    a: float,
    b: float,
    clip_min: float,
    clip_max: float,
    eps: float = 1e-6,
    infty_bound: float = 999.0,
) -> pd.DataFrame:
    nutrients = {
        "acetate": _candidates(regime_cfg, "acetate"),
        "oxygen": _candidates(regime_cfg, "oxygen"),
        "ammonium": _candidates(regime_cfg, "ammonium"),
        "phosphate": _candidates(regime_cfg, "phosphate"),
    }

    rows_out: list[dict] = []
    for _, row in conditions_df.iterrows():
        cid = str(row["condition_id"])
        acetate_mM = float(row.get("acetate_mM", 0.0) or 0.0)
        atpm_eff = _compute_atpm_eff(a=a, b=b, acetate_mM=acetate_mM, clip_min=clip_min, clip_max=clip_max)

        model = load_sbml_model(model_path)
        apply_condition_to_model(model, row.to_dict(), medium_cfg)
        _apply_fixed_flux(model, atpm_rid, atpm_eff)

        sol = model.optimize()
        if sol.status != "optimal":
            rec: dict = {"condition_id": cid, "objective_value": np.nan, "atpm_eff": atpm_eff}
            for nutrient in nutrients.keys():
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
            rows_out.append(rec)
            continue

        rec = {"condition_id": cid, "objective_value": float(sol.objective_value), "atpm_eff": atpm_eff}
        for nutrient, cand in nutrients.items():
            rid_used = pick_first_existing_reaction_id(model, cand) if cand else None
            if rid_used is None:
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

            sat = compute_saturation_for_reaction(
                rid=rid_used,
                model=model,
                solution=sol,
                eps=float(eps),
                infty_bound=float(infty_bound),
            )
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

        rows_out.append(rec)

    out_df = pd.DataFrame(rows_out)
    if out_df["condition_id"].duplicated().any():
        raise RuntimeError("Duplicate condition_id in regime table output")
    return out_df


def _run_xai_report(*, python_exe: str, features_path: Path, outdir: Path, conditions: str, medium: str, regime_fba: Path) -> None:
    import subprocess

    script = Path("acetate_xai") / "scripts" / "xai_report.py"
    cmd = [
        python_exe,
        str(script),
        "--features",
        str(features_path),
        "--outdir",
        str(outdir),
        "--conditions",
        str(conditions),
        "--medium",
        str(medium),
        "--regime-fba",
        str(regime_fba),
    ]
    subprocess.run(cmd, check=True)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    # avoid oversubscription
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    for p in [args.model, args.conditions, args.targets, args.medium, args.atpm_fit, args.regime_config]:
        if not Path(p).exists():
            print(f"[ERROR] input not found: {p}", file=sys.stderr)
            return 2

    try:
        conditions_df = load_conditions_csv(args.conditions)
        medium_cfg = load_config(args.medium)
        regime_cfg = load_config(args.regime_config)
        targets = _load_targets_json(args.targets)
        fit = json.loads(Path(args.atpm_fit).read_text(encoding="utf-8"))
        a = float(fit["a"])
        b = float(fit["b"])
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Failed to load inputs: {e}", file=sys.stderr)
        return 2

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    parts_dir = outdir / "fva_parts"
    failures_csv = outdir / "failed_conditions.csv"

    rows = conditions_df.to_dict(orient="records")
    logging.info(
        "Running calibrated campaign: n_conditions=%d, n_targets=%d, n_jobs=%d, fraction=%.3f",
        len(rows),
        len(targets),
        args.n_jobs,
        float(args.fraction),
    )

    failures: list[Failure] = []
    results = Parallel(n_jobs=int(args.n_jobs), backend=args.backend)(
        delayed(_run_one_condition_fva)(
            model_path=str(args.model),
            medium_cfg=medium_cfg,
            targets=targets,
            atpm_rid=str(args.atpm_rid),
            a=a,
            b=b,
            clip_min=float(args.clip_min),
            clip_max=float(args.clip_max),
            row=row,
            out_parts_dir=parts_dir,
            fraction=float(args.fraction),
        )
        for row in rows
    )
    for r in results:
        if r is not None:
            failures.append(r)

    pd.DataFrame([asdict(f) for f in failures], columns=["condition_id", "error_type", "error_message"]).to_csv(
        failures_csv, index=False
    )
    logging.info("Wrote failures log: %s (failed=%d)", failures_csv, len(failures))

    # Collect parts -> fva_all + features
    try:
        fva_all, features = collect_and_build_features(parts_dir=parts_dir, conditions_df=conditions_df)
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Collect/build features failed: {e}", file=sys.stderr)
        return 2

    fva_all_path = outdir / "fva_all.parquet"
    features_path = outdir / "features.parquet"
    save_table(fva_all, fva_all_path, fmt="parquet")
    save_table(features, features_path, fmt="parquet")

    # Regime table (FBA-based saturation) with per-condition ATPM_eff
    regime_fba = _run_regime_table_per_condition(
        model_path=str(args.model),
        conditions_df=conditions_df,
        medium_cfg=medium_cfg,
        regime_cfg=regime_cfg,
        atpm_rid=str(args.atpm_rid),
        a=a,
        b=b,
        clip_min=float(args.clip_min),
        clip_max=float(args.clip_max),
        eps=1e-6,
        infty_bound=999.0,
    )
    regime_fba_path = outdir / "regime_fba.parquet"
    save_table(regime_fba, regime_fba_path, fmt="parquet")

    # XAI report (baseline) under outdir/xai/
    xai_outdir = outdir / "xai"
    try:
        _run_xai_report(
            python_exe=sys.executable,
            features_path=features_path,
            outdir=xai_outdir,
            conditions=str(args.conditions),
            medium=str(args.medium),
            regime_fba=regime_fba_path,
        )
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] xai_report failed: {e}", file=sys.stderr)
        return 2

    print(f"[OK] Outputs in: {outdir}")
    print(f"[OK] Wrote: {features_path}")
    print(f"[OK] Wrote: {regime_fba_path}")
    print(f"[OK] Wrote: {xai_outdir / 'regime_table.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

