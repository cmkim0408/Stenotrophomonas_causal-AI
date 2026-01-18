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
    sample_id: str
    error_type: str
    error_message: str


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Large LHS sampling campaign (4D): always FBA+regime labels; optional targeted FVA."
    )
    p.add_argument("--model", required=True, help="SBML model path")
    p.add_argument("--outdir", default="results/campaigns/C_random_LHS", help="Output directory")
    p.add_argument("--medium", default="acetate_xai/configs/medium_base.yaml", help="Base medium YAML")
    p.add_argument("--regime-config", default="acetate_xai/configs/regime_exchanges.yaml", help="Regime exchanges YAML")

    p.add_argument("--n-samples", type=int, default=2000, help="Number of LHS samples (default: 2000)")
    p.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")

    p.add_argument("--n-jobs", type=int, default=16, help="Parallel workers (default: 16). Use -1 for all cores.")
    p.add_argument("--backend", default="loky", choices=["loky", "threading"], help="joblib backend")

    # knobs / mapping choices
    p.add_argument(
        "--acetate-mode",
        default="one_to_one",
        choices=["one_to_one", "k_ac"],
        help="How to map acetate_mM to EX_ac_e uptake cap: one_to_one sets lb=-acetate_mM; k_ac uses medium scaling.",
    )
    p.add_argument("--o2-rxn", default="EX_o2_e", help="O2 exchange rxn id (default: EX_o2_e)")
    p.add_argument("--nh4-rxn", default="EX_nh4_e", help="NH4 exchange rxn id (default: EX_nh4_e)")
    p.add_argument("--atpm-rxn", default="ATPM", help="ATPM rxn id (default: ATPM)")

    # FVA controls
    p.add_argument("--do-fva", type=int, default=0, choices=[0, 1], help="Whether to run targeted FVA (default: 0)")
    p.add_argument("--targets", default="acetate_xai/configs/targets_120.json", help="Targets JSON for FVA")
    p.add_argument("--fraction", type=float, default=0.95, help="FVA fraction_of_optimum (default: 0.95)")
    p.add_argument("--fva-subsample", type=int, default=None, help="Run FVA for N random samples only")
    p.add_argument(
        "--fva-phase-boundary",
        type=int,
        default=None,
        help="Auto-select N samples near phase boundaries (recommended) and run FVA only for them",
    )
    p.add_argument("--fva-boundary-k", type=int, default=12, help="kNN neighborhood size for boundary scoring (default: 12)")
    p.add_argument("--keep-parts", action="store_true", help="Keep fva_parts parquet files (default: off)")
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


def _lhs_design(n: int, *, seed: int) -> pd.DataFrame:
    """
    4D LatinHypercube sampling:
      - acetate_mM: 0..200
      - o2_uptake_max: 0..20   (EX_o2_e.lb = -o2_uptake_max)
      - nh4_uptake_max: 0..10  (EX_nh4_e.lb = -nh4_uptake_max)
      - atpm: 0..25            (ATPM lb=ub=atpm)
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
    lo = np.array([0.0, 0.0, 0.0, 0.0], dtype=float)
    hi = np.array([200.0, 20.0, 10.0, 25.0], dtype=float)
    x = qmc.scale(u, lo, hi)
    return pd.DataFrame(x, columns=["acetate_mM", "o2_uptake_max", "nh4_uptake_max", "atpm"])


def _set_exchange_uptake_cap(model, rid: str, uptake_max: float) -> None:
    # uptake is negative flux
    rxn = model.reactions.get_by_id(rid)
    u = max(0.0, float(uptake_max))
    if rxn.upper_bound < 0.0:
        rxn.upper_bound = 0.0
    rxn.lower_bound = -u


def _fix_flux(model, rid: str, value: float) -> None:
    rxn = model.reactions.get_by_id(rid)
    v = float(value)
    # Set ub first to avoid transient lb>ub when raising ub
    rxn.upper_bound = v
    rxn.lower_bound = v


def _primary_regime_label(rec: dict[str, Any]) -> str:
    # Same precedence as build_regime_dataset.py label logic
    if bool(rec.get("acetate_sat", False)):
        return "Ac_limited"
    if bool(rec.get("nh4_sat", False)):
        return "N_limited"
    if bool(rec.get("pi_sat", False)):
        return "Pi_limited"
    if bool(rec.get("o2_sat", False)):
        return "O2_limited"
    return "Unconstrained"


def _run_fba_one(
    *,
    model_path: str,
    medium_cfg: dict[str, Any],
    regime_cfg: dict[str, Any],
    acetate_mode: str,
    o2_rxn: str,
    nh4_rxn: str,
    atpm_rxn: str,
    row: dict[str, Any],
) -> tuple[dict[str, Any] | None, Failure | None]:
    sid = str(row["sample_id"])
    try:
        model = load_sbml_model(model_path)

        # Apply base medium bounds + (optionally) acetate scaling through existing helper.
        # We keep yeast disabled and nh4cl_gL=0 here; we override O2/NH4 by direct uptake caps after.
        cond_row = {
            "condition_id": sid,
            "pH0": 7.0,
            "yeast_extract_gL": 0.0,
            "nh4cl_gL": 0.0,
            "acetate_mM": float(row["acetate_mM"]),
        }
        apply_condition_to_model(model, cond_row, medium_cfg)

        # Override acetate mapping if requested (one_to_one: lb=-acetate_mM)
        if acetate_mode == "one_to_one":
            try:
                _set_exchange_uptake_cap(model, str(medium_cfg["exchanges"]["acetate"]), uptake_max=float(row["acetate_mM"]))
            except Exception:
                # fallback: try EX_ac_e
                _set_exchange_uptake_cap(model, "EX_ac_e", uptake_max=float(row["acetate_mM"]))

        # Override O2/NH4 uptake caps
        _set_exchange_uptake_cap(model, str(o2_rxn), uptake_max=float(row["o2_uptake_max"]))
        _set_exchange_uptake_cap(model, str(nh4_rxn), uptake_max=float(row["nh4_uptake_max"]))

        # ATPM rxnfix
        try:
            _fix_flux(model, str(atpm_rxn), float(row["atpm"]))
        except KeyError:
            logger.warning("ATPM rxn not found in model (skipped): %s", atpm_rxn)

        sol = model.optimize()
        status = str(sol.status)
        if status != "optimal":
            rec = {
                "sample_id": sid,
                "acetate_mM": float(row["acetate_mM"]),
                "o2_uptake_max": float(row["o2_uptake_max"]),
                "nh4_uptake_max": float(row["nh4_uptake_max"]),
                "atpm": float(row["atpm"]),
                "objective_value": np.nan,
                "status": status,
            }
            # still include sat flags as False/NA
            rec.update({"acetate_sat": False, "o2_sat": False, "nh4_sat": False, "pi_sat": False, "primary_regime": "Unconstrained"})
            return rec, None

        objective_value = float(sol.objective_value)

        # Saturation flags (use regime cfg candidates; fall back to canonical IDs)
        def _pick(key: str, fallback: str) -> str | None:
            cand = regime_cfg.get(key, [])
            if isinstance(cand, list) and cand:
                rid = pick_first_existing_reaction_id(model, [str(x) for x in cand])
                if rid:
                    return rid
            try:
                model.reactions.get_by_id(fallback)
                return fallback
            except KeyError:
                return None

        rid_ac = _pick("acetate", "EX_ac_e")
        rid_o2 = _pick("oxygen", str(o2_rxn))
        rid_nh4 = _pick("ammonium", str(nh4_rxn))
        rid_pi = _pick("phosphate", "EX_pi_e")

        def _sat(rid: str | None, prefix: str) -> dict[str, Any]:
            if not rid:
                return {f"{prefix}_rid": "", f"{prefix}_sat": False, f"{prefix}_sat_side": "missing", f"{prefix}_flux": np.nan}
            s = compute_saturation_for_reaction(rid=rid, model=model, solution=sol)
            return {f"{prefix}_rid": s.rid, f"{prefix}_sat": bool(s.saturated), f"{prefix}_sat_side": s.sat_side, f"{prefix}_flux": s.flux}

        rec: dict[str, Any] = {
            "sample_id": sid,
            "acetate_mM": float(row["acetate_mM"]),
            "o2_uptake_max": float(row["o2_uptake_max"]),
            "nh4_uptake_max": float(row["nh4_uptake_max"]),
            "atpm": float(row["atpm"]),
            "objective_value": objective_value,
            "status": status,
        }
        rec.update(_sat(rid_ac, "acetate"))
        rec.update(_sat(rid_o2, "o2"))
        rec.update(_sat(rid_nh4, "nh4"))
        rec.update(_sat(rid_pi, "pi"))
        rec["primary_regime"] = _primary_regime_label(rec)
        return rec, None
    except (FileNotFoundError, MediumConfigError, ValueError, MediumConfigError, FVAError) as e:
        return None, Failure(sample_id=sid, error_type=type(e).__name__, error_message=str(e))
    except Exception as e:  # noqa: BLE001
        return None, Failure(sample_id=sid, error_type=type(e).__name__, error_message=str(e))


def _boundary_subsample(design: pd.DataFrame, labels: pd.Series, *, k: int, n_pick: int, seed: int) -> list[str]:
    """
    Pick points near phase boundaries by kNN label disagreement.
    This keeps the sampling dense around transitions without mixing runs.
    """
    from sklearn.neighbors import NearestNeighbors

    X = design[["acetate_mM", "o2_uptake_max", "nh4_uptake_max", "atpm"]].to_numpy(dtype=float)
    # normalize to [0,1] ranges
    mins = X.min(axis=0)
    maxs = X.max(axis=0)
    denom = np.maximum(1e-12, maxs - mins)
    Xn = (X - mins) / denom

    nn = NearestNeighbors(n_neighbors=min(int(k), len(Xn)), algorithm="auto")
    nn.fit(Xn)
    idx = nn.kneighbors(Xn, return_distance=False)
    y = labels.to_numpy()
    # disagreement score: fraction of neighbors with different label
    scores = []
    for i in range(len(Xn)):
        neigh = idx[i]
        if len(neigh) <= 1:
            scores.append(0.0)
            continue
        scores.append(float(np.mean(y[neigh] != y[i])))
    s = np.asarray(scores)
    # stable tie-break by tiny noise
    rng = np.random.default_rng(int(seed))
    s = s + 1e-9 * rng.standard_normal(size=len(s))
    pick_idx = np.argsort(s)[::-1][: int(min(n_pick, len(s)))]
    return design.iloc[pick_idx]["sample_id"].astype(str).tolist()


def _save_metadata(outdir: Path, payload: dict[str, Any]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "run_metadata.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    model_path = str(Path(args.model))
    if not Path(model_path).exists():
        print(f"[ERROR] model not found: {model_path}", file=sys.stderr)
        return 2

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # load configs
    try:
        medium_cfg = load_config(args.medium)
        regime_cfg = load_config(args.regime_config)
        targets = _load_targets_json(args.targets) if int(args.do_fva) == 1 else []
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] failed to load config/targets: {e}", file=sys.stderr)
        return 2

    # metadata
    _save_metadata(
        outdir,
        {
            "n_samples": int(args.n_samples),
            "seed": int(args.seed),
            "ranges": {
                "acetate_mM": [0, 200],
                "o2_uptake_max": [0, 20],
                "nh4_uptake_max": [0, 10],
                "atpm": [0, 25],
            },
            "model": str(args.model),
            "medium": str(args.medium),
            "regime_config": str(args.regime_config),
            "do_fva": int(args.do_fva),
            "targets": str(args.targets) if int(args.do_fva) == 1 else "",
            "fraction_of_optimum": float(args.fraction),
            "solver": "cobra-default",
        },
    )

    # design
    design = _lhs_design(int(args.n_samples), seed=int(args.seed))
    design.insert(0, "sample_id", [f"LHS_{i:05d}" for i in range(len(design))])
    save_table(design, outdir / "design.parquet", fmt="parquet")

    # oversubscription guardrails
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    if int(args.n_jobs) == -1:
        logger.warning("[WARN] --n-jobs=-1 can cause memory spikes (2000 model loads). Consider --n-jobs 16.")

    rows = design.to_dict(orient="records")
    logger.info("Running FBA+regime for n=%d samples (n_jobs=%d)", len(rows), int(args.n_jobs))

    results = Parallel(n_jobs=int(args.n_jobs), backend=args.backend)(
        delayed(_run_fba_one)(
            model_path=model_path,
            medium_cfg=medium_cfg,
            regime_cfg=regime_cfg,
            acetate_mode=str(args.acetate_mode),
            o2_rxn=str(args.o2_rxn),
            nh4_rxn=str(args.nh4_rxn),
            atpm_rxn=str(args.atpm_rxn),
            row=row,
        )
        for row in rows
    )

    records: list[dict[str, Any]] = []
    failures: list[Failure] = []
    for rec, fail in results:
        if rec is not None:
            records.append(rec)
        if fail is not None:
            failures.append(fail)

    if failures:
        pd.DataFrame([asdict(f) for f in failures]).to_csv(outdir / "failed_samples.csv", index=False)

    labels = pd.DataFrame(records)
    # ensure stable order (one row per sample_id)
    labels = design[["sample_id", "acetate_mM", "o2_uptake_max", "nh4_uptake_max", "atpm"]].merge(
        labels,
        on=["sample_id", "acetate_mM", "o2_uptake_max", "nh4_uptake_max", "atpm"],
        how="left",
    )

    # minimal required columns
    keep = [
        "sample_id",
        "acetate_mM",
        "o2_uptake_max",
        "nh4_uptake_max",
        "atpm",
        "objective_value",
        "primary_regime",
        "acetate_sat",
        "o2_sat",
        "nh4_sat",
        "pi_sat",
        "status",
    ]
    for c in keep:
        if c not in labels.columns:
            labels[c] = pd.NA
    save_table(labels[keep], outdir / "regime_labels.parquet", fmt="parquet")

    # features.parquet always exists (FBA-only baseline)
    base_features = labels[["sample_id", "objective_value", "acetate_mM", "o2_uptake_max", "nh4_uptake_max", "atpm"]].copy()
    save_table(base_features, outdir / "features.parquet", fmt="parquet")

    # Optional FVA
    if int(args.do_fva) != 1:
        print(f"[OK] Wrote: {outdir / 'design.parquet'}")
        print(f"[OK] Wrote: {outdir / 'regime_labels.parquet'}")
        print(f"[OK] Wrote: {outdir / 'features.parquet'} (FBA-only)")
        return 0

    # Choose which samples get FVA
    sample_ids = design["sample_id"].astype(str).tolist()
    chosen: list[str]
    if args.fva_phase_boundary is not None:
        chosen = _boundary_subsample(
            design=design,
            labels=labels["primary_regime"].fillna("Unconstrained").astype(str),
            k=int(args.fva_boundary_k),
            n_pick=int(args.fva_phase_boundary),
            seed=int(args.seed),
        )
        logger.info("FVA selection: phase-boundary n=%d", len(chosen))
    elif args.fva_subsample is not None:
        rng = np.random.default_rng(int(args.seed))
        n_pick = int(min(int(args.fva_subsample), len(sample_ids)))
        chosen = [sample_ids[i] for i in rng.choice(len(sample_ids), size=n_pick, replace=False)]
        logger.info("FVA selection: random subsample n=%d", len(chosen))
    else:
        chosen = sample_ids
        logger.info("FVA selection: all samples n=%d", len(chosen))

    parts_dir = outdir / "fva_parts"
    parts_dir.mkdir(parents=True, exist_ok=True)

    # FVA runner reuses the same worker logic but only for chosen samples
    chosen_rows = design.loc[design["sample_id"].astype(str).isin(set(chosen))].to_dict(orient="records")

    def _run_fva_for_one(row: dict[str, Any]) -> Failure | None:
        sid = str(row["sample_id"])
        out_path = parts_dir / f"sample_id={sid}.parquet"
        try:
            model = load_sbml_model(model_path)
            cond_row = {
                "condition_id": sid,
                "pH0": 7.0,
                "yeast_extract_gL": 0.0,
                "nh4cl_gL": 0.0,
                "acetate_mM": float(row["acetate_mM"]),
            }
            apply_condition_to_model(model, cond_row, medium_cfg)

            if str(args.acetate_mode) == "one_to_one":
                _set_exchange_uptake_cap(model, str(medium_cfg["exchanges"]["acetate"]), uptake_max=float(row["acetate_mM"]))

            _set_exchange_uptake_cap(model, str(args.o2_rxn), uptake_max=float(row["o2_uptake_max"]))
            _set_exchange_uptake_cap(model, str(args.nh4_rxn), uptake_max=float(row["nh4_uptake_max"]))

            try:
                _fix_flux(model, str(args.atpm_rxn), float(row["atpm"]))
            except KeyError:
                pass

            sol = model.optimize()
            if str(sol.status) != "optimal":
                raise FVAError(f"FBA failed: status={sol.status}")
            obj = float(sol.objective_value)

            fva_df = run_targeted_fva(model, targets=targets, fraction_of_optimum=float(args.fraction))
            out = fva_df.copy()
            out.insert(0, "sample_id", sid)
            out.insert(1, "objective_value", obj)
            out = out[["sample_id", "objective_value", "reaction_id", "fva_min", "fva_max"]]
            save_table(out, out_path, fmt="parquet")
            return None
        except Exception as e:  # noqa: BLE001
            return Failure(sample_id=sid, error_type=type(e).__name__, error_message=str(e))

    logger.info("Running targeted FVA: n_samples=%d, n_targets=%d", len(chosen_rows), len(targets))
    fva_failures = Parallel(n_jobs=int(args.n_jobs), backend=args.backend)(delayed(_run_fva_for_one)(r) for r in chosen_rows)
    fva_failures = [f for f in fva_failures if f is not None]
    if fva_failures:
        pd.DataFrame([asdict(f) for f in fva_failures]).to_csv(outdir / "failed_fva_samples.csv", index=False)

    # Build wide features and left-join onto base_features
    from acetate_xai.collect import build_fva_long_features, build_wide_feature_matrix

    files = sorted(parts_dir.glob("*.parquet"))
    if files:
        fva_all = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
        long_feat = build_fva_long_features(
            fva_all.rename(columns={"sample_id": "condition_id"})  # reuse existing helper schema
        ).rename(columns={"condition_id": "sample_id"})
        # adapt wide builder by temporarily renaming
        wide = build_wide_feature_matrix(
            long_feat.rename(columns={"sample_id": "condition_id"})
        ).rename(columns={"condition_id": "sample_id"})

        features = base_features.merge(wide, on="sample_id", how="left")
        save_table(features, outdir / "features.parquet", fmt="parquet")

    if not args.keep_parts:
        try:
            for f in parts_dir.glob("*.parquet"):
                f.unlink()
            parts_dir.rmdir()
        except Exception:
            pass

    print(f"[OK] Wrote: {outdir / 'design.parquet'}")
    print(f"[OK] Wrote: {outdir / 'regime_labels.parquet'}")
    print(f"[OK] Wrote: {outdir / 'features.parquet'} (FBA + optional FVA columns)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

