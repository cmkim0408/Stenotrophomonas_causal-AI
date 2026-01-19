from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from acetate_xai.config import load_config
from acetate_xai.medium import apply_condition_to_model
from acetate_xai.io import load_sbml_model, save_table

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Failure:
    sample_id: str
    error_type: str
    error_message: str


MUST_TRY_TARGETS: list[list[str]] = [
    # Exchanges
    ["EX_ac_e"],
    ["EX_o2_e"],
    ["EX_nh4_e"],
    ["EX_h_e"],
    ["EX_pi_e"],
    ["EX_co2_e"],
    ["EX_h2o_e"],
    # Maintenance/Energy
    ["ATPM"],
    ["ATPS4rpp"],
    ["NADH16pp", "NADH16"],
    ["CYTBO3_4pp", "CYTBO3"],
    # Central carbon
    ["ACONT"],
    ["CS"],
    ["ICL"],
    ["PPC"],
    ["PFK"],
    ["ENO"],
    ["PYK"],
    ["MDH"],
    ["ICDHyr"],
]


def _utc_now_iso() -> str:
    return datetime.now(tz=timezone.utc).isoformat()


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Tier A (Light) FVA for all feasible LHS samples: 30-target panel, chunked + restartable."
    )
    p.add_argument("--model", required=True, help="SBML model path (e.g., acetate_xai/models/model.xml)")
    p.add_argument(
        "--labels",
        required=True,
        help="regime_labels.parquet (must include status=='optimal', sample_id, acetate_mM, o2_uptake_max, nh4_uptake_max, atpm)",
    )
    p.add_argument(
        "--outdir",
        required=True,
        help="Output directory (e.g., results/campaigns/C_random_LHS_FVA_LIGHT)",
    )
    p.add_argument(
        "--medium",
        default="acetate_xai/configs/medium_base.yaml",
        help="Base medium YAML (default: acetate_xai/configs/medium_base.yaml)",
    )
    p.add_argument("--fraction", type=float, default=0.95, help="fraction_of_optimum for FVA (default: 0.95)")
    p.add_argument("--chunk-size", type=int, default=50, help="Chunk size (default: 50)")
    p.add_argument("--n-jobs", type=int, default=8, help="Parallel workers (default: 8)")
    p.add_argument("--backend", default="loky", choices=["loky", "threading"], help="joblib backend")
    p.add_argument("--merge-only", action="store_true", help="Only merge existing fva_parts into features_dense.parquet")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def _resolve_targets(model_path: str, outdir: Path) -> list[str]:
    model = load_sbml_model(model_path)
    rxn_ids = {r.id for r in model.reactions}

    resolved: list[str] = []
    for group in MUST_TRY_TARGETS:
        chosen = None
        for rid in group:
            if rid in rxn_ids:
                chosen = rid
                break
        if chosen is not None and chosen not in resolved:
            resolved.append(chosen)

    # Save the resolved target list for reproducibility
    out_json = outdir / "targets_30_resolved.json"
    out_json.write_text(json.dumps(resolved, indent=2), encoding="utf-8")
    logger.info("Resolved targets: requested_groups=%d, resolved=%d -> %s", len(MUST_TRY_TARGETS), len(resolved), out_json)
    return resolved


def _set_exchange_uptake_cap(model, rid: str, uptake_max: float) -> None:
    rxn = model.reactions.get_by_id(rid)
    u = max(0.0, float(uptake_max))
    if rxn.upper_bound < 0.0:
        rxn.upper_bound = 0.0
    rxn.lower_bound = -u


def _fix_flux(model, rid: str, value: float) -> None:
    rxn = model.reactions.get_by_id(rid)
    v = float(value)
    # Set ub first to avoid transient lb>ub
    rxn.upper_bound = v
    rxn.lower_bound = v


def _safe_filename(s: str) -> str:
    return str(s).strip().replace("/", "_").replace("\\", "_")


def _fva_one_sample(
    *,
    model_path: str,
    medium_cfg: dict[str, Any],
    targets: list[str],
    fraction: float,
    sample: dict[str, Any],
) -> tuple[pd.DataFrame | None, Failure | None]:
    sid = str(sample["sample_id"])
    try:
        model = load_sbml_model(model_path)

        # Apply base medium bounds (minimal). We keep yeast disabled and nh4cl_gL=0 here.
        # Then override caps using the LHS knobs.
        cond_row = {
            "condition_id": sid,
            "pH0": 7.0,
            "yeast_extract_gL": 0.0,
            "nh4cl_gL": 0.0,
            "acetate_mM": float(sample["acetate_mM"]),
        }
        apply_condition_to_model(model, cond_row, medium_cfg)

        # Override LHS knob caps (1:1 mapping)
        _set_exchange_uptake_cap(model, "EX_ac_e", uptake_max=float(sample["acetate_mM"]))
        _set_exchange_uptake_cap(model, "EX_o2_e", uptake_max=float(sample["o2_uptake_max"]))
        _set_exchange_uptake_cap(model, "EX_nh4_e", uptake_max=float(sample["nh4_uptake_max"]))

        # ATPM rxnfix (ub open handled by setting ub=v here)
        try:
            _fix_flux(model, "ATPM", float(sample["atpm"]))
        except KeyError:
            # If missing, we still run; but warn once per worker.
            logger.warning("ATPM not found in model; skipping rxnfix for sample_id=%s", sid)

        sol = model.optimize()
        if str(sol.status) != "optimal":
            raise RuntimeError(f"FBA status={sol.status}")
        objective_value = float(sol.objective_value)

        try:
            from cobra.flux_analysis import flux_variability_analysis
        except Exception as e:  # noqa: BLE001
            raise RuntimeError(f"Missing cobra FVA: {e}") from e

        fva = flux_variability_analysis(model, reaction_list=targets, fraction_of_optimum=float(fraction))
        fva = fva.rename(columns={"minimum": "fva_min", "maximum": "fva_max"}).reset_index().rename(columns={"index": "reaction_id"})
        fva["reaction_id"] = fva["reaction_id"].astype(str)
        fva["fva_min"] = pd.to_numeric(fva["fva_min"], errors="coerce")
        fva["fva_max"] = pd.to_numeric(fva["fva_max"], errors="coerce")
        fva["width"] = fva["fva_max"] - fva["fva_min"]
        fva["mid"] = (fva["fva_max"] + fva["fva_min"]) / 2.0

        base = {
            "sample_id": sid,
            "acetate_mM": float(sample["acetate_mM"]),
            "o2_uptake_max": float(sample["o2_uptake_max"]),
            "nh4_uptake_max": float(sample["nh4_uptake_max"]),
            "atpm": float(sample["atpm"]),
            "objective_value": objective_value,
        }
        if "primary_regime" in sample:
            base["primary_regime"] = str(sample["primary_regime"])

        out = fva.copy()
        for k, v in base.items():
            out.insert(0, k, v)
        out = out[
            [
                "sample_id",
                "acetate_mM",
                "o2_uptake_max",
                "nh4_uptake_max",
                "atpm",
                "objective_value",
                "primary_regime",
                "reaction_id",
                "fva_min",
                "fva_max",
                "width",
                "mid",
            ]
        ]
        return out, None
    except Exception as e:  # noqa: BLE001
        return None, Failure(sample_id=sid, error_type=type(e).__name__, error_message=str(e))


def _write_failed(out_csv: Path, failures: list[Failure]) -> None:
    if not failures:
        return
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    exists = out_csv.exists()
    with out_csv.open("a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["sample_id", "error_type", "error_message"])
        if not exists:
            w.writeheader()
        for fa in failures:
            w.writerow(asdict(fa))


def _merge_parts(parts_dir: Path, out_features: Path) -> tuple[int, int]:
    parts = sorted(parts_dir.glob("part_*.parquet"))
    if not parts:
        raise FileNotFoundError(f"No part_*.parquet found under: {parts_dir}")

    df = pd.concat((pd.read_parquet(p) for p in parts), ignore_index=True)
    # Long -> wide features
    df["reaction_id"] = df["reaction_id"].astype(str)
    idx_cols = ["sample_id", "acetate_mM", "o2_uptake_max", "nh4_uptake_max", "atpm", "objective_value", "primary_regime"]
    base = df[idx_cols].drop_duplicates(subset=["sample_id"]).set_index("sample_id")

    w = df.pivot(index="sample_id", columns="reaction_id", values="width")
    w.columns = [f"width__{c}" for c in w.columns.astype(str)]
    m = df.pivot(index="sample_id", columns="reaction_id", values="mid")
    m.columns = [f"mid__{c}" for c in m.columns.astype(str)]

    features = base.join(w, how="left").join(m, how="left").reset_index()
    out_features.parent.mkdir(parents=True, exist_ok=True)
    save_table(features, out_features, fmt="parquet")
    return int(features.shape[0]), int(df["reaction_id"].nunique())


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    model_path = str(Path(args.model))
    labels_path = Path(args.labels)
    outdir = Path(args.outdir)
    parts_dir = outdir / "fva_parts"
    out_features = outdir / "features_dense.parquet"
    out_failed = outdir / "failed_samples.csv"
    out_meta = outdir / "run_metadata.json"

    outdir.mkdir(parents=True, exist_ok=True)
    parts_dir.mkdir(parents=True, exist_ok=True)

    # Avoid oversubscription on servers
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    if args.merge_only:
        n_rows, n_targets = _merge_parts(parts_dir, out_features)
        print(f"[OK] Merge-only complete: rows={n_rows}, targets={n_targets} -> {out_features}")
        return 0

    if not Path(model_path).exists():
        print(f"[ERROR] Model not found: {model_path}", file=sys.stderr)
        return 2
    if not labels_path.exists():
        print(f"[ERROR] Labels parquet not found: {labels_path}", file=sys.stderr)
        return 2

    medium_cfg = load_config(args.medium)
    targets = _resolve_targets(model_path, outdir)
    if not targets:
        print("[ERROR] Resolved targets list is empty; cannot run FVA.", file=sys.stderr)
        return 2

    df = pd.read_parquet(labels_path)
    if "status" not in df.columns or "sample_id" not in df.columns:
        raise ValueError("labels parquet must include columns: status, sample_id")

    # Feasible set filter
    df = df[df["status"].astype(str) == "optimal"].copy()
    need_cols = ["sample_id", "acetate_mM", "o2_uptake_max", "nh4_uptake_max", "atpm"]
    missing = [c for c in need_cols if c not in df.columns]
    if missing:
        raise ValueError(f"labels parquet missing required columns: {missing}")

    # Keep primary_regime if present
    keep_cols = need_cols + (["primary_regime"] if "primary_regime" in df.columns else [])
    df = df[keep_cols].copy()

    rows = df.to_dict(orient="records")
    n_total = len(rows)
    logger.info("Feasible samples: %d (status==optimal)", n_total)
    logger.info("Targets resolved: %d", len(targets))
    logger.info("Output dir: %s", outdir)

    # Metadata (written early)
    out_meta.write_text(
        json.dumps(
            {
                "timestamp_utc": _utc_now_iso(),
                "model": model_path,
                "labels": str(labels_path),
                "medium": str(args.medium),
                "fraction": float(args.fraction),
                "chunk_size": int(args.chunk_size),
                "n_jobs": int(args.n_jobs),
                "backend": str(args.backend),
                "n_feasible": int(n_total),
                "targets_resolved_path": str((outdir / "targets_30_resolved.json").as_posix()),
                "targets_resolved": targets,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    # tqdm optional
    try:
        from tqdm import tqdm  # type: ignore
    except Exception:  # noqa: BLE001
        tqdm = None

    chunk = int(args.chunk_size)
    if chunk <= 0:
        raise ValueError("--chunk-size must be positive")

    failures_total = 0
    parts_written = 0

    idxs = list(range(0, n_total, chunk))
    iterator = idxs
    if tqdm is not None:
        iterator = tqdm(idxs, desc="Chunks", unit="chunk")

    for start in iterator:
        end = min(n_total, start + chunk)
        part_path = parts_dir / f"part_{start}_{end}.parquet"
        if part_path.exists():
            continue  # restartable: skip completed chunk

        chunk_rows = rows[start:end]
        results = Parallel(n_jobs=int(args.n_jobs), backend=str(args.backend))(
            delayed(_fva_one_sample)(
                model_path=model_path,
                medium_cfg=medium_cfg,
                targets=targets,
                fraction=float(args.fraction),
                sample=s,
            )
            for s in chunk_rows
        )

        outs: list[pd.DataFrame] = []
        fails: list[Failure] = []
        for out_df, fail in results:
            if fail is not None:
                fails.append(fail)
            elif out_df is not None:
                outs.append(out_df)

        if outs:
            part_df = pd.concat(outs, ignore_index=True)
            save_table(part_df, part_path, fmt="parquet")
            parts_written += 1
        if fails:
            _write_failed(out_failed, fails)
            failures_total += len(fails)

    n_rows, n_targets = _merge_parts(parts_dir, out_features)
    print(f"[OK] Completed: feasible={n_total}, failed={failures_total}, parts_written={parts_written}")
    print(f"[OK] features_dense rows={n_rows}, targets={n_targets} -> {out_features}")
    print(f"[OK] targets -> {outdir / 'targets_30_resolved.json'}")
    print(f"[OK] metadata -> {out_meta}")
    if out_failed.exists():
        print(f"[OK] failures -> {out_failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

