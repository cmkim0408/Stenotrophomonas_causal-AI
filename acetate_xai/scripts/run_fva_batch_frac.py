from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed

from acetate_xai.config import load_config
from acetate_xai.fva import FVAError, run_targeted_fva
from acetate_xai.io import load_conditions_csv, load_sbml_model, save_table
from acetate_xai.medium import MediumConfigError, apply_condition_to_model


@dataclass(frozen=True)
class Failure:
    condition_id: str
    error_type: str
    error_message: str


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run FBA + targeted FVA for all conditions (parallel).")
    p.add_argument("--model", required=True, help="SBML model path (e.g., models/model.xml)")
    p.add_argument("--conditions", required=True, help="Conditions CSV (long format)")
    p.add_argument("--targets", required=True, help="targets JSON (array of reaction ids)")
    p.add_argument("--medium", required=True, help="medium YAML config (e.g., configs/medium_base.yaml)")
    p.add_argument("--outdir", required=True, help="Output directory (e.g., results/fva_parts)")
    p.add_argument("--n-jobs", type=int, default=1, help="Parallel workers (joblib). Use 1 to disable.")
    p.add_argument("--backend", default="loky", choices=["loky", "threading"], help="joblib backend.")
    p.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Run only the first N conditions from the CSV (for small desktop tests).",
    )
    p.add_argument(
        "--condition-ids",
        nargs="+",
        default=None,
        help="Run only specified condition_id values (space-separated). Overrides --limit.",
    )
    p.add_argument(
        "--fraction",
        type=float,
        default=0.95,
        help="fraction_of_optimum for targeted FVA (default: 0.95)",
    )
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
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


def _run_one_condition(
    *,
    model_path: str,
    medium_cfg: dict,
    targets: list[str],
    row: dict,
    outdir: str,
    fraction: float,
) -> Failure | None:
    cid = str(row.get("condition_id", ""))
    if not cid:
        return Failure(condition_id="", error_type="ValueError", error_message="Missing condition_id in row")

    out_dir = Path(outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"condition_id={_safe_filename_component(cid)}.parquet"

    try:
        model = load_sbml_model(model_path)
        apply_condition_to_model(model, row, medium_cfg)
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


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(message)s",
    )

    model_path = str(Path(args.model))
    if not Path(model_path).exists():
        print(f"[ERROR] Model file not found: {model_path}", file=sys.stderr)
        return 2

    if not (0.0 < float(args.fraction) <= 1.0):
        print("[ERROR] --fraction must be in (0, 1].", file=sys.stderr)
        return 2

    try:
        medium_cfg = load_config(args.medium)
        targets = _load_targets_json(args.targets)
        cond_df = load_conditions_csv(args.conditions)
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Failed to load inputs: {e}", file=sys.stderr)
        return 2

    if args.condition_ids:
        wanted = set(str(x) for x in args.condition_ids)
        before = len(cond_df)
        cond_df = cond_df.loc[cond_df["condition_id"].astype(str).isin(wanted)].copy()
        after = len(cond_df)
        if after == 0:
            example = ", ".join(load_conditions_csv(args.conditions)["condition_id"].astype(str).head(10).tolist())
            print(
                f"[ERROR] None of the requested --condition-ids were found. Example available ids: {example}",
                file=sys.stderr,
            )
            return 2
        logging.info("Subset by --condition-ids: %d -> %d", before, after)
    elif args.limit is not None:
        if args.limit <= 0:
            print("[ERROR] --limit must be a positive integer", file=sys.stderr)
            return 2
        before = len(cond_df)
        cond_df = cond_df.head(int(args.limit)).copy()
        logging.info("Subset by --limit: %d -> %d", before, len(cond_df))

    outdir = str(Path(args.outdir))
    Path(outdir).mkdir(parents=True, exist_ok=True)

    rows = cond_df.to_dict(orient="records")
    logging.info(
        "Running batch: n_conditions=%d, n_targets=%d, n_jobs=%d, backend=%s, fraction=%.3f",
        len(rows),
        len(targets),
        args.n_jobs,
        args.backend,
        float(args.fraction),
    )

    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    failures: list[Failure] = []
    results = Parallel(n_jobs=int(args.n_jobs), backend=args.backend)(
        delayed(_run_one_condition)(
            model_path=model_path,
            medium_cfg=medium_cfg,
            targets=targets,
            row=row,
            outdir=outdir,
            fraction=float(args.fraction),
        )
        for row in rows
    )
    for r in results:
        if r is not None:
            failures.append(r)

    failed_path = Path("results") / "failed_conditions.csv"
    failed_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [asdict(f) for f in failures],
        columns=["condition_id", "error_type", "error_message"],
    ).to_csv(failed_path, index=False)

    print(f"[OK] Completed conditions={len(rows)}; failed={len(failures)}")
    print(f"[OK] Outputs in: {outdir}")
    print(f"[OK] Failed conditions log: {failed_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

