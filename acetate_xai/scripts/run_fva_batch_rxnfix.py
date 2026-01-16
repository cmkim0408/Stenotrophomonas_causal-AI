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
    p = argparse.ArgumentParser(
        description="Run FBA + targeted FVA for all conditions (parallel) with fixed flux constraints."
    )
    # Same args as scripts/run_fva_batch.py
    p.add_argument("--model", required=True, help="SBML model path (e.g., models/model.xml)")
    p.add_argument("--conditions", required=True, help="Conditions CSV (long format)")
    p.add_argument("--targets", required=True, help="targets JSON (array of reaction ids)")
    p.add_argument("--medium", required=True, help="medium YAML config (e.g., configs/medium_base.yaml)")
    p.add_argument("--outdir", required=True, help="Output directory (e.g., results/fva_parts)")
    p.add_argument("--n-jobs", type=int, default=1, help="Parallel workers (joblib). Use 1 to disable.")
    p.add_argument("--backend", default="loky", choices=["loky", "threading"], help="joblib backend.")
    p.add_argument("--limit", type=int, default=None, help="Run only the first N conditions from the CSV.")
    p.add_argument(
        "--condition-ids",
        nargs="+",
        default=None,
        help="Run only specified condition_id values (space-separated). Overrides --limit.",
    )
    # New: reaction fixed flux constraints
    p.add_argument(
        "--rxn-fix",
        action="append",
        default=[],
        help='Fix reaction flux: "RXNID=value". Repeatable. Example: --rxn-fix ATPM=20',
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


def _parse_rxn_fix(items: list[str]) -> list[tuple[str, float]]:
    """
    Parse ["ATPM=20", "EX_o2_e=-10"] into [(rid, value), ...]
    """
    out: list[tuple[str, float]] = []
    for raw in items:
        if raw is None:
            continue
        s = str(raw).strip()
        if not s:
            continue
        if "=" not in s:
            raise ValueError(f"--rxn-fix must be in form RXNID=value, got: {s!r}")
        rid, val = s.split("=", 1)
        rid = rid.strip()
        if not rid:
            raise ValueError(f"--rxn-fix has empty reaction id: {s!r}")
        try:
            v = float(val.strip())
        except Exception as e:  # noqa: BLE001
            raise ValueError(f"--rxn-fix has non-numeric value: {s!r}") from e
        out.append((rid, v))
    return out


def _apply_rxn_fix(model, rxn_fix: list[tuple[str, float]]) -> None:
    """
    For each (rxn_id, value):
      - if present: set ub=1000 (first), then lb=ub=value (fixed flux)
      - if missing: warn
    """
    present = set(r.id for r in model.reactions)
    for rid, value in rxn_fix:
        if rid not in present:
            logging.warning("Reaction not found for --rxn-fix (skipped): %s", rid)
            continue
        rxn = model.reactions.get_by_id(rid)
        # Make sure we don't hit transient bound invalidity
        rxn.upper_bound = 1000.0
        rxn.lower_bound = float(value)
        rxn.upper_bound = float(value)
        logging.info("Fixed flux %s = %s", rid, value)


def _run_one_condition(
    *,
    model_path: str,
    medium_cfg: dict,
    targets: list[str],
    rxn_fix: list[tuple[str, float]],
    row: dict,
    outdir: str,
    fraction: float = 0.95,
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
        _apply_rxn_fix(model, rxn_fix)

        sol = model.optimize()
        if sol.status != "optimal":
            raise FVAError(f"FBA failed: status={sol.status}")

        objective_value = float(sol.objective_value)
        fva_df = run_targeted_fva(model, targets=targets, fraction_of_optimum=fraction)

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

    try:
        medium_cfg = load_config(args.medium)
        targets = _load_targets_json(args.targets)
        cond_df = load_conditions_csv(args.conditions)
        rxn_fix = _parse_rxn_fix(list(args.rxn_fix or []))
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Failed to load inputs: {e}", file=sys.stderr)
        return 2

    # Subset selection for desktop tests
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
        "Running batch (rxnfix): n_conditions=%d, n_targets=%d, n_jobs=%d, backend=%s, n_rxn_fix=%d",
        len(rows),
        len(targets),
        args.n_jobs,
        args.backend,
        len(rxn_fix),
    )

    # Avoid oversubscription
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    failures: list[Failure] = []
    results = Parallel(n_jobs=int(args.n_jobs), backend=args.backend)(
        delayed(_run_one_condition)(
            model_path=model_path,
            medium_cfg=medium_cfg,
            targets=targets,
            rxn_fix=rxn_fix,
            row=row,
            outdir=outdir,
            fraction=0.95,
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

