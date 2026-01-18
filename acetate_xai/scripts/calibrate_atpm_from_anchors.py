from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from acetate_xai.config import load_config
from acetate_xai.io import load_conditions_csv, load_sbml_model, save_table
from acetate_xai.medium import apply_condition_to_model


@dataclass(frozen=True)
class ScanRow:
    anchor_id: str
    atpm: float
    objective_value: float
    status: str  # "optimal" | "infeasible" | "error"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Calibrate ATPM_eff from anchor conditions (grid scan + f(Ac) fit).")
    p.add_argument("--model", required=True, help="SBML model path (e.g., acetate_xai/models/model.xml)")
    p.add_argument("--conditions", required=True, help="Conditions CSV (long format)")
    p.add_argument("--medium", required=True, help="Medium YAML config (e.g., acetate_xai/configs/medium_base.yaml)")
    p.add_argument(
        "--anchor-ids",
        nargs="+",
        default=["AC_25", "AC_100", "AC_150"],
        help="Anchor condition_id values (default: AC_25 AC_100 AC_150)",
    )
    p.add_argument(
        "--atpm-grid",
        nargs="+",
        default=None,
        help="Explicit ATPM grid values, e.g. --atpm-grid 0 5 10 ... 200. If omitted uses --atpm-min/--atpm-max/--atpm-step.",
    )
    p.add_argument("--atpm-min", type=float, default=0.0, help="ATPM grid min (default: 0)")
    p.add_argument("--atpm-max", type=float, default=200.0, help="ATPM grid max (default: 200)")
    p.add_argument("--atpm-step", type=float, default=5.0, help="ATPM grid step (default: 5)")
    p.add_argument(
        "--mode",
        choices=["norm", "rank"],
        default="norm",
        help="How to pick atpm_best per anchor. norm=OD_norm vs mu_norm matching (default). rank=rank-matching (small anchor sets only).",
    )
    p.add_argument("--out", required=True, help="Output directory (e.g., results/campaigns/C6_atpm_calibrated/calibration/)")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def _parse_grid(args) -> np.ndarray:
    if args.atpm_grid:
        vals: list[float] = []
        for x in args.atpm_grid:
            try:
                vals.append(float(str(x).strip()))
            except Exception as e:  # noqa: BLE001
                raise ValueError(f"Invalid --atpm-grid value: {x!r}") from e
        arr = np.array(vals, dtype=float)
        arr = np.unique(arr)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            raise ValueError("--atpm-grid produced empty grid")
        return np.sort(arr)

    vmin = float(args.atpm_min)
    vmax = float(args.atpm_max)
    step = float(args.atpm_step)
    if step <= 0:
        raise ValueError("--atpm-step must be > 0")
    if vmax < vmin:
        raise ValueError("--atpm-max must be >= --atpm-min")
    # inclusive grid
    n = int(np.floor((vmax - vmin) / step)) + 1
    grid = vmin + step * np.arange(n, dtype=float)
    # ensure vmax included if within tolerance
    if grid.size == 0 or grid[-1] < vmax - 1e-9:
        grid = np.append(grid, vmax)
    return np.unique(np.round(grid, 10))


def _apply_fixed_atpm(model, atpm_value: float, atpm_rid: str = "ATPM") -> None:
    present = set(r.id for r in model.reactions)
    if atpm_rid not in present:
        raise ValueError(f"ATPM reaction id not found in model: {atpm_rid}")
    rxn = model.reactions.get_by_id(atpm_rid)
    rxn.upper_bound = 1000.0
    rxn.lower_bound = float(atpm_value)
    rxn.upper_bound = float(atpm_value)


def _scan_anchor(
    *,
    model_path: str,
    medium_cfg: dict,
    anchor_row: dict,
    grid: np.ndarray,
) -> list[ScanRow]:
    out: list[ScanRow] = []
    cid = str(anchor_row["condition_id"])
    for atpm in grid.tolist():
        try:
            model = load_sbml_model(model_path)
            apply_condition_to_model(model, anchor_row, medium_cfg)
            _apply_fixed_atpm(model, float(atpm), atpm_rid="ATPM")
            sol = model.optimize()
            if sol.status != "optimal":
                out.append(ScanRow(anchor_id=cid, atpm=float(atpm), objective_value=float("nan"), status=str(sol.status)))
                continue
            out.append(ScanRow(anchor_id=cid, atpm=float(atpm), objective_value=float(sol.objective_value), status="optimal"))
        except Exception as e:  # noqa: BLE001
            out.append(ScanRow(anchor_id=cid, atpm=float(atpm), objective_value=float("nan"), status=f"error:{type(e).__name__}"))
    return out


def _pick_best_norm(scan_df: pd.DataFrame, anchors_df: pd.DataFrame) -> pd.DataFrame:
    """
    For each anchor:
      OD_norm = (OD - min)/(max-min) across anchors
      mu_norm(atpm) = objective(atpm) / objective(atpm=0)  [if objective0>0]
      pick atpm minimizing |mu_norm - OD_norm| among feasible rows.
    """
    a = anchors_df[["condition_id", "acetate_mM", "measured_OD"]].copy()
    a["measured_OD"] = pd.to_numeric(a["measured_OD"], errors="coerce")
    if a["measured_OD"].isna().any():
        missing = a.loc[a["measured_OD"].isna(), "condition_id"].astype(str).tolist()
        raise ValueError(f"Anchor measured_OD is missing for: {missing}")

    od = a["measured_OD"].astype(float).values
    od_min = float(np.min(od))
    od_max = float(np.max(od))
    if abs(od_max - od_min) < 1e-12:
        # Degenerate: all ODs identical; pick smallest feasible atpm for all anchors
        out_rows: list[dict] = []
        for cid in a["condition_id"].astype(str).tolist():
            sub = scan_df[(scan_df["anchor_id"] == cid) & (scan_df["status"] == "optimal")].copy()
            if sub.empty:
                out_rows.append({"anchor_id": cid, "acetate_mM": float(a.loc[a["condition_id"] == cid, "acetate_mM"].iloc[0]), "atpm_best": float("nan")})
                continue
            best = float(sub.sort_values("atpm", ascending=True)["atpm"].iloc[0])
            out_rows.append({"anchor_id": cid, "acetate_mM": float(a.loc[a["condition_id"] == cid, "acetate_mM"].iloc[0]), "atpm_best": best})
        return pd.DataFrame(out_rows)

    a["OD_norm"] = (a["measured_OD"].astype(float) - od_min) / (od_max - od_min)

    # objective at atpm=0 per anchor (must be optimal)
    obj0 = (
        scan_df[(scan_df["status"] == "optimal") & (np.isclose(scan_df["atpm"].astype(float), 0.0))]
        .groupby("anchor_id")["objective_value"]
        .first()
    )

    out_rows = []
    for _, row in a.iterrows():
        cid = str(row["condition_id"])
        od_norm = float(row["OD_norm"])

        sub = scan_df[(scan_df["anchor_id"] == cid) & (scan_df["status"] == "optimal")].copy()
        if sub.empty:
            out_rows.append({"anchor_id": cid, "acetate_mM": float(row["acetate_mM"]), "atpm_best": float("nan")})
            continue

        obj0_v = float(obj0.get(cid, np.nan))
        if not np.isfinite(obj0_v) or obj0_v <= 0:
            # fallback: normalize by max feasible objective
            obj0_v = float(np.nanmax(pd.to_numeric(sub["objective_value"], errors="coerce").values))

        sub["mu_norm"] = pd.to_numeric(sub["objective_value"], errors="coerce") / float(obj0_v)
        sub["score"] = (sub["mu_norm"] - od_norm).abs()
        sub = sub.dropna(subset=["score"]).sort_values(["score", "atpm"], ascending=[True, True])
        best = float(sub["atpm"].iloc[0]) if not sub.empty else float("nan")
        out_rows.append({"anchor_id": cid, "acetate_mM": float(row["acetate_mM"]), "atpm_best": best})

    return pd.DataFrame(out_rows)


def _pick_best_rank(scan_df: pd.DataFrame, anchors_df: pd.DataFrame) -> pd.DataFrame:
    """
    Brute-force rank matching for small anchor sets (<=4 recommended).
    Find atpm choices that maximize agreement between rank(objective) and rank(measured_OD).
    """
    a = anchors_df[["condition_id", "acetate_mM", "measured_OD"]].copy()
    a["measured_OD"] = pd.to_numeric(a["measured_OD"], errors="coerce")
    if a["measured_OD"].isna().any():
        missing = a.loc[a["measured_OD"].isna(), "condition_id"].astype(str).tolist()
        raise ValueError(f"Anchor measured_OD is missing for: {missing}")

    ids = a["condition_id"].astype(str).tolist()
    if len(ids) < 2:
        raise ValueError("--mode rank requires >=2 anchors")
    if len(ids) > 4:
        raise ValueError("--mode rank supports up to 4 anchors (avoid combinatorial explosion)")

    # desired order by measured OD
    desired = a.sort_values("measured_OD", ascending=True)["condition_id"].astype(str).tolist()
    desired_rank = {cid: i for i, cid in enumerate(desired)}

    choices: list[list[tuple[float, float]]] = []  # per anchor: list of (atpm, objective)
    for cid in ids:
        sub = scan_df[(scan_df["anchor_id"] == cid) & (scan_df["status"] == "optimal")].copy()
        sub["objective_value"] = pd.to_numeric(sub["objective_value"], errors="coerce")
        sub = sub.dropna(subset=["objective_value"])
        if sub.empty:
            raise ValueError(f"No feasible (optimal) points for anchor: {cid}")
        choices.append(list(zip(sub["atpm"].astype(float).tolist(), sub["objective_value"].astype(float).tolist())))

    # brute force
    best = None
    best_score = -1
    best_tie = None  # minimize sum atpm

    def _agreement_score(order: list[str]) -> int:
        # count pairwise agreements with desired rank
        score = 0
        for i in range(len(order)):
            for j in range(i + 1, len(order)):
                ci, cj = order[i], order[j]
                score += int((desired_rank[ci] < desired_rank[cj]) == (i < j))
        return score

    # iterate nested loops up to 4 anchors
    import itertools

    for combo in itertools.product(*choices):
        # combo is ((atpm,obj), ...)
        obj_map = {cid: combo[i][1] for i, cid in enumerate(ids)}
        # derive order by objective
        order = sorted(ids, key=lambda c: (obj_map[c], c))
        score = _agreement_score(order)
        tie = sum(combo[i][0] for i in range(len(ids)))
        if score > best_score or (score == best_score and (best_tie is None or tie < best_tie)):
            best_score = score
            best_tie = tie
            best = [combo[i][0] for i in range(len(ids))]

    if best is None:
        raise RuntimeError("rank calibration failed unexpectedly")

    out_rows = []
    for cid, atpm in zip(ids, best):
        ac = float(a.loc[a["condition_id"].astype(str) == cid, "acetate_mM"].iloc[0])
        out_rows.append({"anchor_id": cid, "acetate_mM": ac, "atpm_best": float(atpm)})
    return pd.DataFrame(out_rows)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    if not Path(args.model).exists():
        print(f"[ERROR] model not found: {args.model}", file=sys.stderr)
        return 2
    if not Path(args.conditions).exists():
        print(f"[ERROR] conditions not found: {args.conditions}", file=sys.stderr)
        return 2
    if not Path(args.medium).exists():
        print(f"[ERROR] medium yaml not found: {args.medium}", file=sys.stderr)
        return 2

    try:
        conditions_df = load_conditions_csv(args.conditions)
        medium_cfg = load_config(args.medium)
        grid = _parse_grid(args)
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] Failed to load inputs: {e}", file=sys.stderr)
        return 2

    anchors = [str(x) for x in (args.anchor_ids or [])]
    if len(anchors) < 2:
        print("[ERROR] Need at least 2 anchors to fit linear ATPM=f(acetate_mM)", file=sys.stderr)
        return 2

    df = conditions_df.copy()
    df["condition_id"] = df["condition_id"].astype(str)

    anchors_df = df.loc[df["condition_id"].isin(set(anchors))].copy()
    missing = sorted(set(anchors) - set(anchors_df["condition_id"].tolist()))
    if missing:
        print(f"[ERROR] Anchor IDs not found in conditions CSV: {missing}", file=sys.stderr)
        return 2

    # scan
    scan_rows: list[ScanRow] = []
    for _, row in anchors_df.iterrows():
        scan_rows.extend(
            _scan_anchor(
                model_path=str(args.model),
                medium_cfg=medium_cfg,
                anchor_row=row.to_dict(),
                grid=grid,
            )
        )
    scan_df = pd.DataFrame([asdict(r) for r in scan_rows])

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) anchor scan parquet
    scan_path = outdir / "anchor_scan.parquet"
    save_table(scan_df, scan_path, fmt="parquet")

    # 2) pick atpm_best per anchor
    if args.mode == "norm":
        best_df = _pick_best_norm(scan_df, anchors_df)
    else:
        best_df = _pick_best_rank(scan_df, anchors_df)
    best_path = outdir / "anchor_best.csv"
    best_df = best_df.sort_values("acetate_mM", ascending=True)
    best_df.to_csv(best_path, index=False)

    # 3) fit linear atpm = a + b*acetate_mM
    fit_df = best_df.dropna(subset=["atpm_best", "acetate_mM"]).copy()
    if len(fit_df) < 2:
        print("[ERROR] Not enough valid anchors to fit linear model (need >=2).", file=sys.stderr)
        return 2

    x = fit_df["acetate_mM"].astype(float).values
    y = fit_df["atpm_best"].astype(float).values
    # polyfit returns slope, intercept for deg=1 (y = b*x + a)
    b, a = np.polyfit(x, y, 1)

    fit = {
        "fit_type": "linear",
        "mode": str(args.mode),
        "a": float(a),
        "b": float(b),
        "clip_min": 0.0,
        "clip_max": 200.0,
        "anchors_used": [str(x) for x in best_df["anchor_id"].astype(str).tolist()],
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "notes": "ATPM_eff = clip(a + b*acetate_mM, 0, 200)",
    }
    fit_path = outdir / "atpm_fit.json"
    fit_path.write_text(json.dumps(fit, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    # console report (results only; no causal claims)
    print(f"[OK] Wrote: {scan_path}")
    print(f"[OK] Wrote: {best_path}")
    print(f"[OK] Wrote: {fit_path} (a={fit['a']:.4g}, b={fit['b']:.4g})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

