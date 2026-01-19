from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class CandidateScore:
    features_path: Path
    run_id: str
    n_points: int
    rho_overall: float | None
    rho_acetate: float | None
    rho_mean_setwise: float | None


def _utc_now_iso() -> str:
    return datetime.now(tz=timezone.utc).isoformat()


def _is_random_lhs_path(p: Path) -> bool:
    s = str(p).lower()
    return "c_random_lhs" in s or "lhs" in s


def _iter_results_roots() -> list[Path]:
    roots = [Path("results")]
    for p in [
        Path("extracted") / "results_all_v3" / "results",
        Path("extracted") / "xai_package_v2" / "results",
        Path("results_all_v2_extracted") / "results",
    ]:
        if p.exists() and p.is_dir():
            roots.append(p)
    out: list[Path] = []
    seen: set[str] = set()
    for r in roots:
        rp = str(r.resolve())
        if rp not in seen:
            out.append(r)
            seen.add(rp)
    return out


def _spearman_rho(x: np.ndarray, y: np.ndarray) -> float | None:
    # small, dependency-free Spearman (ties -> average ranks) using pandas
    if len(x) < 3:
        return None
    xs = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    ys = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    if np.nanstd(xs) == 0 or np.nanstd(ys) == 0:
        return None
    return float(np.corrcoef(xs, ys)[0, 1])


def _score_candidate(df: pd.DataFrame, *, run_id: str) -> CandidateScore | None:
    # Normalize expected columns
    if "objective_value" not in df.columns or "measured_OD" not in df.columns:
        return None
    dd = df.copy()
    dd["objective_value"] = pd.to_numeric(dd["objective_value"], errors="coerce")
    dd["measured_OD"] = pd.to_numeric(dd["measured_OD"], errors="coerce")
    if "set_name" in dd.columns:
        dd["set_name"] = dd["set_name"].astype(str)
    else:
        dd["set_name"] = "unknown"

    dd = dd.dropna(subset=["objective_value", "measured_OD"]).copy()
    if "condition_id" not in dd.columns:
        dd["condition_id"] = np.arange(len(dd))

    n = len(dd)
    if n < 8:
        return None

    rho_all = _spearman_rho(dd["objective_value"].to_numpy(), dd["measured_OD"].to_numpy())

    # acetate subset score if present
    rho_ac = None
    ac = dd[dd["set_name"] == "acetate_gradient"].copy()
    if len(ac) >= 3:
        rho_ac = _spearman_rho(ac["objective_value"].to_numpy(), ac["measured_OD"].to_numpy())

    # mean of per-set rhos (ignore sets with <3 or undefined)
    rhos = []
    for set_name, g in dd.groupby("set_name"):
        if len(g) < 3:
            continue
        r = _spearman_rho(g["objective_value"].to_numpy(), g["measured_OD"].to_numpy())
        if r is not None and np.isfinite(r):
            rhos.append(float(r))
    rho_setwise = float(np.mean(rhos)) if rhos else None

    return CandidateScore(
        features_path=Path(""),
        run_id=run_id,
        n_points=n,
        rho_overall=rho_all,
        rho_acetate=rho_ac,
        rho_mean_setwise=rho_setwise,
    )


def _load_features_with_objective(features_path: Path) -> pd.DataFrame:
    """
    Our pipeline sometimes stores objective_value in sibling regime_fba.parquet
    (not inside features.parquet). This mirrors the behavior used in figure scripts.
    """
    df = pd.read_parquet(features_path)
    if "objective_value" in df.columns:
        return df

    if "condition_id" not in df.columns:
        # can't merge without condition_id
        return df

    run_dir = features_path.parent
    candidates = [
        run_dir / "regime_fba.parquet",
        run_dir / "xai" / "regime_fba.parquet",
        run_dir / "xai" / "regime_fba_rxnfix.parquet",
        run_dir / "xai" / "regime_fba.parquet",
    ]
    reg = None
    for p in candidates:
        if p.exists():
            reg = p
            break
    if reg is None:
        return df

    try:
        reg_df = pd.read_parquet(reg)
    except Exception:
        return df

    if "objective_value" not in reg_df.columns or "condition_id" not in reg_df.columns:
        return df

    merged = df.merge(reg_df[["condition_id", "objective_value"]], on="condition_id", how="left")
    return merged


def _run_id_from_features_path(results_root: Path, features_path: Path) -> str:
    try:
        rel = features_path.relative_to(results_root / "campaigns")
        parts = rel.parts
        if len(parts) >= 2:
            return f"{parts[0]}__{parts[1]}"
        return rel.as_posix().replace("/", "__")
    except Exception:
        return features_path.as_posix().replace("/", "__")


def _find_candidate_features(results_root: Path) -> list[Path]:
    base = results_root / "campaigns"
    if not base.exists():
        return []
    feats = list(base.rglob("features.parquet"))
    # Exclude random LHS (validation fig must be single 22-condition run)
    feats = [p for p in feats if not _is_random_lhs_path(p)]
    return feats


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Pick the single best 22-condition run and recreate Fig02E scatter.")
    p.add_argument(
        "--results-root",
        default=None,
        help="Results root that contains campaigns/. If omitted, searches results/ then extracted results roots.",
    )
    p.add_argument(
        "--score",
        choices=["overall", "acetate", "mean_setwise"],
        default="acetate",
        help="Best-run selection criterion (default: acetate).",
    )
    p.add_argument(
        "--outdir",
        default=str(Path("results") / "figures_draft"),
        help="Output directory (default: results/figures_draft).",
    )
    p.add_argument("--min-n", type=int, default=8, help="Minimum non-NaN points required (default: 8).")
    p.add_argument("--seed", type=int, default=42, help="Seed for any randomness (default: 42).")
    p.add_argument("--color-by-set", action="store_true", help="Color points by set_name.")
    return p


def _pick_best(scores: list[CandidateScore], criterion: str) -> CandidateScore:
    def key(s: CandidateScore) -> float:
        v = None
        if criterion == "overall":
            v = s.rho_overall
        elif criterion == "acetate":
            v = s.rho_acetate
        else:
            v = s.rho_mean_setwise
        if v is None or not np.isfinite(v):
            return -np.inf
        return float(v)

    scores_sorted = sorted(scores, key=key, reverse=True)
    if not scores_sorted or key(scores_sorted[0]) == -np.inf:
        raise RuntimeError(f"No valid candidate had a finite score for criterion={criterion}.")
    return scores_sorted[0]


def _plot_scatter(
    df: pd.DataFrame,
    *,
    out_png: Path,
    title: str,
    rho: float | None,
    color_by_set: bool,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_png.parent.mkdir(parents=True, exist_ok=True)

    dd = df.copy()
    dd = dd.dropna(subset=["objective_value", "measured_OD"]).copy()

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.grid(True, alpha=0.25)

    if color_by_set and "set_name" in dd.columns:
        for set_name, g in dd.groupby("set_name"):
            ax.scatter(
                g["objective_value"].to_numpy(),
                g["measured_OD"].to_numpy(),
                s=55,
                alpha=0.85,
                label=str(set_name),
            )
        ax.legend(title="set_name", fontsize=8, title_fontsize=9, loc="best", frameon=True)
    else:
        ax.scatter(dd["objective_value"].to_numpy(), dd["measured_OD"].to_numpy(), s=55, alpha=0.85)

    ax.set_xlabel("Predicted growth (FBA objective)")
    ax.set_ylabel("Measured OD600 at 32 h")
    ax.set_title(title)
    if rho is not None and np.isfinite(rho):
        ax.text(0.98, 0.98, f"Spearman ρ = {rho:.2f}", ha="right", va="top", transform=ax.transAxes, fontsize=11, fontweight="bold")

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    rng = np.random.default_rng(args.seed)
    _ = rng  # reserved for future (e.g., subsampling)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    roots = [Path(args.results_root)] if args.results_root else _iter_results_roots()

    all_scores: list[CandidateScore] = []
    best_inputs: dict[str, str] = {}

    for results_root in roots:
        feats = _find_candidate_features(results_root)
        for fp in feats:
            try:
                df = _load_features_with_objective(fp)
            except Exception:
                continue
            run_id = _run_id_from_features_path(results_root, fp)
            sc = _score_candidate(df, run_id=run_id)
            if sc is None or sc.n_points < int(args.min_n):
                continue
            sc = CandidateScore(
                features_path=fp,
                run_id=sc.run_id,
                n_points=sc.n_points,
                rho_overall=sc.rho_overall,
                rho_acetate=sc.rho_acetate,
                rho_mean_setwise=sc.rho_mean_setwise,
            )
            all_scores.append(sc)
            best_inputs[str(fp)] = str(results_root)

    if not all_scores:
        raise RuntimeError("No candidate 22-condition runs found with objective_value & measured_OD.")

    best = _pick_best(all_scores, args.score)

    # Load best df and write data CSV
    df_best = _load_features_with_objective(best.features_path)
    keep = ["condition_id", "set_name", "objective_value", "measured_OD"]
    for c in keep:
        if c not in df_best.columns:
            if c == "set_name":
                df_best["set_name"] = "unknown"
            elif c == "condition_id":
                df_best["condition_id"] = np.arange(len(df_best))
            else:
                raise RuntimeError(
                    f"Best run missing required column: {c}. "
                    f"(features={best.features_path}) "
                    "If objective_value is missing, ensure regime_fba.parquet exists next to features.parquet."
                )
    out_csv = outdir / "Fig02E_best_run_recreated_data.csv"
    df_best[keep].to_csv(out_csv, index=False)

    # Compute rho for annotation (based on criterion)
    rho_annot = best.rho_acetate if args.score == "acetate" else (best.rho_mean_setwise if args.score == "mean_setwise" else best.rho_overall)

    out_png = outdir / "Fig02E_best_run_recreated.png"
    title = f"Model validation (best run; score={args.score})"
    _plot_scatter(df_best[keep], out_png=out_png, title=title, rho=rho_annot, color_by_set=bool(args.color_by_set))

    # Save ranking table
    rank_csv = outdir / "Fig02E_best_run_candidates.csv"
    with rank_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["run_id", "features_path", "n_points", "rho_overall", "rho_acetate", "rho_mean_setwise"])
        for s in sorted(all_scores, key=lambda x: (x.rho_acetate or -9e9), reverse=True)[:200]:
            w.writerow([s.run_id, s.features_path.as_posix(), s.n_points, s.rho_overall, s.rho_acetate, s.rho_mean_setwise])

    meta = {
        "timestamp_utc": _utc_now_iso(),
        "score": args.score,
        "min_n": int(args.min_n),
        "picked": {
            "run_id": best.run_id,
            "features_path": best.features_path.as_posix(),
            "n_points": best.n_points,
            "rho_overall": best.rho_overall,
            "rho_acetate": best.rho_acetate,
            "rho_mean_setwise": best.rho_mean_setwise,
        },
        "out_png": out_png.as_posix(),
        "out_csv": out_csv.as_posix(),
    }
    (outdir / "Fig02E_best_run_recreated_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"[OK] Best run: {best.run_id}")
    print(f"[OK] Wrote: {out_png}")
    print(f"[OK] Wrote: {out_csv}")
    print(f"[OK] Wrote: {rank_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

