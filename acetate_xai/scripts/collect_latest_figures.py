from __future__ import annotations

import argparse
import csv
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class Picked:
    figure_id: str
    role: str
    chosen_path: Path
    source_rule: str


def _mtime_iso(p: Path) -> str:
    ts = p.stat().st_mtime
    return datetime.fromtimestamp(ts, tz=timezone.utc).isoformat()


def _size_bytes(p: Path) -> int:
    return int(p.stat().st_size)


def _exists(p: Path) -> bool:
    return p.exists()


def _first_existing(paths: Iterable[Path]) -> Path | None:
    for p in paths:
        if p.exists():
            return p
    return None


def _all_results_roots() -> list[Path]:
    """
    Search order:
      1) local results/
      2) extracted results (results_all_v3, xai_package_v2) if present
    """
    roots = [Path("results")]
    for p in [
        Path("extracted") / "results_all_v3" / "results",
        Path("extracted") / "xai_package_v2" / "results",
        Path("results_all_v2_extracted") / "results",
    ]:
        if p.exists() and p.is_dir():
            roots.append(p)
    # de-dup
    out = []
    seen = set()
    for r in roots:
        rp = str(r.resolve())
        if rp not in seen:
            out.append(r)
            seen.add(rp)
    return out


def _latest_by_mtime(candidates: Iterable[Path]) -> Path | None:
    c = [p for p in candidates if p.exists() and p.is_file()]
    if not c:
        return None
    c.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return c[0]


def _glob_all(roots: list[Path], pattern: str) -> list[Path]:
    out: list[Path] = []
    for r in roots:
        out.extend(list(r.glob(pattern)))
    # keep only files
    out = [p for p in out if p.exists() and p.is_file()]
    return out


def _rel_to_workspace(p: Path) -> str:
    # Prefer workspace-relative display
    try:
        return str(p.relative_to(Path.cwd())).replace("\\", "/")
    except Exception:
        return str(p).replace("\\", "/")


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _copy(src: Path, dst: Path) -> None:
    _ensure_dir(dst.parent)
    shutil.copy2(src, dst)

def _to_markdown_table(rows: list[dict]) -> str:
    if not rows:
        return ""
    cols = list(rows[0].keys())
    # stringify
    srows = []
    for r in rows:
        srows.append({k: str(r.get(k, "")) for k in cols})
    widths = {c: max(len(c), max(len(rr[c]) for rr in srows)) for c in cols}
    header = "| " + " | ".join(c.ljust(widths[c]) for c in cols) + " |"
    sep = "| " + " | ".join("-" * widths[c] for c in cols) + " |"
    body = ["| " + " | ".join(rr[c].ljust(widths[c]) for c in cols) + " |" for rr in srows]
    return "\n".join([header, sep] + body)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Collect 'latest confirmed figures' into results/figures_final with manifest.")
    p.add_argument("--outdir", default="results/figures_final", help="Output figure folder (default: results/figures_final)")
    p.add_argument("--summary-dir", default="results/summary", help="Summary folder (default: results/summary)")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    roots = _all_results_roots()
    outdir = Path(args.outdir)
    summary_dir = Path(args.summary_dir)
    _ensure_dir(outdir)
    _ensure_dir(summary_dir)

    picked: list[Picked] = []

    # A) Regime SHAP (fixed folder)
    fixed_rel = Path("campaigns") / "C_random_LHS_FVA" / "xai_xgb_regime_FVAonly_bigshap_v2"
    fixed_dir = _first_existing([r / fixed_rel for r in roots])
    if fixed_dir is None:
        raise FileNotFoundError(f"Could not find fixed Regime SHAP dir under any results root: {fixed_rel}")

    def _must(p: Path) -> Path:
        if not p.exists():
            raise FileNotFoundError(f"Missing required file: {p}")
        return p

    regime_files = {
        "Fig03_shap_regime_beeswarm.png": _must(fixed_dir / "shap_beeswarm.png"),
        "Fig03_shap_regime_bar.png": _must(fixed_dir / "shap_bar.png"),
        "Fig03_shap_regime_class_O2.png": _must(fixed_dir / "shap_beeswarm_class_O2_limited.png"),
        "Fig03_shap_regime_class_N.png": _must(fixed_dir / "shap_beeswarm_class_N_limited.png"),
        "Fig03_shap_regime_class_Ac.png": _must(fixed_dir / "shap_beeswarm_class_Ac_limited.png"),
    }
    # evidence check only
    _must(fixed_dir / "shap_richness_report.json")

    for dst_name, src in regime_files.items():
        picked.append(Picked(figure_id=dst_name, role="Fig03", chosen_path=src, source_rule="A: fixed bigshap_v2 folder"))

    # B) Severity SHAP (latest shap_beeswarm.png among xai_xgb*severity*)
    severity_candidates = _glob_all(roots, "campaigns/**/xai_xgb*severity*/shap_beeswarm.png")
    if severity_candidates:
        # prefer ones under C_random_LHS_FVA if present
        preferred = [p for p in severity_candidates if "C_random_LHS_FVA" in str(p)]
        pool = preferred if preferred else severity_candidates
        best = _latest_by_mtime(pool)
        if best is not None:
            picked.append(Picked("Fig04_shap_severity_beeswarm.png", "Fig04", best, "B: latest severity beeswarm (prefer C_random_LHS_FVA)"))
            bar = best.parent / "shap_bar.png"
            if bar.exists():
                picked.append(Picked("Fig04_shap_severity_bar.png", "Fig04", bar, "B: same folder shap_bar.png"))

    # C) Fig01/Fig02 core figures from results/figures_draft or results/summary
    # We select latest by mtime among matching filenames across local roots (prefer local results/ folders).
    def _pick_latest_named(names: list[str]) -> Path | None:
        cands: list[Path] = []
        for r in roots:
            for base in [r / "figures_draft", r / "summary", r / "summary" / "paper_assets" / "03_figures_existing"]:
                for nm in names:
                    cands.append(base / nm)
        return _latest_by_mtime(cands)

    fig01 = _pick_latest_named(["Fig01_workflow.png"])
    if fig01:
        picked.append(Picked("Fig01_workflow.png", "Fig01", fig01, "C: latest Fig01* under figures_draft/summary"))

    fig02a = _pick_latest_named(["Fig02A_experiment_anchors_4panel.png"])
    if fig02a:
        picked.append(
            Picked("Fig02A_experiment_anchors_4panel.png", "Fig02", fig02a, "C: latest Fig02A experimental anchors 4-panel")
        )

    fig02_map = _pick_latest_named(["Fig02_regime_map.png", "regime_map.png"])
    if fig02_map:
        picked.append(Picked("Fig02_regime_map.png", "Fig02", fig02_map, "C: latest Fig02_regime_map/regime_map"))

    fig02b = _pick_latest_named(["Fig02B_regime_vs_OD_v2.png", "Fig02B_regime_vs_OD.png"])
    if fig02b:
        picked.append(Picked("Fig02B_regime_vs_OD.png", "Fig02", fig02b, "C: latest Fig02B"))

    fig02d = _pick_latest_named(["Fig02D_OD_vs_severity.png"])
    if fig02d:
        picked.append(Picked("Fig02D_OD_vs_severity.png", "Fig02", fig02d, "C: latest Fig02D"))

    # Fig02E anchor: pick latest among known Fig02E candidates (after calib preferred by recency)
    fig02e_anchor = _pick_latest_named(
        [
            "Fig02E_best_run_recreated.png",
            "Fig02E_after_calibObjective_vs_OD.png",
            "Fig02E_best_run_correlation.png",
            "Fig02E_anchor_scatter.png",
            "Fig02E_before_objective_vs_OD.png",
        ]
    )
    if fig02e_anchor:
        picked.append(Picked("Fig02E_anchor.png", "Fig02", fig02e_anchor, "C: latest Fig02E anchor-like plot"))

    # D) Causal DAG (prefer results/causal_dag; else figures_draft; choose latest)
    causal_candidates: list[Path] = []
    for r in roots:
        causal_candidates.extend(
            [
                r / "causal_dag" / "causal_dag.png",
                r / "figures_draft" / "Fig05_causal_dag.png",
            ]
        )
    fig05 = _latest_by_mtime(causal_candidates)
    if fig05:
        picked.append(Picked("Fig05_causal_dag.png", "Fig05", fig05, "D: latest causal DAG"))

    fig06 = _pick_latest_named(["Fig06_top_edges_bar.png"])
    if fig06:
        picked.append(Picked("Fig06_top_edges_bar.png", "Fig06", fig06, "D: latest top edges bar"))

    # Copy to outdir with standardized names (figure_id)
    for rec in picked:
        _copy(rec.chosen_path, outdir / rec.figure_id)

    # Manifest CSV/MD
    man_csv = summary_dir / "latest_figures_manifest.csv"
    man_md = summary_dir / "latest_figures_manifest.md"
    rows = []
    for rec in picked:
        p = rec.chosen_path
        rows.append(
            {
                "figure_id": rec.figure_id,
                "role": rec.role,
                "chosen_path": _rel_to_workspace(p),
                "source_rule": rec.source_rule,
                "mtime_iso": _mtime_iso(p),
                "size_bytes": _size_bytes(p),
            }
        )
    rows.sort(key=lambda r: r["figure_id"])
    with man_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["figure_id"])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    md = []
    md.append("## Latest figures manifest (auto-selected)")
    md.append("")
    md.append("Selection rules were applied exactly as specified in the conversation.")
    md.append("")
    if rows:
        md.append(_to_markdown_table(rows))
    else:
        md.append("(No figures selected.)")
    man_md.write_text("\n".join(md) + "\n", encoding="utf-8")

    # List file sizes after copy
    list_txt = summary_dir / "figures_final_list.txt"
    lines = []
    for p in sorted(outdir.glob("*.png")):
        lines.append(f"{p.as_posix()}\t{p.stat().st_size}")
    list_txt.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[OK] Wrote: {man_csv}")
    print(f"[OK] Wrote: {man_md}")
    print(f"[OK] Copied figures to: {outdir}")
    print(f"[OK] Wrote: {list_txt}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

