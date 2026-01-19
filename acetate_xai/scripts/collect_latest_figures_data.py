from __future__ import annotations

import argparse
import csv
import json
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class DataPick:
    figure_id: str
    data_id: str
    chosen_path: Path
    source_rule: str


def _mtime_iso(p: Path) -> str:
    ts = p.stat().st_mtime
    return datetime.fromtimestamp(ts, tz=timezone.utc).isoformat()


def _size_bytes(p: Path) -> int:
    return int(p.stat().st_size)


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _copy(src: Path, dst: Path) -> None:
    _ensure_dir(dst.parent)
    shutil.copy2(src, dst)


def _all_results_roots() -> list[Path]:
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


def _glob_all(roots: list[Path], pattern: str) -> list[Path]:
    out: list[Path] = []
    for r in roots:
        out.extend(list(r.glob(pattern)))
    return [p for p in out if p.exists() and p.is_file()]


def _latest_by_mtime(paths: Iterable[Path]) -> Path | None:
    c = [p for p in paths if p.exists() and p.is_file()]
    if not c:
        return None
    c.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return c[0]


def _rel_to_workspace(p: Path) -> str:
    try:
        return str(p.relative_to(Path.cwd())).replace("\\", "/")
    except Exception:
        return str(p).replace("\\", "/")


def _to_markdown_table(rows: list[dict]) -> str:
    if not rows:
        return ""
    cols = list(rows[0].keys())
    srows = [{k: str(r.get(k, "")) for k in cols} for r in rows]
    widths = {c: max(len(c), max(len(rr[c]) for rr in srows)) for c in cols}
    header = "| " + " | ".join(c.ljust(widths[c]) for c in cols) + " |"
    sep = "| " + " | ".join("-" * widths[c] for c in cols) + " |"
    body = ["| " + " | ".join(rr[c].ljust(widths[c]) for c in cols) + " |" for rr in srows]
    return "\n".join([header, sep] + body)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Collect data files backing the selected latest figures.")
    p.add_argument(
        "--figures-manifest",
        default=str(Path("results") / "summary" / "latest_figures_manifest.csv"),
        help="CSV produced by collect_latest_figures.py (default: results/summary/latest_figures_manifest.csv)",
    )
    p.add_argument(
        "--outdir",
        default=str(Path("results") / "figures_final" / "data"),
        help="Data output folder (default: results/figures_final/data)",
    )
    p.add_argument(
        "--summary-dir",
        default=str(Path("results") / "summary"),
        help="Summary dir for manifest outputs (default: results/summary)",
    )
    return p


def _read_figures_manifest(path: Path) -> list[dict]:
    if not path.exists():
        raise FileNotFoundError(f"Figures manifest not found: {path}")
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _safe_add(picks: list[DataPick], *, fig: str, data_id: str, p: Path | None, rule: str) -> None:
    if p is None:
        return
    if not p.exists():
        return
    picks.append(DataPick(fig, data_id, p, rule))


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    roots = _all_results_roots()
    fig_manifest = Path(args.figures_manifest)
    outdir = Path(args.outdir)
    summary_dir = Path(args.summary_dir)
    _ensure_dir(outdir)
    _ensure_dir(summary_dir)

    figs = _read_figures_manifest(fig_manifest)
    by_id = {r["figure_id"]: Path(r["chosen_path"]) for r in figs if "figure_id" in r and "chosen_path" in r}

    picks: list[DataPick] = []

    # --- Fig02A
    fig02a = by_id.get("Fig02A_experiment_anchors_4panel.png")
    if fig02a:
        _safe_add(
            picks,
            fig="Fig02A_experiment_anchors_4panel.png",
            data_id="Fig02A_experiment_anchors_4panel.csv",
            p=fig02a.with_suffix(".csv"),
            rule="Same basename .csv in figures_draft",
        )

    # --- Fig02B / Fig02C summaries
    fig02b = by_id.get("Fig02B_regime_vs_OD.png")
    if fig02b:
        # Prefer v2 summary CSV if present
        cand = [
            fig02b.parent / "Fig02B_regime_vs_OD_summary.csv",
            fig02b.parent / "Fig02B_regime_vs_OD.csv",
        ]
        _safe_add(picks, fig="Fig02B_regime_vs_OD.png", data_id=cand[0].name, p=_latest_by_mtime(cand), rule="Prefer summary csv, else csv")

    # Severity-by-regime (if exists in figures_draft)
    for cand_name in ["Fig02C_regime_vs_severity_summary.csv", "Fig02C_regime_vs_severity.csv"]:
        p = Path("results") / "figures_draft" / cand_name
        if p.exists():
            _safe_add(picks, fig="Fig02C_regime_vs_severity.png", data_id=cand_name, p=p, rule="Figures_draft summary/csv if exists")
            break

    # --- Fig02E anchor-like plots: collect all known CSV companions if present
    for nm in [
        "Fig02E_best_run_recreated_data.csv",
        "Fig02E_best_run_recreated_meta.json",
        "Fig02E_best_run_candidates.csv",
        "Fig02E_anchor_scatter_data.csv",
        "Fig02E_before_objective_vs_OD_data.csv",
        "Fig02E_after_calibObjective_vs_OD_data.csv",
        "Fig02E_before_after_join.csv",
        "Fig02E_correlation_summary.csv",
        "Fig02E_acetate_train_test_metrics.txt",
    ]:
        p = Path("results") / "figures_draft" / nm
        if p.exists():
            _safe_add(picks, fig="Fig02E_anchor.png", data_id=nm, p=p, rule="Known Fig02E companion outputs under results/figures_draft")

    # --- Fig02 regime map: back it with regime_map.csv if available
    fig02map = by_id.get("Fig02_regime_map.png")
    if fig02map:
        # Most robust: use any results_root/platform_summary/regime_map.csv if present
        candidates = []
        for r in roots:
            candidates.append(r / "platform_summary" / "regime_map.csv")
        _safe_add(
            picks,
            fig="Fig02_regime_map.png",
            data_id="regime_map.csv",
            p=_latest_by_mtime(candidates),
            rule="Latest platform_summary/regime_map.csv across results roots",
        )

    # --- Fig03 SHAP regime: pick support files from same folder as the chosen beeswarm if possible
    fig03_bee = by_id.get("Fig03_shap_regime_beeswarm.png")
    if fig03_bee:
        d = fig03_bee.parent
        for nm in ["shap_importance.csv", "shap_by_class_importance.csv", "shap_richness_report.json", "confusion_matrix.csv", "classification_report.txt", "label_mapping.csv"]:
            _safe_add(picks, fig="Fig03_shap_regime_beeswarm.png", data_id=nm, p=d / nm, rule="Same xai_xgb_regime folder support files")

    # --- Fig04 SHAP severity: pick support files from same folder as chosen beeswarm
    fig04_bee = by_id.get("Fig04_shap_severity_beeswarm.png")
    if fig04_bee:
        d = fig04_bee.parent
        for nm in ["shap_importance.csv", "shap_richness_report.json", "regression_metrics.csv", "model.json"]:
            _safe_add(picks, fig="Fig04_shap_severity_beeswarm.png", data_id=nm, p=d / nm, rule="Same xai_xgb_severity folder support files")

    # --- Fig05 causal DAG: edges + stability + metadata from latest causal_dag folder
    candidates = []
    for r in roots:
        candidates.append(r / "causal_dag" / "dag_edges.csv")
    dag_edges = _latest_by_mtime(candidates)
    if dag_edges:
        d = dag_edges.parent
        for nm in ["dag_edges.csv", "edge_stability.csv", "dag_metadata.json"]:
            _safe_add(picks, fig="Fig05_causal_dag.png", data_id=nm, p=d / nm, rule="Latest causal_dag support files")

    # --- Fig06 top edges bar: if there is a matching CSV near figures_draft, collect it
    for nm in ["Fig06_top_edges_bar.csv", "top_edges.csv", "edge_stability.csv"]:
        p = Path("results") / "figures_draft" / nm
        if p.exists():
            _safe_add(picks, fig="Fig06_top_edges_bar.png", data_id=nm, p=p, rule="Nearby csv for Fig06 if present")
            break

    # Copy picks to outdir with stable naming: <figureStem>__<data_id>
    copied: list[dict] = []
    for rec in picks:
        src = rec.chosen_path
        figure_stem = Path(rec.figure_id).stem
        dst_name = f"{figure_stem}__{Path(rec.data_id).name}"
        dst = outdir / dst_name
        _copy(src, dst)
        copied.append(
            {
                "figure_id": rec.figure_id,
                "data_id": rec.data_id,
                "chosen_path": _rel_to_workspace(src),
                "copied_to": _rel_to_workspace(dst),
                "source_rule": rec.source_rule,
                "mtime_iso": _mtime_iso(src),
                "size_bytes": _size_bytes(src),
            }
        )

    # Write manifests
    man_csv = summary_dir / "latest_figures_data_manifest.csv"
    man_md = summary_dir / "latest_figures_data_manifest.md"
    copied.sort(key=lambda r: (r["figure_id"], r["data_id"]))
    with man_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(copied[0].keys()) if copied else ["figure_id"])
        w.writeheader()
        for r in copied:
            w.writerow(r)

    md_lines = [
        "## Latest figures data manifest (auto-collected)",
        "",
        f"- figures manifest: `{_rel_to_workspace(fig_manifest)}`",
        f"- data output dir: `{_rel_to_workspace(outdir)}`",
        "",
    ]
    if copied:
        md_lines.append(_to_markdown_table(copied))
    else:
        md_lines.append("(No data files collected.)")
    man_md.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    # Quick inventory text
    list_txt = summary_dir / "figures_final_data_list.txt"
    lines = []
    for p in sorted(outdir.glob("*")):
        if p.is_file():
            lines.append(f"{p.as_posix()}\t{p.stat().st_size}")
    list_txt.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # Optional: also emit a machine-readable index JSON
    idx_json = summary_dir / "latest_figures_data_manifest.json"
    idx_json.write_text(json.dumps(copied, indent=2), encoding="utf-8")

    print(f"[OK] Copied {len(copied)} data files to: {outdir}")
    print(f"[OK] Wrote: {man_csv}")
    print(f"[OK] Wrote: {man_md}")
    print(f"[OK] Wrote: {list_txt}")
    print(f"[OK] Wrote: {idx_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

