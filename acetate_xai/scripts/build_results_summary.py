from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _rel(p: Path, root: Path) -> str:
    try:
        return str(p.relative_to(root)).replace("\\", "/")
    except Exception:
        return str(p).replace("\\", "/")


def _is_dir(p: Path) -> bool:
    return p.exists() and p.is_dir()


def _is_file(p: Path) -> bool:
    return p.exists() and p.is_file()


def _count_parquet_files(p: Path) -> int:
    if not _is_dir(p):
        return 0
    return sum(1 for _ in p.glob("*.parquet"))


def _load_json_safe(p: Path) -> dict[str, Any] | None:
    if not _is_file(p):
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None


def _find_xai_dirs(run_dir: Path) -> list[Path]:
    out: list[Path] = []
    for name in ["xai", "xai_rxnfix", "xai_xgb", "xai_xgb_maintenance", "xai_xgb_regime_FVAonly_bigshap_v2"]:
        p = run_dir / name
        if _is_dir(p):
            out.append(p)
    # also include any xai_xgb* folders
    for p in run_dir.glob("xai_xgb*"):
        if _is_dir(p):
            out.append(p)
    # de-dup
    seen = set()
    uniq = []
    for p in out:
        s = str(p.resolve())
        if s not in seen:
            uniq.append(p)
            seen.add(s)
    return uniq


def _pick_best_shap_dir(run_dir: Path) -> tuple[Path | None, dict[str, Any] | None]:
    """
    Find a directory under run_dir that contains shap_richness_report.json.
    Prefer the one with largest shap_n_used_actual (or n_used fallback).
    """
    candidates = []
    for xd in _find_xai_dirs(run_dir):
        rep = xd / "shap_richness_report.json"
        j = _load_json_safe(rep)
        if j is None:
            continue
        n = j.get("shap_n_used_actual", None)
        if n is None:
            n = j.get("n_used", None)
        try:
            n_int = int(n)
        except Exception:
            n_int = -1
        candidates.append((n_int, xd, j))
    if not candidates:
        return None, None
    candidates.sort(key=lambda t: t[0], reverse=True)
    return candidates[0][1], candidates[0][2]


@dataclass(frozen=True)
class RunRecord:
    results_root_id: str
    results_root: Path
    campaign: str
    run: str
    run_dir: Path
    design_parquet: bool
    regime_labels_or_fba: str  # "regime_labels.parquet" | "regime_fba.parquet" | ""
    features_parquet: bool
    fva_all_parquet: bool
    fva_parts_count: int
    has_xai: bool
    has_xai_xgb: bool
    shap_beeswarm: bool
    shap_bar: bool
    shap_richness_json: bool
    shap_n_used: int | None
    domination_index: float | None


def _scan_campaign_runs(results_dir: Path, *, results_root_id: str) -> list[RunRecord]:
    campaigns_root = results_dir / "campaigns"
    if not _is_dir(campaigns_root):
        return []

    records: list[RunRecord] = []
    for camp_dir in sorted([p for p in campaigns_root.iterdir() if p.is_dir()]):
        campaign = camp_dir.name

        # Determine run directories:
        # - If camp_dir contains features.parquet/regime_labels.parquet directly, treat as a "run".
        # - Also include each subdir under camp_dir as a run if it contains at least one key artifact.
        run_dirs: list[tuple[str, Path]] = []

        def _looks_like_run(d: Path) -> bool:
            keys = [
                d / "features.parquet",
                d / "regime_labels.parquet",
                d / "regime_fba.parquet",
                d / "fva_all.parquet",
                d / "xai",
            ]
            return any(p.exists() for p in keys)

        if _looks_like_run(camp_dir):
            run_dirs.append((camp_dir.name, camp_dir))

        for sub in sorted([p for p in camp_dir.iterdir() if p.is_dir()]):
            if _looks_like_run(sub):
                run_dirs.append((sub.name, sub))

        # de-dup by path
        seen = set()
        uniq: list[tuple[str, Path]] = []
        for run_name, d in run_dirs:
            r = str(d.resolve())
            if r in seen:
                continue
            uniq.append((run_name, d))
            seen.add(r)

        for run_name, run_dir in uniq:
            design = _is_file(run_dir / "design.parquet")
            regime_file = ""
            if _is_file(run_dir / "regime_labels.parquet"):
                regime_file = "regime_labels.parquet"
            elif _is_file(run_dir / "regime_fba.parquet"):
                regime_file = "regime_fba.parquet"

            features = _is_file(run_dir / "features.parquet")
            fva_all = _is_file(run_dir / "fva_all.parquet")
            fva_parts_count = _count_parquet_files(run_dir / "fva_parts")

            xai_dirs = _find_xai_dirs(run_dir)
            has_xai = any(p.name in {"xai", "xai_rxnfix"} for p in xai_dirs)
            has_xai_xgb = any(p.name.startswith("xai_xgb") for p in xai_dirs)

            shap_dir, shap_json = _pick_best_shap_dir(run_dir)
            shap_beeswarm = bool(shap_dir and _is_file(shap_dir / "shap_beeswarm.png"))
            shap_bar = bool(shap_dir and _is_file(shap_dir / "shap_bar.png"))
            shap_richness = bool(shap_dir and _is_file(shap_dir / "shap_richness_report.json"))
            shap_n_used = None
            dom = None
            if shap_json:
                n = shap_json.get("shap_n_used_actual", shap_json.get("n_used", None))
                try:
                    shap_n_used = int(n) if n is not None else None
                except Exception:
                    shap_n_used = None
                di = shap_json.get("domination_index", None)
                try:
                    dom = float(di) if di is not None else None
                except Exception:
                    dom = None

            records.append(
                RunRecord(
                    results_root_id=results_root_id,
                    results_root=results_dir,
                    campaign=campaign,
                    run=run_name,
                    run_dir=run_dir,
                    design_parquet=design,
                    regime_labels_or_fba=regime_file,
                    features_parquet=features,
                    fva_all_parquet=fva_all,
                    fva_parts_count=fva_parts_count,
                    has_xai=has_xai,
                    has_xai_xgb=has_xai_xgb,
                    shap_beeswarm=shap_beeswarm,
                    shap_bar=shap_bar,
                    shap_richness_json=shap_richness,
                    shap_n_used=shap_n_used,
                    domination_index=dom,
                )
            )

    return records


def _write_inventory(records: list[RunRecord], *, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for r in records:
        rows.append(
            {
                "results_root": r.results_root_id,
                "campaign": r.campaign,
                "run": r.run,
                "run_path": _rel(r.run_dir, r.results_root),
                "design_parquet": int(r.design_parquet),
                "regime_file": r.regime_labels_or_fba,
                "features_parquet": int(r.features_parquet),
                "fva_all_parquet": int(r.fva_all_parquet),
                "fva_parts_count": int(r.fva_parts_count),
                "has_xai_dir": int(r.has_xai),
                "has_xai_xgb_dir": int(r.has_xai_xgb),
                "shap_beeswarm": int(r.shap_beeswarm),
                "shap_bar": int(r.shap_bar),
                "shap_richness_report": int(r.shap_richness_json),
                "shap_n_used": r.shap_n_used if r.shap_n_used is not None else "",
                "domination_index": f"{r.domination_index:.4g}" if r.domination_index is not None else "",
            }
        )
    pd.DataFrame(rows).sort_values(["results_root", "campaign", "run"]).to_csv(outdir / "inventory_runs.csv", index=False)

    # Markdown summary
    md = []
    md.append("## Results inventory (campaign/run)")
    md.append("")
    roots = sorted(set(r.results_root_id for r in records))
    md.append(f"- **results_roots_scanned**: {', '.join(f'`{x}`' for x in roots) if roots else '(none)'}")
    md.append(f"- **n_runs_detected**: {len(records)}")
    md.append("")
    if not records:
        md.append("No runs detected under `results/campaigns/`.")
    else:
        df = pd.DataFrame(rows)
        n_shap = int(df["shap_richness_report"].astype(int).sum())
        md.append(f"- **runs_with_shap_richness_report**: {n_shap}")
        md.append("")
        md.append("### Quick table (top 20 by shap_n_used)")
        md.append("")
        df2 = df.copy()
        df2["shap_n_used_num"] = pd.to_numeric(df2["shap_n_used"], errors="coerce")
        df2 = df2.sort_values(["shap_n_used_num"], ascending=False).head(20)
        md.append(
            df2[["results_root", "campaign", "run", "shap_n_used", "domination_index", "run_path"]].to_markdown(
                index=False
            )
        )

    (outdir / "inventory_runs.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def _select_runs(records: list[RunRecord]) -> dict[str, Any]:
    """
    Main run:
      - has shap_richness_report.json
      - maximize shap_n_used
      - domination_index < 0.4

    Comparison run:
      - prefer campaign named C_random_LHS (FBA-only), else any run with design+regime_labels but no shap.
    """
    # main candidates
    cands = []
    for r in records:
        if not r.shap_richness_json:
            continue
        if r.shap_n_used is None:
            continue
        if r.domination_index is None or r.domination_index >= 0.4:
            continue
        cands.append(r)
    cands.sort(key=lambda rr: (rr.shap_n_used or -1), reverse=True)
    main = cands[0] if cands else None

    # comparison: prefer FBA-only LHS run (no SHAP richness) if available
    comp = None
    if main is not None:
        # Try: any C_random_LHS run without SHAP richness and different from main (across all roots)
        lhs_no_shap = [
            r
            for r in records
            if (r.campaign == "C_random_LHS" or r.run == "C_random_LHS")
            and (not r.shap_richness_json)
            and (r.run_dir.resolve() != main.run_dir.resolve())
        ]
        if lhs_no_shap:
            comp = lhs_no_shap[0]
        else:
            # Fallback: different run in same campaign without SHAP richness
            same_campaign = [
                r
                for r in records
                if (r.campaign == main.campaign)
                and (r.run_dir.resolve() != main.run_dir.resolve())
                and (not r.shap_richness_json)
            ]
            if same_campaign:
                comp = same_campaign[0]

    if comp is None:
        for r in records:
            if r.campaign == "C_random_LHS" or r.run == "C_random_LHS":
                # if only one exists, accept it
                comp = r
                break
    if comp is None:
        # fallback: any run that looks FBA-only (design+regime but no fva_all and no shap)
        fba_only = [r for r in records if r.design_parquet and r.regime_labels_or_fba and (not r.fva_all_parquet) and (not r.shap_richness_json)]
        if fba_only:
            comp = fba_only[0]

    out: dict[str, Any] = {}
    out["main"] = {
        "results_root": main.results_root_id if main else "",
        "campaign": main.campaign if main else "",
        "run": main.run if main else "",
        "run_path": _rel(main.run_dir, main.results_root) if main else "",
        "shap_n_used": int(main.shap_n_used) if (main and main.shap_n_used is not None) else None,
        "domination_index": float(main.domination_index) if (main and main.domination_index is not None) else None,
    }
    out["comparison"] = {
        "results_root": comp.results_root_id if comp else "",
        "campaign": comp.campaign if comp else "",
        "run": comp.run if comp else "",
        "run_path": _rel(comp.run_dir, comp.results_root) if comp else "",
    }
    return out


def _write_selected(selected: dict[str, Any], *, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "selected_runs.json").write_text(json.dumps(selected, indent=2), encoding="utf-8")

    md = []
    md.append("## Selected runs (paper main + comparison)")
    md.append("")
    main = selected.get("main", {})
    comp = selected.get("comparison", {})
    md.append("### Main run (final)")
    md.append(f"- **results_root**: `{main.get('results_root','')}`")
    md.append(f"- **campaign**: `{main.get('campaign','')}`")
    md.append(f"- **run**: `{main.get('run','')}`")
    md.append(f"- **shap_n_used**: `{main.get('shap_n_used', '')}`")
    md.append(f"- **domination_index**: `{main.get('domination_index', '')}` (must be < 0.4)")
    md.append("")
    md.append("Selection rationale (auto):")
    md.append("- Has `shap_richness_report.json`")
    md.append("- Maximizes `shap_n_used_actual` among eligible runs")
    md.append("- Satisfies `domination_index < 0.4` (avoid SHAP collapse)")
    md.append("")
    md.append("### Comparison run")
    md.append(f"- **results_root**: `{comp.get('results_root','')}`")
    md.append(f"- **campaign**: `{comp.get('campaign','')}`")
    md.append(f"- **run**: `{comp.get('run','')}`")
    md.append("")
    md.append("Selection rationale (auto):")
    md.append("- Prefer campaign `C_random_LHS` (FBA-only LHS)")
    md.append("- Otherwise pick a design+regime run with no FVA/SHAP")
    if (
        main.get("results_root") == comp.get("results_root")
        and main.get("run_path") == comp.get("run_path")
        and main.get("run_path")
    ):
        md.append("")
        md.append("**Note**: No distinct FBA-only LHS comparison run was detected; comparison currently points to the same run.")
    (outdir / "selected_runs.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def _manifest_row(kind: str, label: str, path: Path, results_dir: Path) -> dict[str, Any]:
    return {
        "kind": kind,
        "label": label,
        "path": _rel(path, results_dir),
        "exists": int(path.exists()),
    }


def _write_manifests(selected: dict[str, Any], *, results_roots: dict[str, Path], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    figs: list[dict[str, Any]] = []
    tabs: list[dict[str, Any]] = []

    # Global figures (if present) – scan each results root
    for rid, rdir in results_roots.items():
        figs_root = rdir / "figures_draft"
        for name, label in [
            ("Fig01_workflow.png", "Workflow"),
            ("Fig02_regime_map.png", "Regime map"),
            ("Fig02A_experiment_anchors_4panel.png", "Experimental anchors (4-panel)"),
            ("Fig02B_regime_vs_OD_v2.png", "OD vs regime (paper v2)"),
            ("Fig02C_regime_vs_severity.png", "Severity vs regime"),
            ("Fig02D_OD_vs_severity.png", "OD vs severity scatter"),
            ("Fig02E_best_run_correlation.png", "Best-run correlation panels"),
            ("Fig03_shap_regime_beeswarm_v2.png", "SHAP beeswarm regime (v2)"),
            ("Fig04_shap_severity_beeswarm_v3.png", "SHAP beeswarm severity (v3)"),
            ("Fig04_shap_severity_beeswarm_v4_jittered.png", "SHAP beeswarm severity (v4 jittered)"),
            ("Fig05_causal_dag.png", "Causal DAG"),
        ]:
            figs.append(
                {
                    "kind": "figure",
                    "label": label,
                    "results_root": rid,
                    "path": _rel(figs_root / name, rdir),
                    "exists": int((figs_root / name).exists()),
                }
            )

    # Main run XAI artifacts
    main = selected.get("main", {})
    main_root_id = str(main.get("results_root", ""))
    main_rel = str(main.get("run_path", ""))
    if main_root_id and main_rel and main_root_id in results_roots:
        root = results_roots[main_root_id]
        run_dir = root / main_rel
        shap_dir, _j = _pick_best_shap_dir(run_dir)
        if shap_dir:
            figs.append({"kind": "figure", "label": "Main run SHAP beeswarm", "results_root": main_root_id, "path": _rel(shap_dir / "shap_beeswarm.png", root), "exists": int((shap_dir / "shap_beeswarm.png").exists())})
            figs.append({"kind": "figure", "label": "Main run SHAP bar", "results_root": main_root_id, "path": _rel(shap_dir / "shap_bar.png", root), "exists": int((shap_dir / "shap_bar.png").exists())})
            for c in ["Ac_limited", "N_limited", "O2_limited"]:
                figs.append({"kind": "figure", "label": f"Main run SHAP beeswarm class {c}", "results_root": main_root_id, "path": _rel(shap_dir / f"shap_beeswarm_class_{c}.png", root), "exists": int((shap_dir / f"shap_beeswarm_class_{c}.png").exists())})
            tabs.append({"kind": "table", "label": "Main run shap_importance.csv", "results_root": main_root_id, "path": _rel(shap_dir / "shap_importance.csv", root), "exists": int((shap_dir / "shap_importance.csv").exists())})
            tabs.append({"kind": "table", "label": "Main run shap_by_class_importance.csv", "results_root": main_root_id, "path": _rel(shap_dir / "shap_by_class_importance.csv", root), "exists": int((shap_dir / "shap_by_class_importance.csv").exists())})
            tabs.append({"kind": "table", "label": "Main run shap_richness_report.json", "results_root": main_root_id, "path": _rel(shap_dir / "shap_richness_report.json", root), "exists": int((shap_dir / "shap_richness_report.json").exists())})

    pd.DataFrame(figs).to_csv(outdir / "figure_manifest.csv", index=False)
    pd.DataFrame(tabs).to_csv(outdir / "table_manifest.csv", index=False)


def _tar_gz_create(archive_path: Path, root: Path, rel_paths: Iterable[str]) -> None:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, "w:gz") as tf:
        for rp in rel_paths:
            src = root / rp
            if not src.exists():
                continue
            tf.add(src, arcname=rp)


def _collect_si_paths(selected: dict[str, Any], *, results_dir: Path) -> dict[str, list[str]]:
    """
    Build relative path lists for SI_min / SI_recommended / SI_full.
    """
    main_path = selected.get("main", {}).get("run_path", "")
    if not main_path:
        return {"SI_min": [], "SI_recommended": [], "SI_full": []}
    run_dir = results_dir / str(main_path)

    shap_dir, _j = _pick_best_shap_dir(run_dir)
    if shap_dir is None:
        shap_dir = run_dir

    # Minimal: xai artifacts + shap reports
    si_min: list[str] = []
    for f in [
        "model.json",
        "run_metadata.json",
        "confusion_matrix.csv",
        "classification_report.txt",
        "label_mapping.csv",
        "shap_beeswarm.png",
        "shap_bar.png",
        "shap_importance.csv",
        "shap_by_class_importance.csv",
        "shap_richness_report.json",
        "shap_richness_report.csv",
    ]:
        si_min.append(_rel(shap_dir / f, results_dir))
    for c in ["Ac_limited", "N_limited", "O2_limited"]:
        si_min.append(_rel(shap_dir / f"shap_beeswarm_class_{c}.png", results_dir))

    # Recommended: include design/labels/features/fva_all if present
    si_rec = list(si_min)
    for f in ["design.parquet", "regime_labels.parquet", "regime_fba.parquet", "features.parquet", "fva_all.parquet"]:
        p = run_dir / f
        if p.exists():
            si_rec.append(_rel(p, results_dir))

    # Full: include fva_parts if present (all parquet under directory)
    si_full = list(si_rec)
    parts_dir = run_dir / "fva_parts"
    if parts_dir.exists() and parts_dir.is_dir():
        for f in sorted(parts_dir.glob("*.parquet")):
            si_full.append(_rel(f, results_dir))

    # De-dup and keep only existing
    def _dedup(xs: list[str]) -> list[str]:
        seen = set()
        out = []
        for x in xs:
            if x in seen:
                continue
            seen.add(x)
            if (results_dir / x).exists():
                out.append(x)
        return out

    return {"SI_min": _dedup(si_min), "SI_recommended": _dedup(si_rec), "SI_full": _dedup(si_full)}


def _write_si_packages(selected: dict[str, Any], *, results_dir: Path, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    bundles = _collect_si_paths(selected, results_dir=results_dir)
    checks: list[str] = []
    for name, rels in bundles.items():
        archive = outdir / f"{name}.tar.gz"
        _tar_gz_create(archive, results_dir, rels)
        checks.append(f"{archive.name}\t{_sha256(archive)}")
    (outdir / "SI_checksums.txt").write_text("\n".join(checks) + "\n", encoding="utf-8")


def _write_readme(selected: dict[str, Any], *, results_dir: Path, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    md = []
    md.append("## Results summary (paper + SI packaging)")
    md.append("")
    md.append("This folder was auto-generated. **Do not edit `results/` in-place**; only `results/summary/` is created/updated.")
    md.append("")
    md.append("### Key outputs produced here")
    md.append("- `inventory_runs.csv` / `inventory_runs.md`: campaign/run inventory")
    md.append("- `selected_runs.json` / `selected_runs.md`: recommended main vs comparison run")
    md.append("- `figure_manifest.csv`: figure file list (exists flags)")
    md.append("- `table_manifest.csv`: table file list (exists flags)")
    md.append("- `SI_min.tar.gz`, `SI_recommended.tar.gz`, `SI_full.tar.gz`: SI bundles")
    md.append("- `SI_checksums.txt`: sha256 for SI bundles")
    md.append("")
    main = selected.get("main", {})
    comp = selected.get("comparison", {})
    md.append("### Selected runs")
    md.append(f"- **Main**: `{main.get('campaign','')}/{main.get('run','')}`")
    md.append(f"  - shap_n_used={main.get('shap_n_used','')}, domination_index={main.get('domination_index','')}")
    md.append(f"- **Comparison**: `{comp.get('campaign','')}/{comp.get('run','')}`")
    md.append("")
    md.append("### SI bundle tiers")
    md.append("- **SI_min**: SHAP/XAI results (beeswarm/bar/importance/richness + model + reports)")
    md.append("- **SI_recommended**: SI_min + design/labels/features/fva_all (if available)")
    md.append("- **SI_full**: SI_recommended + fva_parts (if available)")
    md.append("")
    md.append("### Notes")
    md.append("- Inventory scans `results/campaigns/**` and detects both `<campaign>/<run>` and `<campaign>/` as runs.")
    md.append("- Representative run selection prefers runs with `shap_richness_report.json` and low SHAP domination.")
    (outdir / "README_results.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Build results/summary inventory, run selection, manifests, and SI bundles.")
    p.add_argument(
        "--results-dir",
        nargs="+",
        default=["results"],
        help="One or more results directories to scan (default: results).",
    )
    p.add_argument("--outdir", default="results/summary", help="Output directory (default: results/summary)")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    results_roots: dict[str, Path] = {}
    for p in args.results_dir:
        rp = Path(p)
        # Make a stable, human-friendly id:
        # - for extracted/.../results, use the parent folder name (e.g., results_all_v3, xai_package_v2)
        rid = rp.parent.name if rp.name.lower() == "results" else rp.name
        # keep id stable but unique
        if rid in results_roots:
            rid = f"{rid}__{len(results_roots)+1}"
        results_roots[rid] = rp

    records: list[RunRecord] = []
    for rid, rdir in results_roots.items():
        records.extend(_scan_campaign_runs(rdir, results_root_id=rid))

    _write_inventory(records, outdir=outdir)
    selected = _select_runs(records)
    _write_selected(selected, outdir=outdir)
    _write_manifests(selected, results_roots=results_roots, outdir=outdir)

    # SI packaging uses the main selected run's root
    main_root_id = str(selected.get("main", {}).get("results_root", ""))
    if main_root_id and main_root_id in results_roots:
        root = results_roots[main_root_id]
        _write_si_packages(selected, results_dir=root, outdir=outdir)
        _write_readme(selected, results_dir=root, outdir=outdir)
    else:
        # still write README with empty selection
        _write_readme(selected, results_dir=Path(args.results_dir[0]), outdir=outdir)

    print(f"[OK] Wrote summary to: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

