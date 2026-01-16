from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Summarize regime platform results into compact tables + report.")
    p.add_argument("--dataset", required=True, help="Input dataset parquet (e.g., results/regime_dataset.parquet)")
    p.add_argument("--xai-dir", required=True, help="Directory with xai platform outputs (e.g., results/xai_platform)")
    p.add_argument(
        "--maintenance-dir",
        required=True,
        help="Directory with maintenance regression outputs (e.g., results/xai_platform_maintenance)",
    )
    p.add_argument("--outdir", required=True, help="Output directory (e.g., results/platform_summary)")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def _mode(series: pd.Series) -> Any:
    s = series.dropna()
    if s.empty:
        return np.nan
    vc = s.value_counts()
    top = vc[vc == vc.max()].index.astype(str).tolist()
    return sorted(top)[0]


def _ensure_severity(df: pd.DataFrame) -> pd.DataFrame:
    if "severity" in df.columns and df["severity"].notna().any():
        return df
    if "objective_value" not in df.columns:
        return df
    tmp = df.copy()
    tmp["objective_value"] = pd.to_numeric(tmp["objective_value"], errors="coerce")
    run_max = tmp.groupby("run_id")["objective_value"].max().rename("objective_value_max")
    tmp = tmp.merge(run_max.reset_index(), on="run_id", how="left")
    den = tmp["objective_value_max"].replace({0.0: np.nan})
    tmp["severity"] = tmp["objective_value"] / den
    return tmp


def _feature_to_module(feature_name: str) -> str:
    """
    Heuristic mapping from a feature name (e.g., width__EX_ac_e) to a module label.
    This is intentionally simple and robust across models.
    """
    s = str(feature_name)
    if "__" in s:
        rxn = s.split("__", 1)[1]
    else:
        rxn = s
    rxn_u = rxn.upper()

    if rxn_u.startswith("EX_"):
        return "Exchange"
    if "ATPM" in rxn_u or rxn_u == "ATPM":
        return "Maintenance"
    if "NH4" in rxn_u or "NH3" in rxn_u:
        return "Nitrogen"
    if "O2" in rxn_u:
        return "Oxygen"
    if rxn_u.endswith("_PI_E") or "EX_PI" in rxn_u or "PI_" in rxn_u:
        return "Phosphate"

    tca_like = ("CS", "GLTA", "ACONT", "ICD", "ICDH", "AKGDH", "FUM", "MDH", "SDH", "ICL", "MALS", "ACEA", "ACEB")
    if any(tok in rxn_u for tok in tca_like):
        return "TCA_Glyoxylate"

    ppp_like = ("G6PD", "6PGDH", "TKT", "TALA", "RPI", "RPE", "PRPP")
    if any(tok in rxn_u for tok in ppp_like):
        return "PPP_NADPH"

    gng_like = ("PPS", "PCK", "PYK", "ENO", "FBP", "PFK", "PGK", "GAPD", "TPI", "FBA")
    if any(tok in rxn_u for tok in gng_like):
        return "GNG_EMP"

    return "Other"


def _load_linear_features_csv(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        df = pd.read_csv(path)
    except Exception:
        return None
    if "feature" not in df.columns:
        return None
    return df


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dataset = pd.read_parquet(args.dataset)
    if "run_id" not in dataset.columns:
        raise ValueError("dataset must contain run_id")
    if "label" not in dataset.columns:
        raise ValueError("dataset must contain label")
    if "objective_value" not in dataset.columns:
        raise ValueError("dataset must contain objective_value")

    dataset["objective_value"] = pd.to_numeric(dataset["objective_value"], errors="coerce")
    dataset = _ensure_severity(dataset)

    # (1) platform_summary.csv: run-level stats
    run_stats = (
        dataset.groupby("run_id")
        .agg(
            primary_regime=("label", _mode),
            objective_mean=("objective_value", "mean"),
            objective_min=("objective_value", "min"),
            maintenance_severity_mean=("severity", "mean"),
            n_conditions=("condition_id", "nunique") if "condition_id" in dataset.columns else ("label", "size"),
        )
        .reset_index()
    )
    run_stats.to_csv(outdir / "platform_summary.csv", index=False)

    # (2) regime_coverage.csv: overall label distribution
    cov = dataset["label"].astype(str).value_counts().rename_axis("label").reset_index(name="count")
    cov["fraction"] = cov["count"] / cov["count"].sum()
    cov.to_csv(outdir / "regime_coverage.csv", index=False)

    # (4) signature_modules.csv: module frequency from linear top features (classification + maintenance)
    xai_linear = _load_linear_features_csv(Path(args.xai_dir) / "linear_top_features.csv")
    maint_linear = _load_linear_features_csv(Path(args.maintenance_dir) / "linear_top_features.csv")
    mod_frames: list[pd.DataFrame] = []
    for source, df in [("xai_platform", xai_linear), ("maintenance", maint_linear)]:
        if df is None:
            continue
        mods = df["feature"].astype(str).map(_feature_to_module)
        m = mods.value_counts().rename_axis("module").reset_index(name="count")
        m["fraction"] = m["count"] / m["count"].sum() if m["count"].sum() else 0.0
        m.insert(0, "source", source)
        mod_frames.append(m)
    if mod_frames:
        pd.concat(mod_frames, ignore_index=True).to_csv(outdir / "signature_modules.csv", index=False)
    else:
        pd.DataFrame(columns=["source", "module", "count", "fraction"]).to_csv(outdir / "signature_modules.csv", index=False)

    # (3) robustness_summary.txt: fixed narrative + inserted numbers
    n_rows = len(dataset)
    n_runs = int(dataset["run_id"].nunique())
    labels = cov.set_index("label")["count"].to_dict()
    obj_mean = float(dataset["objective_value"].mean())
    obj_min = float(dataset["objective_value"].min())
    sev_mean = float(dataset["severity"].mean()) if "severity" in dataset.columns else float("nan")

    lines = [
        "Platform Robustness Summary",
        "===========================",
        "",
        f"- Total rows: {n_rows}",
        f"- Total runs: {n_runs}",
        f"- Objective (mean/min): {obj_mean:.6g} / {obj_min:.6g}",
        f"- Maintenance severity (mean): {sev_mean:.6g}",
        "",
        "Regime coverage (counts):",
    ]
    for lab, cnt in labels.items():
        lines.append(f"- {lab}: {cnt}")
    lines.extend(
        [
            "",
            "Interpretation note:",
            "- 'primary_regime' is the mode of labels within each run_id.",
            "- 'severity' is normalized within run_id by objective_value / max(objective_value).",
            "- 'signature_modules' summarizes module frequencies derived from top linear features.",
            "",
        ]
    )
    (outdir / "robustness_summary.txt").write_text("\n".join(lines), encoding="utf-8")

    print(f"[OK] Wrote: {outdir / 'platform_summary.csv'}")
    print(f"[OK] Wrote: {outdir / 'regime_coverage.csv'}")
    print(f"[OK] Wrote: {outdir / 'signature_modules.csv'}")
    print(f"[OK] Wrote: {outdir / 'robustness_summary.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

