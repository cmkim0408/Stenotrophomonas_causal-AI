from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Plot regime map: primary_regime vs maintenance_severity colored by run_id.")
    p.add_argument("--dataset", required=True, help="Input dataset parquet (e.g., results/regime_dataset.parquet)")
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
    """
    Ensure `severity` exists:
      severity = objective_value / max(objective_value) within run_id
    """
    if "severity" in df.columns and df["severity"].notna().any():
        return df
    tmp = df.copy()
    tmp["objective_value"] = pd.to_numeric(tmp["objective_value"], errors="coerce")
    run_max = tmp.groupby("run_id")["objective_value"].max().rename("objective_value_max")
    tmp = tmp.merge(run_max.reset_index(), on="run_id", how="left")
    den = tmp["objective_value_max"].replace({0.0: np.nan})
    tmp["severity"] = tmp["objective_value"] / den
    return tmp


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(args.dataset)
    required = {"run_id", "condition_id", "label", "objective_value"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"dataset missing required columns: {sorted(missing)}")

    df = _ensure_severity(df)

    # 1) primary_regime per run_id (mode of label)
    run_primary = df.groupby("run_id")["label"].apply(_mode).rename("primary_regime").reset_index()

    # 2) condition-level map table
    map_df = df[["run_id", "condition_id", "label", "severity"]].copy()
    map_df = map_df.merge(run_primary, on="run_id", how="left")
    map_df = map_df.rename(columns={"severity": "maintenance_severity"})
    map_df.to_csv(outdir / "regime_map.csv", index=False)

    # 3) plot: x=categorical primary_regime, y=maintenance_severity, color=run_id
    try:
        import matplotlib.pyplot as plt
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(f"matplotlib is required to plot. Install it and retry. Error: {e}") from e

    # Categorical x positions
    x_categories = sorted(map_df["primary_regime"].astype(str).unique().tolist())
    x_map = {cat: i for i, cat in enumerate(x_categories)}
    map_df["x"] = map_df["primary_regime"].astype(str).map(x_map)

    # Color by run_id (hash to colormap)
    run_ids = sorted(map_df["run_id"].astype(str).unique().tolist())
    cmap = plt.get_cmap("tab20")
    colors = {rid: cmap(i % 20) for i, rid in enumerate(run_ids)}

    fig, ax = plt.subplots(figsize=(10, 4))
    for rid in run_ids:
        sub = map_df.loc[map_df["run_id"].astype(str) == rid]
        ax.scatter(
            sub["x"].values,
            sub["maintenance_severity"].values,
            s=18,
            alpha=0.8,
            color=colors[rid],
            label=rid,
        )

    ax.set_xticks(list(x_map.values()))
    ax.set_xticklabels(x_categories, rotation=0)
    ax.set_xlabel("primary_regime (mode per run_id)")
    ax.set_ylabel("maintenance_severity (objective / run max)")
    ax.set_title("Regime map (colored by run_id)")
    ax.grid(True, axis="y", linestyle="--", alpha=0.3)

    # Legend: if too many runs, place outside and limit size
    if len(run_ids) <= 12:
        ax.legend(fontsize=8, loc="best")
    else:
        ax.legend(fontsize=6, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.0, ncol=1)
        fig.tight_layout(rect=[0, 0, 0.80, 1])

    out_png = outdir / "regime_map.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] Wrote: {outdir / 'regime_map.csv'}")
    print(f"[OK] Wrote: {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

