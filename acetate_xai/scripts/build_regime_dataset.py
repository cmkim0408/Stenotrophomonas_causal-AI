from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import pandas as pd


def _normpath(p: Path) -> str:
    # normalize separators for substring checks
    return str(p).replace("\\", "/")


def _is_true(x) -> bool:
    if isinstance(x, bool):
        return x
    if x is None:
        return False
    s = str(x).strip().lower()
    return s in {"true", "1", "yes", "y"}


def _label_row(r: pd.Series) -> str:
    # priority: acetate > nh4 > pi > o2 > unconstrained
    if _is_true(r.get("acetate_sat")):
        return "Ac_limited"
    if _is_true(r.get("nh4_sat")):
        return "N_limited"
    if _is_true(r.get("pi_sat")):
        return "Pi_limited"
    if _is_true(r.get("o2_sat")):
        return "O2_limited"
    return "Unconstrained"


def _find_regime_tables(campaign_root: Path) -> list[Path]:
    return sorted(campaign_root.rglob("regime_table.csv"))


def _choose_tables(regime_tables: list[Path]) -> list[tuple[Path, str, str, Path]]:
    """
    Apply selection rules and return tuples:
      (regime_table_csv, campaign_name, run_folder_name, mode) where mode in {"xai","xai_rxnfix"}
    """
    chosen: dict[tuple[str, str], tuple[Path, str]] = {}
    # key: (campaign, run_folder) -> (csv_path, mode)

    for csv_path in regime_tables:
        pnorm = _normpath(csv_path)

        # a) exclude C4_atpm_sweep (exact substring)
        if "/C4_atpm_sweep/" in pnorm:
            continue

        # Expect .../<campaign>/<run>/<xai or xai_rxnfix>/regime_table.csv
        xai_dir = csv_path.parent.name
        if xai_dir not in {"xai", "xai_rxnfix"}:
            continue
        run_dir = csv_path.parent.parent
        campaign_dir = run_dir.parent
        campaign_name = campaign_dir.name
        run_name = run_dir.name

        # c) C4_atpm_sweep_v2: only xai_rxnfix
        if campaign_name == "C4_atpm_sweep_v2" and xai_dir != "xai_rxnfix":
            continue

        key = (campaign_name, run_name)
        prev = chosen.get(key)
        if prev is None:
            chosen[key] = (csv_path, xai_dir)
            continue

        prev_path, prev_mode = prev
        # b) if xai_rxnfix exists, prefer it over xai
        if prev_mode == "xai" and xai_dir == "xai_rxnfix":
            chosen[key] = (csv_path, xai_dir)
        # otherwise keep existing

    out: list[tuple[Path, str, str, Path]] = []
    for (campaign_name, run_name), (csv_path, mode) in sorted(chosen.items()):
        out.append((csv_path, campaign_name, run_name, Path(mode)))
    return out


def _load_run_features(run_dir: Path) -> pd.DataFrame | None:
    fp = run_dir / "features.parquet"
    if fp.exists():
        return pd.read_parquet(fp)
    return None


def build_regime_dataset(*, campaign_root: Path) -> pd.DataFrame:
    all_csv = _find_regime_tables(campaign_root)
    if not all_csv:
        raise FileNotFoundError(f"No regime_table.csv found under: {campaign_root}")

    selected = _choose_tables(all_csv)
    if not selected:
        raise FileNotFoundError(f"No regime_table.csv left after filtering under: {campaign_root}")

    frames: list[pd.DataFrame] = []
    for csv_path, campaign_name, run_name, mode in selected:
        run_dir = csv_path.parent.parent
        run_id = f"{campaign_name}__{run_name}__{mode}"

        reg = pd.read_csv(csv_path)
        if "condition_id" not in reg.columns:
            logging.warning("Skip (no condition_id): %s", csv_path)
            continue

        reg["run_id"] = run_id
        reg["label"] = reg.apply(_label_row, axis=1)

        feat = _load_run_features(run_dir)
        if feat is None:
            logging.warning("Skip (features.parquet missing in run dir): %s", run_dir)
            continue
        if "condition_id" not in feat.columns:
            logging.warning("Skip (features.parquet missing condition_id): %s", run_dir / "features.parquet")
            continue

        merged = reg.merge(feat, on="condition_id", how="inner")
        merged["campaign"] = campaign_name
        merged["run_folder"] = run_name
        merged["xai_mode"] = str(mode)

        frames.append(merged)

    if not frames:
        raise RuntimeError("No runs produced merged data (all skipped due to missing columns/files).")
    return pd.concat(frames, ignore_index=True)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Build regime_dataset.parquet from campaigns regime_table.csv + features.parquet.")
    p.add_argument("--campaign-root", required=True, help="Root directory containing campaign runs (e.g., results/campaigns)")
    p.add_argument("--out", required=True, help="Output parquet path (e.g., results/regime_dataset.parquet)")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    campaign_root = Path(args.campaign_root)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # help reduce oversubscription when run on servers
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")

    df = build_regime_dataset(campaign_root=campaign_root)
    df.to_parquet(out_path, index=False)
    print(f"[OK] Wrote {out_path} (rows={len(df)}, runs={df['run_id'].nunique()})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

