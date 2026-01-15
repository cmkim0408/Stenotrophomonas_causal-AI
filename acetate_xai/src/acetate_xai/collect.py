from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


class CollectError(RuntimeError):
    """Raised when collection/feature building fails."""


def list_part_files(parts_dir: str | Path) -> list[Path]:
    p = Path(parts_dir)
    if not p.exists():
        raise FileNotFoundError(f"parts-dir not found: {p}")
    files = sorted(p.glob("*.parquet"))
    if not files:
        raise CollectError(f"No parquet files found under parts-dir: {p}")
    return files


def load_fva_parts(parts_dir: str | Path) -> pd.DataFrame:
    files = list_part_files(parts_dir)
    logger.info("Loading %d parquet parts from %s", len(files), parts_dir)
    dfs: list[pd.DataFrame] = []
    for f in files:
        dfs.append(pd.read_parquet(f))
    df = pd.concat(dfs, ignore_index=True)
    # minimal schema check
    required = {"condition_id", "objective_value", "reaction_id", "fva_min", "fva_max"}
    missing = required - set(df.columns)
    if missing:
        raise CollectError(f"Missing required columns in parts concat: {sorted(missing)}")
    return df


def build_fva_long_features(fva_all: pd.DataFrame) -> pd.DataFrame:
    """
    Add per-(condition_id, reaction_id) features to the long table:
    - fva_width = fva_max - fva_min
    - fva_mid = (fva_max + fva_min)/2
    - fva_abs_width = abs(fva_width)
    - sign_change = (fva_min < 0) & (fva_max > 0)
    """
    df = fva_all.copy()
    df["fva_width"] = df["fva_max"] - df["fva_min"]
    df["fva_mid"] = (df["fva_max"] + df["fva_min"]) / 2.0
    df["fva_abs_width"] = df["fva_width"].abs()
    df["sign_change"] = (df["fva_min"] < 0) & (df["fva_max"] > 0)
    return df


def build_wide_feature_matrix(long_df: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot long features into wide format: 1 row per condition_id.

    Output columns:
    - width__RXNID
    - mid__RXNID
    - signchange__RXNID
    """
    needed = {"condition_id", "reaction_id", "fva_width", "fva_mid", "sign_change"}
    missing = needed - set(long_df.columns)
    if missing:
        raise CollectError(f"Missing columns to build wide features: {sorted(missing)}")

    df = long_df[["condition_id", "reaction_id", "fva_width", "fva_mid", "sign_change"]].copy()
    # Ensure uniqueness (should be exactly one row per pair)
    dup = df.duplicated(subset=["condition_id", "reaction_id"]).any()
    if dup:
        raise CollectError("Duplicate rows found for (condition_id, reaction_id).")

    w = df.pivot(index="condition_id", columns="reaction_id", values="fva_width")
    m = df.pivot(index="condition_id", columns="reaction_id", values="fva_mid")
    s = df.pivot(index="condition_id", columns="reaction_id", values="sign_change")

    w.columns = [f"width__{c}" for c in w.columns.astype(str)]
    m.columns = [f"mid__{c}" for c in m.columns.astype(str)]
    s.columns = [f"signchange__{c}" for c in s.columns.astype(str)]

    wide = pd.concat([w, m, s], axis=1).reset_index()
    return wide


def join_conditions_features(
    *,
    features_wide: pd.DataFrame,
    conditions_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Join condition metadata + measured_OD onto wide features.
    """
    keep_cols = [
        "condition_id",
        "set_name",
        "measured_OD",
        "pH0",
        "yeast_extract_gL",
        "nh4cl_gL",
        "acetate_mM",
    ]
    missing = [c for c in keep_cols if c not in conditions_df.columns]
    if missing:
        raise CollectError(f"conditions CSV missing required columns for join: {missing}")

    meta = conditions_df[keep_cols].copy()
    out = meta.merge(features_wide, on="condition_id", how="inner", validate="one_to_one")
    return out


def collect_and_build_features(
    *,
    parts_dir: str | Path,
    conditions_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns (fva_all, features_joined).
    """
    fva_all = load_fva_parts(parts_dir)
    long_feat = build_fva_long_features(fva_all)
    wide = build_wide_feature_matrix(long_feat)
    features = join_conditions_features(features_wide=wide, conditions_df=conditions_df)
    return fva_all, features

