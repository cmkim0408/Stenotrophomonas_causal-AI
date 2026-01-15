from __future__ import annotations

import logging
from pathlib import Path
from typing import Literal

import pandas as pd

logger = logging.getLogger(__name__)

ALLOWED_SET_NAMES: tuple[str, ...] = (
    "yeast_gradient",
    "pH_YE_toggle",
    "nh4_gradient",
    "acetate_gradient",
)

REQUIRED_CONDITIONS_COLUMNS: tuple[str, ...] = (
    "condition_id",
    "set_name",
    "pH0",
    "yeast_extract_gL",
    "nh4cl_gL",
    "acetate_mM",
    "notes",
    "measured_OD",
)


def load_sbml_model(sbml_path: str | Path):
    """
    Load an SBML model using cobra.

    Returns
    -------
    cobra.Model
    """
    from cobra.io import read_sbml_model

    p = Path(sbml_path)
    if not p.exists():
        raise FileNotFoundError(f"SBML file not found: {p}")
    logger.info("Loading SBML model: %s", p)
    return read_sbml_model(str(p))


def load_conditions_csv(path: str | Path) -> pd.DataFrame:
    """
    Load the experiment conditions table (one row = one condition).

    Required columns
    ----------------
    - condition_id
    - set_name: one of {"yeast_gradient","pH_YE_toggle","nh4_gradient","acetate_gradient"}
    - pH0
    - yeast_extract_gL
    - nh4cl_gL
    - acetate_mM
    - notes
    - measured_OD (if missing in file, will be created with NaN)
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Conditions CSV not found: {p}")

    df = pd.read_csv(p)

    # measured_OD: allow missing in early drafts, but always expose it downstream.
    if "measured_OD" not in df.columns:
        df["measured_OD"] = pd.NA

    missing = [c for c in REQUIRED_CONDITIONS_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            "Missing required columns in conditions CSV: "
            + ", ".join(missing)
            + f". Required columns: {', '.join(REQUIRED_CONDITIONS_COLUMNS)}"
        )

    # Normalize/validate set_name
    df["set_name"] = df["set_name"].astype(str).str.strip()
    bad = sorted(set(df.loc[~df["set_name"].isin(ALLOWED_SET_NAMES), "set_name"].tolist()))
    if bad:
        raise ValueError(
            "Invalid set_name values: "
            + ", ".join(bad)
            + f". Allowed: {', '.join(ALLOWED_SET_NAMES)}"
        )

    # Light type coercions (keep it forgiving; errors become NaN)
    for col in ("pH0", "yeast_extract_gL", "nh4cl_gL", "acetate_mM", "measured_OD"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["condition_id"] = df["condition_id"].astype(str).str.strip()
    df["notes"] = df["notes"].fillna("").astype(str)

    return df


def save_table(
    df: pd.DataFrame,
    out_path: str | Path,
    *,
    fmt: Literal["parquet", "csv"] | None = None,
) -> Path:
    """
    Save a table to parquet or CSV, inferred by extension unless fmt is provided.
    """
    p = Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)

    if fmt is None:
        suffix = p.suffix.lower()
        if suffix == ".parquet":
            fmt = "parquet"
        elif suffix == ".csv":
            fmt = "csv"
        else:
            raise ValueError(f"Cannot infer format from extension: {p.suffix} (use .parquet or .csv)")

    if fmt == "parquet":
        df.to_parquet(p, index=False)
    elif fmt == "csv":
        df.to_csv(p, index=False)
    else:
        raise ValueError(f"Unsupported fmt: {fmt}")

    logger.info("Saved table: %s (rows=%d, cols=%d)", p, len(df), len(df.columns))
    return p

