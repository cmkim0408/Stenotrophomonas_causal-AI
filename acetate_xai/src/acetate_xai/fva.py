from __future__ import annotations

import logging
from typing import Iterable

import pandas as pd

logger = logging.getLogger(__name__)


class FVAError(RuntimeError):
    """Raised when FBA/FVA cannot be computed."""


def run_targeted_fva(
    model,
    targets: Iterable[str],
    fraction_of_optimum: float = 0.95,
) -> pd.DataFrame:
    """
    Run targeted Flux Variability Analysis (FVA) for a set of reaction IDs.

    Parameters
    ----------
    model:
        cobra.Model (already configured with medium/bounds).
    targets:
        Reaction IDs to run FVA on.
    fraction_of_optimum:
        Fraction of optimal objective to enforce during FVA (e.g., 0.95).

    Returns
    -------
    DataFrame with columns: reaction_id, fva_min, fva_max
    """
    from cobra.flux_analysis import flux_variability_analysis

    target_list = list(dict.fromkeys([str(t) for t in targets]))
    if not target_list:
        raise ValueError("targets is empty")
    if not (0.0 < float(fraction_of_optimum) <= 1.0):
        raise ValueError("fraction_of_optimum must be in (0, 1].")

    logger.info("Running targeted FVA: n_targets=%d, fraction_of_optimum=%.3f", len(target_list), fraction_of_optimum)
    try:
        fva_df = flux_variability_analysis(
            model,
            reaction_list=target_list,
            fraction_of_optimum=float(fraction_of_optimum),
            loopless=False,
        )
    except Exception as e:  # noqa: BLE001
        raise FVAError(f"FVA failed: {e}") from e

    # cobra returns index=reaction_id, columns=["minimum","maximum"]
    out = (
        fva_df.rename(columns={"minimum": "fva_min", "maximum": "fva_max"})
        .reset_index()
        .rename(columns={"index": "reaction_id"})
    )
    if "reaction_id" not in out.columns:
        # depending on cobra/pandas, reset_index may name it differently
        out = out.rename(columns={out.columns[0]: "reaction_id"})

    out["reaction_id"] = out["reaction_id"].astype(str)
    return out[["reaction_id", "fva_min", "fva_max"]]

