from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class RegimeError(RuntimeError):
    """Raised when FBA regime computation fails unexpectedly."""


@dataclass(frozen=True)
class SatResult:
    rid: str
    flux: float
    lb: float
    ub: float
    is_constrained: bool
    saturated: bool
    sat_side: str  # "lb" | "ub" | "fixed" | "none" | "open" | "missing"


def pick_first_existing_reaction_id(model, candidates: list[str]) -> str | None:
    present = set(r.id for r in model.reactions)
    for rid in candidates:
        if rid in present:
            return rid
    return None


def _is_open_bound(lb: float, ub: float, infty_bound: float) -> bool:
    return (lb <= -infty_bound) and (ub >= infty_bound)


def compute_saturation_for_reaction(
    *,
    rid: str,
    model,
    solution,
    eps: float = 1e-6,
    infty_bound: float = 999.0,
) -> SatResult:
    rxn = model.reactions.get_by_id(rid)
    lb = float(rxn.lower_bound)
    ub = float(rxn.upper_bound)
    flux = float(solution.fluxes[rid]) if rid in solution.fluxes.index else float("nan")

    if _is_open_bound(lb, ub, infty_bound=infty_bound):
        return SatResult(rid=rid, flux=flux, lb=lb, ub=ub, is_constrained=False, saturated=False, sat_side="open")

    # If either side is finite, treat as constrained (meaningful for limitation)
    is_constrained = True

    # Fixed
    if abs(lb - ub) <= eps and abs(flux - lb) <= eps:
        return SatResult(rid=rid, flux=flux, lb=lb, ub=ub, is_constrained=is_constrained, saturated=True, sat_side="fixed")

    # Lower bound saturation
    if abs(flux - lb) <= eps:
        return SatResult(rid=rid, flux=flux, lb=lb, ub=ub, is_constrained=is_constrained, saturated=True, sat_side="lb")

    # Upper bound saturation
    if abs(flux - ub) <= eps:
        return SatResult(rid=rid, flux=flux, lb=lb, ub=ub, is_constrained=is_constrained, saturated=True, sat_side="ub")

    return SatResult(rid=rid, flux=flux, lb=lb, ub=ub, is_constrained=is_constrained, saturated=False, sat_side="none")


def run_fba_regime_table(
    *,
    model_path: str,
    conditions_df: pd.DataFrame,
    medium_cfg: dict[str, Any],
    regime_cfg: dict[str, Any],
    limit: int | None = None,
    condition_ids: list[str] | None = None,
    eps: float = 1e-6,
    infty_bound: float = 999.0,
) -> pd.DataFrame:
    """
    Run FBA once per condition and compute saturation flags using FBA flux + bounds.

    Output is 1 row per condition_id (wide).
    """
    from acetate_xai.io import load_sbml_model
    from acetate_xai.medium import apply_condition_to_model

    df = conditions_df.copy()
    if condition_ids:
        wanted = set(str(x) for x in condition_ids)
        df = df.loc[df["condition_id"].astype(str).isin(wanted)].copy()
    elif limit is not None:
        df = df.head(int(limit)).copy()

    # regime_cfg schema: nutrient -> list of candidate rxn_ids
    def _candidates(key: str) -> list[str]:
        c = regime_cfg.get(key, [])
        if isinstance(c, list):
            return [str(x) for x in c]
        return []

    nutrients = {
        "acetate": _candidates("acetate"),
        "oxygen": _candidates("oxygen"),
        "ammonium": _candidates("ammonium"),
        "phosphate": _candidates("phosphate"),
    }

    out_rows: list[dict[str, Any]] = []
    for _, row in df.iterrows():
        cid = str(row["condition_id"])

        model = load_sbml_model(model_path)
        apply_condition_to_model(model, row.to_dict(), medium_cfg)

        sol = model.optimize()
        if sol.status != "optimal":
            # record row with NaNs but keep table shape stable
            out_rows.append({"condition_id": cid, "objective_value": np.nan})
            continue

        rec: dict[str, Any] = {"condition_id": cid, "objective_value": float(sol.objective_value)}

        for nutrient, cand in nutrients.items():
            rid_used = pick_first_existing_reaction_id(model, cand) if cand else None
            if rid_used is None:
                rec.update(
                    {
                        f"{nutrient}_rid": "",
                        f"{nutrient}_flux": np.nan,
                        f"{nutrient}_lb": np.nan,
                        f"{nutrient}_ub": np.nan,
                        f"{nutrient}_is_constrained": False,
                        f"{nutrient}_sat": False,
                        f"{nutrient}_sat_side": "missing",
                    }
                )
                continue

            sat = compute_saturation_for_reaction(
                rid=rid_used,
                model=model,
                solution=sol,
                eps=eps,
                infty_bound=infty_bound,
            )
            rec.update(
                {
                    f"{nutrient}_rid": sat.rid,
                    f"{nutrient}_flux": sat.flux,
                    f"{nutrient}_lb": sat.lb,
                    f"{nutrient}_ub": sat.ub,
                    f"{nutrient}_is_constrained": sat.is_constrained,
                    f"{nutrient}_sat": sat.saturated,
                    f"{nutrient}_sat_side": sat.sat_side,
                }
            )

        out_rows.append(rec)

    out = pd.DataFrame(out_rows)
    # Ensure one row per condition
    if out["condition_id"].duplicated().any():
        raise RegimeError("Duplicate condition_id rows produced in FBA regime table.")
    return out

