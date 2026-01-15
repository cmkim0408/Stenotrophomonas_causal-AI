from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

logger = logging.getLogger(__name__)


class MediumConfigError(ValueError):
    """Raised when medium config is invalid or inconsistent with the model."""


@dataclass(frozen=True)
class MediumApplyResult:
    """
    Records what we changed on the model for traceability.

    Notes
    -----
    We keep pH0 only as metadata in this initial version (no direct constraints).
    """

    condition_id: str | None
    pH0: float | None
    yeast_enabled: bool
    changed_bounds: list[tuple[str, float | None, float | None, float | None, float | None]]
    # (rxn_id, old_lb, old_ub, new_lb, new_ub)


def _get_rxn(model, rxn_id: str):
    try:
        return model.reactions.get_by_id(rxn_id)
    except KeyError as e:
        raise MediumConfigError(f"Reaction not found in model: {rxn_id}") from e


def _set_bounds_with_log(model, rxn_id: str, *, lb: float | None, ub: float | None, changes: list):
    rxn = _get_rxn(model, rxn_id)
    old_lb, old_ub = rxn.lower_bound, rxn.upper_bound
    new_lb = old_lb if lb is None else float(lb)
    new_ub = old_ub if ub is None else float(ub)
    if new_lb != old_lb or new_ub != old_ub:
        # Set ub first to avoid transient invalid bounds when widening ub from negative values.
        rxn.upper_bound = new_ub
        rxn.lower_bound = new_lb
        changes.append((rxn_id, old_lb, old_ub, new_lb, new_ub))
        logger.info("Bound update %s: lb %.6g -> %.6g, ub %.6g -> %.6g", rxn_id, old_lb, new_lb, old_ub, new_ub)

def _set_exchange_uptake_max_with_log(model, rxn_id: str, *, uptake_max: float, changes: list) -> None:
    """
    Set a maximum uptake for an exchange reaction using the standard COBRA convention:
    uptake is negative flux, so we set lower_bound = -uptake_max.

    We also sanitize cases where the model has a negative upper bound (which would make
    it impossible to set a less-negative lower bound without violating lb <= ub).
    In v0 we clamp upper_bound to 0.0 in that case (no secretion; uptake-only).
    """
    rxn = _get_rxn(model, rxn_id)
    old_lb, old_ub = rxn.lower_bound, rxn.upper_bound

    u = max(0.0, float(uptake_max))
    new_lb = -u
    new_ub = old_ub
    if new_ub < 0.0:
        new_ub = 0.0

    if new_lb > new_ub:
        # Safety: ensure bounds are consistent.
        new_ub = 0.0

    if new_lb != old_lb or new_ub != old_ub:
        # Set ub first to avoid transient invalid bounds when widening ub from negative values.
        rxn.upper_bound = new_ub
        rxn.lower_bound = new_lb
        changes.append((rxn_id, old_lb, old_ub, new_lb, new_ub))
        logger.info(
            "Exchange uptake cap %s: uptake_max=%.6g => lb %.6g -> %.6g, ub %.6g -> %.6g",
            rxn_id,
            u,
            old_lb,
            new_lb,
            old_ub,
            new_ub,
        )


def _try_set_bounds_with_log(
    model, rxn_id: str, *, lb: float | None, ub: float | None, changes: list
) -> None:
    """
    Like _set_bounds_with_log but does not hard-fail if reaction is missing.
    Intended for optional base-medium bounds or optional YE-open exchanges.
    """
    try:
        _set_bounds_with_log(model, rxn_id, lb=lb, ub=ub, changes=changes)
    except MediumConfigError:
        logger.warning("Reaction not found in model (skipped): %s", rxn_id)


def apply_condition_to_model(
    model,
    condition_row: dict[str, Any],
    medium_config: dict[str, Any],
) -> MediumApplyResult:
    """
    Apply one experimental condition row to a COBRA model by updating exchange bounds.

    Expected condition_row keys (subset used here)
    ----------------------------------------------
    - condition_id (optional)
    - pH0 (stored only; not applied to constraints yet)
    - acetate_mM
    - nh4cl_gL
    - yeast_extract_gL

    Scaling rules (initial)
    -----------------------
    - acetate_uptake_max = k_ac * acetate_mM
      -> set EX_ac_e.lower_bound = -acetate_uptake_max
    - nh4_uptake_max = k_nh4 * nh4cl_gL
      -> set EX_nh4_e.lower_bound = -nh4_uptake_max
    - yeast_extract_gL:
      if > enabled_if_gL_gt, open configured "vitamin/cofactor" exchanges to a small uptake.
    """
    if not isinstance(medium_config, dict):
        raise MediumConfigError("medium_config must be a dict (loaded from YAML/JSON).")

    exchanges = medium_config.get("exchanges", {})
    scaling = medium_config.get("scaling", {})
    base_bounds = medium_config.get("base_bounds", {})
    yeast_cfg = medium_config.get("yeast_extract", {})

    if not isinstance(exchanges, dict) or not exchanges:
        raise MediumConfigError("medium_config.exchanges must be a non-empty mapping.")

    changes: list[tuple[str, float | None, float | None, float | None, float | None]] = []

    # 1) Apply base medium bounds (if present)
    if isinstance(base_bounds, dict):
        for rxn_id, b in base_bounds.items():
            if not isinstance(b, dict):
                continue
            _try_set_bounds_with_log(
                model,
                str(rxn_id),
                lb=b.get("lb", None),
                ub=b.get("ub", None),
                changes=changes,
            )

    # Pull out condition values (forgiving parsing)
    def _get_float(key: str) -> float | None:
        v = condition_row.get(key, None)
        if v is None:
            return None
        try:
            return float(v)
        except Exception:  # noqa: BLE001
            return None

    acetate_mM = _get_float("acetate_mM")
    nh4cl_gL = _get_float("nh4cl_gL")
    yeast_gL = _get_float("yeast_extract_gL")
    pH0 = _get_float("pH0")
    condition_id = condition_row.get("condition_id", None)
    if condition_id is not None:
        condition_id = str(condition_id)

    # 2) Condition-specific scaling updates
    k_ac = float(scaling.get("k_ac", 0.0))
    k_nh4 = float(scaling.get("k_nh4", 0.0))

    ac_ex = exchanges.get("acetate", None)
    nh4_ex = exchanges.get("ammonium", None)
    if ac_ex is None or nh4_ex is None:
        raise MediumConfigError("medium_config.exchanges must include keys: acetate, ammonium")

    if acetate_mM is not None:
        uptake_max = max(0.0, k_ac * acetate_mM)
        _set_exchange_uptake_max_with_log(model, str(ac_ex), uptake_max=uptake_max, changes=changes)

    if nh4cl_gL is not None:
        uptake_max = max(0.0, k_nh4 * nh4cl_gL)
        _set_exchange_uptake_max_with_log(model, str(nh4_ex), uptake_max=uptake_max, changes=changes)

    # 3) Yeast extract handling: toggle + open exchanges (optional)
    enabled_if = float(yeast_cfg.get("enabled_if_gL_gt", 0.0))
    yeast_enabled = (yeast_gL is not None) and (yeast_gL > enabled_if)
    if yeast_enabled:
        open_lb = float(yeast_cfg.get("open_uptake_lb", -1.0))
        open_list = yeast_cfg.get("open_exchanges_when_enabled", [])
        if isinstance(open_list, list):
            for rid in open_list:
                rid_s = str(rid)
                _try_set_bounds_with_log(model, rid_s, lb=open_lb, ub=None, changes=changes)

    # pH0: intentionally not applied to constraints in initial version.
    if pH0 is not None:
        logger.info("Condition metadata pH0=%.3f (not applied to constraints in v0)", pH0)

    return MediumApplyResult(
        condition_id=condition_id,
        pH0=pH0,
        yeast_enabled=yeast_enabled,
        changed_bounds=changes,
    )

