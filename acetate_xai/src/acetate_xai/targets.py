from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import yaml

logger = logging.getLogger(__name__)


class TargetSelectionError(RuntimeError):
    """Raised when target selection cannot be completed."""


@dataclass(frozen=True)
class Anchor:
    name: str
    keywords: tuple[str, ...]


def load_anchors_yaml(path: str | Path) -> list[Anchor]:
    """
    Load anchors YAML.

    Expected schema
    ---------------
    anchors:
      - name: "acetate exchange"
        keywords: ["EX_ac", "acetate exchange"]
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Anchors YAML not found: {p}")

    with p.open("r", encoding="utf-8") as f:
        raw = yaml.safe_load(f) or {}

    if not isinstance(raw, dict) or "anchors" not in raw:
        raise ValueError("anchors.yaml must be a mapping with top-level key 'anchors'.")
    if not isinstance(raw["anchors"], list):
        raise ValueError("anchors.yaml['anchors'] must be a list.")

    anchors: list[Anchor] = []
    for i, item in enumerate(raw["anchors"]):
        if not isinstance(item, dict):
            raise ValueError(f"Anchor #{i} must be a mapping with keys: name, keywords.")
        name = str(item.get("name", "")).strip()
        kws = item.get("keywords", [])
        if not name:
            raise ValueError(f"Anchor #{i} has empty 'name'.")
        if not isinstance(kws, list) or not kws:
            raise ValueError(f"Anchor '{name}' must have a non-empty keywords list.")
        keywords = tuple(str(k).strip() for k in kws if str(k).strip())
        if not keywords:
            raise ValueError(f"Anchor '{name}' has no usable keywords.")
        anchors.append(Anchor(name=name, keywords=keywords))

    return anchors


def _kw_match(text: str, keywords: Iterable[str]) -> bool:
    t = text.lower()
    for kw in keywords:
        if kw.lower() in t:
            return True
    return False


def find_anchor_matches(model, anchors: list[Anchor]) -> list[str]:
    """
    Return reaction IDs matched by any anchor keyword against reaction.id or reaction.name.
    """
    selected: list[str] = []
    for rxn in model.reactions:
        rid = str(rxn.id)
        rname = str(rxn.name or "")
        for a in anchors:
            if _kw_match(rid, a.keywords) or _kw_match(rname, a.keywords):
                selected.append(rid)
                break
    # preserve order but remove duplicates
    return list(dict.fromkeys(selected))


def find_blocked_reaction_ids(model) -> set[str]:
    """
    Identify blocked reactions under the current model bounds/medium.

    Note: This uses cobra's blocked-reaction analysis (fast/standard). It does not
    require per-reaction FVA and is appropriate for large GSMs.
    """
    from cobra.flux_analysis import find_blocked_reactions

    blocked = find_blocked_reactions(model)
    return set(str(rid) for rid in blocked)


def rank_reactions_by_pfba_flux(model) -> list[str]:
    """
    Run pFBA once and rank reactions by |flux| descending.
    """
    from cobra.flux_analysis import pfba

    sol = pfba(model)
    if sol.status != "optimal":
        raise TargetSelectionError(f"pFBA failed: status={sol.status}")

    s = sol.fluxes.abs().sort_values(ascending=False)
    return [str(rid) for rid in s.index.tolist()]


def select_targets_anchored(
    *,
    model,
    anchors: list[Anchor],
    target_count: int = 120,
    blocked_eps_note: str | None = None,
) -> list[str]:
    """
    Select targets using anchor matching + blocked filtering + pFBA-based auto-fill.

    Returns
    -------
    List of reaction IDs of length == target_count.
    """
    if target_count <= 0:
        raise ValueError("target_count must be > 0")

    anchor_ids = find_anchor_matches(model, anchors)
    logger.info("Anchor matched %d reactions (pre-dedup).", len(anchor_ids))

    blocked = find_blocked_reaction_ids(model)
    logger.info("Blocked reactions detected: %d", len(blocked))
    if blocked_eps_note:
        logger.info("Blocked note: %s", blocked_eps_note)

    targets: list[str] = [rid for rid in anchor_ids if rid not in blocked]
    targets = list(dict.fromkeys(targets))
    logger.info("After removing blocked + dedup: %d targets", len(targets))

    if len(targets) >= target_count:
        return targets[:target_count]

    ranked = rank_reactions_by_pfba_flux(model)
    for rid in ranked:
        if rid in blocked:
            continue
        if rid in targets:
            continue
        targets.append(rid)
        if len(targets) >= target_count:
            break

    if len(targets) != target_count:
        raise TargetSelectionError(
            f"Could not reach target_count={target_count}. Got {len(targets)}. "
            "This can happen if the model has too many blocked reactions under the base medium."
        )

    return targets


def save_targets_json(target_ids: list[str], out_path: str | Path) -> Path:
    """
    Save targets as a JSON array (length must be exact for downstream checks).
    """
    p = Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as f:
        json.dump(target_ids, f, indent=2)
    return p

