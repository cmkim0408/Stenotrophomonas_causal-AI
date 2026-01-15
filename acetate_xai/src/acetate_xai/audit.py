from __future__ import annotations

import csv
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class AuditRow:
    requested_id: str
    status: str  # "present" | "missing"
    suggestion_1: str
    suggestion_2: str
    suggestion_3: str


def collect_requested_exchange_ids(medium_config: dict[str, Any]) -> list[str]:
    """
    Collect exchange reaction IDs referenced by medium_base.yaml.

    Sources
    -------
    - exchanges: mapping values
    - base_bounds: mapping keys
    - yeast_extract.open_exchanges_when_enabled: list values
    """
    ids: list[str] = []

    exchanges = medium_config.get("exchanges", {})
    if isinstance(exchanges, dict):
        ids.extend([str(v) for v in exchanges.values() if v is not None and str(v).strip()])

    base_bounds = medium_config.get("base_bounds", {})
    if isinstance(base_bounds, dict):
        ids.extend([str(k) for k in base_bounds.keys() if k is not None and str(k).strip()])

    yeast_cfg = medium_config.get("yeast_extract", {})
    if isinstance(yeast_cfg, dict):
        open_list = yeast_cfg.get("open_exchanges_when_enabled", [])
        if isinstance(open_list, list):
            ids.extend([str(x) for x in open_list if x is not None and str(x).strip()])

    # de-dup while preserving order
    return list(dict.fromkeys(ids))


def _normalize_text(s: str) -> str:
    return re.sub(r"\s+", " ", (s or "").lower()).strip()


_SYNONYMS: dict[str, list[str]] = {
    "btn": ["biotin", "vitamin b7"],
    "thm": ["thiamine", "vitamin b1"],
    "ribflv": ["riboflavin", "vitamin b2"],
    "fol": ["folate", "folic", "vitamin b9", "pteroyl"],
    "nad": ["nicotinamide", "niacin", "nad", "nadh"],
    "pydx": ["pyridoxine", "vitamin b6"],
    "pan4p": ["pantothenate", "vitamin b5"],
}


def _keywords_from_requested_id(requested_id: str) -> list[str]:
    """
    Heuristic keyword extraction from an exchange reaction id like 'EX_btn_e'.
    """
    rid = _normalize_text(requested_id)
    rid = rid.replace("exchange", "")

    # Common patterns
    rid = rid.replace("ex_", "")
    rid = rid.replace("ex", "")
    rid = rid.replace("(e)", "_e")

    tokens = re.split(r"[^a-z0-9]+", rid)
    tokens = [t for t in tokens if t and t not in {"e", "lp", "rev"}]

    kws: list[str] = []
    for t in tokens:
        kws.append(t)
        kws.extend(_SYNONYMS.get(t, []))
    # de-dup
    return list(dict.fromkeys(kws))


def _is_exchange_reaction(rxn) -> bool:
    """
    Try to recognize exchange reactions robustly.
    """
    rid = str(getattr(rxn, "id", ""))
    if rid.startswith("EX_") or rid.startswith("DM_") or rid.startswith("SK_"):
        return True
    # cobra has rxn.boundary for exchange/demand/sink in many models
    if bool(getattr(rxn, "boundary", False)):
        return True
    return False


def _reaction_search_text(rxn) -> str:
    parts: list[str] = []
    parts.append(str(getattr(rxn, "id", "")))
    parts.append(str(getattr(rxn, "name", "")))
    # metabolite ids/names on exchange often contain the clue
    for met in getattr(rxn, "metabolites", {}).keys():
        parts.append(str(getattr(met, "id", "")))
        parts.append(str(getattr(met, "name", "")))
    return _normalize_text(" ".join(parts))


def suggest_exchange_replacements(model, requested_id: str, top_k: int = 3) -> list[str]:
    """
    Suggest alternative exchange reaction IDs by keyword matching against:
    - reaction.id / reaction.name
    - metabolite id/name participating in the exchange

    Returns up to top_k reaction IDs, ranked by simple match score.
    """
    keywords = _keywords_from_requested_id(requested_id)
    if not keywords:
        return []

    candidates = [rxn for rxn in model.reactions if _is_exchange_reaction(rxn)]

    scored: list[tuple[int, str]] = []
    for rxn in candidates:
        text = _reaction_search_text(rxn)
        score = 0
        for kw in keywords:
            kw_n = _normalize_text(kw)
            if not kw_n or len(kw_n) < 2:
                continue
            if kw_n in text:
                score += 1
        if score > 0:
            scored.append((score, str(rxn.id)))

    scored.sort(key=lambda x: (-x[0], x[1]))
    out: list[str] = []
    for _, rid in scored:
        if rid not in out:
            out.append(rid)
        if len(out) >= top_k:
            break
    return out


def audit_exchange_ids(model, medium_config: dict[str, Any]) -> list[AuditRow]:
    requested = collect_requested_exchange_ids(medium_config)
    present_ids = set(str(r.id) for r in model.reactions)

    rows: list[AuditRow] = []
    for rid in requested:
        status = "present" if rid in present_ids else "missing"
        suggestions = suggest_exchange_replacements(model, rid, top_k=3) if status == "missing" else []
        s1 = suggestions[0] if len(suggestions) > 0 else ""
        s2 = suggestions[1] if len(suggestions) > 1 else ""
        s3 = suggestions[2] if len(suggestions) > 2 else ""
        rows.append(AuditRow(requested_id=rid, status=status, suggestion_1=s1, suggestion_2=s2, suggestion_3=s3))
    return rows


def write_audit_csv(rows: Iterable[AuditRow], out_path: str | Path) -> Path:
    p = Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["requested_id", "status", "suggestion_1", "suggestion_2", "suggestion_3"],
        )
        w.writeheader()
        for r in rows:
            w.writerow(
                {
                    "requested_id": r.requested_id,
                    "status": r.status,
                    "suggestion_1": r.suggestion_1,
                    "suggestion_2": r.suggestion_2,
                    "suggestion_3": r.suggestion_3,
                }
            )
    return p

