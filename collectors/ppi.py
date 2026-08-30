"""Protein-protein interaction (PPI) partner collector.

Partners for a gene are pooled from up to three public sources:

* **SIGNOR**  — curated signalling interactions (no key required).
* **STRING**  — functional association network (no key required).
* **BioGRID** — curated physical/genetic interactions (API key required).

Every source is best-effort: if one is unreachable the others still produce a
result, and a gene with no partners anywhere simply scores zero downstream.

Symbols are normalised to upper case throughout so that the self-gene is
reliably excluded and so that the overlap comparison in
:mod:`nw_overlap` (which upper-cases both sides) behaves consistently.
"""
from __future__ import annotations

import os
import threading
from typing import Any, Iterable, Optional

from ._http import SESSION

SIGNOR_URL = os.environ.get("SIGNOR_URL", "https://signor.uniroma2.it/getData.php")
STRING_URL = os.environ.get("STRING_URL", "https://string-db.org/api/json/network")
BIOGRID_URL = os.environ.get("BIOGRID_URL", "https://webservice.thebiogrid.org/interactions/")

TIMEOUT = (5, 20)
# SIGNOR is fetched as one bulk download, so it gets a longer read budget.
SIGNOR_TIMEOUT = (5, 120)

HUMAN_TAX_ID = 9606

# SIGNOR has no per-gene endpoint, so its human interactome is downloaded once
# per process and indexed in memory. Without this the bulk file would be
# re-fetched for every gene in every request.
_signor_lock = threading.Lock()
_signor_index: Optional[dict[str, set[str]]] = None


def _enabled(env_name: str, default: bool = True) -> bool:
    raw = os.environ.get(env_name)
    if raw is None:
        return default
    return raw.strip().lower() not in {"0", "false", "no", "off", ""}


def reset_signor_cache() -> None:
    """Drop the cached SIGNOR index (used by tests and by long-lived servers)."""
    global _signor_index
    with _signor_lock:
        _signor_index = None


def _iter_rows(payload: Any) -> Iterable[dict]:
    """Yield interaction rows from a list-of-rows or dict-of-rows payload."""
    if isinstance(payload, list):
        rows = payload
    elif isinstance(payload, dict):
        rows = payload.values()
    else:
        return
    for row in rows:
        if isinstance(row, dict):
            yield row


def _load_signor_index() -> dict[str, set[str]]:
    """Fetch and index the human SIGNOR interactome, caching the result.

    A failed download is cached as an empty index so that one unreachable
    source does not stall every gene in the request.
    """
    global _signor_index
    if _signor_index is not None:
        return _signor_index

    with _signor_lock:
        if _signor_index is not None:  # another thread won the race
            return _signor_index

        index: dict[str, set[str]] = {}
        try:
            response = SESSION.get(
                SIGNOR_URL,
                params={"organism": str(HUMAN_TAX_ID), "format": "json"},
                timeout=SIGNOR_TIMEOUT,
            )
            response.raise_for_status()
            for row in _iter_rows(response.json()):
                a = str(row.get("ENTITYA") or "").strip().upper()
                b = str(row.get("ENTITYB") or "").strip().upper()
                if not a or not b or a == b:
                    continue
                index.setdefault(a, set()).add(b)
                index.setdefault(b, set()).add(a)
        except Exception:
            index = {}

        _signor_index = index
        return _signor_index


def _signor_partners(gene: str) -> set[str]:
    if not _enabled("NW_ENABLE_SIGNOR"):
        return set()
    return set(_load_signor_index().get(gene.upper(), ()))


def _string_partners(gene: str, top_n: int) -> dict[str, float]:
    """Return ``{partner: combined_score}`` from the STRING network endpoint."""
    if not _enabled("NW_ENABLE_STRING"):
        return {}

    partners: dict[str, float] = {}
    try:
        response = SESSION.get(
            STRING_URL,
            params={
                "identifiers": gene,
                "species": HUMAN_TAX_ID,
                "limit": top_n,
                "caller_identity": "nw_overlap_app",
            },
            timeout=TIMEOUT,
        )
        response.raise_for_status()
        for row in _iter_rows(response.json()):
            a = str(row.get("preferredName_A") or "").strip().upper()
            b = str(row.get("preferredName_B") or "").strip().upper()
            try:
                score = float(row.get("score") or 0.0)
            except (TypeError, ValueError):
                score = 0.0
            target = gene.upper()
            partner = b if a == target else a if b == target else ""
            if partner and partner != target:
                partners[partner] = max(partners.get(partner, 0.0), score)
    except Exception:
        return partners
    return partners


def _biogrid_partners(gene: str, biogrid_key: str, top_n: int) -> set[str]:
    if not biogrid_key:
        return set()

    partners: set[str] = set()
    try:
        response = SESSION.get(
            BIOGRID_URL,
            params={
                "geneList": gene,
                "taxId": HUMAN_TAX_ID,
                "accesskey": biogrid_key,
                "format": "json",
                "max": top_n,
                "interSpeciesExcluded": "true",
                "searchNames": "true",
            },
            timeout=TIMEOUT,
        )
        response.raise_for_status()
        for row in _iter_rows(response.json()):
            a = str(row.get("OFFICIAL_SYMBOL_A") or "").strip().upper()
            b = str(row.get("OFFICIAL_SYMBOL_B") or "").strip().upper()
            target = gene.upper()
            partner = b if a == target else a if b == target else ""
            if partner and partner != target:
                partners.add(partner)
    except Exception:
        return partners
    return partners


def get_ppi_partners(gene: str, biogrid_key: str = "", top_n: int = 30) -> list[str]:
    """Collect PPI partners of ``gene`` and return the top ``top_n`` symbols.

    When more partners are found than ``top_n``, they are ranked by how many of
    the three sources support the interaction, then by STRING combined score,
    then alphabetically. Ranking (rather than an arbitrary slice of a set)
    keeps the result deterministic and keeps the best-supported interactions
    when the list has to be cut.
    """
    gene = (gene or "").strip()
    if not gene:
        return []

    top_n = max(1, int(top_n))
    target = gene.upper()

    signor = _signor_partners(gene)
    string_scores = _string_partners(gene, top_n)
    biogrid = _biogrid_partners(gene, biogrid_key, top_n)

    support: dict[str, int] = {}
    for source in (signor, set(string_scores), biogrid):
        for partner in source:
            support[partner] = support.get(partner, 0) + 1

    support.pop(target, None)

    ranked = sorted(
        support,
        key=lambda p: (-support[p], -string_scores.get(p, 0.0), p),
    )
    return ranked[:top_n]
