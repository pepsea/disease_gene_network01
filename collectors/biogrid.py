"""BioGRID PPI collector.

Note: BioGRID's terms restrict use to non-commercial academic research. A free
API key is available at https://webservice.thebiogrid.org/ and is required —
without one this collector returns nothing and the other sources carry the
analysis.
"""
from __future__ import annotations

import logging
import os
from typing import Any, Optional

from ._http import SESSION

log = logging.getLogger(__name__)

BIOGRID_API = os.environ.get(
    "BIOGRID_URL", "https://webservice.thebiogrid.org/interactions/"
)
TIMEOUT = (5, 20)


def get_interactions(
    gene_symbol: str, api_key: Optional[str] = None
) -> list[dict[str, Any]]:
    """Return BioGRID PPIs for ``gene_symbol``.

    Falls back to the ``BIOGRID_API_KEY`` environment variable when no key is
    passed. BioGRID returns one record per experiment/publication, so records
    are collapsed to one entry per partner: the highest score when scores are
    present, otherwise the first record seen.
    """
    key = api_key or os.environ.get("BIOGRID_API_KEY", "")
    if not key:
        log.info("[BioGRID] no API key (BIOGRID_API_KEY); skipping")
        return []

    params = {
        "accessKey": key,
        "geneList": gene_symbol,
        "searchNames": "true",
        "includeHeader": "true",
        "taxId": "9606",
        "interSpeciesExcluded": "true",
        "selfInteractionsExcluded": "true",
        "format": "json",
        "max": 200,
        "start": 0,
    }

    try:
        response = SESSION.get(BIOGRID_API, params=params, timeout=TIMEOUT)
        response.raise_for_status()
        data = response.json()
    except Exception as exc:
        log.warning("[BioGRID] %s: %s", gene_symbol, exc)
        return []

    if not isinstance(data, dict):
        return []

    best_by_partner: dict[str, dict[str, Any]] = {}
    gene_upper = gene_symbol.upper()

    for item in data.values():
        if not isinstance(item, dict):
            continue
        sym_a = str(item.get("OFFICIAL_SYMBOL_A") or "")
        sym_b = str(item.get("OFFICIAL_SYMBOL_B") or "")
        partner = sym_b if sym_a.upper() == gene_upper else sym_a
        if not partner or partner.upper() == gene_upper:
            continue

        score_raw = item.get("SCORE", None)
        try:
            score = float(score_raw) if score_raw not in (None, "", "-") else None
        except (TypeError, ValueError):
            score = None

        key_p = partner.upper()
        existing = best_by_partner.get(key_p)
        if existing is not None:
            existing_score = existing.get("score")
            if existing_score is not None and (score is None or score <= existing_score):
                continue
            if existing_score is None and score is None:
                continue

        best_by_partner[key_p] = {
            "source": sym_a,
            "target": sym_b,
            "partner": partner,
            "partner_type": "gene",
            "direction": "—",  # BioGRID is undirected
            "effect": "physical association",
            "mechanism": str(item.get("EXPERIMENTAL_SYSTEM") or ""),
            "exp_type": str(item.get("EXPERIMENTAL_SYSTEM_TYPE") or ""),
            "pmid": str(item.get("PUBMED_ID", "")),
            "score": score,
            "db": "BioGRID",
        }

    return sorted(
        best_by_partner.values(),
        key=lambda x: x.get("score") if x.get("score") is not None else -1,
        reverse=True,
    )
