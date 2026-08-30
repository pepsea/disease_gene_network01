"""STRING protein interaction database (CC BY 4.0).

Confidence-scored functional and physical associations. ``required_score`` is
the 0-1000 threshold STRING itself applies (400 = medium, 700 = high).
Per-gene results are cached on disk for three days.
"""
from __future__ import annotations

import json
import logging
import os
import time
from typing import Any

from ._cache import cache_enabled, cache_path
from ._http import SESSION

log = logging.getLogger(__name__)

BASE = os.environ.get("STRING_URL", "https://string-db.org/api/json/network")
CACHE_TTL = 3 * 24 * 3600  # 3 days
TIMEOUT = (5, 25)


def _cache_file(gene_symbol: str, required_score: int):
    return cache_path("string", f"{gene_symbol.upper()}_{required_score}.json")


def _load_cache(gene_symbol: str, required_score: int):
    if not cache_enabled():
        return None
    path = _cache_file(gene_symbol, required_score)
    try:
        if path.exists() and time.time() - path.stat().st_mtime < CACHE_TTL:
            return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return None


def _save_cache(gene_symbol: str, required_score: int, data: list[dict]) -> None:
    if not cache_enabled():
        return
    try:
        _cache_file(gene_symbol, required_score).write_text(
            json.dumps(data, ensure_ascii=False), encoding="utf-8"
        )
    except Exception:
        log.warning("[STRING] could not write cache for %s", gene_symbol)


def get_interactions(
    gene_symbol: str,
    required_score: int = 400,
    species: int = 9606,
    max_results: int = 100,
) -> list[dict[str, Any]]:
    """Return STRING interactions for a gene, best score per partner.

    STRING's network endpoint also returns edges *between* partners, not just
    edges touching the queried gene; those rows are dropped. When a partner
    appears on several rows only its highest-scoring row is kept.
    """
    cached = _load_cache(gene_symbol, required_score)
    if cached is not None:
        return cached

    try:
        response = SESSION.get(
            BASE,
            params={
                "identifiers": gene_symbol,
                "species": species,
                "required_score": required_score,
                "limit": max_results,
                "caller_identity": "nw_overlap_app",
            },
            timeout=TIMEOUT,
        )
        if response.status_code == 404:
            _save_cache(gene_symbol, required_score, [])
            return []
        response.raise_for_status()
        rows = response.json()
    except Exception as exc:
        log.warning("[STRING] %s: %s", gene_symbol, exc)
        return []

    if not isinstance(rows, list):
        return []

    center = gene_symbol.upper()
    best_by_partner: dict[str, dict[str, Any]] = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        a = (row.get("preferredName_A") or "").strip()
        b = (row.get("preferredName_B") or "").strip()

        if center == a.upper():
            partner = b
        elif center == b.upper():
            partner = a
        else:
            continue  # partner-to-partner edge, not ours

        if not partner or partner.upper() == center:
            continue

        try:
            score = float(row.get("score") or 0.0)
        except (TypeError, ValueError):
            score = 0.0

        key = partner.upper()
        existing = best_by_partner.get(key)
        if existing is not None and (existing.get("score") or 0) >= score:
            continue

        best_by_partner[key] = {
            "source": gene_symbol,
            "target": partner,
            "partner": partner,
            "partner_type": "gene",  # STRING carries genes only
            "effect": "",
            "mechanism": "functional/physical association",
            "direction": "—",
            "score": score,  # STRING score is 0-1 here
            "subscores": {
                "experimental": row.get("escore"),
                "database": row.get("dscore"),
                "textmining": row.get("tscore"),
                "coexpression": row.get("ascore"),
                "neighborhood": row.get("nscore"),
                "fusion": row.get("fscore"),
                "cooccurrence": row.get("pscore"),
            },
            "pmid": "",
            "db": "STRING",
        }

    interactions = sorted(
        best_by_partner.values(), key=lambda x: x.get("score") or 0, reverse=True
    )
    _save_cache(gene_symbol, required_score, interactions)
    return interactions
