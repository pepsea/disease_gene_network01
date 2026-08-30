"""g:Profiler — pathway/GO enrichment for a gene list (BSD 2-Clause, no key).

Ported from the original project's ``gprofiler.py``. KEGG and other sources
with commercial licence restrictions are deliberately excluded; only Reactome,
GO:BP and WikiPathways are queried by default.

Results are cached on disk per (gene list, sources) so that re-running the same
analysis, or two genes with the same neighbourhood, does not re-query.
"""
from __future__ import annotations

import hashlib
import json
import logging
import os
import time
from typing import Any, Optional

from ._cache import cache_enabled, cache_path
from ._http import SESSION

log = logging.getLogger(__name__)

GPROFILER_API = os.environ.get(
    "GPROFILER_URL", "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
)
TIMEOUT = (5, 60)
CACHE_TTL = 3 * 24 * 3600  # 3 days

# KEGG is excluded: its licence does not allow commercial use.
DEFAULT_SOURCES = ("REAC", "GO:BP", "WP")


def _cache_key(genes: list[str], sources: tuple[str, ...], threshold: float) -> str:
    payload = json.dumps(
        {"g": sorted(g.upper() for g in genes), "s": list(sources), "t": threshold},
        sort_keys=True,
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:32]


def _load_cache(key: str) -> Optional[list[dict]]:
    if not cache_enabled():
        return None
    path = cache_path("gprofiler", f"{key}.json")
    try:
        if path.exists() and time.time() - path.stat().st_mtime < CACHE_TTL:
            return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return None


def _save_cache(key: str, data: list[dict]) -> None:
    if not cache_enabled():
        return
    try:
        cache_path("gprofiler", f"{key}.json").write_text(
            json.dumps(data, ensure_ascii=False), encoding="utf-8"
        )
    except Exception:
        log.warning("[g:Profiler] could not write cache")


def enrich_gene_list(
    gene_symbols: list[str],
    organism: str = "hsapiens",
    sources: tuple[str, ...] | list[str] | None = None,
    significance_threshold: float = 0.05,
    max_results: int = 50,
) -> list[dict[str, Any]]:
    """Run g:Profiler enrichment on a gene list, best p-value first.

    Returns ``[{"source", "term_id", "name", "p_value", "intersection_size",
    "term_size", "genes"}, ...]`` where ``genes`` are the query genes annotated
    to that term.
    """
    genes = [g.strip() for g in (gene_symbols or []) if g and g.strip()]
    if not genes:
        return []

    sources = tuple(sources) if sources else DEFAULT_SOURCES

    key = _cache_key(genes, sources, significance_threshold)
    cached = _load_cache(key)
    if cached is not None:
        return cached[:max_results]

    try:
        response = SESSION.post(
            GPROFILER_API,
            json={
                "organism": organism,
                "query": genes,
                "sources": list(sources),
                "user_threshold": significance_threshold,
                "ordered": False,
                "all_results": False,
                # `no_evidences` is omitted on purpose: without it g:Profiler
                # returns `intersections`, which is how the hit genes are found.
            },
            timeout=TIMEOUT,
        )
        response.raise_for_status()
        data = response.json()
    except Exception as exc:
        log.warning("[g:Profiler] %s", exc)
        return []

    if not isinstance(data, dict):
        return []

    raw_results = data.get("result") or []
    meta = data.get("meta") or {}
    # The query gene order that `intersections` indexes into.
    queries = (meta.get("query_metadata") or {}).get("queries") or {}
    ordered_genes = queries.get("query_1") or genes

    results: list[dict[str, Any]] = []
    for item in raw_results:
        if not isinstance(item, dict):
            continue
        # intersections[i] non-empty => the i-th query gene is annotated here.
        raw_ix = item.get("intersections") or []
        hit_genes = [
            ordered_genes[i]
            for i, annotation in enumerate(raw_ix)
            if annotation and i < len(ordered_genes)
        ]
        try:
            p_value = float(item.get("p_value", 1.0))
        except (TypeError, ValueError):
            p_value = 1.0

        results.append(
            {
                "source": item.get("source", ""),
                "term_id": item.get("native", ""),
                "name": item.get("name", ""),
                "p_value": p_value,
                "intersection_size": item.get("intersection_size", 0),
                "term_size": item.get("term_size", 0),
                "genes": hit_genes,
            }
        )

    results.sort(key=lambda x: x["p_value"])
    _save_cache(key, results)
    return results[:max_results]
