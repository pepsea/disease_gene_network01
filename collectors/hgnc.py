"""HGNC (HUGO Gene Nomenclature Committee) symbol validation.

Input gene symbols are checked against the HGNC REST API so that the analysis
runs on approved symbols. Three outcomes matter:

* the symbol is an **approved** HGNC symbol — used as-is;
* it is an **alias** or a **previous** symbol — mapped to the approved one, so
  that e.g. ``PS1`` is analysed as ``PSEN1``;
* it matches nothing — reported to the user, but still analysed, since a very
  new symbol may simply not be in the cached view yet.

Validation is best-effort: if HGNC is unreachable the symbol is marked
``unverified`` and the analysis proceeds with what the user typed, rather than
failing the whole request.
"""
from __future__ import annotations

import os
import threading
from typing import Any, Optional
from urllib.parse import quote

from ._http import SESSION

HGNC_API = os.environ.get("HGNC_API_URL", "https://rest.genenames.org")
TIMEOUT = (5, 20)

# Symbols repeat heavily across requests, so resolutions are memoised for the
# life of the process. Only definitive answers are cached; a network failure is
# left uncached so the next request can retry.
_cache_lock = threading.Lock()
_cache: dict[str, dict[str, Any]] = {}

# Which HGNC field matched, in the order they are tried.
_LOOKUP_ORDER = (
    ("symbol", "approved"),
    ("alias_symbol", "alias"),
    ("prev_symbol", "previous"),
)


def reset_cache() -> None:
    """Clear the memoised symbol resolutions (used by tests)."""
    with _cache_lock:
        _cache.clear()


def _fetch(field: str, value: str) -> Optional[list[dict[str, Any]]]:
    """Return HGNC docs for ``field=value``, or ``None`` if HGNC is unreachable."""
    url = f"{HGNC_API}/fetch/{field}/{quote(value, safe='')}"
    try:
        response = SESSION.get(url, headers={"Accept": "application/json"}, timeout=TIMEOUT)
        response.raise_for_status()
        payload = response.json()
    except Exception:
        return None
    if not isinstance(payload, dict):
        return None
    docs = (payload.get("response") or {}).get("docs")
    return docs if isinstance(docs, list) else []


def _entry(doc: dict[str, Any], raw: str, status: str) -> dict[str, Any]:
    approved = str(doc.get("symbol") or "").strip() or raw
    withdrawn = str(doc.get("status") or "").lower().startswith("entry withdrawn")
    return {
        "input": raw,
        "symbol": approved,
        "status": "withdrawn" if withdrawn else status,
        "hgnc_id": doc.get("hgnc_id", ""),
        "name": doc.get("name", ""),
        "ensembl_gene_id": doc.get("ensembl_gene_id", ""),
    }


def resolve_gene_symbol(symbol: str) -> dict[str, Any]:
    """Resolve one symbol against HGNC.

    Returns a dict with ``input``, ``symbol`` (the approved symbol to analyse),
    ``status`` (``approved`` / ``alias`` / ``previous`` / ``withdrawn`` /
    ``unknown`` / ``unverified``), ``hgnc_id``, ``name`` and, when an alias is
    ambiguous, ``candidates``.
    """
    raw = (symbol or "").strip()
    if not raw:
        return {"input": raw, "symbol": raw, "status": "unknown", "hgnc_id": "", "name": ""}

    key = raw.upper()
    with _cache_lock:
        cached = _cache.get(key)
    if cached is not None:
        return {**cached, "input": raw}

    unreachable = False
    for field, status in _LOOKUP_ORDER:
        docs = _fetch(field, raw)
        if docs is None:
            unreachable = True
            continue
        if not docs:
            continue

        result = _entry(docs[0], raw, status)
        if len(docs) > 1:
            # An alias can point at several genes; surface the ambiguity rather
            # than silently picking one.
            result["candidates"] = [
                str(d.get("symbol") or "") for d in docs if d.get("symbol")
            ]
        with _cache_lock:
            _cache[key] = result
        return result

    if unreachable:
        # HGNC could not be reached; do not cache, and do not block the analysis.
        return {"input": raw, "symbol": raw, "status": "unverified", "hgnc_id": "", "name": ""}

    result = {"input": raw, "symbol": raw, "status": "unknown", "hgnc_id": "", "name": ""}
    with _cache_lock:
        _cache[key] = result
    return result


def resolve_gene_symbols(symbols: list[str]) -> list[dict[str, Any]]:
    """Resolve a list of symbols, dropping entries that map to the same gene.

    Input order is preserved, and the first spelling of a duplicate wins.
    """
    resolved: list[dict[str, Any]] = []
    seen: set[str] = set()
    for symbol in symbols:
        entry = resolve_gene_symbol(symbol)
        key = entry["symbol"].upper()
        if key in seen:
            continue
        seen.add(key)
        resolved.append(entry)
    return resolved
