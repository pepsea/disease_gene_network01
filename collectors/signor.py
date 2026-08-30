"""SIGNOR — Signaling Network Open Resource (CC BY 4.0).

Causal signalling interactions (phosphorylation, activation, inhibition, ...).
No API key. The full human dataset is downloaded as TSV and cached on disk,
then filtered for the gene of interest — the same approach as the original
project.
"""
from __future__ import annotations

import logging
import os
import threading
import time
from typing import Any, Optional

from ._cache import cache_enabled, cache_path
from ._http import SESSION

log = logging.getLogger(__name__)

SIGNOR_TSV = os.environ.get(
    "SIGNOR_TSV_URL", "https://signor.uniroma2.it/getData.php?organism=9606&format=tsv"
)
CACHE_TTL = 7 * 24 * 3600  # 7 days
TIMEOUT = (5, 120)

# TSV column order, as documented by SIGNOR.
COLS = [
    "entityA", "typeA", "idA", "dbA",
    "entityB", "typeB", "idB", "dbB",
    "effect", "mechanism", "residue", "sequence",
    "taxId", "cellData", "tissueData", "modA", "modB",
    "pmid", "direct", "sentence_id", "annotated_by",
    "notes", "signor_id", "score",
]

# The original re-scanned the whole TSV for every gene. The rows are indexed by
# gene once per process instead, which is the same filtering with one pass.
_index_lock = threading.Lock()
_index: Optional[dict[str, list[dict[str, Any]]]] = None


def reset_cache() -> None:
    """Drop the in-memory index (used by tests)."""
    global _index
    with _index_lock:
        _index = None


def _get_signor_tsv() -> str:
    """Return the human SIGNOR TSV, from disk cache when it is fresh."""
    if cache_enabled():
        path = cache_path("signor_9606.tsv")
        if path.exists() and time.time() - path.stat().st_mtime < CACHE_TTL:
            return path.read_text(encoding="utf-8")

    log.info("[SIGNOR] downloading TSV (first run or cache older than 7 days)")
    response = SESSION.get(SIGNOR_TSV, timeout=TIMEOUT)
    response.raise_for_status()
    text = response.text
    if cache_enabled():
        cache_path("signor_9606.tsv").write_text(text, encoding="utf-8")
    return text


def _parse(text: str) -> dict[str, list[dict[str, Any]]]:
    """Index SIGNOR rows by the upper-cased symbol of each participant."""
    index: dict[str, list[dict[str, Any]]] = {}

    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 9:
            continue

        entity_a_raw, entity_b_raw = parts[0].strip(), parts[4].strip()
        entity_a, entity_b = entity_a_raw.upper(), entity_b_raw.upper()
        if not entity_a or not entity_b:
            continue

        type_a, type_b = parts[1].strip(), parts[5].strip()
        effect = parts[8].strip()
        mechanism = parts[9].strip() if len(parts) > 9 else ""
        residue = parts[10].strip() if len(parts) > 10 else ""
        pmid = parts[17].strip() if len(parts) > 17 else ""
        score_raw = parts[23].strip() if len(parts) > 23 else ""
        try:
            score = float(score_raw) if score_raw else None
        except ValueError:
            score = None

        # One row yields an entry for whichever side the queried gene is on.
        for gene, gene_raw, partner_raw, partner_type, direction in (
            (entity_a, entity_a_raw, entity_b_raw, type_b, "→"),
            (entity_b, entity_b_raw, entity_a_raw, type_a, "←"),
        ):
            # Keep only protein/complex partners, so that chemicals and
            # phenotypes never enter the analysis.
            if partner_type not in ("protein", "complex"):
                continue
            index.setdefault(gene, []).append(
                {
                    "source": gene_raw if direction == "→" else partner_raw,
                    "target": partner_raw if direction == "→" else gene_raw,
                    "partner": partner_raw,
                    "partner_type": "gene" if partner_type == "protein" else partner_type,
                    "direction": direction,
                    "effect": effect,
                    "mechanism": mechanism,
                    "residue": residue,
                    "pmid": pmid,
                    "score": score,
                    "db": "SIGNOR",
                }
            )

    return index


def _get_index() -> dict[str, list[dict[str, Any]]]:
    global _index
    if _index is not None:
        return _index
    with _index_lock:
        if _index is None:
            _index = _parse(_get_signor_tsv())
        return _index


def get_interactions(gene_symbol: str) -> list[dict[str, Any]]:
    """Return SIGNOR causal interactions involving the gene."""
    return list(_get_index().get((gene_symbol or "").strip().upper(), ()))
