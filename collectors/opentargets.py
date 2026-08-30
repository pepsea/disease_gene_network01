"""Open Targets Platform collector.

Two things are needed from Open Targets:

* ``resolve_disease_id`` — turn a free-text disease name into an EFO/MONDO id.
* ``get_disease_top_genes`` — the top scoring target genes for that disease.

Both go through the public Open Targets Platform GraphQL API (v4), which needs
no API key.
"""
from __future__ import annotations

import os
import re
from typing import Any, Optional

from ._http import SESSION

OT_API = os.environ.get(
    "OT_API_URL", "https://api.platform.opentargets.org/api/v4/graphql"
)

# Timeouts (connect, read) in seconds.
TIMEOUT = (5, 30)

# Open Targets caps a single associations page; keep requests inside that.
MAX_PAGE_SIZE = 500

# An ontology id the user may paste directly instead of a disease name,
# e.g. "EFO_0000249" or "MONDO_0004975".
_ONTOLOGY_ID_RE = re.compile(r"^(EFO|MONDO|HP|DOID|Orphanet|NCIT|GO|MP|OTAR)_\d+$", re.I)

SEARCH_QUERY = """
query SearchEntity($q: String!, $entity: [String!]) {
  search(queryString: $q, entityNames: $entity) {
    hits {
      id
      name
      entity
      score
    }
  }
}
"""

DISEASE_TOP_GENES_QUERY = """
query DiseaseTopGenes($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    name
    associatedTargets(page: { index: 0, size: $size }) {
      count
      rows {
        target {
          id
          approvedSymbol
        }
        score
      }
    }
  }
}
"""


class OpenTargetsError(RuntimeError):
    """Raised when the Open Targets API cannot be queried or returns errors."""


def _post(query: str, variables: dict[str, Any]) -> dict[str, Any]:
    """Run a GraphQL query and return its ``data`` payload."""
    try:
        response = SESSION.post(
            OT_API, json={"query": query, "variables": variables}, timeout=TIMEOUT
        )
        response.raise_for_status()
        payload = response.json()
    except ValueError as exc:  # malformed JSON
        raise OpenTargetsError(f"Open Targets returned a non-JSON response: {exc}") from exc
    except Exception as exc:  # network / HTTP error
        raise OpenTargetsError(f"Open Targets request failed: {exc}") from exc

    if payload.get("errors"):
        messages = "; ".join(
            str(e.get("message", e)) for e in payload["errors"] if isinstance(e, dict)
        )
        raise OpenTargetsError(f"Open Targets GraphQL error: {messages or payload['errors']}")

    return payload.get("data") or {}


def resolve_disease_id(disease_name: str) -> tuple[Optional[str], str]:
    """Resolve a disease name to ``(ontology_id, label)``.

    Returns ``(None, disease_name)`` when nothing matches. An exact
    case-insensitive name match wins over the search engine's own ranking; a
    string that already looks like an ontology id is passed straight through.
    """
    disease_name = (disease_name or "").strip()
    if not disease_name:
        return None, disease_name

    if _ONTOLOGY_ID_RE.match(disease_name):
        ontology_id = disease_name.upper().replace("ORPHANET", "Orphanet")
        label = get_disease_label(ontology_id)
        return ontology_id, label or disease_name

    data = _post(SEARCH_QUERY, {"q": disease_name, "entity": ["disease"]})
    hits = [
        hit
        for hit in (data.get("search") or {}).get("hits") or []
        if hit.get("entity") == "disease" and hit.get("id")
    ]
    if not hits:
        return None, disease_name

    wanted = disease_name.casefold()
    exact = [h for h in hits if (h.get("name") or "").casefold() == wanted]
    best = exact[0] if exact else hits[0]
    return best["id"], best.get("name") or disease_name


def get_disease_label(disease_id: str) -> Optional[str]:
    """Return the display label for an ontology id, or ``None`` if unknown."""
    data = _post(
        "query DiseaseLabel($efoId: String!) { disease(efoId: $efoId) { id name } }",
        {"efoId": disease_id},
    )
    disease = data.get("disease")
    return disease.get("name") if disease else None


def get_disease_top_genes(disease_id: str, top_n: int = 100) -> list[dict[str, Any]]:
    """Return the top ``top_n`` Open Targets associated genes for a disease.

    Each entry is ``{"symbol": str, "score": float, "target_id": str}``, ordered
    by descending association score (the order Open Targets returns them in).
    """
    size = max(1, min(int(top_n), MAX_PAGE_SIZE))
    data = _post(DISEASE_TOP_GENES_QUERY, {"efoId": disease_id, "size": size})

    disease = data.get("disease")
    if not disease:
        raise OpenTargetsError(f"Unknown disease id: {disease_id}")

    rows = (disease.get("associatedTargets") or {}).get("rows") or []

    genes: list[dict[str, Any]] = []
    for row in rows:
        target = row.get("target") or {}
        symbol = (target.get("approvedSymbol") or "").strip()
        if not symbol:
            continue
        genes.append(
            {
                "symbol": symbol,
                "score": float(row.get("score") or 0.0),
                "target_id": target.get("id", ""),
            }
        )
    return genes[:size]
