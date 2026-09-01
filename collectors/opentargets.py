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
      description
    }
  }
}
"""

# Same query without `description`, used if that field is ever unavailable so
# that disease search keeps working against a changed schema.
SEARCH_QUERY_MINIMAL = """
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

# Set once the richer query is known to fail, so the fallback is not re-probed.
_search_query_supported = True

# Phenotype queries, richest first. `phenotypes` and the shape of its evidence
# rows cannot be verified from this sandbox, so a rejected field falls back
# rather than failing the request.
PHENOTYPE_QUERIES = [
    """
query DiseasePhenotypes($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    name
    phenotypes(page: { index: 0, size: $size }) {
      count
      rows {
        phenotypeHPO { id name description }
        evidence { aspect frequency qualifierNot resource }
      }
    }
  }
}
""",
    """
query DiseasePhenotypes($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    name
    phenotypes(page: { index: 0, size: $size }) {
      count
      rows { phenotypeHPO { id name description } }
    }
  }
}
""",
    """
query DiseasePhenotypes($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    name
    phenotypes(page: { index: 0, size: $size }) {
      count
      rows { phenotypeHPO { id name } }
    }
  }
}
""",
]

# Set once a variant is known to work, so later calls skip the rejected ones.
_phenotype_query_index: Optional[int] = None

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


def is_ontology_id(value: str) -> bool:
    """True when the string already looks like an EFO/MONDO-style ontology id."""
    return bool(_ONTOLOGY_ID_RE.match((value or "").strip()))


def normalise_ontology_id(value: str) -> str:
    """Canonicalise the casing of an ontology id (``efo_1`` -> ``EFO_1``)."""
    return (value or "").strip().upper().replace("ORPHANET", "Orphanet")


def search_diseases(query: str, limit: int = 10) -> list[dict[str, Any]]:
    """Return candidate diseases for a free-text query, best match first.

    Each candidate is ``{"id", "name", "description", "exact"}``, where
    ``exact`` marks a case-insensitive full-name match. Exact matches are
    promoted above the search engine's own ranking, so that "Alzheimer disease"
    does not lose to "Alzheimer disease 2".
    """
    global _search_query_supported

    query = (query or "").strip()
    if not query:
        return []

    if _search_query_supported:
        try:
            data = _post(SEARCH_QUERY, {"q": query, "entity": ["disease"]})
        except OpenTargetsError:
            # Retry once without `description`; if that works, stop trying the
            # richer query for the rest of the process.
            data = _post(SEARCH_QUERY_MINIMAL, {"q": query, "entity": ["disease"]})
            _search_query_supported = False
    else:
        data = _post(SEARCH_QUERY_MINIMAL, {"q": query, "entity": ["disease"]})

    wanted = query.casefold()
    candidates: list[dict[str, Any]] = []
    for hit in (data.get("search") or {}).get("hits") or []:
        if hit.get("entity") != "disease" or not hit.get("id"):
            continue
        name = hit.get("name") or ""
        candidates.append(
            {
                "id": hit["id"],
                "name": name,
                "description": (hit.get("description") or "")[:300],
                "exact": name.casefold() == wanted,
            }
        )

    candidates.sort(key=lambda c: not c["exact"])  # stable: exact matches first
    return candidates[: max(1, int(limit))]


def resolve_disease_id(disease_name: str) -> tuple[Optional[str], str]:
    """Resolve a disease name to ``(ontology_id, label)``.

    Returns ``(None, disease_name)`` when nothing matches. A string that already
    looks like an ontology id is passed straight through; otherwise the best
    search candidate is taken.
    """
    disease_name = (disease_name or "").strip()
    if not disease_name:
        return None, disease_name

    if is_ontology_id(disease_name):
        ontology_id = normalise_ontology_id(disease_name)
        label = get_disease_label(ontology_id)
        return ontology_id, label or disease_name

    candidates = search_diseases(disease_name, limit=1)
    if not candidates:
        return None, disease_name
    return candidates[0]["id"], candidates[0]["name"] or disease_name


def get_disease_label(disease_id: str) -> Optional[str]:
    """Return the display label for an ontology id, or ``None`` if unknown."""
    data = _post(
        "query DiseaseLabel($efoId: String!) { disease(efoId: $efoId) { id name } }",
        {"efoId": disease_id},
    )
    disease = data.get("disease")
    return disease.get("name") if disease else None


def get_disease_with_top_genes(
    disease_id: str, top_n: int = 100
) -> tuple[str, list[dict[str, Any]]]:
    """Return ``(disease_label, top_genes)`` in a single API round trip.

    The associations query already carries the disease name, so callers that
    need both do not have to look the label up separately.
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
    return disease.get("name") or disease_id, genes[:size]


def get_disease_top_genes(disease_id: str, top_n: int = 100) -> list[dict[str, Any]]:
    """Return the top ``top_n`` Open Targets associated genes for a disease.

    Each entry is ``{"symbol": str, "score": float, "target_id": str}``, ordered
    by descending association score (the order Open Targets returns them in).
    """
    return get_disease_with_top_genes(disease_id, top_n)[1]


def hpo_to_ontology_id(hpo_id: str) -> str:
    """``HP:0002354`` -> ``HP_0002354``, the form Open Targets indexes it under."""
    return str(hpo_id or "").strip().replace(":", "_").upper()


def get_disease_phenotypes(disease_id: str, limit: int = 50) -> list[dict[str, Any]]:
    """Return the disease's HPO phenotypes (its symptoms), most relevant first.

    Each entry is ``{"hpo_id", "ontology_id", "name", "description", "frequency",
    "aspect", "resources", "resource", "excluded"}``. ``resources`` lists every
    source annotating the phenotype — Open Targets carries HPO's merge of OMIM,
    Orphanet and DECIPHER, so a symptom can be supported by several. Phenotypes every source marks as
    *absent* in this disease (``qualifierNot``) are flagged ``excluded`` — they
    describe what the disease does not cause, so callers should not use them to
    seed genes.

    Returns ``[]`` when the API exposes no phenotypes for the disease, or when
    the field is unavailable; symptom analysis is optional, never fatal.
    """
    global _phenotype_query_index

    size = max(1, min(int(limit), MAX_PAGE_SIZE))
    order = (
        [_phenotype_query_index]
        if _phenotype_query_index is not None
        else range(len(PHENOTYPE_QUERIES))
    )

    disease = None
    for index in order:
        try:
            data = _post(PHENOTYPE_QUERIES[index], {"efoId": disease_id, "size": size})
        except OpenTargetsError:
            continue
        disease = data.get("disease")
        if disease is None:
            return []
        _phenotype_query_index = index
        break

    if disease is None:
        return []

    rows = (disease.get("phenotypes") or {}).get("rows") or []

    phenotypes: list[dict[str, Any]] = []
    for row in rows:
        hpo = row.get("phenotypeHPO") or {}
        hpo_id = str(hpo.get("id") or "").strip()
        if not hpo_id:
            continue

        evidence = [e for e in (row.get("evidence") or []) if isinstance(e, dict)]
        # Only treat a phenotype as excluded when every source says so.
        excluded = bool(evidence) and all(e.get("qualifierNot") for e in evidence)
        frequency = next((e.get("frequency") for e in evidence if e.get("frequency")), "")
        aspect = next((e.get("aspect") for e in evidence if e.get("aspect")), "")
        # A phenotype is often annotated by several sources (Open Targets
        # carries HPO's merge of OMIM, Orphanet and DECIPHER), so keep them all
        # rather than the first — Orphanet in particular is where most of the
        # frequency annotation comes from.
        resources = list(
            dict.fromkeys(
                str(e.get("resource")).strip()
                for e in evidence
                if e.get("resource")
            )
        )

        phenotypes.append(
            {
                "hpo_id": hpo_id,
                "ontology_id": hpo_to_ontology_id(hpo_id),
                "name": hpo.get("name") or hpo_id,
                "description": (hpo.get("description") or "")[:300],
                "frequency": frequency or "",
                "aspect": aspect or "",
                "resources": resources,
                # First source, kept for callers that want a single label.
                "resource": resources[0] if resources else "",
                "excluded": excluded,
            }
        )
    return phenotypes
