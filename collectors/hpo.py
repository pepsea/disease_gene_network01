"""HPO (Human Phenotype Ontology) direct collector.

Open Targets carries HPO's disease-phenotype annotations, but it does not
always expose them, and its symptom-to-gene step depends on each HP term also
being indexed as a disease. This reads HPO's own release files instead, which
removes both dependencies:

* ``phenotype.hpoa``          — disease -> phenotype, with frequency, aspect and
                                the NOT qualifier. This is the file Open Targets
                                itself ingests (OMIM + Orphanet + DECIPHER).
* ``phenotype_to_genes.txt``  — phenotype -> genes, curated, so a symptom's genes
                                come straight from HPO.

Both are downloaded once and indexed in memory, cached on disk like the SIGNOR
bulk file. Parsing is driven by the header line, so a column reordering in a
future release does not silently shift fields.

HPO keys diseases on OMIM / ORPHA / DECIPHER ids, so callers pass the cross
references for the disease (see ``opentargets.get_disease_xrefs``).
"""
from __future__ import annotations

import logging
import os
import re
import threading
import time
from typing import Any, Iterable, Optional

from ._cache import cache_enabled, cache_path
from ._http import SESSION

log = logging.getLogger(__name__)

HPOA_URL = os.environ.get(
    "HPO_HPOA_URL", "https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa"
)
PHENOTYPE_TO_GENES_URL = os.environ.get(
    "HPO_PHENOTYPE_TO_GENES_URL",
    "https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt",
)

CACHE_TTL = 7 * 24 * 3600  # 7 days; HPO releases monthly
TIMEOUT = (5, 180)

# database_id prefix -> the label used elsewhere in the app.
SOURCE_LABELS = {"OMIM": "OMIM", "ORPHA": "ORPHANET", "ORPHANET": "ORPHANET",
                 "DECIPHER": "DECIPHER"}

_lock = threading.Lock()
_hpoa_index: Optional[dict[str, list[dict[str, Any]]]] = None
_gene_index: Optional[dict[str, list[str]]] = None
_term_names: Optional[dict[str, str]] = None


def reset_cache() -> None:
    """Drop the in-memory indexes (used by tests)."""
    global _hpoa_index, _gene_index, _term_names
    with _lock:
        _hpoa_index = None
        _gene_index = None
        _term_names = None


def _fetch(url: str, cache_name: str) -> str:
    """Download a release file, serving a fresh disk cache when there is one."""
    if cache_enabled():
        path = cache_path("hpo", cache_name)
        if path.exists() and time.time() - path.stat().st_mtime < CACHE_TTL:
            return path.read_text(encoding="utf-8")

    log.info("[HPO] downloading %s", cache_name)
    response = SESSION.get(url, timeout=TIMEOUT)
    response.raise_for_status()
    text = response.text
    if cache_enabled():
        try:
            cache_path("hpo", cache_name).write_text(text, encoding="utf-8")
        except Exception:
            log.warning("[HPO] could not cache %s", cache_name)
    return text


# HPO has used several spellings for the same columns across releases, and the
# gene file's header is prefixed "#Format:".
_COLUMN_ALIASES = {
    "hpo_label": "hpo_name",
    "hpo_term_id": "hpo_id",
    "hpo_term_name": "hpo_name",
    "entrez_gene_id": "ncbi_gene_id",
    "entrez_gene_symbol": "gene_symbol",
    "gene_id": "ncbi_gene_id",
    "gene_symbol_": "gene_symbol",
    "disease_id_for_link": "disease_id",
    "diseaseid": "disease_id",
    "databaseid": "database_id",
    "additional_info_from_g_d_source": "disease_source_info",
    "g_d_source": "source",
    "frequency_hpo": "frequency",
    "frequency_raw": "frequency_raw",
}


def _normalise_column(name: str) -> str:
    """Canonicalise one header cell."""
    cleaned = name.strip().lstrip("#").strip().lower()
    # "#Format: HPO-id" -> "hpo-id"
    if cleaned.startswith("format:"):
        cleaned = cleaned.split(":", 1)[1]
    cleaned = re.sub(r"[^a-z0-9]+", "_", cleaned).strip("_")
    return _COLUMN_ALIASES.get(cleaned, cleaned)


def _rows(text: str, expected: Iterable[str]) -> Iterable[dict[str, str]]:
    """Yield rows as dicts, using the header line to name the columns.

    HPO files start with ``#``-prefixed metadata; the header is the last such
    line before the data, or a plain first line. Columns are matched by name so
    that a reordered release still parses; if no usable header is found the
    caller's expected order is used positionally.
    """
    expected = list(expected)
    header: Optional[list[str]] = None

    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            candidate = [_normalise_column(c) for c in line.lstrip("#").split("\t")]
            # A metadata line has one field; a header names several columns.
            if len(candidate) > 1 and any(c in candidate for c in expected):
                header = candidate
            continue

        parts = line.split("\t")
        if header is None:
            candidate = [_normalise_column(c) for c in parts]
            # A header line names columns rather than carrying ids.
            if any(c in candidate for c in expected):
                header = candidate
                continue
            header = expected

        row = {name: (parts[i].strip() if i < len(parts) else "")
               for i, name in enumerate(header)}
        yield row


def _load_hpoa() -> dict[str, list[dict[str, Any]]]:
    """Index ``phenotype.hpoa`` by disease id."""
    global _hpoa_index
    if _hpoa_index is not None:
        return _hpoa_index

    with _lock:
        if _hpoa_index is not None:
            return _hpoa_index

        index: dict[str, list[dict[str, Any]]] = {}
        try:
            text = _fetch(HPOA_URL, "phenotype.hpoa")
        except Exception as exc:
            log.warning("[HPO] phenotype.hpoa unavailable: %s", exc)
            _hpoa_index = {}
            return _hpoa_index

        expected = ["database_id", "disease_name", "qualifier", "hpo_id", "reference",
                    "evidence", "onset", "frequency", "sex", "modifier", "aspect",
                    "biocuration"]
        for row in _rows(text, expected):
            disease_id = (row.get("database_id") or "").strip()
            hpo_id = (row.get("hpo_id") or "").strip()
            if not disease_id or not hpo_id.upper().startswith("HP:"):
                continue
            index.setdefault(disease_id.upper(), []).append(
                {
                    "hpo_id": hpo_id,
                    "disease_name": row.get("disease_name", ""),
                    "qualifier": (row.get("qualifier") or "").strip().upper(),
                    "frequency": (row.get("frequency") or "").strip(),
                    "aspect": (row.get("aspect") or "").strip(),
                    "reference": row.get("reference", ""),
                }
            )

        _hpoa_index = index
        return _hpoa_index


def _load_genes() -> tuple[dict[str, list[str]], dict[str, str]]:
    """Index ``phenotype_to_genes.txt`` by HPO id, and collect term labels."""
    global _gene_index, _term_names
    if _gene_index is not None and _term_names is not None:
        return _gene_index, _term_names

    with _lock:
        if _gene_index is not None and _term_names is not None:
            return _gene_index, _term_names

        genes: dict[str, list[str]] = {}
        names: dict[str, str] = {}
        try:
            text = _fetch(PHENOTYPE_TO_GENES_URL, "phenotype_to_genes.txt")
        except Exception as exc:
            log.warning("[HPO] phenotype_to_genes.txt unavailable: %s", exc)
            _gene_index, _term_names = {}, {}
            return _gene_index, _term_names

        expected = ["hpo_id", "hpo_name", "ncbi_gene_id", "gene_symbol",
                    "disease_source_info", "source", "disease_id"]
        seen: dict[str, set[str]] = {}
        for row in _rows(text, expected):
            hpo_id = (row.get("hpo_id") or "").strip().upper()
            symbol = (row.get("gene_symbol") or "").strip()
            if not hpo_id.startswith("HP:") or not symbol:
                continue
            name = (row.get("hpo_name") or "").strip()
            if name and hpo_id not in names:
                names[hpo_id] = name
            bucket = seen.setdefault(hpo_id, set())
            if symbol.upper() in bucket:
                continue
            bucket.add(symbol.upper())
            genes.setdefault(hpo_id, []).append(symbol)

        _gene_index, _term_names = genes, names
        return _gene_index, _term_names


_ORPHANET_PREFIX = re.compile(r"^(ORPHANET|ORPHA|ORPHACODE)[:_]", re.I)


def normalise_disease_id(disease_id: str) -> str:
    """Canonicalise a disease id to the form ``phenotype.hpoa`` uses.

    Orphanet ids appear as ``Orphanet:1020``, ``ORPHA:1020`` or ``ORPHAcode``
    depending on the source; HPO writes them as ``ORPHA:1020``.
    """
    value = str(disease_id or "").strip()
    if not value:
        return ""
    value = _ORPHANET_PREFIX.sub("ORPHA:", value)
    return value.upper()


def _source_of(disease_id: str) -> str:
    prefix = disease_id.split(":", 1)[0].upper()
    return SOURCE_LABELS.get(prefix, prefix)


def get_disease_phenotypes(disease_ids: list[str], limit: int = 50) -> list[dict[str, Any]]:
    """Return HPO phenotypes for one or more OMIM / ORPHA / DECIPHER ids.

    The result matches the shape of
    :func:`collectors.opentargets.get_disease_phenotypes`, so it is a drop-in
    replacement. When several ids annotate the same phenotype their sources are
    merged; a phenotype is ``excluded`` only when every annotation carries the
    NOT qualifier.

    Only phenotypic abnormalities (aspect ``P``) are returned — inheritance and
    clinical-course terms are not symptoms.
    """
    index = _load_hpoa()
    if not index:
        return []

    _, term_names = _load_genes()

    merged: dict[str, dict[str, Any]] = {}
    for disease_id in disease_ids or []:
        key = normalise_disease_id(disease_id)
        if not key:
            continue
        for row in index.get(key, ()):
            hpo_id = row["hpo_id"]
            aspect = row.get("aspect", "")
            if aspect and aspect.upper() != "P":
                continue

            entry = merged.setdefault(
                hpo_id,
                {
                    "hpo_id": hpo_id,
                    "ontology_id": hpo_id.replace(":", "_").upper(),
                    "name": term_names.get(hpo_id.upper(), hpo_id),
                    "description": "",
                    "frequency": "",
                    "frequency_source": "",
                    "aspect": aspect or "P",
                    "resources": [],
                    "resource": "",
                    "excluded": True,      # tightened below by any non-NOT row
                    "_annotations": 0,
                },
            )
            entry["_annotations"] += 1
            source = _source_of(key)
            if source and source not in entry["resources"]:
                entry["resources"].append(source)
            if row.get("qualifier") != "NOT":
                entry["excluded"] = False

            frequency = (row.get("frequency") or "").strip()
            if frequency:
                # Orphanet is where the HPO frequency terms come from, so its
                # value wins over one picked up from another source.
                if not entry["frequency"] or (
                    source == "ORPHANET" and entry["frequency_source"] != "ORPHANET"
                ):
                    entry["frequency"] = frequency
                    entry["frequency_source"] = source

    phenotypes = list(merged.values())
    # Most-annotated first, as a stand-in for Open Targets' relevance ordering.
    phenotypes.sort(key=lambda p: (-p["_annotations"], p["hpo_id"]))
    for entry in phenotypes:
        entry.pop("_annotations", None)
        entry["resource"] = entry["resources"][0] if entry["resources"] else ""

    return phenotypes[: max(1, int(limit))]


def get_phenotype_genes(hpo_id: str, top_n: int = 50) -> list[dict[str, Any]]:
    """Return the genes HPO curates for a phenotype.

    HPO's phenotype-gene links are curated and unweighted, so every gene is
    given a score of 1.0: for an HPO-sourced symptom the weighted overlap
    reduces to the unweighted one, which is the honest reading of a set with no
    confidence values.
    """
    genes, _ = _load_genes()
    key = str(hpo_id or "").strip().upper().replace("_", ":")
    return [
        {"symbol": symbol, "score": 1.0, "target_id": ""}
        for symbol in genes.get(key, ())[: max(1, int(top_n))]
    ]


def list_diseases() -> list[dict[str, Any]]:
    """Every disease HPO annotates, as ``{"id", "name", "source", "phenotype_count"}``.

    HPO keeps its own disease registry — OMIM, Orphanet and DECIPHER entries with
    their own names — so a disease can be picked here directly, without going
    through Open Targets and its cross references.
    """
    index = _load_hpoa()
    diseases: list[dict[str, Any]] = []
    for disease_id, rows in index.items():
        if not rows:
            continue
        name = (rows[0].get("disease_name") or "").strip()
        diseases.append(
            {
                "id": disease_id,
                "name": name or disease_id,
                "source": _source_of(disease_id),
                # Phenotypic abnormalities only, matching what a search result
                # would actually yield as symptoms.
                "phenotype_count": sum(
                    1 for r in rows
                    if (not r.get("aspect") or r["aspect"].upper() == "P")
                    and r.get("qualifier") != "NOT"
                ),
            }
        )
    return diseases


def search_diseases(query: str, limit: int = 10) -> list[dict[str, Any]]:
    """Search HPO's disease registry by name.

    Exact matches rank first, then names starting with the query, then names
    containing it; within a tier the better-annotated disease comes first. The
    tiering matters: across thousands of rare diseases a bare substring match
    surfaces the wrong one, but ranked and shown as candidates for the user to
    choose from, partial matching is exactly what is wanted.
    """
    wanted = (query or "").strip().casefold()
    if not wanted:
        return []

    scored: list[tuple[int, int, str, dict[str, Any]]] = []
    for disease in list_diseases():
        name = disease["name"].casefold()
        if name == wanted:
            tier = 0
        elif name.startswith(wanted):
            tier = 1
        elif wanted in name:
            tier = 2
        else:
            continue
        scored.append((tier, -disease["phenotype_count"], disease["name"], disease))

    scored.sort(key=lambda row: row[:3])
    return [
        {**disease, "exact": tier == 0}
        for tier, _, _, disease in scored[: max(1, int(limit))]
    ]


def find_disease_ids_by_name(name: str, limit: int = 5) -> list[str]:
    """Look up OMIM / ORPHA ids by an exact disease name.

    Used for the automatic fallback, where nobody is present to choose between
    candidates; :func:`search_diseases` is the interactive counterpart.
    """
    wanted = (name or "").strip().casefold()
    if not wanted:
        return []

    index = _load_hpoa()
    hits: list[str] = []
    for disease_id, rows in index.items():
        if rows and (rows[0].get("disease_name") or "").strip().casefold() == wanted:
            hits.append(disease_id)
            if len(hits) >= limit:
                break
    return hits
