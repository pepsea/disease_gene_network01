"""Build a gene set from a disease's symptoms.

The disease-level gene list in table 1 answers "which genes are associated with
this disease?". This answers a narrower question: "which genes drive the
symptoms this disease presents with?"

Each HPO phenotype Open Targets holds for the disease is itself an entry in the
Open Targets disease index (``HP:0002354`` is indexed as ``HP_0002354``), so its
associated targets are read with the same associations query used for the
disease. The per-phenotype gene lists are then pooled into one symptom-derived
gene set, which is scored against a gene's PPI neighbourhood exactly like the
disease gene list.

A gene's score is the **sum** of its association scores across the symptoms it
is linked to, so both strength and breadth count: a gene behind five of the
disease's symptoms outweighs one behind a single symptom at the same strength.
Because the overlap is a ratio of sums, this stays bounded in [0, 1].
"""
from __future__ import annotations

import concurrent.futures
import logging
from typing import Any, Optional

from collectors.opentargets import OpenTargetsError, get_disease_top_genes

log = logging.getLogger(__name__)

DEFAULT_MAX_PHENOTYPES = 20
DEFAULT_GENES_PER_PHENOTYPE = 50
DEFAULT_MAX_WORKERS = 5


def build_symptom_gene_set(
    phenotypes: list[dict[str, Any]],
    genes_per_phenotype: int = DEFAULT_GENES_PER_PHENOTYPE,
    max_phenotypes: int = DEFAULT_MAX_PHENOTYPES,
    max_workers: int = DEFAULT_MAX_WORKERS,
) -> dict[str, Any]:
    """Expand symptoms into a pooled gene set.

    Args:
        phenotypes: entries from ``get_disease_phenotypes``.
        genes_per_phenotype: top associated targets to take per symptom.
        max_phenotypes: how many symptoms to expand (they arrive most relevant
            first, and each one costs an API call).

    Returns:
        ``{"genes", "per_phenotype", "expanded", "empty", "failed",
        "phenotype_count"}``. ``genes`` is the pooled set,
        ``[{"symbol", "score", "max_score", "phenotype_count", "phenotypes"}]``
        sorted by descending score — the shape
        :func:`nw_overlap.calc_network_overlap` consumes. ``per_phenotype``
        keeps each symptom's own gene list, so a gene can also be scored
        symptom by symptom.
    """
    # Phenotypes marked as absent in this disease describe what it does not
    # cause, so they must not seed genes.
    usable = [p for p in (phenotypes or []) if isinstance(p, dict) and not p.get("excluded")]
    usable = [p for p in usable if p.get("ontology_id")][: max(1, int(max_phenotypes))]

    if not usable:
        return {"genes": [], "per_phenotype": [], "expanded": [], "empty": [],
                "failed": [], "phenotype_count": 0}

    def fetch(phenotype: dict[str, Any]) -> tuple[dict[str, Any], Optional[list[dict]]]:
        try:
            genes = get_disease_top_genes(
                phenotype["ontology_id"], top_n=genes_per_phenotype
            )
            return phenotype, genes
        except OpenTargetsError as exc:
            # A phenotype Open Targets does not index as a disease is expected,
            # not an error worth failing the analysis over.
            log.info("[symptoms] %s: %s", phenotype.get("hpo_id"), exc)
            return phenotype, None
        except Exception as exc:
            log.warning("[symptoms] %s: %s", phenotype.get("hpo_id"), exc)
            return phenotype, None

    pooled: dict[str, dict[str, Any]] = {}
    per_phenotype: list[dict[str, Any]] = []
    expanded: list[dict[str, Any]] = []
    empty: list[dict[str, Any]] = []
    failed: list[dict[str, Any]] = []

    workers = max(1, min(int(max_workers), len(usable)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        for phenotype, genes in executor.map(fetch, usable):
            label = {"hpo_id": phenotype.get("hpo_id", ""),
                     "name": phenotype.get("name", ""),
                     "frequency": phenotype.get("frequency", "")}
            if genes is None:
                failed.append(label)
                continue
            if not genes:
                empty.append(label)
                continue

            label = {**label, "gene_count": len(genes),
                     "resources": phenotype.get("resources") or [],
                     "ontology_id": phenotype.get("ontology_id", "")}
            # Keep the per-symptom list so each symptom can be scored on its own.
            per_phenotype.append({**label, "genes": [dict(g) for g in genes]})
            expanded.append(label)
            for gene in genes:
                symbol = str(gene.get("symbol") or "").strip()
                if not symbol:
                    continue
                try:
                    score = float(gene.get("score") or 0.0)
                except (TypeError, ValueError):
                    score = 0.0

                entry = pooled.setdefault(
                    symbol.upper(),
                    {"symbol": symbol, "score": 0.0, "max_score": 0.0,
                     "phenotype_count": 0, "phenotypes": []},
                )
                entry["score"] += score
                entry["max_score"] = max(entry["max_score"], score)
                entry["phenotype_count"] += 1
                if len(entry["phenotypes"]) < 5:
                    entry["phenotypes"].append(phenotype.get("name", ""))

    genes_out = sorted(pooled.values(), key=lambda g: g["score"], reverse=True)
    for gene in genes_out:
        gene["score"] = round(gene["score"], 4)
        gene["max_score"] = round(gene["max_score"], 4)

    # Preserve the input ordering (most relevant symptom first), which the
    # thread pool does not guarantee.
    order = {p["ontology_id"]: i for i, p in enumerate(usable)}
    per_phenotype.sort(key=lambda p: order.get(p.get("ontology_id", ""), 1_000_000))
    expanded.sort(key=lambda p: order.get(p.get("ontology_id", ""), 1_000_000))

    return {
        "genes": genes_out,
        "per_phenotype": per_phenotype,
        "expanded": expanded,
        "empty": empty,
        "failed": failed,
        "phenotype_count": len(usable),
    }
