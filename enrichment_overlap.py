"""Overlap between a gene's pathway enrichment and the disease's.

The network overlap in :mod:`nw_overlap` asks "how many of the disease's genes
does this gene reach?". This asks the same question one level up: "how much of
the disease's enriched biology does this gene's neighbourhood reproduce?"

Two genes can sit in the same pathways without interacting directly, so this
catches functional proximity that the gene-level overlap misses.

Also carried through is the original project's ``pathway_fit`` signal: whether
the target gene is itself annotated to the disease's enriched pathways.
"""
from __future__ import annotations

import math
from typing import Any

# Interpretation bands for the weighted pathway coverage.
STRONG_THRESHOLD = 0.3
MODERATE_THRESHOLD = 0.1

# A p-value floor, so that a single astronomically significant term cannot
# dominate the weighting.
_MIN_P = 1e-300


def interpret_score(score: float) -> str:
    """Return ``"strong"`` / ``"moderate"`` / ``"weak"`` for a coverage score."""
    if score >= STRONG_THRESHOLD:
        return "strong"
    if score >= MODERATE_THRESHOLD:
        return "moderate"
    return "weak"


def _weight(pathway: dict[str, Any]) -> float:
    """Weight a disease pathway by its significance: -log10(p)."""
    try:
        p = float(pathway.get("p_value", 1.0))
    except (TypeError, ValueError):
        return 0.0
    if p <= 0:
        p = _MIN_P
    if p >= 1:
        return 0.0
    return -math.log10(p)


def _term_id(pathway: dict[str, Any]) -> str:
    return str(pathway.get("term_id") or "").strip().upper()


def calc_enrichment_overlap(
    gene: str,
    gene_pathways: list[dict[str, Any]],
    disease_pathways: list[dict[str, Any]],
) -> dict[str, Any]:
    """Score how much of the disease's enriched biology a gene reproduces.

    Args:
        gene: the target gene symbol.
        gene_pathways: enrichment of the gene plus its PPI partners.
        disease_pathways: enrichment of the disease's top genes.

    Returns:
        ``pathway_weighted_percent`` — share of the disease's total pathway
        significance (-log10 p) that the gene covers, counting both the
        network's enrichment and the target's own pathway membership;
        ``pathway_overlap_percent`` — the same unweighted, by term count;
        ``pathway_matched_count`` — pathways covered in total;
        ``pathway_overlap_count`` — of those, the ones the network reached;
        ``target_fit_percent`` — the original project's ``pathway_fit``: the
        share of disease pathways the target gene is itself annotated to;
        plus the overlapping terms, each tagged with how it was reached.
    """
    target = (gene or "").strip().upper()

    disease_pathways = [
        p for p in (disease_pathways or []) if isinstance(p, dict) and _term_id(p)
    ]
    gene_pathways = [
        p for p in (gene_pathways or []) if isinstance(p, dict) and _term_id(p)
    ]

    disease_count = len(disease_pathways)
    gene_term_ids = {_term_id(p) for p in gene_pathways}

    # Disease pathways reached through the network's enrichment.
    network_hits = {_term_id(p) for p in disease_pathways if _term_id(p) in gene_term_ids}

    # Disease pathways the target gene itself is annotated to (the original
    # project's pathway_fit). The target counts on its own, exactly as its OT
    # score does in the gene-level overlap: a single gene rarely produces a
    # significant enrichment term, so relying on the enrichment alone would
    # silently drop the target's own membership.
    target_in = [
        p
        for p in disease_pathways
        if target in {str(g).strip().upper() for g in (p.get("genes") or [])}
    ]
    target_hits = {_term_id(p) for p in target_in}
    target_fit = round(len(target_in) / disease_count, 3) if disease_count else 0.0

    # Covered = reached through the network OR through the target itself.
    covered = network_hits | target_hits
    overlap = [p for p in disease_pathways if _term_id(p) in covered]

    total_weight = sum(_weight(p) for p in disease_pathways)
    overlap_weight = sum(_weight(p) for p in overlap)
    weighted = round(overlap_weight / total_weight, 3) if total_weight > 0 else 0.0
    simple = round(len(overlap) / disease_count, 3) if disease_count else 0.0

    return {
        "pathway_weighted_score": weighted,
        "pathway_weighted_percent": round(weighted * 100, 1),
        "pathway_simple_ratio": simple,
        "pathway_overlap_percent": round(simple * 100, 1),
        # Total disease pathways covered, target included.
        "pathway_matched_count": len(overlap),
        # Of those, the ones the network's enrichment reached.
        "pathway_overlap_count": len(network_hits),
        "disease_pathway_count": disease_count,
        "gene_pathway_count": len(gene_pathways),
        "target_fit_score": target_fit,
        "target_fit_percent": round(target_fit * 100, 1),
        "target_in_pathway_count": len(target_in),
        "target_in_pathways": [
            {
                "term_id": p.get("term_id", ""),
                "name": p.get("name", ""),
                "source": p.get("source", ""),
                "p_value": p.get("p_value", 1.0),
            }
            for p in sorted(target_in, key=lambda p: p.get("p_value", 1.0))
        ],
        "overlapping_pathways": [
            {
                "term_id": p.get("term_id", ""),
                "name": p.get("name", ""),
                "source": p.get("source", ""),
                "p_value": p.get("p_value", 1.0),
                "term_size": p.get("term_size", 0),
                # How this pathway was reached: the network's enrichment, the
                # target's own membership, or both.
                "via": (
                    "both"
                    if _term_id(p) in network_hits and _term_id(p) in target_hits
                    else "network"
                    if _term_id(p) in network_hits
                    else "target"
                ),
            }
            for p in sorted(overlap, key=lambda p: p.get("p_value", 1.0))
        ],
        "pathway_interpretation": interpret_score(weighted),
    }
