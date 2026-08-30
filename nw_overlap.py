"""Network-disease gene overlap scoring.

For one target gene, this measures how much of a disease's Open Targets gene
signal is reachable through the gene's PPI neighbourhood::

    weighted_score = ( Σ OT score of overlapping partners + OT score of the
                       target itself, if it is a disease gene )
                     / Σ OT score of all disease genes

A gene whose interaction partners are themselves strongly disease-associated
scores high; an isolated gene, or one whose partners are unrelated to the
disease, scores near zero.
"""
from __future__ import annotations

from typing import Any, Optional

# Interpretation bands for the weighted score, used by the UI legend.
STRONG_THRESHOLD = 0.3
MODERATE_THRESHOLD = 0.1


def interpret_score(weighted_score: float) -> str:
    """Return ``"strong"`` / ``"moderate"`` / ``"weak"`` for a weighted score."""
    if weighted_score >= STRONG_THRESHOLD:
        return "strong"
    if weighted_score >= MODERATE_THRESHOLD:
        return "moderate"
    return "weak"


def calc_network_overlap(
    gene: str,
    ppi_partners: list[str],
    disease_genes: list[dict[str, Any]],
) -> dict[str, Any]:
    """Score the overlap between a gene's PPI neighbourhood and a disease.

    Args:
        gene: The target gene symbol.
        ppi_partners: PPI partner symbols of ``gene``.
        disease_genes: ``[{"symbol": str, "score": float}, ...]`` from Open
            Targets, ordered by descending association score.

    Returns:
        A dict with ``weighted_score`` and ``weighted_percent`` (the share of
        the disease network's total OT score that the neighbourhood covers),
        ``simple_ratio`` and ``overlap_percent`` (the share of disease genes
        matched, unweighted), ``overlap_count``, ``matched_count`` (overlap plus
        the target itself), ``disease_gene_count``, ``ppi_partner_count``,
        ``target_self`` (the symbol when the target is itself a disease gene,
        else ``None``), ``target_self_score``, ``overlapping_genes``
        (score-sorted) and ``interpretation``.
    """
    target = (gene or "").strip().upper()

    def _score(entry: dict[str, Any]) -> float:
        try:
            return float(entry.get("score") or 0.0)
        except (TypeError, ValueError):
            return 0.0

    def _symbol(entry: dict[str, Any]) -> str:
        return str(entry.get("symbol") or "").strip().upper()

    # Entries without a symbol can never be matched, so they are dropped rather
    # than left to dilute the denominator.
    disease_genes = [
        g for g in (disease_genes or []) if isinstance(g, dict) and _symbol(g)
    ]
    total_ot_score = sum(_score(g) for g in disease_genes)
    ppi_set = {str(p).strip().upper() for p in (ppi_partners or []) if str(p).strip()}

    # Is the target itself among the disease's top genes?
    self_entry = next((g for g in disease_genes if _symbol(g) == target), None)
    self_score: Optional[float] = _score(self_entry) if self_entry else None

    # Disease genes reachable through the PPI neighbourhood (self excluded, so
    # that it is never counted twice).
    overlap = [g for g in disease_genes if _symbol(g) in ppi_set and _symbol(g) != target]

    weighted = sum(_score(g) for g in overlap)
    if self_score is not None:
        weighted += self_score

    simple_n = len(overlap) + (1 if self_entry else 0)
    simple_ratio = round(simple_n / max(1, len(disease_genes)), 3)
    weighted_score = round(weighted / total_ot_score, 3) if total_ot_score > 0 else 0.0

    return {
        "weighted_score": weighted_score,
        # Share of the disease network covered, as a percentage.
        "weighted_percent": round(weighted_score * 100, 1),
        "simple_ratio": simple_ratio,
        "overlap_percent": round(simple_ratio * 100, 1),
        "overlap_count": len(overlap),
        "matched_count": simple_n,
        "disease_gene_count": len(disease_genes),
        "ppi_partner_count": len(ppi_set),
        "target_self": gene.strip() if self_entry else None,
        "target_self_score": self_score,
        "overlapping_genes": sorted(overlap, key=_score, reverse=True),
        "interpretation": interpret_score(weighted_score),
    }
