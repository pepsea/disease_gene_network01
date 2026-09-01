"""Significance testing for the symptom matrix.

A gene that covers one symptom strongly and nothing else is a different claim
from a gene that covers many symptoms. The matrix mean blends the two; this
separates them and puts a number on the second.

For one symptom, the question is whether the gene's PPI neighbourhood overlaps
that symptom's gene set more than a neighbourhood of the same size drawn at
random would. That is a hypergeometric tail probability. Across symptoms the
per-symptom p-values are combined with Fisher's method, which is exactly the
behaviour wanted here: consistent moderate overlap across many symptoms beats
one strong hit and nothing else.

Genes are then corrected for multiple testing with Benjamini-Hochberg, since
every input gene is tested against the same symptom set.

Implemented in log space with the standard library only — no SciPy — so the
container stays small and the numbers stay exact enough for tail probabilities
far below float resolution.
"""
from __future__ import annotations

import math
from typing import Any, Iterable

# Protein-coding genes, the population the neighbourhood is drawn from.
DEFAULT_BACKGROUND = 20000

# p-values below this are floored so that Fisher's method never takes log(0).
_P_FLOOR = 1e-300


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k: int, background: int, set_size: int, draws: int) -> float:
    """``P(X >= k)`` for a hypergeometric draw.

    Args:
        k: observed overlap.
        background: total genes in the population.
        set_size: genes in the symptom's set.
        draws: genes in the gene's neighbourhood (itself included).
    """
    k, background, set_size, draws = int(k), int(background), int(set_size), int(draws)
    if k <= 0 or draws <= 0 or set_size <= 0:
        return 1.0
    if background <= 0 or set_size > background or draws > background:
        return 1.0
    # More overlap than is possible cannot happen.
    upper = min(draws, set_size)
    if k > upper:
        return 0.0

    log_denominator = _log_comb(background, draws)
    total = 0.0
    for i in range(k, upper + 1):
        log_term = (
            _log_comb(set_size, i)
            + _log_comb(background - set_size, draws - i)
            - log_denominator
        )
        if log_term > -745:  # below this exp() underflows to zero anyway
            total += math.exp(log_term)
    return min(1.0, max(0.0, total))


def fisher_combine(p_values: Iterable[float]) -> float:
    """Combine independent p-values with Fisher's method.

    The combined statistic is chi-square distributed with ``2m`` degrees of
    freedom. Because the degrees of freedom are always even, the survival
    function has a closed form and needs no special functions.
    """
    values = [min(1.0, max(_P_FLOOR, float(p))) for p in p_values]
    if not values:
        return 1.0

    statistic = -2.0 * sum(math.log(p) for p in values)
    m = len(values)

    # P(chi2_{2m} >= x) = exp(-x/2) * sum_{i=0}^{m-1} (x/2)^i / i!
    half = statistic / 2.0
    if half <= 0:
        return 1.0

    log_total_terms = []
    for i in range(m):
        log_total_terms.append(i * math.log(half) - math.lgamma(i + 1) if half > 0 else 0.0)
    max_log = max(log_total_terms)
    summed = sum(math.exp(t - max_log) for t in log_total_terms)
    log_p = -half + max_log + math.log(summed)
    return min(1.0, max(0.0, math.exp(log_p)))


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    """Return BH-adjusted q-values, in the order the p-values were given."""
    n = len(p_values)
    if n == 0:
        return []

    order = sorted(range(n), key=lambda i: p_values[i])
    adjusted = [1.0] * n
    previous = 1.0
    # Walk from the largest p-value down, keeping the running minimum.
    for rank in range(n, 0, -1):
        index = order[rank - 1]
        value = min(previous, p_values[index] * n / rank)
        adjusted[index] = min(1.0, max(0.0, value))
        previous = adjusted[index]
    return adjusted


def score_symptom_breadth(
    cells: list[dict[str, Any]],
    ppi_partner_count: int,
    background: int = DEFAULT_BACKGROUND,
    target_counts_as_draw: bool = True,
) -> dict[str, Any]:
    """Test a gene's overlap across every symptom.

    Args:
        cells: one entry per symptom, carrying ``matched_count`` and
            ``gene_count`` (the symptom's set size).
        ppi_partner_count: the neighbourhood size the overlap was drawn from.

    Returns per-symptom p-values plus the Fisher-combined p-value, the number of
    symptoms actually covered, and that as a share of the symptoms tested.
    """
    draws = int(ppi_partner_count) + (1 if target_counts_as_draw else 0)

    p_values: list[float] = []
    covered = 0
    for cell in cells or []:
        matched = int(cell.get("matched_count") or 0)
        set_size = int(cell.get("gene_count") or 0)
        p = hypergeom_sf(matched, background, set_size, draws)
        cell["p_value"] = p
        p_values.append(p)
        if matched > 0:
            covered += 1

    total = len(p_values)
    return {
        "symptom_p_value": fisher_combine(p_values) if p_values else 1.0,
        "symptom_breadth_count": covered,
        "symptom_breadth_percent": round(covered / total * 100, 1) if total else 0.0,
        "symptom_tested_count": total,
    }


def apply_fdr(results: list[dict[str, Any]], key: str = "symptom_p_value",
              q_key: str = "symptom_q_value") -> None:
    """Add BH q-values across the tested genes, in place."""
    testable = [r for r in results if isinstance(r.get(key), (int, float))]
    if not testable:
        return
    q_values = benjamini_hochberg([float(r[key]) for r in testable])
    for result, q in zip(testable, q_values):
        result[q_key] = q
