import math

import pytest

from symptom_stats import (
    apply_fdr,
    benjamini_hochberg,
    fisher_combine,
    hypergeom_sf,
    score_symptom_breadth,
)


# --- hypergeometric --------------------------------------------------------

def test_no_overlap_is_not_significant():
    assert hypergeom_sf(0, 20000, 50, 30) == 1.0


def test_more_overlap_is_more_significant():
    p1 = hypergeom_sf(1, 20000, 50, 30)
    p3 = hypergeom_sf(3, 20000, 50, 30)
    p5 = hypergeom_sf(5, 20000, 50, 30)
    assert 1.0 > p1 > p3 > p5 > 0.0


def test_a_bigger_neighbourhood_makes_the_same_overlap_less_surprising():
    """Finding 3 hits among 100 partners is weaker than among 10."""
    assert hypergeom_sf(3, 20000, 50, 10) < hypergeom_sf(3, 20000, 50, 100)


def test_a_bigger_symptom_set_makes_the_same_overlap_less_surprising():
    assert hypergeom_sf(3, 20000, 20, 30) < hypergeom_sf(3, 20000, 500, 30)


def test_overlap_beyond_what_is_possible_is_zero():
    assert hypergeom_sf(40, 20000, 50, 30) == 0.0   # cannot exceed the draws


def test_complete_overlap_of_a_small_set():
    # Every one of 3 partners is in a 3-gene symptom set.
    p = hypergeom_sf(3, 20000, 3, 3)
    assert 0 < p < 1e-9


@pytest.mark.parametrize("args", [
    (5, 0, 50, 30), (5, 20000, 0, 30), (5, 20000, 50, 0),
    (5, 20000, 50, 30000), (5, 20000, 30000, 30),
])
def test_degenerate_inputs_are_not_significant(args):
    assert hypergeom_sf(*args) == 1.0


def test_the_tail_is_a_probability():
    for k in range(0, 10):
        assert 0.0 <= hypergeom_sf(k, 20000, 50, 30) <= 1.0


# --- Fisher's method -------------------------------------------------------

def test_a_single_p_value_passes_through():
    assert fisher_combine([0.001]) == pytest.approx(0.001, rel=1e-6)


def test_breadth_beats_a_single_strong_hit():
    """The point of the test: many moderate hits outrank one strong one."""
    broad = fisher_combine([0.05] * 5)
    narrow = fisher_combine([0.001, 1.0, 1.0, 1.0, 1.0])
    assert broad < narrow


def test_consistent_coverage_compounds():
    assert fisher_combine([0.05] * 5) < fisher_combine([0.05] * 2) < fisher_combine([0.05])


def test_uninformative_p_values_combine_to_one():
    assert fisher_combine([1.0, 1.0, 1.0]) == pytest.approx(1.0, abs=1e-9)


def test_no_p_values_is_not_significant():
    assert fisher_combine([]) == 1.0


def test_zero_p_values_do_not_produce_a_nan():
    result = fisher_combine([0.0, 0.5])
    assert not math.isnan(result)
    assert 0.0 <= result <= 1.0


def test_the_result_is_a_probability():
    assert 0.0 <= fisher_combine([0.2, 0.3, 0.9]) <= 1.0


# --- Benjamini-Hochberg ----------------------------------------------------

def test_bh_adjusts_upward():
    q = benjamini_hochberg([0.001, 0.01, 0.05, 0.5])
    assert all(qi >= pi for qi, pi in zip(q, [0.001, 0.01, 0.05, 0.5]))


def test_bh_preserves_input_order():
    q = benjamini_hochberg([0.5, 0.001])
    assert q[0] > q[1]


def test_bh_is_monotone_in_rank():
    p = [0.001, 0.02, 0.03, 0.9]
    q = benjamini_hochberg(p)
    ranked = [q[i] for i in sorted(range(len(p)), key=lambda i: p[i])]
    assert ranked == sorted(ranked)


def test_bh_never_exceeds_one():
    assert all(v <= 1.0 for v in benjamini_hochberg([0.9, 0.95, 0.99]))


def test_bh_of_nothing():
    assert benjamini_hochberg([]) == []


# --- the per-gene wrapper --------------------------------------------------

def cells(*pairs):
    return [{"matched_count": m, "gene_count": n} for m, n in pairs]


def test_breadth_counts_the_symptoms_actually_covered():
    got = score_symptom_breadth(cells((2, 50), (0, 50), (1, 50), (0, 50)), 30)
    assert got["symptom_breadth_count"] == 2
    assert got["symptom_tested_count"] == 4
    assert got["symptom_breadth_percent"] == 50.0


def test_a_broad_gene_is_more_significant_than_a_narrow_one():
    broad = score_symptom_breadth(cells((2, 50), (2, 50), (2, 50), (2, 50)), 30)
    narrow = score_symptom_breadth(cells((3, 50), (0, 50), (0, 50), (0, 50)), 30)
    assert broad["symptom_breadth_count"] > narrow["symptom_breadth_count"]
    assert broad["symptom_p_value"] < narrow["symptom_p_value"]


def test_per_symptom_p_values_are_recorded_on_the_cells():
    rows = cells((2, 50), (0, 50))
    score_symptom_breadth(rows, 30)
    assert rows[0]["p_value"] < 1.0
    assert rows[1]["p_value"] == 1.0


def test_the_target_itself_counts_as_a_draw():
    with_self = score_symptom_breadth(cells((1, 50)), 0)   # no partners, self only
    assert with_self["symptom_p_value"] < 1.0


def test_a_gene_covering_nothing_is_not_significant():
    got = score_symptom_breadth(cells((0, 50), (0, 50)), 30)
    assert got["symptom_p_value"] == pytest.approx(1.0, abs=1e-9)
    assert got["symptom_breadth_count"] == 0


def test_no_symptoms_is_not_significant():
    got = score_symptom_breadth([], 30)
    assert got["symptom_p_value"] == 1.0
    assert got["symptom_breadth_percent"] == 0.0


def test_fdr_is_applied_across_genes():
    results = [
        {"gene": "A", "symptom_p_value": 0.001},
        {"gene": "B", "symptom_p_value": 0.5},
        {"gene": "C"},                              # untested, e.g. errored
    ]
    apply_fdr(results)
    assert results[0]["symptom_q_value"] > 0.001
    assert results[0]["symptom_q_value"] < results[1]["symptom_q_value"]
    assert "symptom_q_value" not in results[2]


def test_fdr_on_an_empty_list_is_a_no_op():
    results = []
    apply_fdr(results)
    assert results == []
