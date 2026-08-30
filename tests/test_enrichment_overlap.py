import pytest

from enrichment_overlap import calc_enrichment_overlap, interpret_score

DISEASE = [
    {"term_id": "R-HSA-1", "name": "Amyloid fiber formation", "source": "REAC",
     "p_value": 1e-10, "genes": ["APP", "PSEN1"], "term_size": 30},
    {"term_id": "GO:1", "name": "neuron death", "source": "GO:BP",
     "p_value": 1e-5, "genes": ["MAPT", "APP"], "term_size": 200},
    {"term_id": "WP:1", "name": "unrelated", "source": "WP",
     "p_value": 1e-1, "genes": ["ZZZ"], "term_size": 50},
]
# Weights are -log10(p): 10, 5, 1 -> total 16


def test_full_overlap_scores_one_hundred_percent():
    gene = [{"term_id": t["term_id"], "p_value": 0.01} for t in DISEASE]
    r = calc_enrichment_overlap("APP", gene, DISEASE)
    assert r["pathway_overlap_percent"] == 100.0
    assert r["pathway_weighted_percent"] == 100.0
    assert r["pathway_overlap_count"] == 3


def test_weighting_favours_the_most_significant_terms():
    r = calc_enrichment_overlap("XYZ", [{"term_id": "R-HSA-1", "p_value": 0.01}], DISEASE)
    assert r["pathway_overlap_count"] == 1
    assert r["pathway_overlap_percent"] == 33.3        # 1 of 3 terms
    assert r["pathway_weighted_percent"] == 62.5       # but 10 of 16 weight


def test_covering_only_a_weak_term_scores_low():
    r = calc_enrichment_overlap("XYZ", [{"term_id": "WP:1", "p_value": 0.02}], DISEASE)
    assert r["pathway_overlap_percent"] == 33.3
    assert r["pathway_weighted_percent"] == 6.2        # 1 of 16 weight
    assert r["pathway_interpretation"] == "weak"


def test_no_overlap_scores_zero():
    r = calc_enrichment_overlap("XYZ", [{"term_id": "GO:999", "p_value": 0.01}], DISEASE)
    assert r["pathway_overlap_percent"] == 0.0
    assert r["pathway_weighted_percent"] == 0.0
    assert r["overlapping_pathways"] == []


def test_target_fit_counts_pathways_containing_the_gene():
    r = calc_enrichment_overlap("APP", [], DISEASE)
    assert r["target_in_pathway_count"] == 2           # both list APP
    assert r["target_fit_percent"] == 66.7
    assert [p["name"] for p in r["target_in_pathways"]] == [
        "Amyloid fiber formation", "neuron death"]


def test_target_fit_is_case_insensitive():
    assert calc_enrichment_overlap("app", [], DISEASE)["target_in_pathway_count"] == 2


def test_target_fit_is_independent_of_network_overlap():
    """A gene can sit in the disease pathways while its network covers none."""
    r = calc_enrichment_overlap("APP", [], DISEASE)
    assert r["pathway_overlap_percent"] == 0.0
    assert r["target_fit_percent"] == 66.7


def test_term_ids_match_case_insensitively():
    r = calc_enrichment_overlap("X", [{"term_id": "r-hsa-1", "p_value": 0.01}], DISEASE)
    assert r["pathway_overlap_count"] == 1


def test_overlapping_pathways_are_sorted_by_significance():
    gene = [{"term_id": "GO:1", "p_value": 0.01}, {"term_id": "R-HSA-1", "p_value": 0.01}]
    r = calc_enrichment_overlap("X", gene, DISEASE)
    assert [p["term_id"] for p in r["overlapping_pathways"]] == ["R-HSA-1", "GO:1"]


def test_empty_disease_enrichment_does_not_divide_by_zero():
    r = calc_enrichment_overlap("APP", [{"term_id": "R-HSA-1", "p_value": 0.01}], [])
    assert r["pathway_weighted_percent"] == 0.0
    assert r["pathway_overlap_percent"] == 0.0
    assert r["disease_pathway_count"] == 0


def test_empty_gene_enrichment_scores_zero():
    r = calc_enrichment_overlap("APP", [], DISEASE)
    assert r["pathway_overlap_percent"] == 0.0
    assert r["gene_pathway_count"] == 0


def test_all_insignificant_disease_terms_do_not_divide_by_zero():
    flat = [{"term_id": "A", "p_value": 1.0}, {"term_id": "B", "p_value": 1.0}]
    r = calc_enrichment_overlap("X", [{"term_id": "A", "p_value": 1.0}], flat)
    assert r["pathway_weighted_percent"] == 0.0   # zero total weight
    assert r["pathway_overlap_percent"] == 50.0   # counts still work


def test_zero_p_value_does_not_break_the_weighting():
    disease = [{"term_id": "A", "p_value": 0.0}, {"term_id": "B", "p_value": 0.5}]
    r = calc_enrichment_overlap("X", [{"term_id": "A", "p_value": 0.0}], disease)
    assert 0 < r["pathway_weighted_percent"] <= 100.0


def test_malformed_entries_are_tolerated():
    disease = [
        {"term_id": "A", "p_value": "1e-5"},   # numeric string
        {"term_id": "", "p_value": 1e-5},      # no id -> dropped
        {"p_value": 1e-5},                     # no id -> dropped
        "junk",                                # wrong type -> dropped
        {"term_id": "B", "p_value": None},     # unusable p -> zero weight
    ]
    r = calc_enrichment_overlap("X", [{"term_id": "A"}, "junk"], disease)
    assert r["disease_pathway_count"] == 2
    assert r["gene_pathway_count"] == 1
    assert r["pathway_overlap_count"] == 1


@pytest.mark.parametrize("score,expected", [
    (0.9, "strong"), (0.3, "strong"), (0.299, "moderate"),
    (0.1, "moderate"), (0.099, "weak"), (0.0, "weak"),
])
def test_interpretation_bands(score, expected):
    assert interpret_score(score) == expected
