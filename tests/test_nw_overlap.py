import pytest

from nw_overlap import calc_network_overlap, interpret_score

DISEASE = [
    {"symbol": "APP", "score": 0.90},
    {"symbol": "PSEN1", "score": 0.60},
    {"symbol": "MAPT", "score": 0.50},
    {"symbol": "APOE", "score": 0.40},
    {"symbol": "BACE1", "score": 0.10},
]
TOTAL = 2.50


def test_partner_overlap_is_score_weighted():
    r = calc_network_overlap("XYZ", ["PSEN1", "MAPT"], DISEASE)
    assert r["weighted_score"] == round(1.10 / TOTAL, 3)
    assert r["overlap_count"] == 2
    assert r["simple_ratio"] == round(2 / 5, 3)
    assert r["target_self"] is None
    assert r["target_self_score"] is None


def test_target_own_score_is_added_when_it_is_a_disease_gene():
    r = calc_network_overlap("APP", ["PSEN1"], DISEASE)
    assert r["weighted_score"] == round((0.60 + 0.90) / TOTAL, 3)
    assert r["target_self"] == "APP"
    assert r["target_self_score"] == 0.90
    # The target itself is not counted among the overlapping partners.
    assert r["overlap_count"] == 1
    assert r["simple_ratio"] == round(2 / 5, 3)


def test_target_is_never_double_counted_when_listed_as_its_own_partner():
    r = calc_network_overlap("APP", ["APP", "PSEN1"], DISEASE)
    assert r["weighted_score"] == round((0.60 + 0.90) / TOTAL, 3)
    assert [g["symbol"] for g in r["overlapping_genes"]] == ["PSEN1"]


def test_matching_is_case_insensitive():
    upper = calc_network_overlap("app", ["psen1", "mapt"], DISEASE)
    lower = calc_network_overlap("APP", ["PSEN1", "MAPT"], DISEASE)
    assert upper["weighted_score"] == lower["weighted_score"]
    assert upper["target_self"] == "app"


def test_overlapping_genes_are_sorted_by_descending_score():
    r = calc_network_overlap("XYZ", ["BACE1", "APP", "MAPT"], DISEASE)
    assert [g["symbol"] for g in r["overlapping_genes"]] == ["APP", "MAPT", "BACE1"]


def test_no_partners_scores_zero():
    r = calc_network_overlap("XYZ", [], DISEASE)
    assert r["weighted_score"] == 0.0
    assert r["overlap_count"] == 0
    assert r["ppi_partner_count"] == 0
    assert r["interpretation"] == "weak"


def test_empty_disease_gene_list_does_not_divide_by_zero():
    r = calc_network_overlap("APP", ["PSEN1"], [])
    assert r["weighted_score"] == 0.0
    assert r["simple_ratio"] == 0.0
    assert r["disease_gene_count"] == 0


def test_all_zero_disease_scores_do_not_divide_by_zero():
    r = calc_network_overlap("A", ["B"], [{"symbol": "B", "score": 0.0}])
    assert r["weighted_score"] == 0.0
    assert r["overlap_count"] == 1


def test_duplicate_partners_are_counted_once():
    r = calc_network_overlap("XYZ", ["MAPT", "mapt", "MAPT"], DISEASE)
    assert r["ppi_partner_count"] == 1
    assert r["overlap_count"] == 1


def test_full_overlap_scores_one():
    r = calc_network_overlap("APP", ["PSEN1", "MAPT", "APOE", "BACE1"], DISEASE)
    assert r["weighted_score"] == 1.0
    assert r["simple_ratio"] == 1.0


def test_malformed_entries_are_tolerated():
    dirty = [
        {"symbol": "PSEN1", "score": "0.6"},   # numeric string
        {"symbol": "MAPT"},                     # missing score
        {"score": 0.3},                         # missing symbol
        {"symbol": "APOE", "score": None},      # null score
        "not-a-dict",                           # wrong type
    ]
    r = calc_network_overlap("XYZ", ["PSEN1", "MAPT", "APOE"], dirty)
    # The symbol-less row and the non-dict are dropped; three usable rows remain.
    assert r["disease_gene_count"] == 3
    assert r["weighted_score"] == 1.0  # 0.6 overlapping out of a 0.6 total
    assert r["overlap_count"] == 3


@pytest.mark.parametrize(
    "score,expected",
    [(0.9, "strong"), (0.3, "strong"), (0.299, "moderate"), (0.1, "moderate"),
     (0.099, "weak"), (0.0, "weak")],
)
def test_interpretation_bands(score, expected):
    assert interpret_score(score) == expected
