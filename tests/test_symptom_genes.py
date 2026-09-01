import pytest

import symptom_genes as sg
from collectors.opentargets import OpenTargetsError


def pheno(hpo, name, excluded=False, freq=""):
    return {"hpo_id": hpo, "ontology_id": hpo.replace(":", "_"), "name": name,
            "frequency": freq, "aspect": "P", "resource": "HPO", "excluded": excluded}


PHENOTYPES = [
    pheno("HP:0002354", "Memory impairment", freq="HP:0040281"),
    pheno("HP:0000726", "Dementia", freq="HP:0040282"),
]
GENES = {
    "HP_0002354": [{"symbol": "APP", "score": 0.8}, {"symbol": "MAPT", "score": 0.5}],
    "HP_0000726": [{"symbol": "APP", "score": 0.6}, {"symbol": "PSEN1", "score": 0.7}],
}


def patch_genes(monkeypatch, table=None, fail=(), empty=()):
    table = GENES if table is None else table

    def fake(ontology_id, top_n=50):
        if ontology_id in fail:
            raise OpenTargetsError(f"Unknown disease id: {ontology_id}")
        if ontology_id in empty:
            return []
        return [dict(g) for g in table.get(ontology_id, [])]

    monkeypatch.setattr(sg, "get_disease_top_genes", fake)


def by_symbol(result):
    return {g["symbol"]: g for g in result["genes"]}


# --- pooling and scoring ---------------------------------------------------

def test_genes_are_pooled_across_symptoms(monkeypatch):
    patch_genes(monkeypatch)
    genes = by_symbol(sg.build_symptom_gene_set(PHENOTYPES))
    assert set(genes) == {"APP", "PSEN1", "MAPT"}


def test_score_sums_across_the_symptoms_a_gene_drives(monkeypatch):
    patch_genes(monkeypatch)
    genes = by_symbol(sg.build_symptom_gene_set(PHENOTYPES))
    assert genes["APP"]["score"] == 1.4        # 0.8 + 0.6
    assert genes["APP"]["max_score"] == 0.8
    assert genes["APP"]["phenotype_count"] == 2
    assert genes["PSEN1"]["score"] == 0.7
    assert genes["PSEN1"]["phenotype_count"] == 1


def test_breadth_outranks_a_single_strong_link(monkeypatch):
    table = {
        "HP_0002354": [{"symbol": "BROAD", "score": 0.5}, {"symbol": "NARROW", "score": 0.9}],
        "HP_0000726": [{"symbol": "BROAD", "score": 0.5}],
    }
    patch_genes(monkeypatch, table)
    result = sg.build_symptom_gene_set(PHENOTYPES)
    assert [g["symbol"] for g in result["genes"]] == ["BROAD", "NARROW"]  # 1.0 vs 0.9


def test_genes_are_sorted_by_descending_score(monkeypatch):
    patch_genes(monkeypatch)
    scores = [g["score"] for g in sg.build_symptom_gene_set(PHENOTYPES)["genes"]]
    assert scores == sorted(scores, reverse=True)


def test_contributing_symptoms_are_recorded(monkeypatch):
    patch_genes(monkeypatch)
    genes = by_symbol(sg.build_symptom_gene_set(PHENOTYPES))
    assert set(genes["APP"]["phenotypes"]) == {"Memory impairment", "Dementia"}
    assert genes["MAPT"]["phenotypes"] == ["Memory impairment"]


def test_contributing_symptom_list_is_capped(monkeypatch):
    many = [pheno(f"HP:{i:07d}", f"symptom {i}") for i in range(10)]
    patch_genes(monkeypatch, {p["ontology_id"]: [{"symbol": "HUB", "score": 0.1}] for p in many})
    genes = by_symbol(sg.build_symptom_gene_set(many))
    assert genes["HUB"]["phenotype_count"] == 10
    assert len(genes["HUB"]["phenotypes"]) == 5   # only the first few are listed


# --- excluded phenotypes ---------------------------------------------------

def test_excluded_phenotypes_do_not_seed_genes(monkeypatch):
    """A phenotype marked absent describes what the disease does NOT cause."""
    absent = pheno("HP:0009999", "Absent finding", excluded=True)
    patch_genes(monkeypatch, {**GENES, "HP_0009999": [{"symbol": "WRONG", "score": 0.9}]})
    result = sg.build_symptom_gene_set(PHENOTYPES + [absent])
    assert "WRONG" not in by_symbol(result)
    assert result["phenotype_count"] == 2
    assert "Absent finding" not in [e["name"] for e in result["expanded"]]


# --- limits ----------------------------------------------------------------

def test_max_phenotypes_caps_the_expansion(monkeypatch):
    many = [pheno(f"HP:{i:07d}", f"symptom {i}") for i in range(10)]
    patch_genes(monkeypatch, {p["ontology_id"]: [{"symbol": f"G{i}", "score": 0.5}]
                              for i, p in enumerate(many)})
    result = sg.build_symptom_gene_set(many, max_phenotypes=3)
    assert len(result["expanded"]) == 3
    assert result["phenotype_count"] == 3


def test_genes_per_phenotype_is_passed_through(monkeypatch):
    seen = {}

    def fake(ontology_id, top_n=50):
        seen["top_n"] = top_n
        return [{"symbol": "APP", "score": 0.5}]

    monkeypatch.setattr(sg, "get_disease_top_genes", fake)
    sg.build_symptom_gene_set(PHENOTYPES, genes_per_phenotype=17)
    assert seen["top_n"] == 17


# --- degraded inputs -------------------------------------------------------

def test_a_phenotype_open_targets_does_not_index_is_reported_not_fatal(monkeypatch):
    patch_genes(monkeypatch, fail={"HP_0000726"})
    result = sg.build_symptom_gene_set(PHENOTYPES)
    assert [f["name"] for f in result["failed"]] == ["Dementia"]
    assert set(by_symbol(result)) == {"APP", "MAPT"}   # the other symptom still counts


def test_a_phenotype_with_no_genes_is_reported_separately(monkeypatch):
    patch_genes(monkeypatch, empty={"HP_0000726"})
    result = sg.build_symptom_gene_set(PHENOTYPES)
    assert [e["name"] for e in result["empty"]] == ["Dementia"]
    assert result["failed"] == []


def test_no_phenotypes_returns_an_empty_set(monkeypatch):
    patch_genes(monkeypatch)
    result = sg.build_symptom_gene_set([])
    assert result == {"genes": [], "per_phenotype": [], "expanded": [],
                      "empty": [], "failed": [], "phenotype_count": 0}


def test_only_excluded_phenotypes_returns_an_empty_set(monkeypatch):
    patch_genes(monkeypatch)
    result = sg.build_symptom_gene_set([pheno("HP:1", "x", excluded=True)])
    assert result["genes"] == []
    assert result["phenotype_count"] == 0


def test_malformed_phenotypes_are_skipped(monkeypatch):
    patch_genes(monkeypatch)
    messy = [PHENOTYPES[0], {"hpo_id": "HP:0", "name": "no ontology id"}, "junk", None]
    result = sg.build_symptom_gene_set(messy)
    assert [e["name"] for e in result["expanded"]] == ["Memory impairment"]


def test_malformed_gene_rows_are_skipped(monkeypatch):
    patch_genes(monkeypatch, {"HP_0002354": [
        {"symbol": "APP", "score": 0.8}, {"symbol": "", "score": 0.5},
        {"score": 0.4}, {"symbol": "BAD", "score": "nan-ish"},
    ]})
    genes = by_symbol(sg.build_symptom_gene_set([PHENOTYPES[0]]))
    assert set(genes) == {"APP", "BAD"}
    assert genes["BAD"]["score"] == 0.0


def test_symbol_casing_is_pooled_but_display_casing_is_kept(monkeypatch):
    patch_genes(monkeypatch, {
        "HP_0002354": [{"symbol": "App", "score": 0.4}],
        "HP_0000726": [{"symbol": "APP", "score": 0.6}],
    })
    result = sg.build_symptom_gene_set(PHENOTYPES)
    assert len(result["genes"]) == 1
    assert result["genes"][0]["score"] == 1.0
    assert result["genes"][0]["symbol"] in {"App", "APP"}


# --- per-symptom lists (the matrix columns) --------------------------------

def test_each_symptom_keeps_its_own_gene_list(monkeypatch):
    patch_genes(monkeypatch)
    per = sg.build_symptom_gene_set(PHENOTYPES)["per_phenotype"]
    assert [p["name"] for p in per] == ["Memory impairment", "Dementia"]
    assert [g["symbol"] for g in per[0]["genes"]] == ["APP", "MAPT"]
    assert [g["symbol"] for g in per[1]["genes"]] == ["APP", "PSEN1"]


def test_per_symptom_order_follows_the_input_not_the_thread_pool(monkeypatch):
    many = [pheno(f"HP:{i:07d}", f"symptom {i}") for i in range(8)]
    patch_genes(monkeypatch, {p["ontology_id"]: [{"symbol": f"G{i}", "score": 0.5}]
                              for i, p in enumerate(many)})
    per = sg.build_symptom_gene_set(many, max_workers=8)["per_phenotype"]
    assert [p["name"] for p in per] == [f"symptom {i}" for i in range(8)]


def test_per_symptom_entries_carry_the_annotation_sources(monkeypatch):
    p1 = {**pheno("HP:0002354", "Memory impairment"), "resources": ["ORPHANET", "OMIM"]}
    patch_genes(monkeypatch, {"HP_0002354": [{"symbol": "APP", "score": 0.8}]})
    per = sg.build_symptom_gene_set([p1])["per_phenotype"]
    assert per[0]["resources"] == ["ORPHANET", "OMIM"]


def test_symptoms_without_genes_are_not_matrix_columns(monkeypatch):
    patch_genes(monkeypatch, empty={"HP_0000726"})
    per = sg.build_symptom_gene_set(PHENOTYPES)["per_phenotype"]
    assert [p["name"] for p in per] == ["Memory impairment"]


# --- integration with the shared scorer ------------------------------------

def test_the_gene_set_feeds_the_shared_overlap_scorer(monkeypatch):
    from nw_overlap import calc_network_overlap

    patch_genes(monkeypatch)
    genes = sg.build_symptom_gene_set(PHENOTYPES)["genes"]
    r = calc_network_overlap("APP", ["PSEN1"], genes)
    # Total 2.6; APP itself 1.4 plus PSEN1 0.7 = 2.1
    assert r["weighted_percent"] == pytest.approx(80.8, abs=0.1)
    assert r["matched_count"] == 2
    assert r["target_self"] == "APP"
