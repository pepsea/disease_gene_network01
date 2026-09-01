"""Which source supplies the symptoms, and what happens when one fails."""
import pytest

import app as app_module
from app import collect_symptoms
from collectors.opentargets import OpenTargetsError

OT_PHENOTYPES = [{"hpo_id": "HP:0002354", "ontology_id": "HP_0002354",
                  "name": "Memory impairment", "frequency": "", "aspect": "P",
                  "resources": ["HPO"], "resource": "HPO", "excluded": False}]
HPO_PHENOTYPES = [{"hpo_id": "HP:0000726", "ontology_id": "HP_0000726",
                   "name": "Dementia", "frequency": "HP:0040281", "aspect": "P",
                   "resources": ["ORPHANET"], "resource": "ORPHANET", "excluded": False}]


@pytest.fixture
def patched(monkeypatch):
    monkeypatch.setattr(app_module, "get_disease_phenotypes",
                        lambda d, limit=50: list(OT_PHENOTYPES))
    monkeypatch.setattr(app_module, "get_disease_xrefs",
                        lambda d: ["OMIM:104300", "Orphanet:1020"])
    monkeypatch.setattr(app_module.hpo, "get_disease_phenotypes",
                        lambda ids, limit=50: list(HPO_PHENOTYPES))
    monkeypatch.setattr(app_module.hpo, "find_disease_ids_by_name", lambda n, limit=5: [])
    return monkeypatch


# --- source selection ------------------------------------------------------

def test_auto_prefers_open_targets(patched):
    phenotypes, meta = collect_symptoms("EFO_1", "Alzheimer disease", "auto", 50)
    assert meta["phenotype_source"] == "opentargets"
    assert [p["hpo_id"] for p in phenotypes] == ["HP:0002354"]


def test_auto_falls_back_to_hpo_when_open_targets_has_none(patched):
    patched.setattr(app_module, "get_disease_phenotypes", lambda d, limit=50: [])
    phenotypes, meta = collect_symptoms("EFO_1", "Alzheimer disease", "auto", 50)
    assert meta["phenotype_source"] == "hpo"
    assert [p["hpo_id"] for p in phenotypes] == ["HP:0000726"]
    assert meta["xrefs"] == ["OMIM:104300", "Orphanet:1020"]


def test_auto_falls_back_when_the_phenotypes_field_errors(patched):
    def boom(d, limit=50):
        raise OpenTargetsError("Cannot query field phenotypes")
    patched.setattr(app_module, "get_disease_phenotypes", boom)
    _, meta = collect_symptoms("EFO_1", "Alzheimer disease", "auto", 50)
    assert meta["phenotype_source"] == "hpo"


def test_opentargets_only_does_not_fall_back(patched):
    patched.setattr(app_module, "get_disease_phenotypes", lambda d, limit=50: [])
    def should_not_run(ids, limit=50):
        raise AssertionError("HPO must not be consulted for source=opentargets")
    patched.setattr(app_module.hpo, "get_disease_phenotypes", should_not_run)
    phenotypes, meta = collect_symptoms("EFO_1", "AD", "opentargets", 50)
    assert phenotypes == []
    assert meta["phenotype_source"] == ""


def test_hpo_only_skips_open_targets(patched):
    def should_not_run(d, limit=50):
        raise AssertionError("Open Targets must not be consulted for source=hpo")
    patched.setattr(app_module, "get_disease_phenotypes", should_not_run)
    phenotypes, meta = collect_symptoms("EFO_1", "AD", "hpo", 50)
    assert meta["phenotype_source"] == "hpo"
    assert [p["hpo_id"] for p in phenotypes] == ["HP:0000726"]


# --- cross references ------------------------------------------------------

def test_hpo_needs_a_cross_reference(patched):
    patched.setattr(app_module, "get_disease_xrefs", lambda d: [])
    phenotypes, meta = collect_symptoms("EFO_1", "AD", "hpo", 50)
    assert phenotypes == []
    assert meta["xrefs"] == []


def test_a_missing_xref_falls_back_to_an_exact_name_match(patched):
    patched.setattr(app_module, "get_disease_xrefs", lambda d: [])
    patched.setattr(app_module.hpo, "find_disease_ids_by_name",
                    lambda n, limit=5: ["ORPHA:1020"] if n == "Alzheimer disease" else [])
    phenotypes, meta = collect_symptoms("EFO_1", "Alzheimer disease", "hpo", 50)
    assert meta["xrefs"] == ["ORPHA:1020"]
    assert [p["hpo_id"] for p in phenotypes] == ["HP:0000726"]


def test_an_xref_lookup_failure_is_not_fatal(patched):
    def boom(d):
        raise OpenTargetsError("Cannot query field dbXRefs")
    patched.setattr(app_module, "get_disease_xrefs", boom)
    phenotypes, meta = collect_symptoms("EFO_1", "AD", "hpo", 50)
    assert phenotypes == []          # no ids, but no exception either


# --- a disease chosen directly in HPO --------------------------------------

def test_a_chosen_hpo_disease_skips_open_targets_entirely(patched):
    def should_not_run(*a, **k):
        raise AssertionError("neither phenotypes nor xrefs should be consulted")
    patched.setattr(app_module, "get_disease_phenotypes", should_not_run)
    patched.setattr(app_module, "get_disease_xrefs", should_not_run)
    phenotypes, meta = collect_symptoms("EFO_1", "AD", "auto", 50,
                                        hpo_disease_ids=["ORPHA:1020"])
    assert meta["phenotype_source"] == "hpo"
    assert meta["xrefs"] == ["ORPHA:1020"]
    assert meta["xref_origin"] == "selected"
    assert [p["hpo_id"] for p in phenotypes] == ["HP:0000726"]


def test_a_chosen_disease_uses_hpo_genes(patched):
    _, meta = collect_symptoms("EFO_1", "AD", "auto", 50, hpo_disease_ids=["OMIM:104300"])
    assert meta["gene_fetcher"] is app_module._hpo_gene_fetcher


def test_several_registrations_can_be_chosen_together(patched):
    captured = {}
    patched.setattr(app_module.hpo, "get_disease_phenotypes",
                    lambda ids, limit=50: captured.setdefault("ids", ids) and [] or list(HPO_PHENOTYPES))
    collect_symptoms("EFO_1", "AD", "auto", 50,
                     hpo_disease_ids=["OMIM:104300", "ORPHA:1020"])
    assert captured["ids"] == ["OMIM:104300", "ORPHA:1020"]


def test_blank_chosen_ids_are_ignored(patched):
    phenotypes, meta = collect_symptoms("EFO_1", "Alzheimer disease", "auto", 50,
                                        hpo_disease_ids=["", "   "])
    # Falls through to the normal Open Targets path.
    assert meta["phenotype_source"] == "opentargets"
    assert meta["xref_origin"] == ""


def test_the_xref_origin_is_recorded(patched):
    patched.setattr(app_module, "get_disease_phenotypes", lambda d, limit=50: [])
    _, meta = collect_symptoms("EFO_1", "AD", "auto", 50)
    assert meta["xref_origin"] == "dbxrefs"

    patched.setattr(app_module, "get_disease_xrefs", lambda d: [])
    patched.setattr(app_module.hpo, "find_disease_ids_by_name",
                    lambda n, limit=5: ["ORPHA:1020"])
    _, meta = collect_symptoms("EFO_1", "Alzheimer disease", "auto", 50)
    assert meta["xref_origin"] == "name"


# --- symptoms named directly by HP id --------------------------------------

CHOSEN_TERMS = [
    {"hpo_id": "HP:0002354", "ontology_id": "HP_0002354", "name": "Memory impairment",
     "frequency": "", "aspect": "P", "resources": [], "resource": "",
     "excluded": False, "gene_count": 3},
]


def test_chosen_hp_terms_outrank_everything(patched):
    def should_not_run(*a, **k):
        raise AssertionError("no disease lookup should happen for chosen HP terms")
    patched.setattr(app_module, "get_disease_phenotypes", should_not_run)
    patched.setattr(app_module, "get_disease_xrefs", should_not_run)
    patched.setattr(app_module.hpo, "get_disease_phenotypes", should_not_run)
    patched.setattr(app_module.hpo, "get_phenotypes_by_id",
                    lambda ids: [dict(t) for t in CHOSEN_TERMS])

    phenotypes, meta = collect_symptoms("EFO_1", "AD", "auto", 50,
                                        hpo_disease_ids=["ORPHA:1020"],
                                        hpo_phenotype_ids=["HP:0002354"])
    assert meta["phenotype_source"] == "hpo"
    assert meta["xref_origin"] == "phenotypes"
    assert meta["xrefs"] == ["HP:0002354"]
    assert [p["hpo_id"] for p in phenotypes] == ["HP:0002354"]


def test_chosen_terms_use_hpo_genes(patched):
    patched.setattr(app_module.hpo, "get_phenotypes_by_id",
                    lambda ids: [dict(t) for t in CHOSEN_TERMS])
    _, meta = collect_symptoms("EFO_1", "AD", "auto", 50, hpo_phenotype_ids=["HP:0002354"])
    assert meta["gene_fetcher"] is app_module._hpo_gene_fetcher


def test_blank_term_ids_fall_through_to_the_normal_path(patched):
    _, meta = collect_symptoms("EFO_1", "AD", "auto", 50, hpo_phenotype_ids=["", "  "])
    assert meta["phenotype_source"] == "opentargets"


def test_terms_that_resolve_to_nothing_do_not_claim_a_source(patched):
    patched.setattr(app_module.hpo, "get_phenotypes_by_id", lambda ids: [])
    phenotypes, meta = collect_symptoms("EFO_1", "AD", "auto", 50,
                                        hpo_phenotype_ids=["OMIM:104300"])
    assert phenotypes == []
    assert meta["phenotype_source"] == ""


# --- gene fetchers ---------------------------------------------------------

def test_hpo_sourced_symptoms_use_hpo_genes(patched):
    _, meta = collect_symptoms("EFO_1", "AD", "hpo", 50)
    assert meta["gene_fetcher"] is app_module._hpo_gene_fetcher


def test_open_targets_symptoms_use_the_hybrid_fetcher(patched):
    _, meta = collect_symptoms("EFO_1", "AD", "auto", 50)
    assert meta["gene_fetcher"] is app_module._hybrid_gene_fetcher


def test_the_hybrid_fetcher_prefers_open_targets(monkeypatch):
    monkeypatch.setattr(app_module, "get_disease_top_genes",
                        lambda oid, top_n=50: [{"symbol": "APP", "score": 0.8}])
    monkeypatch.setattr(app_module.hpo, "get_phenotype_genes",
                        lambda oid, top_n=50: [{"symbol": "WRONG", "score": 1.0}])
    got = app_module._hybrid_gene_fetcher("HP_0002354", 50)
    assert got["source"] == "opentargets"
    assert [g["symbol"] for g in got["genes"]] == ["APP"]


def test_the_hybrid_fetcher_falls_back_when_open_targets_has_no_genes(monkeypatch):
    """Open Targets only resolves an HP term it indexes as a disease."""
    monkeypatch.setattr(app_module, "get_disease_top_genes", lambda oid, top_n=50: [])
    monkeypatch.setattr(app_module.hpo, "get_phenotype_genes",
                        lambda oid, top_n=50: [{"symbol": "APP", "score": 1.0}])
    got = app_module._hybrid_gene_fetcher("HP_0002354", 50)
    assert got["source"] == "hpo"
    assert [g["symbol"] for g in got["genes"]] == ["APP"]


def test_the_hybrid_fetcher_falls_back_when_open_targets_errors(monkeypatch):
    def boom(oid, top_n=50):
        raise OpenTargetsError("Unknown disease id: HP_0002354")
    monkeypatch.setattr(app_module, "get_disease_top_genes", boom)
    monkeypatch.setattr(app_module.hpo, "get_phenotype_genes",
                        lambda oid, top_n=50: [{"symbol": "APP", "score": 1.0}])
    assert app_module._hybrid_gene_fetcher("HP_0002354", 50)["source"] == "hpo"


def test_both_sources_empty_yields_nothing(monkeypatch):
    monkeypatch.setattr(app_module, "get_disease_top_genes", lambda oid, top_n=50: [])
    monkeypatch.setattr(app_module.hpo, "get_phenotype_genes", lambda oid, top_n=50: [])
    assert app_module._hybrid_gene_fetcher("HP_0002354", 50)["genes"] == []
