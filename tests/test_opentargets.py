import pytest

from collectors import opentargets as ot
from tests.helpers import FakeResponse, FakeSession


def search_payload(hits):
    return {"data": {"search": {"hits": hits}}}


def genes_payload(rows, name="Alzheimer disease", disease_id="EFO_0000249"):
    return {
        "data": {
            "disease": {
                "id": disease_id,
                "name": name,
                "associatedTargets": {"count": len(rows), "rows": rows},
            }
        }
    }


def row(symbol, score, target_id=None):
    return {"target": {"id": target_id or f"ENSG_{symbol}", "approvedSymbol": symbol}, "score": score}


# --- resolve_disease_id ----------------------------------------------------

def test_resolve_prefers_an_exact_name_match_over_search_ranking(monkeypatch):
    hits = [
        {"id": "EFO_9999", "name": "Alzheimer disease 2", "entity": "disease"},
        {"id": "EFO_0000249", "name": "Alzheimer disease", "entity": "disease"},
    ]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(search_payload(hits))))
    assert ot.resolve_disease_id("alzheimer DISEASE") == ("EFO_0000249", "Alzheimer disease")


def test_resolve_falls_back_to_the_top_hit(monkeypatch):
    hits = [{"id": "EFO_1", "name": "Alzheimer disease 2", "entity": "disease"}]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(search_payload(hits))))
    assert ot.resolve_disease_id("alzheimers") == ("EFO_1", "Alzheimer disease 2")


def test_resolve_ignores_non_disease_entities(monkeypatch):
    hits = [
        {"id": "ENSG000", "name": "APP", "entity": "target"},
        {"id": "EFO_1", "name": "amyloidosis", "entity": "disease"},
    ]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(search_payload(hits))))
    assert ot.resolve_disease_id("APP") == ("EFO_1", "amyloidosis")


def test_resolve_returns_none_when_nothing_matches(monkeypatch):
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(search_payload([]))))
    assert ot.resolve_disease_id("not a disease") == (None, "not a disease")


def test_resolve_returns_none_for_blank_input():
    assert ot.resolve_disease_id("   ") == (None, "")


def test_an_ontology_id_skips_search_and_is_labelled(monkeypatch):
    session = FakeSession(
        post=lambda url, **kw: FakeResponse(
            {"data": {"disease": {"id": "EFO_0000249", "name": "Alzheimer disease"}}}
        )
    )
    monkeypatch.setattr(ot, "SESSION", session)
    assert ot.resolve_disease_id("efo_0000249") == ("EFO_0000249", "Alzheimer disease")
    assert "SearchEntity" not in session.post_calls[0]["json"]["query"]


def test_graphql_errors_raise(monkeypatch):
    payload = {"errors": [{"message": "boom"}]}
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(payload)))
    with pytest.raises(ot.OpenTargetsError, match="boom"):
        ot.resolve_disease_id("Alzheimer disease")


def test_network_failure_raises_opentargets_error(monkeypatch):
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=ConnectionError("down")))
    with pytest.raises(ot.OpenTargetsError, match="request failed"):
        ot.resolve_disease_id("Alzheimer disease")


def test_non_json_response_raises_opentargets_error(monkeypatch):
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(text="<html>502</html>")))
    with pytest.raises(ot.OpenTargetsError):
        ot.resolve_disease_id("Alzheimer disease")


def test_http_error_raises_opentargets_error(monkeypatch):
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(status=503)))
    with pytest.raises(ot.OpenTargetsError, match="request failed"):
        ot.resolve_disease_id("Alzheimer disease")


# --- get_disease_top_genes -------------------------------------------------

def test_top_genes_are_parsed_in_order(monkeypatch):
    rows = [row("APP", 0.9), row("PSEN1", 0.6), row("MAPT", 0.5)]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(genes_payload(rows))))
    genes = ot.get_disease_top_genes("EFO_0000249", top_n=100)
    assert [g["symbol"] for g in genes] == ["APP", "PSEN1", "MAPT"]
    assert genes[0]["score"] == 0.9
    assert genes[0]["target_id"] == "ENSG_APP"


def test_top_n_is_passed_and_clamped(monkeypatch):
    session = FakeSession(post=FakeResponse(genes_payload([])))
    monkeypatch.setattr(ot, "SESSION", session)
    ot.get_disease_top_genes("EFO_1", top_n=10_000)
    assert session.post_calls[0]["json"]["variables"]["size"] == ot.MAX_PAGE_SIZE
    ot.get_disease_top_genes("EFO_1", top_n=0)
    assert session.post_calls[1]["json"]["variables"]["size"] == 1


def test_results_are_truncated_to_top_n(monkeypatch):
    rows = [row(f"G{i}", 1.0 / (i + 1)) for i in range(10)]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(genes_payload(rows))))
    assert len(ot.get_disease_top_genes("EFO_1", top_n=3)) == 3


def test_rows_without_a_symbol_are_skipped(monkeypatch):
    rows = [row("APP", 0.9), {"target": {"id": "X"}, "score": 0.5}, {"score": 0.4}]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(genes_payload(rows))))
    assert [g["symbol"] for g in ot.get_disease_top_genes("EFO_1")] == ["APP"]


def test_null_scores_become_zero(monkeypatch):
    rows = [{"target": {"id": "X", "approvedSymbol": "APP"}, "score": None}]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(genes_payload(rows))))
    assert ot.get_disease_top_genes("EFO_1")[0]["score"] == 0.0


def test_unknown_disease_id_raises(monkeypatch):
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse({"data": {"disease": None}})))
    with pytest.raises(ot.OpenTargetsError, match="Unknown disease id"):
        ot.get_disease_top_genes("EFO_nope")


def test_disease_with_no_associations_returns_empty(monkeypatch):
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(genes_payload([]))))
    assert ot.get_disease_top_genes("EFO_1") == []


# --- search_diseases -------------------------------------------------------

def test_search_returns_candidates_with_exact_first(monkeypatch):
    hits = [
        {"id": "EFO_1", "name": "Alzheimer disease 2", "entity": "disease", "description": "d2"},
        {"id": "EFO_0000249", "name": "Alzheimer disease", "entity": "disease", "description": "d1"},
    ]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(search_payload(hits))))
    results = ot.search_diseases("alzheimer disease")
    assert [r["id"] for r in results] == ["EFO_0000249", "EFO_1"]
    assert results[0]["exact"] is True
    assert results[0]["description"] == "d1"


def test_search_preserves_ranking_when_nothing_matches_exactly(monkeypatch):
    hits = [
        {"id": "EFO_1", "name": "Alzheimer disease 2", "entity": "disease"},
        {"id": "EFO_2", "name": "Alzheimer disease 3", "entity": "disease"},
    ]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(search_payload(hits))))
    assert [r["id"] for r in ot.search_diseases("alzheimer")] == ["EFO_1", "EFO_2"]


def test_search_drops_non_disease_entities(monkeypatch):
    hits = [
        {"id": "ENSG1", "name": "APP", "entity": "target"},
        {"id": "EFO_1", "name": "amyloidosis", "entity": "disease"},
    ]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(search_payload(hits))))
    assert [r["id"] for r in ot.search_diseases("APP")] == ["EFO_1"]


def test_search_respects_the_limit(monkeypatch):
    hits = [{"id": f"EFO_{i}", "name": f"d{i}", "entity": "disease"} for i in range(20)]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(search_payload(hits))))
    assert len(ot.search_diseases("d", limit=5)) == 5


def test_search_of_a_blank_query_makes_no_request(monkeypatch):
    session = FakeSession(post=FakeResponse(search_payload([])))
    monkeypatch.setattr(ot, "SESSION", session)
    assert ot.search_diseases("   ") == []
    assert session.post_calls == []


def test_search_falls_back_when_description_is_unsupported(monkeypatch):
    """A schema without `description` must not break disease search."""
    calls = []

    def handler(url, **kwargs):
        query = kwargs["json"]["query"]
        calls.append(query)
        if "description" in query:
            return FakeResponse({"errors": [{"message": "Cannot query field description"}]})
        return FakeResponse(search_payload(
            [{"id": "EFO_1", "name": "Alzheimer disease", "entity": "disease"}]))

    monkeypatch.setattr(ot, "SESSION", FakeSession(post=handler))
    monkeypatch.setattr(ot, "_search_query_index", None)

    assert [r["id"] for r in ot.search_diseases("alzheimer")] == ["EFO_1"]
    assert len(calls) == 2  # rich query, then the one without description

    # The rejected variant is not retried afterwards.
    ot.search_diseases("alzheimer")
    assert len(calls) == 3


def test_search_falls_back_when_paging_is_unsupported(monkeypatch):
    """A schema without `page` on search must not break it either."""
    calls = []

    def handler(url, **kwargs):
        query = kwargs["json"]["query"]
        calls.append(query)
        if "page:" in query.replace(" ", "") or "page: {" in query:
            return FakeResponse({"errors": [{"message": "Unknown argument page"}]})
        return FakeResponse(search_payload(
            [{"id": "EFO_1", "name": "Alzheimer disease", "entity": "disease"}]))

    monkeypatch.setattr(ot, "SESSION", FakeSession(post=handler))
    monkeypatch.setattr(ot, "_search_query_index", None)
    assert [r["id"] for r in ot.search_diseases("alzheimer")] == ["EFO_1"]
    # Both paged variants rejected, then the unpaged one works.
    assert len(calls) == 3


def test_search_raises_when_every_variant_fails(monkeypatch):
    monkeypatch.setattr(ot, "SESSION", FakeSession(
        post=FakeResponse({"errors": [{"message": "boom"}]})))
    monkeypatch.setattr(ot, "_search_query_index", None)
    with pytest.raises(ot.OpenTargetsError):
        ot.search_diseases("alzheimer")


def test_the_requested_size_is_sent_as_the_page_size(monkeypatch):
    session = FakeSession(post=FakeResponse(search_payload([])))
    monkeypatch.setattr(ot, "SESSION", session)
    monkeypatch.setattr(ot, "_search_query_index", None)
    ot.search_diseases("alzheimer", limit=120)
    assert session.post_calls[0]["json"]["variables"]["size"] == 120
    ot.search_diseases("alzheimer", limit=9999)
    assert session.post_calls[1]["json"]["variables"]["size"] == 500


def test_ontology_id_helpers():
    assert ot.is_ontology_id("EFO_0000249")
    assert ot.is_ontology_id("mondo_0004975")
    assert not ot.is_ontology_id("Alzheimer disease")
    assert ot.normalise_ontology_id("efo_0000249") == "EFO_0000249"
    assert ot.normalise_ontology_id("orphanet_123") == "Orphanet_123"


def test_get_disease_with_top_genes_returns_the_label(monkeypatch):
    rows = [row("APP", 0.9)]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(genes_payload(rows))))
    label, genes = ot.get_disease_with_top_genes("EFO_0000249")
    assert label == "Alzheimer disease"
    assert [g["symbol"] for g in genes] == ["APP"]


# --- get_disease_phenotypes ------------------------------------------------

def phenotype_payload(rows, name="Alzheimer disease"):
    return {"data": {"disease": {"id": "EFO_0000249", "name": name,
                                 "phenotypes": {"count": len(rows), "rows": rows}}}}


def hpo_row(hpo_id, name, evidence=None):
    return {"phenotypeHPO": {"id": hpo_id, "name": name, "description": "d"},
            "evidence": evidence if evidence is not None else
            [{"aspect": "P", "frequency": "HP:0040281", "qualifierNot": False,
              "resource": "HPO"}]}


@pytest.fixture(autouse=True)
def reset_phenotype_variant(monkeypatch):
    monkeypatch.setattr(ot, "_phenotype_query_index", None)


def test_phenotypes_are_parsed(monkeypatch):
    rows = [hpo_row("HP:0002354", "Memory impairment")]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(phenotype_payload(rows))))
    p = ot.get_disease_phenotypes("EFO_0000249")[0]
    assert p["hpo_id"] == "HP:0002354"
    assert p["ontology_id"] == "HP_0002354"
    assert p["name"] == "Memory impairment"
    assert p["frequency"] == "HP:0040281"
    assert p["aspect"] == "P"
    assert p["excluded"] is False


def test_hpo_ids_are_converted_to_the_open_targets_form():
    assert ot.hpo_to_ontology_id("HP:0002354") == "HP_0002354"
    assert ot.hpo_to_ontology_id("hp:0002354") == "HP_0002354"
    assert ot.hpo_to_ontology_id("") == ""


def test_a_phenotype_all_sources_call_absent_is_excluded(monkeypatch):
    rows = [hpo_row("HP:1", "Absent finding",
                    [{"qualifierNot": True, "aspect": "P"},
                     {"qualifierNot": True, "aspect": "P"}])]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(phenotype_payload(rows))))
    assert ot.get_disease_phenotypes("EFO_1")[0]["excluded"] is True


def test_a_phenotype_only_some_sources_call_absent_is_kept(monkeypatch):
    rows = [hpo_row("HP:1", "Contested finding",
                    [{"qualifierNot": True}, {"qualifierNot": False, "frequency": "HP:0040282"}])]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(phenotype_payload(rows))))
    p = ot.get_disease_phenotypes("EFO_1")[0]
    assert p["excluded"] is False
    assert p["frequency"] == "HP:0040282"


def test_rows_without_an_hpo_id_are_skipped(monkeypatch):
    rows = [hpo_row("HP:1", "ok"), {"phenotypeHPO": {"name": "no id"}}, {}]
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(phenotype_payload(rows))))
    assert [p["hpo_id"] for p in ot.get_disease_phenotypes("EFO_1")] == ["HP:1"]


def test_a_disease_with_no_phenotypes_returns_empty(monkeypatch):
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse(phenotype_payload([]))))
    assert ot.get_disease_phenotypes("EFO_1") == []


def test_unknown_disease_returns_empty_rather_than_raising(monkeypatch):
    monkeypatch.setattr(ot, "SESSION", FakeSession(post=FakeResponse({"data": {"disease": None}})))
    assert ot.get_disease_phenotypes("EFO_nope") == []


def test_phenotype_query_falls_back_when_evidence_is_unsupported(monkeypatch):
    """A schema without the evidence fields must still yield the symptom list."""
    calls = []

    def handler(url, **kwargs):
        query = kwargs["json"]["query"]
        calls.append(query)
        if "evidence" in query:
            return FakeResponse({"errors": [{"message": "Cannot query field evidence"}]})
        return FakeResponse(phenotype_payload(
            [{"phenotypeHPO": {"id": "HP:1", "name": "Memory impairment"}}]))

    monkeypatch.setattr(ot, "SESSION", FakeSession(post=handler))
    result = ot.get_disease_phenotypes("EFO_1")
    assert [p["name"] for p in result] == ["Memory impairment"]
    assert result[0]["frequency"] == ""      # unavailable, not invented
    assert len(calls) == 2                   # rich query, then the fallback


def test_the_working_phenotype_variant_is_remembered(monkeypatch):
    calls = []

    def handler(url, **kwargs):
        query = kwargs["json"]["query"]
        calls.append(query)
        if "evidence" in query:
            return FakeResponse({"errors": [{"message": "no evidence field"}]})
        return FakeResponse(phenotype_payload([{"phenotypeHPO": {"id": "HP:1", "name": "x"}}]))

    monkeypatch.setattr(ot, "SESSION", FakeSession(post=handler))
    ot.get_disease_phenotypes("EFO_1")
    ot.get_disease_phenotypes("EFO_2")
    assert len(calls) == 3   # 2 for the first call, 1 for the second


def test_phenotypes_unavailable_entirely_returns_empty(monkeypatch):
    monkeypatch.setattr(ot, "SESSION", FakeSession(
        post=FakeResponse({"errors": [{"message": "Cannot query field phenotypes"}]})))
    assert ot.get_disease_phenotypes("EFO_1") == []


def test_phenotype_limit_is_clamped(monkeypatch):
    session = FakeSession(post=FakeResponse(phenotype_payload([])))
    monkeypatch.setattr(ot, "SESSION", session)
    ot.get_disease_phenotypes("EFO_1", limit=10_000)
    assert session.post_calls[0]["json"]["variables"]["size"] == ot.MAX_PAGE_SIZE
