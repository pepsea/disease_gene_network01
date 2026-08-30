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
