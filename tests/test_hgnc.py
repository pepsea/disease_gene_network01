import pytest

from collectors import hgnc
from tests.helpers import FakeResponse, FakeSession


@pytest.fixture(autouse=True)
def clear_cache():
    hgnc.reset_cache()
    yield
    hgnc.reset_cache()


def docs(*entries):
    return {"response": {"numFound": len(entries), "docs": list(entries)}}


def doc(symbol, name="", status="Approved", hgnc_id=None):
    return {"symbol": symbol, "name": name, "status": status,
            "hgnc_id": hgnc_id or f"HGNC:{symbol}", "ensembl_gene_id": f"ENSG_{symbol}"}


def route(**by_field):
    """Serve /fetch/<field>/<value> from a {field: {value: docs}} mapping."""
    def handler(url, **kwargs):
        _, field, value = url.rsplit("/", 2)
        table = by_field.get(field)
        if table is None:
            return FakeResponse(docs())
        if isinstance(table, Exception):
            raise table
        return FakeResponse(table.get(value.upper(), docs()))
    return handler


def test_approved_symbol_is_returned_as_is(monkeypatch):
    monkeypatch.setattr(hgnc, "SESSION", FakeSession(
        get=route(symbol={"TP53": docs(doc("TP53", "tumor protein p53"))})))
    r = hgnc.resolve_gene_symbol("TP53")
    assert r["symbol"] == "TP53"
    assert r["status"] == "approved"
    assert r["hgnc_id"] == "HGNC:TP53"
    assert r["name"] == "tumor protein p53"


def test_alias_is_mapped_to_the_approved_symbol(monkeypatch):
    monkeypatch.setattr(hgnc, "SESSION", FakeSession(
        get=route(alias_symbol={"PS1": docs(doc("PSEN1", "presenilin 1"))})))
    r = hgnc.resolve_gene_symbol("PS1")
    assert (r["input"], r["symbol"], r["status"]) == ("PS1", "PSEN1", "alias")


def test_previous_symbol_is_mapped(monkeypatch):
    monkeypatch.setattr(hgnc, "SESSION", FakeSession(
        get=route(prev_symbol={"MLL": docs(doc("KMT2A"))})))
    r = hgnc.resolve_gene_symbol("MLL")
    assert (r["symbol"], r["status"]) == ("KMT2A", "previous")


def test_approved_wins_over_alias(monkeypatch):
    # A symbol that is approved for one gene and an alias for another.
    monkeypatch.setattr(hgnc, "SESSION", FakeSession(get=route(
        symbol={"ABC": docs(doc("ABC"))},
        alias_symbol={"ABC": docs(doc("XYZ"))},
    )))
    assert hgnc.resolve_gene_symbol("ABC")["status"] == "approved"


def test_unmatched_symbol_is_unknown(monkeypatch):
    monkeypatch.setattr(hgnc, "SESSION", FakeSession(get=route()))
    r = hgnc.resolve_gene_symbol("NOTAGENE")
    assert (r["symbol"], r["status"]) == ("NOTAGENE", "unknown")


def test_withdrawn_entries_are_flagged(monkeypatch):
    monkeypatch.setattr(hgnc, "SESSION", FakeSession(
        get=route(symbol={"OLD1": docs(doc("OLD1", status="Entry Withdrawn"))})))
    assert hgnc.resolve_gene_symbol("OLD1")["status"] == "withdrawn"


def test_ambiguous_alias_lists_candidates(monkeypatch):
    monkeypatch.setattr(hgnc, "SESSION", FakeSession(
        get=route(alias_symbol={"AMB": docs(doc("GENE1"), doc("GENE2"))})))
    r = hgnc.resolve_gene_symbol("AMB")
    assert r["symbol"] == "GENE1"
    assert r["candidates"] == ["GENE1", "GENE2"]


def test_hgnc_outage_does_not_block_the_analysis(monkeypatch):
    monkeypatch.setattr(hgnc, "SESSION", FakeSession(get=ConnectionError("hgnc down")))
    r = hgnc.resolve_gene_symbol("TP53")
    # The symbol passes through untouched so the run can still proceed.
    assert (r["symbol"], r["status"]) == ("TP53", "unverified")


def test_unverified_results_are_not_cached(monkeypatch):
    monkeypatch.setattr(hgnc, "SESSION", FakeSession(get=ConnectionError("down")))
    assert hgnc.resolve_gene_symbol("TP53")["status"] == "unverified"
    monkeypatch.setattr(hgnc, "SESSION", FakeSession(
        get=route(symbol={"TP53": docs(doc("TP53"))})))
    assert hgnc.resolve_gene_symbol("TP53")["status"] == "approved"


def test_resolutions_are_cached(monkeypatch):
    session = FakeSession(get=route(symbol={"TP53": docs(doc("TP53"))}))
    monkeypatch.setattr(hgnc, "SESSION", session)
    hgnc.resolve_gene_symbol("TP53")
    hgnc.resolve_gene_symbol("tp53")
    hgnc.resolve_gene_symbol("TP53")
    assert len(session.get_calls) == 1


def test_lookup_is_case_insensitive(monkeypatch):
    monkeypatch.setattr(hgnc, "SESSION", FakeSession(
        get=route(symbol={"TP53": docs(doc("TP53"))})))
    assert hgnc.resolve_gene_symbol("tp53")["symbol"] == "TP53"


def test_symbols_with_special_characters_are_url_encoded(monkeypatch):
    session = FakeSession(get=route(symbol={"HLA-A": docs(doc("HLA-A"))}))
    monkeypatch.setattr(hgnc, "SESSION", session)
    assert hgnc.resolve_gene_symbol("HLA-A")["symbol"] == "HLA-A"
    assert session.get_calls[0]["url"].endswith("/fetch/symbol/HLA-A")


def test_blank_symbol_makes_no_request(monkeypatch):
    session = FakeSession(get=route())
    monkeypatch.setattr(hgnc, "SESSION", session)
    assert hgnc.resolve_gene_symbol("  ")["status"] == "unknown"
    assert session.get_calls == []


def test_batch_resolution_deduplicates_by_approved_symbol(monkeypatch):
    monkeypatch.setattr(hgnc, "SESSION", FakeSession(get=route(
        symbol={"PSEN1": docs(doc("PSEN1")), "TP53": docs(doc("TP53"))},
        alias_symbol={"PS1": docs(doc("PSEN1"))},
    )))
    resolved = hgnc.resolve_gene_symbols(["PSEN1", "PS1", "TP53"])
    # PS1 resolves to PSEN1, which is already present, so it is dropped.
    assert [e["symbol"] for e in resolved] == ["PSEN1", "TP53"]


def test_batch_preserves_input_order(monkeypatch):
    monkeypatch.setattr(hgnc, "SESSION", FakeSession(get=route(
        symbol={"C": docs(doc("C")), "A": docs(doc("A")), "B": docs(doc("B"))})))
    assert [e["symbol"] for e in hgnc.resolve_gene_symbols(["C", "A", "B"])] == ["C", "A", "B"]


def test_malformed_payloads_are_tolerated(monkeypatch):
    monkeypatch.setattr(hgnc, "SESSION", FakeSession(get=lambda url, **kw: FakeResponse("junk")))
    assert hgnc.resolve_gene_symbol("TP53")["status"] == "unverified"
