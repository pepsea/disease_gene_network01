import pytest

from collectors import ppi
from tests.helpers import FakeResponse, FakeSession


@pytest.fixture(autouse=True)
def clear_signor_cache():
    ppi.reset_signor_cache()
    yield
    ppi.reset_signor_cache()


def signor_rows(*pairs):
    return [{"ENTITYA": a, "ENTITYB": b} for a, b in pairs]


def string_rows(*triples):
    return [{"preferredName_A": a, "preferredName_B": b, "score": s} for a, b, s in triples]


def biogrid_rows(*pairs):
    return {str(i): {"OFFICIAL_SYMBOL_A": a, "OFFICIAL_SYMBOL_B": b}
            for i, (a, b) in enumerate(pairs)}


def route(signor=None, string=None, biogrid=None):
    """Dispatch a fake GET to the right per-source payload."""
    def handler(url, **kwargs):
        if url == ppi.SIGNOR_URL:
            if signor is None:
                raise ConnectionError("signor down")
            return FakeResponse(signor)
        if url == ppi.STRING_URL:
            if string is None:
                raise ConnectionError("string down")
            return FakeResponse(string)
        if url == ppi.BIOGRID_URL:
            if biogrid is None:
                raise ConnectionError("biogrid down")
            return FakeResponse(biogrid)
        raise AssertionError(f"unexpected url {url}")
    return handler


def test_partners_are_pooled_across_all_three_sources(monkeypatch):
    monkeypatch.setattr(ppi, "SESSION", FakeSession(get=route(
        signor=signor_rows(("APP", "SIG1")),
        string=string_rows(("APP", "STR1", 0.9)),
        biogrid=biogrid_rows(("APP", "BIO1")),
    )))
    assert set(ppi.get_ppi_partners("APP", biogrid_key="k")) == {"SIG1", "STR1", "BIO1"}


def test_interactions_are_read_in_both_directions(monkeypatch):
    monkeypatch.setattr(ppi, "SESSION", FakeSession(get=route(
        signor=signor_rows(("PARTNER_A", "APP"), ("APP", "PARTNER_B")),
        string=[],
    )))
    assert set(ppi.get_ppi_partners("APP")) == {"PARTNER_A", "PARTNER_B"}


def test_the_gene_itself_is_never_a_partner(monkeypatch):
    monkeypatch.setattr(ppi, "SESSION", FakeSession(get=route(
        signor=signor_rows(("APP", "app"), ("APP", "PSEN1")),
        string=string_rows(("APP", "App", 0.9)),
    )))
    assert ppi.get_ppi_partners("APP") == ["PSEN1"]


def test_symbols_are_upper_cased(monkeypatch):
    monkeypatch.setattr(ppi, "SESSION", FakeSession(get=route(
        signor=signor_rows(("app", "psen1")), string=[],
    )))
    assert ppi.get_ppi_partners("APP") == ["PSEN1"]


def test_one_dead_source_does_not_break_the_others(monkeypatch):
    monkeypatch.setattr(ppi, "SESSION", FakeSession(get=route(
        signor=None,  # raises
        string=string_rows(("APP", "STR1", 0.8)),
    )))
    assert ppi.get_ppi_partners("APP") == ["STR1"]


def test_all_sources_down_returns_empty(monkeypatch):
    monkeypatch.setattr(ppi, "SESSION", FakeSession(get=route()))
    assert ppi.get_ppi_partners("APP", biogrid_key="k") == []


def test_biogrid_is_skipped_without_a_key(monkeypatch):
    session = FakeSession(get=route(signor=[], string=[], biogrid=biogrid_rows(("APP", "BIO1"))))
    monkeypatch.setattr(ppi, "SESSION", session)
    assert ppi.get_ppi_partners("APP") == []
    assert all(call["url"] != ppi.BIOGRID_URL for call in session.get_calls)


def test_multi_source_partners_outrank_single_source_ones(monkeypatch):
    monkeypatch.setattr(ppi, "SESSION", FakeSession(get=route(
        signor=signor_rows(("APP", "BOTH")),
        string=string_rows(("APP", "BOTH", 0.4), ("APP", "STRING_ONLY", 0.99)),
    )))
    assert ppi.get_ppi_partners("APP", top_n=1) == ["BOTH"]


def test_string_score_breaks_ties_within_a_support_level(monkeypatch):
    monkeypatch.setattr(ppi, "SESSION", FakeSession(get=route(
        signor=[],
        string=string_rows(("APP", "LOW", 0.2), ("APP", "HIGH", 0.95), ("APP", "MID", 0.5)),
    )))
    assert ppi.get_ppi_partners("APP", top_n=2) == ["HIGH", "MID"]


def test_results_are_capped_at_top_n_and_deterministic(monkeypatch):
    rows = signor_rows(*[("APP", f"P{i:02d}") for i in range(50)])
    monkeypatch.setattr(ppi, "SESSION", FakeSession(get=route(signor=rows, string=[])))
    first = ppi.get_ppi_partners("APP", top_n=5)
    ppi.reset_signor_cache()
    second = ppi.get_ppi_partners("APP", top_n=5)
    assert first == second
    assert len(first) == 5


def test_signor_bulk_download_happens_once_across_genes(monkeypatch):
    session = FakeSession(get=route(signor=signor_rows(("APP", "P1"), ("PSEN1", "P2")), string=[]))
    monkeypatch.setattr(ppi, "SESSION", session)
    ppi.get_ppi_partners("APP")
    ppi.get_ppi_partners("PSEN1")
    ppi.get_ppi_partners("MAPT")
    signor_calls = [c for c in session.get_calls if c["url"] == ppi.SIGNOR_URL]
    assert len(signor_calls) == 1


def test_signor_index_is_shared_between_threads(monkeypatch):
    import concurrent.futures

    session = FakeSession(get=route(signor=signor_rows(("APP", "P1")), string=[]))
    monkeypatch.setattr(ppi, "SESSION", session)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        list(ex.map(lambda g: ppi.get_ppi_partners(g), ["APP"] * 8))
    assert len([c for c in session.get_calls if c["url"] == ppi.SIGNOR_URL]) == 1


def test_signor_can_be_disabled_by_env(monkeypatch):
    monkeypatch.setenv("NW_ENABLE_SIGNOR", "0")
    session = FakeSession(get=route(string=string_rows(("APP", "STR1", 0.5))))
    monkeypatch.setattr(ppi, "SESSION", session)
    assert ppi.get_ppi_partners("APP") == ["STR1"]
    assert all(c["url"] != ppi.SIGNOR_URL for c in session.get_calls)


def test_blank_gene_makes_no_requests(monkeypatch):
    session = FakeSession(get=route())
    monkeypatch.setattr(ppi, "SESSION", session)
    assert ppi.get_ppi_partners("   ") == []
    assert session.get_calls == []


def test_self_interaction_rows_are_ignored(monkeypatch):
    monkeypatch.setattr(ppi, "SESSION", FakeSession(get=route(
        signor=signor_rows(("APP", "APP")), string=[],
    )))
    assert ppi.get_ppi_partners("APP") == []


def test_malformed_rows_are_tolerated(monkeypatch):
    monkeypatch.setattr(ppi, "SESSION", FakeSession(get=route(
        signor=[{"ENTITYA": "APP"}, {}, "junk", {"ENTITYA": "APP", "ENTITYB": "OK"}],
        string=[{"preferredName_A": "APP", "preferredName_B": "S1", "score": "bad"}],
    )))
    assert set(ppi.get_ppi_partners("APP")) == {"OK", "S1"}


def test_human_taxon_is_requested(monkeypatch):
    session = FakeSession(get=route(signor=[], string=[], biogrid={}))
    monkeypatch.setattr(ppi, "SESSION", session)
    ppi.get_ppi_partners("APP", biogrid_key="k")
    by_url = {c["url"]: c["params"] for c in session.get_calls}
    assert by_url[ppi.SIGNOR_URL]["organism"] == "9606"
    assert by_url[ppi.STRING_URL]["species"] == 9606
    assert by_url[ppi.BIOGRID_URL]["taxId"] == 9606
    assert by_url[ppi.BIOGRID_URL]["accesskey"] == "k"
