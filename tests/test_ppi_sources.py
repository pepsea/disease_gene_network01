"""Tests for the individual PPI source collectors."""
import pytest

from collectors import biogrid, signor, string_db
from tests.helpers import FakeResponse, FakeSession


@pytest.fixture(autouse=True)
def no_disk_cache(monkeypatch):
    monkeypatch.setenv("PPI_CACHE_DISABLED", "1")
    signor.reset_cache()
    yield
    signor.reset_cache()


# --- SIGNOR ----------------------------------------------------------------

def tsv(*rows):
    """Build a SIGNOR TSV. Each row: (entityA, typeA, entityB, typeB, effect, score)."""
    lines = []
    for a, ta, b, tb, effect, score in rows:
        cols = [""] * 24
        cols[0], cols[1], cols[4], cols[5], cols[8] = a, ta, b, tb, effect
        cols[9], cols[17], cols[23] = "phosphorylation", "12345", str(score)
        lines.append("\t".join(cols))
    return "\n".join(lines)


def signor_session(tsv_text):
    """A session whose GET returns the given TSV as the response body."""
    return FakeSession(get=lambda url, **kw: FakeResponse(body=tsv_text))


def test_signor_finds_the_gene_on_either_side(monkeypatch):
    monkeypatch.setattr(signor, "SESSION", signor_session(
        tsv(("APP", "protein", "PSEN1", "protein", "up-regulates", 0.7),
            ("MAPT", "protein", "APP", "protein", "binds", 0.5))))
    partners = {i["partner"].upper() for i in signor.get_interactions("APP")}
    assert partners == {"PSEN1", "MAPT"}


def test_signor_records_direction(monkeypatch):
    monkeypatch.setattr(signor, "SESSION", signor_session(
        tsv(("APP", "protein", "PSEN1", "protein", "up-regulates", 0.7))))
    assert signor.get_interactions("APP")[0]["direction"] == "→"
    assert signor.get_interactions("PSEN1")[0]["direction"] == "←"


def test_signor_excludes_non_protein_partners(monkeypatch):
    monkeypatch.setattr(signor, "SESSION", signor_session(
        tsv(("APP", "protein", "ASPIRIN", "chemical", "down-regulates", 0.4),
            ("APP", "protein", "PSEN1", "protein", "binds", 0.7))))
    assert [i["partner"] for i in signor.get_interactions("APP")] == ["PSEN1"]


def test_signor_keeps_complexes(monkeypatch):
    monkeypatch.setattr(signor, "SESSION", signor_session(
        tsv(("APP", "protein", "GAMMA_SEC", "complex", "binds", 0.6))))
    assert signor.get_interactions("APP")[0]["partner_type"] == "complex"


def test_signor_parses_scores(monkeypatch):
    monkeypatch.setattr(signor, "SESSION", signor_session(
        tsv(("APP", "protein", "PSEN1", "protein", "binds", 0.73))))
    assert signor.get_interactions("APP")[0]["score"] == 0.73


def test_signor_tolerates_short_and_malformed_rows(monkeypatch):
    text = "too\tshort\n" + tsv(("APP", "protein", "PSEN1", "protein", "binds", 0.7))
    monkeypatch.setattr(signor, "SESSION", signor_session(text))
    assert [i["partner"] for i in signor.get_interactions("APP")] == ["PSEN1"]


def test_signor_downloads_once_for_many_genes(monkeypatch):
    session = signor_session(tsv(("APP", "protein", "PSEN1", "protein", "binds", 0.7)))
    monkeypatch.setattr(signor, "SESSION", session)
    signor.get_interactions("APP")
    signor.get_interactions("PSEN1")
    signor.get_interactions("MAPT")
    assert len(session.get_calls) == 1


def test_signor_lookup_is_case_insensitive(monkeypatch):
    monkeypatch.setattr(signor, "SESSION", signor_session(
        tsv(("app", "protein", "psen1", "protein", "binds", 0.7))))
    assert len(signor.get_interactions("APP")) == 1


# --- STRING ----------------------------------------------------------------

def string_row(a, b, score, **extra):
    return {"preferredName_A": a, "preferredName_B": b, "score": score, **extra}


def test_string_keeps_only_edges_touching_the_query_gene(monkeypatch):
    monkeypatch.setattr(string_db, "SESSION", FakeSession(get=FakeResponse([
        string_row("APP", "PSEN1", 0.9),
        string_row("PSEN1", "MAPT", 0.8),  # partner-to-partner edge
    ])))
    assert [i["partner"] for i in string_db.get_interactions("APP")] == ["PSEN1"]


def test_string_keeps_the_best_score_per_partner(monkeypatch):
    monkeypatch.setattr(string_db, "SESSION", FakeSession(get=FakeResponse([
        string_row("APP", "PSEN1", 0.4), string_row("PSEN1", "APP", 0.92),
    ])))
    interactions = string_db.get_interactions("APP")
    assert len(interactions) == 1
    assert interactions[0]["score"] == 0.92


def test_string_results_are_sorted_by_score(monkeypatch):
    monkeypatch.setattr(string_db, "SESSION", FakeSession(get=FakeResponse([
        string_row("APP", "LOW", 0.2), string_row("APP", "HIGH", 0.95),
    ])))
    assert [i["partner"] for i in string_db.get_interactions("APP")] == ["HIGH", "LOW"]


def test_string_required_score_is_sent(monkeypatch):
    session = FakeSession(get=FakeResponse([]))
    monkeypatch.setattr(string_db, "SESSION", session)
    string_db.get_interactions("APP", required_score=900)
    assert session.get_calls[0]["params"]["required_score"] == 900
    assert session.get_calls[0]["params"]["species"] == 9606


def test_string_captures_subscores(monkeypatch):
    monkeypatch.setattr(string_db, "SESSION", FakeSession(get=FakeResponse([
        string_row("APP", "PSEN1", 0.9, escore=0.4, dscore=0.5, tscore=0.1),
    ])))
    subs = string_db.get_interactions("APP")[0]["subscores"]
    assert subs["experimental"] == 0.4 and subs["database"] == 0.5


def test_string_404_returns_empty(monkeypatch):
    monkeypatch.setattr(string_db, "SESSION", FakeSession(get=FakeResponse([], status=404)))
    assert string_db.get_interactions("NOSUCHGENE") == []


def test_string_outage_returns_empty(monkeypatch):
    monkeypatch.setattr(string_db, "SESSION", FakeSession(get=ConnectionError("down")))
    assert string_db.get_interactions("APP") == []


def test_string_self_edges_are_dropped(monkeypatch):
    monkeypatch.setattr(string_db, "SESSION", FakeSession(get=FakeResponse([
        string_row("APP", "APP", 0.99), string_row("APP", "PSEN1", 0.5),
    ])))
    assert [i["partner"] for i in string_db.get_interactions("APP")] == ["PSEN1"]


# --- BioGRID ---------------------------------------------------------------

def biogrid_payload(*entries):
    return {str(i): e for i, e in enumerate(entries)}


def bg(a, b, score=None, system="Affinity Capture-MS", pmid="1"):
    return {"OFFICIAL_SYMBOL_A": a, "OFFICIAL_SYMBOL_B": b, "SCORE": score,
            "EXPERIMENTAL_SYSTEM": system, "EXPERIMENTAL_SYSTEM_TYPE": "physical",
            "PUBMED_ID": pmid}


def test_biogrid_needs_a_key(monkeypatch):
    monkeypatch.delenv("BIOGRID_API_KEY", raising=False)
    session = FakeSession(get=FakeResponse({}))
    monkeypatch.setattr(biogrid, "SESSION", session)
    assert biogrid.get_interactions("APP") == []
    assert session.get_calls == []


def test_biogrid_sends_the_key_and_human_taxon(monkeypatch):
    session = FakeSession(get=FakeResponse({}))
    monkeypatch.setattr(biogrid, "SESSION", session)
    biogrid.get_interactions("APP", api_key="secret")
    params = session.get_calls[0]["params"]
    assert params["accessKey"] == "secret"
    assert params["taxId"] == "9606"
    assert params["selfInteractionsExcluded"] == "true"


def test_biogrid_falls_back_to_the_environment_key(monkeypatch):
    monkeypatch.setenv("BIOGRID_API_KEY", "from-env")
    session = FakeSession(get=FakeResponse({}))
    monkeypatch.setattr(biogrid, "SESSION", session)
    biogrid.get_interactions("APP")
    assert session.get_calls[0]["params"]["accessKey"] == "from-env"


def test_biogrid_collapses_duplicate_records_per_partner(monkeypatch):
    monkeypatch.setattr(biogrid, "SESSION", FakeSession(get=FakeResponse(biogrid_payload(
        bg("APP", "PSEN1", score=0.2, pmid="1"),
        bg("APP", "PSEN1", score=0.8, pmid="2"),
        bg("APP", "PSEN1", score=None, pmid="3"),
    ))))
    interactions = biogrid.get_interactions("APP", api_key="k")
    assert len(interactions) == 1
    assert interactions[0]["score"] == 0.8


def test_biogrid_reads_both_orientations(monkeypatch):
    monkeypatch.setattr(biogrid, "SESSION", FakeSession(get=FakeResponse(biogrid_payload(
        bg("PSEN1", "APP"), bg("APP", "MAPT"),
    ))))
    assert {i["partner"].upper() for i in biogrid.get_interactions("APP", api_key="k")} == {"PSEN1", "MAPT"}


def test_biogrid_outage_returns_empty(monkeypatch):
    monkeypatch.setattr(biogrid, "SESSION", FakeSession(get=ConnectionError("down")))
    assert biogrid.get_interactions("APP", api_key="k") == []


def test_biogrid_handles_missing_scores(monkeypatch):
    monkeypatch.setattr(biogrid, "SESSION", FakeSession(get=FakeResponse(biogrid_payload(
        bg("APP", "PSEN1", score="-"), bg("APP", "MAPT", score=""),
    ))))
    assert all(i["score"] is None for i in biogrid.get_interactions("APP", api_key="k"))
