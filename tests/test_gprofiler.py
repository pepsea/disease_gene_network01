import pytest

from collectors import gprofiler
from tests.helpers import FakeResponse, FakeSession


@pytest.fixture(autouse=True)
def no_disk_cache(monkeypatch):
    monkeypatch.setenv("PPI_CACHE_DISABLED", "1")


def payload(*terms, query_genes=None):
    """Build a g:Profiler response. Each term: (native, name, source, p, intersections)."""
    result = [
        {"native": n, "name": nm, "source": src, "p_value": p,
         "intersections": ix, "intersection_size": sum(1 for a in ix if a),
         "term_size": 100}
        for n, nm, src, p, ix in terms
    ]
    meta = {"query_metadata": {"queries": {"query_1": query_genes or []}}}
    return {"result": result, "meta": meta}


def test_terms_are_parsed_and_sorted_by_p_value(monkeypatch):
    monkeypatch.setattr(gprofiler, "SESSION", FakeSession(post=FakeResponse(payload(
        ("GO:1", "weak term", "GO:BP", 0.04, [[], ["x"]]),
        ("R-HSA-1", "strong term", "REAC", 1e-10, [["x"], ["x"]]),
        query_genes=["APP", "PSEN1"],
    ))))
    results = gprofiler.enrich_gene_list(["APP", "PSEN1"])
    assert [r["term_id"] for r in results] == ["R-HSA-1", "GO:1"]
    assert results[0]["name"] == "strong term"
    assert results[0]["source"] == "REAC"
    assert results[0]["p_value"] == 1e-10


def test_hit_genes_come_from_the_intersections_index(monkeypatch):
    monkeypatch.setattr(gprofiler, "SESSION", FakeSession(post=FakeResponse(payload(
        ("R-HSA-1", "term", "REAC", 1e-9, [["IEA"], [], ["IDA"]]),
        query_genes=["APP", "PSEN1", "MAPT"],
    ))))
    assert gprofiler.enrich_gene_list(["APP", "PSEN1", "MAPT"])[0]["genes"] == ["APP", "MAPT"]


def test_query_falls_back_to_the_input_order(monkeypatch):
    """When meta carries no query list, the submitted order is used."""
    monkeypatch.setattr(gprofiler, "SESSION", FakeSession(post=FakeResponse(
        {"result": [{"native": "R-1", "name": "t", "source": "REAC", "p_value": 0.01,
                     "intersections": [["x"], []]}], "meta": {}})))
    assert gprofiler.enrich_gene_list(["APP", "PSEN1"])[0]["genes"] == ["APP"]


def test_commercial_sources_are_not_requested(monkeypatch):
    session = FakeSession(post=FakeResponse(payload()))
    monkeypatch.setattr(gprofiler, "SESSION", session)
    gprofiler.enrich_gene_list(["APP"])
    sources = session.post_calls[0]["json"]["sources"]
    assert sources == ["REAC", "GO:BP", "WP"]
    assert "KEGG" not in sources


def test_request_shape(monkeypatch):
    session = FakeSession(post=FakeResponse(payload()))
    monkeypatch.setattr(gprofiler, "SESSION", session)
    gprofiler.enrich_gene_list(["APP", "PSEN1"], significance_threshold=0.01)
    body = session.post_calls[0]["json"]
    assert body["organism"] == "hsapiens"
    assert body["query"] == ["APP", "PSEN1"]
    assert body["user_threshold"] == 0.01
    # no_evidences must stay absent, or intersections are not returned.
    assert "no_evidences" not in body


def test_max_results_caps_the_output(monkeypatch):
    terms = tuple((f"GO:{i}", f"t{i}", "GO:BP", i / 1000, [["x"]]) for i in range(1, 30))
    monkeypatch.setattr(gprofiler, "SESSION",
                        FakeSession(post=FakeResponse(payload(*terms, query_genes=["APP"]))))
    assert len(gprofiler.enrich_gene_list(["APP"], max_results=5)) == 5


def test_empty_gene_list_makes_no_request(monkeypatch):
    session = FakeSession(post=FakeResponse(payload()))
    monkeypatch.setattr(gprofiler, "SESSION", session)
    assert gprofiler.enrich_gene_list([]) == []
    assert gprofiler.enrich_gene_list(["", "  "]) == []
    assert session.post_calls == []


def test_outage_returns_empty(monkeypatch):
    monkeypatch.setattr(gprofiler, "SESSION", FakeSession(post=ConnectionError("down")))
    assert gprofiler.enrich_gene_list(["APP"]) == []


def test_http_error_returns_empty(monkeypatch):
    monkeypatch.setattr(gprofiler, "SESSION", FakeSession(post=FakeResponse(status=500)))
    assert gprofiler.enrich_gene_list(["APP"]) == []


def test_malformed_payload_returns_empty(monkeypatch):
    monkeypatch.setattr(gprofiler, "SESSION", FakeSession(post=FakeResponse("not-a-dict")))
    assert gprofiler.enrich_gene_list(["APP"]) == []


def test_non_numeric_p_value_defaults_to_one(monkeypatch):
    monkeypatch.setattr(gprofiler, "SESSION", FakeSession(post=FakeResponse(
        {"result": [{"native": "R-1", "name": "t", "source": "REAC",
                     "p_value": "bad", "intersections": []}], "meta": {}})))
    assert gprofiler.enrich_gene_list(["APP"])[0]["p_value"] == 1.0


def test_results_are_cached_on_disk(monkeypatch, tmp_path):
    monkeypatch.delenv("PPI_CACHE_DISABLED", raising=False)
    monkeypatch.setattr(gprofiler, "cache_enabled", lambda: True)
    monkeypatch.setattr(gprofiler, "cache_path",
                        lambda *parts: tmp_path.joinpath(*parts))
    (tmp_path / "gprofiler").mkdir(parents=True, exist_ok=True)

    session = FakeSession(post=FakeResponse(payload(
        ("R-1", "t", "REAC", 0.01, [["x"]]), query_genes=["APP"])))
    monkeypatch.setattr(gprofiler, "SESSION", session)

    first = gprofiler.enrich_gene_list(["APP"])
    second = gprofiler.enrich_gene_list(["APP"])
    assert first == second
    assert len(session.post_calls) == 1  # served from cache the second time


def test_cache_key_ignores_gene_order(monkeypatch, tmp_path):
    monkeypatch.delenv("PPI_CACHE_DISABLED", raising=False)
    monkeypatch.setattr(gprofiler, "cache_enabled", lambda: True)
    monkeypatch.setattr(gprofiler, "cache_path", lambda *parts: tmp_path.joinpath(*parts))
    (tmp_path / "gprofiler").mkdir(parents=True, exist_ok=True)

    session = FakeSession(post=FakeResponse(payload(
        ("R-1", "t", "REAC", 0.01, [["x"]]), query_genes=["APP"])))
    monkeypatch.setattr(gprofiler, "SESSION", session)
    gprofiler.enrich_gene_list(["APP", "PSEN1"])
    gprofiler.enrich_gene_list(["PSEN1", "APP"])
    assert len(session.post_calls) == 1
