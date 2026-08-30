import pytest

import app as app_module
from app import create_app, parse_genes, parse_ppi_options
from collectors.opentargets import OpenTargetsError

DISEASE_GENES = [
    {"symbol": "APP", "score": 0.9},
    {"symbol": "PSEN1", "score": 0.6},
    {"symbol": "MAPT", "score": 0.5},
]
PPI = {"APP": ["PSEN1", "MAPT"], "PSEN1": ["APP"], "LONELY": []}


def fake_collect(gene, **kwargs):
    return {
        "partners": PPI.get(gene.upper(), []),
        "excluded_hubs": [],
        "source_counts": {"SIGNOR": len(PPI.get(gene.upper(), []))},
        "graph_size": len(PPI.get(gene.upper(), [])) + 1,
    }


def fake_resolve_symbols(genes):
    return [
        {"input": g, "symbol": g.upper(), "status": "approved", "hgnc_id": f"HGNC:{g}", "name": ""}
        for g in genes
    ]


@pytest.fixture
def patched(monkeypatch):
    monkeypatch.setattr(app_module, "resolve_disease_id",
                        lambda name: ("EFO_0000249", "Alzheimer disease"))
    monkeypatch.setattr(app_module, "get_disease_with_top_genes",
                        lambda disease_id, top_n=100: ("Alzheimer disease", list(DISEASE_GENES)))
    monkeypatch.setattr(app_module, "collect_ppi_partners", fake_collect)
    monkeypatch.setattr(app_module, "resolve_gene_symbols", fake_resolve_symbols)
    return monkeypatch


@pytest.fixture
def client(patched):
    return create_app({"TESTING": True}).test_client()


# --- parse_genes -----------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("APP\nPSEN1\nAPOE", ["APP", "PSEN1", "APOE"]),
    ("APP, PSEN1, APOE", ["APP", "PSEN1", "APOE"]),
    ("APP; PSEN1\tAPOE", ["APP", "PSEN1", "APOE"]),
    (["APP", " PSEN1 ", ""], ["APP", "PSEN1"]),
    ("APP\napp\nAPP", ["APP"]),
    ("", []), (None, []), (42, []),
])
def test_parse_genes(raw, expected):
    assert parse_genes(raw) == expected


# --- parse_ppi_options -----------------------------------------------------

DEFAULTS = {"sources": ["signor", "string"], "string_score": 700,
            "hub_threshold": 1000, "max_nodes": 100, "exclude_hubs": True}


def test_ppi_defaults_apply_when_absent():
    assert parse_ppi_options(None, DEFAULTS) == {
        "sources": ["signor", "string"], "string_score": 700, "min_score": None,
        "hub_threshold": 1000, "max_nodes": 100, "exclude_hubs": True,
    }


def test_ppi_sources_are_selectable():
    opts = parse_ppi_options({"sources": ["biogrid", "signor"]}, DEFAULTS)
    assert opts["sources"] == ["signor", "biogrid"]  # normalised to a fixed order


def test_unknown_ppi_sources_are_dropped():
    assert parse_ppi_options({"sources": ["signor", "bogus"]}, DEFAULTS)["sources"] == ["signor"]


def test_empty_source_selection_falls_back_to_defaults():
    assert parse_ppi_options({"sources": []}, DEFAULTS)["sources"] == ["signor", "string"]


@pytest.mark.parametrize("value,expected", [
    (900, 900), (0, 0), (1000, 1000), (5000, 1000), (-10, 0), ("700", 700), ("bad", 700),
])
def test_string_score_is_clamped(value, expected):
    assert parse_ppi_options({"string_score": value}, DEFAULTS)["string_score"] == expected


@pytest.mark.parametrize("value,expected", [
    (0.4, 0.4), ("0.7", 0.7), (2.0, 1.0), (-1, 0.0),
    (None, None), ("", None), ("null", None), ("bad", None),
])
def test_min_score_is_parsed_and_clamped(value, expected):
    assert parse_ppi_options({"min_score": value}, DEFAULTS)["min_score"] == expected


def test_hub_threshold_and_max_nodes_are_clamped():
    opts = parse_ppi_options({"hub_threshold": 50, "max_nodes": 9999}, DEFAULTS)
    assert opts["hub_threshold"] == 50
    assert opts["max_nodes"] == 500


def test_exclude_hubs_can_be_turned_off():
    assert parse_ppi_options({"exclude_hubs": False}, DEFAULTS)["exclude_hubs"] is False


# --- disease search --------------------------------------------------------

def test_disease_search_returns_candidates(monkeypatch, client):
    monkeypatch.setattr(app_module, "search_diseases", lambda q, limit=10: [
        {"id": "EFO_0000249", "name": "Alzheimer disease", "description": "d", "exact": True},
        {"id": "EFO_1", "name": "Alzheimer disease 2", "description": "", "exact": False},
    ])
    body = client.get("/api/diseases?q=alzheimer").get_json()
    assert body["query"] == "alzheimer"
    assert [r["id"] for r in body["results"]] == ["EFO_0000249", "EFO_1"]


def test_disease_search_requires_a_query(client):
    assert client.get("/api/diseases?q=  ").status_code == 400


def test_disease_search_accepts_an_ontology_id(monkeypatch, client):
    monkeypatch.setattr(app_module, "get_disease_label", lambda d: "Alzheimer disease")
    body = client.get("/api/diseases?q=efo_0000249").get_json()
    assert body["results"] == [
        {"id": "EFO_0000249", "name": "Alzheimer disease", "description": "", "exact": True}
    ]


def test_disease_search_on_an_unknown_id_returns_no_results(monkeypatch, client):
    monkeypatch.setattr(app_module, "get_disease_label", lambda d: None)
    assert client.get("/api/diseases?q=EFO_9999999").get_json()["results"] == []


def test_disease_search_reports_upstream_failure(monkeypatch, client):
    def boom(q, limit=10):
        raise OpenTargetsError("Open Targets request failed")
    monkeypatch.setattr(app_module, "search_diseases", boom)
    assert client.get("/api/diseases?q=x").status_code == 502


# --- gene validation -------------------------------------------------------

def test_gene_validation_endpoint(monkeypatch, client):
    monkeypatch.setattr(app_module, "resolve_gene_symbols", lambda genes: [
        {"input": "PS1", "symbol": "PSEN1", "status": "alias", "hgnc_id": "HGNC:8828", "name": ""},
        {"input": "FOO", "symbol": "FOO", "status": "unknown", "hgnc_id": "", "name": ""},
    ])
    body = client.post("/api/genes/validate", json={"genes": "PS1, FOO"}).get_json()
    assert body["summary"] == {"alias": 1, "unknown": 1}
    assert body["genes"][0]["symbol"] == "PSEN1"


def test_gene_validation_requires_genes(client):
    assert client.post("/api/genes/validate", json={}).status_code == 400


def test_gene_validation_respects_the_limit(patched):
    c = create_app({"TESTING": True, "MAX_GENES": 2}).test_client()
    res = c.post("/api/genes/validate", json={"genes": ["A", "B", "C"]})
    assert res.status_code == 400


# --- analyze: validation ---------------------------------------------------

@pytest.mark.parametrize("payload", [
    {}, {"disease": "AD"}, {"genes": ["APP"]}, {"disease": "  ", "genes": ["APP"]},
    {"disease": "AD", "genes": []},
])
def test_missing_input_is_rejected(client, payload):
    res = client.post("/api/analyze", json=payload)
    assert res.status_code == 400
    assert "required" in res.get_json()["error"]


def test_gene_limit_is_enforced(client):
    res = client.post("/api/analyze",
                      json={"disease": "AD", "genes": [f"G{i}" for i in range(31)]})
    assert res.status_code == 400
    assert "Too many genes" in res.get_json()["error"]


def test_unknown_disease_returns_404(monkeypatch, client):
    monkeypatch.setattr(app_module, "resolve_disease_id", lambda name: (None, name))
    res = client.post("/api/analyze", json={"disease": "nope", "genes": ["APP"]})
    assert res.status_code == 404


def test_opentargets_outage_returns_502(monkeypatch, client):
    def boom(name):
        raise OpenTargetsError("timeout")
    monkeypatch.setattr(app_module, "resolve_disease_id", boom)
    assert client.post("/api/analyze", json={"disease": "AD", "genes": ["APP"]}).status_code == 502


# --- analyze: disease selection -------------------------------------------

def test_disease_id_skips_the_name_search(monkeypatch, client):
    def should_not_run(name):
        raise AssertionError("resolve_disease_id must not be called when disease_id is given")
    monkeypatch.setattr(app_module, "resolve_disease_id", should_not_run)
    body = client.post("/api/analyze",
                       json={"disease_id": "EFO_0000249", "genes": ["APP"]}).get_json()
    assert body["disease_id"] == "EFO_0000249"
    assert body["disease_label"] == "Alzheimer disease"


def test_disease_id_casing_is_normalised(monkeypatch, client):
    seen = {}
    def capture(disease_id, top_n=100):
        seen["id"] = disease_id
        return "Alzheimer disease", list(DISEASE_GENES)
    monkeypatch.setattr(app_module, "get_disease_with_top_genes", capture)
    client.post("/api/analyze", json={"disease_id": "efo_0000249", "genes": ["APP"]})
    assert seen["id"] == "EFO_0000249"


# --- analyze: gene symbols -------------------------------------------------

def test_aliases_are_analysed_under_the_approved_symbol(monkeypatch, client):
    monkeypatch.setattr(app_module, "resolve_gene_symbols", lambda genes: [
        {"input": "PS1", "symbol": "PSEN1", "status": "alias", "hgnc_id": "HGNC:8828", "name": "presenilin 1"},
    ])
    body = client.post("/api/analyze", json={"disease": "AD", "genes": ["PS1"]}).get_json()
    result = body["results"][0]
    assert result["gene"] == "PSEN1"
    assert result["input_gene"] == "PS1"
    assert result["symbol_status"] == "alias"
    assert result["hgnc_id"] == "HGNC:8828"
    assert body["genes"]["summary"] == {"alias": 1}


def test_unknown_symbols_are_flagged_but_still_analysed(monkeypatch, client):
    monkeypatch.setattr(app_module, "resolve_gene_symbols", lambda genes: [
        {"input": "NOTAGENE", "symbol": "NOTAGENE", "status": "unknown", "hgnc_id": "", "name": ""},
    ])
    body = client.post("/api/analyze", json={"disease": "AD", "genes": ["NOTAGENE"]}).get_json()
    assert body["results"][0]["symbol_status"] == "unknown"
    assert body["results"][0]["weighted_score"] == 0.0


def test_symbol_validation_can_be_disabled(patched):
    c = create_app({"TESTING": True, "VALIDATE_SYMBOLS": False}).test_client()
    body = c.post("/api/analyze", json={"disease": "AD", "genes": ["app"]}).get_json()
    assert body["results"][0]["gene"] == "app"
    assert body["results"][0]["symbol_status"] == "unchecked"


# --- analyze: PPI options --------------------------------------------------

def test_ppi_options_reach_the_collector(monkeypatch, client):
    seen = {}
    def capture(gene, **kwargs):
        seen.update(kwargs)
        return fake_collect(gene, **kwargs)
    monkeypatch.setattr(app_module, "collect_ppi_partners", capture)
    client.post("/api/analyze", json={
        "disease": "AD", "genes": ["APP"],
        "ppi": {"sources": ["string"], "string_score": 900, "min_score": 0.5,
                "hub_threshold": 500, "max_nodes": 20, "exclude_hubs": False},
    })
    assert seen["sources"] == ["string"]
    assert seen["string_required_score"] == 900
    assert seen["min_score"] == 0.5
    assert seen["hub_threshold"] == 500
    assert seen["max_nodes"] == 20
    assert seen["exclude_hubs"] is False


def test_biogrid_is_dropped_without_a_key(client):
    body = client.post("/api/analyze", json={
        "disease": "AD", "genes": ["APP"], "ppi": {"sources": ["signor", "biogrid"]},
    }).get_json()
    assert body["ppi"]["sources"] == ["signor"]
    assert body["ppi"]["biogrid_skipped"] is True


def test_biogrid_is_kept_when_a_key_is_configured(patched):
    c = create_app({"TESTING": True, "BIOGRID_KEY": "k"}).test_client()
    body = c.post("/api/analyze", json={
        "disease": "AD", "genes": ["APP"], "ppi": {"sources": ["biogrid"]},
    }).get_json()
    assert body["ppi"]["sources"] == ["biogrid"]
    assert "biogrid_skipped" not in body["ppi"]


def test_the_options_used_are_echoed_back(client):
    body = client.post("/api/analyze", json={"disease": "AD", "genes": ["APP"]}).get_json()
    assert body["ppi"]["sources"] == ["signor", "string"]
    assert body["ppi"]["string_score"] == 700


# --- analyze: results ------------------------------------------------------

def test_analyze_returns_percentages_against_the_disease_network(client):
    body = client.post("/api/analyze", json={"disease": "AD", "genes": ["APP"]}).get_json()
    r = body["results"][0]
    assert body["disease_gene_count"] == 3
    # APP overlaps PSEN1 + MAPT and is itself a disease gene: all 3 of 3.
    assert r["overlap_percent"] == 100.0
    assert r["weighted_percent"] == 100.0
    assert r["matched_count"] == 3


def test_partial_overlap_percentages(monkeypatch, client):
    monkeypatch.setattr(app_module, "collect_ppi_partners",
                        lambda gene, **kw: {"partners": ["PSEN1"], "excluded_hubs": [],
                                            "source_counts": {}, "graph_size": 2})
    body = client.post("/api/analyze", json={"disease": "AD", "genes": ["XYZ"]}).get_json()
    r = body["results"][0]
    assert r["matched_count"] == 1
    assert r["overlap_percent"] == 33.3          # 1 of 3 disease genes
    assert r["weighted_percent"] == 30.0         # 0.6 of a 2.0 total score


def test_results_are_sorted_by_descending_score(client):
    body = client.post("/api/analyze",
                       json={"disease": "AD", "genes": ["LONELY", "PSEN1", "APP"]}).get_json()
    scores = [r["weighted_score"] for r in body["results"]]
    assert scores == sorted(scores, reverse=True)


def test_hub_exclusions_are_reported(monkeypatch, client):
    monkeypatch.setattr(app_module, "collect_ppi_partners",
                        lambda gene, **kw: {"partners": ["PSEN1"], "excluded_hubs": ["UBC", "ACTB"],
                                            "source_counts": {"SIGNOR": 3}, "graph_size": 4})
    r = client.post("/api/analyze", json={"disease": "AD", "genes": ["APP"]}).get_json()["results"][0]
    assert r["excluded_hub_count"] == 2
    assert r["source_counts"] == {"SIGNOR": 3}


def test_one_failing_gene_does_not_sink_the_request(monkeypatch, client):
    def flaky(gene, **kwargs):
        if gene == "BAD":
            raise RuntimeError("ppi exploded")
        return fake_collect(gene, **kwargs)
    monkeypatch.setattr(app_module, "collect_ppi_partners", flaky)
    body = client.post("/api/analyze", json={"disease": "AD", "genes": ["APP", "BAD"]}).get_json()
    by_gene = {r["gene"]: r for r in body["results"]}
    assert by_gene["BAD"]["error"] == "ppi exploded"
    assert by_gene["APP"]["weighted_score"] > 0


# --- pages -----------------------------------------------------------------

def test_index_page_renders(client):
    body = client.get("/").get_data(as_text=True)
    assert "NW重複評価" in body
    assert 'id="disease"' in body and 'id="genes"' in body


def test_healthz(client):
    body = client.get("/healthz").get_json()
    assert body["status"] == "ok"
    assert body["ppi_defaults"]["string_score"] == 700
