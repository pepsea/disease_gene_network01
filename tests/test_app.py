import pytest

import app as app_module
from app import create_app, parse_genes
from collectors.opentargets import OpenTargetsError

DISEASE_GENES = [
    {"symbol": "APP", "score": 0.9},
    {"symbol": "PSEN1", "score": 0.6},
    {"symbol": "MAPT", "score": 0.5},
]
PPI = {"APP": ["PSEN1", "MAPT"], "PSEN1": ["APP"], "LONELY": []}


@pytest.fixture
def client(monkeypatch):
    monkeypatch.setattr(app_module, "resolve_disease_id",
                        lambda name: ("EFO_0000249", "Alzheimer disease"))
    monkeypatch.setattr(app_module, "get_disease_top_genes",
                        lambda disease_id, top_n=100: list(DISEASE_GENES))
    monkeypatch.setattr(app_module, "get_ppi_partners",
                        lambda gene, biogrid_key="", top_n=30: PPI.get(gene.upper(), []))
    flask_app = create_app({"TESTING": True})
    return flask_app.test_client()


# --- parse_genes -----------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("APP\nPSEN1\nAPOE", ["APP", "PSEN1", "APOE"]),
    ("APP, PSEN1, APOE", ["APP", "PSEN1", "APOE"]),
    ("APP; PSEN1\tAPOE", ["APP", "PSEN1", "APOE"]),
    ("APP\n\n  \nPSEN1", ["APP", "PSEN1"]),
    (["APP", " PSEN1 ", ""], ["APP", "PSEN1"]),
    ("APP\napp\nAPP", ["APP"]),           # case-insensitive de-duplication
    (["APP,PSEN1"], ["APP", "PSEN1"]),    # separators inside a list item
    ("", []), (None, []), (42, []), ({"a": 1}, []),
])
def test_parse_genes(raw, expected):
    assert parse_genes(raw) == expected


def test_parse_genes_keeps_the_first_spelling():
    assert parse_genes("app\nAPP") == ["app"]


# --- validation ------------------------------------------------------------

@pytest.mark.parametrize("payload", [
    {}, {"disease": "Alzheimer disease"}, {"genes": ["APP"]},
    {"disease": "  ", "genes": ["APP"]}, {"disease": "AD", "genes": []},
    {"disease": "AD", "genes": ["  "]},
])
def test_missing_input_is_rejected(client, payload):
    res = client.post("/api/analyze", json=payload)
    assert res.status_code == 400
    assert "required" in res.get_json()["error"]


def test_non_json_body_is_rejected(client):
    res = client.post("/api/analyze", data="not json", content_type="text/plain")
    assert res.status_code == 400


def test_gene_limit_is_enforced(client):
    genes = [f"GENE{i}" for i in range(31)]
    res = client.post("/api/analyze", json={"disease": "AD", "genes": genes})
    assert res.status_code == 400
    assert "Too many genes" in res.get_json()["error"]


def test_gene_limit_is_configurable(monkeypatch):
    monkeypatch.setattr(app_module, "resolve_disease_id", lambda n: ("EFO_1", "AD"))
    monkeypatch.setattr(app_module, "get_disease_top_genes",
                        lambda d, top_n=100: list(DISEASE_GENES))
    monkeypatch.setattr(app_module, "get_ppi_partners",
                        lambda gene, biogrid_key="", top_n=30: [])
    c = create_app({"TESTING": True, "MAX_GENES": 2}).test_client()
    assert c.post("/api/analyze", json={"disease": "AD", "genes": ["A", "B"]}).status_code == 200
    assert c.post("/api/analyze", json={"disease": "AD", "genes": ["A", "B", "C"]}).status_code == 400


def test_duplicates_do_not_count_towards_the_limit(monkeypatch):
    monkeypatch.setattr(app_module, "resolve_disease_id", lambda n: ("EFO_1", "AD"))
    monkeypatch.setattr(app_module, "get_disease_top_genes",
                        lambda d, top_n=100: list(DISEASE_GENES))
    monkeypatch.setattr(app_module, "get_ppi_partners",
                        lambda gene, biogrid_key="", top_n=30: [])
    c = create_app({"TESTING": True, "MAX_GENES": 2}).test_client()
    res = c.post("/api/analyze", json={"disease": "AD", "genes": ["APP", "app", "APP"]})
    assert res.status_code == 200
    assert len(res.get_json()["results"]) == 1


# --- upstream failures -----------------------------------------------------

def test_unknown_disease_returns_404(monkeypatch, client):
    monkeypatch.setattr(app_module, "resolve_disease_id", lambda name: (None, name))
    res = client.post("/api/analyze", json={"disease": "nope", "genes": ["APP"]})
    assert res.status_code == 404
    assert "Disease not found: nope" in res.get_json()["error"]


def test_disease_without_genes_returns_404(monkeypatch, client):
    monkeypatch.setattr(app_module, "get_disease_top_genes", lambda d, top_n=100: [])
    res = client.post("/api/analyze", json={"disease": "AD", "genes": ["APP"]})
    assert res.status_code == 404
    assert "No associated genes" in res.get_json()["error"]


def test_opentargets_outage_returns_502(monkeypatch, client):
    def boom(name):
        raise OpenTargetsError("Open Targets request failed: timeout")
    monkeypatch.setattr(app_module, "resolve_disease_id", boom)
    res = client.post("/api/analyze", json={"disease": "AD", "genes": ["APP"]})
    assert res.status_code == 502
    assert "timeout" in res.get_json()["error"]


def test_one_failing_gene_does_not_sink_the_request(monkeypatch, client):
    def flaky(gene, biogrid_key="", top_n=30):
        if gene == "BAD":
            raise RuntimeError("ppi exploded")
        return PPI.get(gene.upper(), [])
    monkeypatch.setattr(app_module, "get_ppi_partners", flaky)
    res = client.post("/api/analyze", json={"disease": "AD", "genes": ["APP", "BAD"]})
    assert res.status_code == 200
    by_gene = {r["gene"]: r for r in res.get_json()["results"]}
    assert by_gene["BAD"]["error"] == "ppi exploded"
    assert by_gene["APP"]["weighted_score"] > 0


# --- successful analysis ---------------------------------------------------

def test_analyze_returns_scored_results(client):
    res = client.post("/api/analyze", json={"disease": "Alzheimer disease",
                                            "genes": ["APP", "LONELY"]})
    assert res.status_code == 200
    body = res.get_json()
    assert body["disease_id"] == "EFO_0000249"
    assert body["disease_label"] == "Alzheimer disease"
    assert body["disease_gene_count"] == 3
    assert body["ppi_top_n"] == 30

    app_result = next(r for r in body["results"] if r["gene"] == "APP")
    assert app_result["weighted_score"] == 1.0
    assert app_result["target_self"] == "APP"
    assert app_result["overlap_count"] == 2
    assert app_result["interpretation"] == "strong"

    lonely = next(r for r in body["results"] if r["gene"] == "LONELY")
    assert lonely["weighted_score"] == 0.0


def test_results_are_sorted_by_descending_score(client):
    res = client.post("/api/analyze", json={"disease": "AD",
                                            "genes": ["LONELY", "PSEN1", "APP"]})
    scores = [r["weighted_score"] for r in res.get_json()["results"]]
    assert scores == sorted(scores, reverse=True)


def test_genes_accepts_a_raw_text_blob(client):
    res = client.post("/api/analyze", json={"disease": "AD", "genes": "APP, PSEN1"})
    assert res.status_code == 200
    assert {r["gene"] for r in res.get_json()["results"]} == {"APP", "PSEN1"}


def test_top_n_override_is_forwarded(monkeypatch, client):
    seen = {}
    def capture(disease_id, top_n=100):
        seen["top_n"] = top_n
        return list(DISEASE_GENES)
    monkeypatch.setattr(app_module, "get_disease_top_genes", capture)
    client.post("/api/analyze", json={"disease": "AD", "genes": ["APP"], "top_n": 50})
    assert seen["top_n"] == 50
    client.post("/api/analyze", json={"disease": "AD", "genes": ["APP"], "top_n": 9999})
    assert seen["top_n"] == 500
    client.post("/api/analyze", json={"disease": "AD", "genes": ["APP"]})
    assert seen["top_n"] == 100


def test_biogrid_key_is_passed_to_the_collector(monkeypatch):
    monkeypatch.setattr(app_module, "resolve_disease_id", lambda n: ("EFO_1", "AD"))
    monkeypatch.setattr(app_module, "get_disease_top_genes",
                        lambda d, top_n=100: list(DISEASE_GENES))
    seen = {}
    def capture(gene, biogrid_key="", top_n=30):
        seen["key"] = biogrid_key
        seen["top_n"] = top_n
        return []
    monkeypatch.setattr(app_module, "get_ppi_partners", capture)
    c = create_app({"TESTING": True, "BIOGRID_KEY": "secret", "PPI_TOP_N": 15}).test_client()
    c.post("/api/analyze", json={"disease": "AD", "genes": ["APP"]})
    assert seen == {"key": "secret", "top_n": 15}


# --- pages -----------------------------------------------------------------

def test_index_page_renders(client):
    res = client.get("/")
    assert res.status_code == 200
    body = res.get_data(as_text=True)
    assert "NW重複評価" in body
    assert 'id="form"' in body and 'id="disease"' in body and 'id="genes"' in body


def test_healthz(client):
    body = client.get("/healthz").get_json()
    assert body["status"] == "ok"
    assert body["biogrid_enabled"] is False


def test_config_comes_from_the_environment(monkeypatch):
    monkeypatch.setenv("NW_MAX_GENES", "7")
    monkeypatch.setenv("NW_PPI_TOP_N", "12")
    monkeypatch.setenv("BIOGRID_KEY", "abc")
    cfg = create_app().config
    assert cfg["MAX_GENES"] == 7
    assert cfg["PPI_TOP_N"] == 12
    assert cfg["BIOGRID_KEY"] == "abc"


@pytest.mark.parametrize("bad", ["0", "-5", "abc", ""])
def test_invalid_env_values_fall_back_to_defaults(monkeypatch, bad):
    monkeypatch.setenv("NW_MAX_GENES", bad)
    assert create_app().config["MAX_GENES"] == app_module.DEFAULT_MAX_GENES
