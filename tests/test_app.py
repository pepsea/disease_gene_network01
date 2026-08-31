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


DISEASE_PATHWAYS = [
    {"term_id": "R-HSA-1", "name": "Amyloid fiber formation", "source": "REAC",
     "p_value": 1e-10, "genes": ["APP", "PSEN1"]},
    {"term_id": "GO:1", "name": "neuron death", "source": "GO:BP",
     "p_value": 1e-5, "genes": ["MAPT"]},
    {"term_id": "WP:1", "name": "unrelated", "source": "WP",
     "p_value": 1e-2, "genes": ["ZZZ"]},
]
# Enrichment of each gene's neighbourhood, keyed by the first query gene.
GENE_PATHWAYS = {
    "APP": [{"term_id": "R-HSA-1", "p_value": 1e-9}, {"term_id": "GO:1", "p_value": 1e-4}],
    "PSEN1": [{"term_id": "R-HSA-1", "p_value": 1e-8}],
    "LONELY": [],
}


def make_fake_enrich():
    """Stub for enrich_gene_list.

    The disease signature is enriched once, before any gene work, so the first
    call is the disease list; later calls are per-gene neighbourhoods. Matching
    on the gene list alone would be ambiguous, since a gene's neighbourhood can
    be exactly the disease gene set.
    """
    state = {"disease_done": False}

    def fake_enrich(genes, max_results=50, **kwargs):
        if not state["disease_done"]:
            state["disease_done"] = True
            return list(DISEASE_PATHWAYS)
        return list(GENE_PATHWAYS.get((genes[0] if genes else "").upper(), []))

    return fake_enrich


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
    monkeypatch.setattr(app_module, "enrich_gene_list", make_fake_enrich())
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


# --- enrichment overlap (second table) -------------------------------------

def test_enrichment_overlap_is_returned_per_gene(client):
    body = client.post("/api/analyze", json={"disease": "AD", "genes": ["APP"]}).get_json()
    assert body["enrichment"]["enabled"] is True
    assert body["enrichment"]["disease_pathway_count"] == 3
    assert [p["name"] for p in body["enrichment"]["top_pathways"]][0] == "Amyloid fiber formation"

    r = body["results"][0]
    # APP's neighbourhood reproduces 2 of the 3 disease pathways...
    assert r["pathway_overlap_count"] == 2
    assert r["disease_pathway_count"] == 3
    assert r["pathway_overlap_percent"] == 66.7
    # ...and those two carry almost all the significance weight.
    assert r["pathway_weighted_percent"] == 88.2
    assert r["pathway_interpretation"] == "strong"
    assert [p["name"] for p in r["overlapping_pathways"]] == [
        "Amyloid fiber formation", "neuron death"]


def test_target_fit_reports_pathways_containing_the_gene(client):
    r = client.post("/api/analyze", json={"disease": "AD", "genes": ["APP"]}).get_json()["results"][0]
    # APP is annotated to 1 of the 3 disease pathways.
    assert r["target_in_pathway_count"] == 1
    assert r["target_fit_percent"] == 33.3
    assert r["target_in_pathways"][0]["name"] == "Amyloid fiber formation"


def test_the_target_itself_counts_in_the_pathway_overlap(client):
    """MAPT's neighbourhood enriches to nothing, but MAPT is in a disease pathway."""
    body = client.post("/api/analyze", json={"disease": "AD", "genes": ["MAPT"]}).get_json()
    r = body["results"][0]
    assert r["gene_pathway_count"] == 0        # nothing came from the network
    assert r["pathway_overlap_count"] == 0
    assert r["target_in_pathway_count"] == 1   # MAPT is listed in "neuron death"
    assert r["pathway_matched_count"] == 1     # so it is counted as covered
    assert r["pathway_overlap_percent"] == 33.3
    assert r["overlapping_pathways"][0]["via"] == "target"


def test_pathway_matched_count_separates_network_and_target_hits(client):
    r = client.post("/api/analyze", json={"disease": "AD", "genes": ["APP"]}).get_json()["results"][0]
    # APP's network reaches two pathways; APP itself is listed in one of them.
    assert r["pathway_overlap_count"] == 2
    assert r["target_in_pathway_count"] == 1
    assert r["pathway_matched_count"] == 2     # union, not 3
    via = {p["term_id"]: p["via"] for p in r["overlapping_pathways"]}
    assert via == {"R-HSA-1": "both", "GO:1": "network"}


def test_gene_with_no_enriched_pathways_scores_zero(client):
    """LONELY has no network enrichment and is in no disease pathway."""
    body = client.post("/api/analyze", json={"disease": "AD", "genes": ["LONELY"]}).get_json()
    r = body["results"][0]
    assert r["pathway_overlap_percent"] == 0.0
    assert r["pathway_weighted_percent"] == 0.0
    assert r["target_in_pathway_count"] == 0
    assert r["pathway_matched_count"] == 0


def test_enrichment_can_be_disabled_per_request(client):
    body = client.post("/api/analyze",
                       json={"disease": "AD", "genes": ["APP"], "enrichment": False}).get_json()
    assert body["enrichment"]["enabled"] is False
    assert body["enrichment"]["disease_pathway_count"] == 0
    assert "pathway_overlap_percent" not in body["results"][0]
    # The gene-level overlap is unaffected.
    assert body["results"][0]["weighted_percent"] == 100.0


def test_enrichment_can_be_disabled_by_config(patched):
    c = create_app({"TESTING": True, "ENRICHMENT": False}).test_client()
    body = c.post("/api/analyze", json={"disease": "AD", "genes": ["APP"]}).get_json()
    assert body["enrichment"]["enabled"] is False
    assert "pathway_overlap_percent" not in body["results"][0]


def test_a_request_cannot_re_enable_disabled_enrichment(patched):
    c = create_app({"TESTING": True, "ENRICHMENT": False}).test_client()
    body = c.post("/api/analyze",
                  json={"disease": "AD", "genes": ["APP"], "enrichment": True}).get_json()
    assert body["enrichment"]["enabled"] is False


def test_enrichment_outage_leaves_the_gene_table_intact(monkeypatch, client):
    monkeypatch.setattr(app_module, "enrich_gene_list", lambda genes, **kw: [])
    body = client.post("/api/analyze", json={"disease": "AD", "genes": ["APP"]}).get_json()
    assert body["enrichment"]["disease_pathway_count"] == 0
    assert "pathway_overlap_percent" not in body["results"][0]
    assert body["results"][0]["weighted_percent"] == 100.0


def test_disease_signature_is_enriched_once_not_per_gene(monkeypatch, client):
    calls = []
    fake = make_fake_enrich()
    def counting(genes, **kwargs):
        calls.append(list(genes))
        return fake(genes, **kwargs)
    monkeypatch.setattr(app_module, "enrich_gene_list", counting)
    client.post("/api/analyze", json={"disease": "AD", "genes": ["APP", "PSEN1", "LONELY"]})
    # One call for the disease signature, then one per gene.
    assert len(calls) == 4
    assert calls[0] == ["APP", "PSEN1", "MAPT"]


def test_enrichment_uses_the_same_disease_genes_as_the_overlap(monkeypatch, client):
    seen = {}
    fake = make_fake_enrich()
    def capture(genes, **kwargs):
        seen.setdefault("first", list(genes))
        return fake(genes, **kwargs)
    monkeypatch.setattr(app_module, "enrich_gene_list", capture)
    body = client.post("/api/analyze", json={"disease": "AD", "genes": ["APP"]}).get_json()
    assert seen["first"] == [g["symbol"] for g in DISEASE_GENES]
    assert body["enrichment"]["disease_gene_n"] == body["disease_gene_count"]


# --- pages -----------------------------------------------------------------

def test_index_page_renders(client):
    body = client.get("/").get_data(as_text=True)
    assert "NW重複評価" in body
    assert 'id="disease"' in body and 'id="genes"' in body


def test_healthz(client):
    body = client.get("/healthz").get_json()
    assert body["status"] == "ok"
    assert body["ppi_defaults"]["string_score"] == 700
