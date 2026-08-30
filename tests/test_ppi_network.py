import pytest

import ppi_network as net


@pytest.fixture(autouse=True)
def no_external_calls(monkeypatch):
    """Hub lookups hit IntAct; default them to "unknown" unless a test says otherwise."""
    net.reset_hub_cache()
    monkeypatch.setattr(net, "global_interactor_count", lambda gene: -1)
    yield
    net.reset_hub_cache()


def signor_rows(*pairs, score=None, partner_type="gene"):
    return [{"source": a, "target": b, "partner": b, "partner_type": partner_type,
             "score": score, "db": "SIGNOR"} for a, b in pairs]


def string_rows(gene, *pairs):
    return [{"source": gene, "target": p, "partner": p, "partner_type": "gene",
             "score": s, "db": "STRING"} for p, s in pairs]


def patch_sources(monkeypatch, signor=(), string=(), biogrid=()):
    monkeypatch.setattr(net.signor, "get_interactions", lambda g: list(signor))
    monkeypatch.setattr(net.string_db, "get_interactions",
                        lambda g, required_score=400: list(string))
    monkeypatch.setattr(net.biogrid, "get_interactions",
                        lambda g, api_key=None: list(biogrid))


# --- source selection ------------------------------------------------------

def test_only_selected_sources_are_queried(monkeypatch):
    called = []
    monkeypatch.setattr(net.signor, "get_interactions",
                        lambda g: called.append("signor") or [])
    monkeypatch.setattr(net.string_db, "get_interactions",
                        lambda g, required_score=400: called.append("string") or [])
    monkeypatch.setattr(net.biogrid, "get_interactions",
                        lambda g, api_key=None: called.append("biogrid") or [])

    net.collect_ppi_partners("APP", sources=["string"])
    assert called == ["string"]

    called.clear()
    net.collect_ppi_partners("APP", sources=["signor", "biogrid"], biogrid_api_key="k")
    assert called == ["signor", "biogrid"]


def test_biogrid_is_skipped_without_a_key(monkeypatch):
    called = []
    monkeypatch.setattr(net.signor, "get_interactions", lambda g: [])
    monkeypatch.setattr(net.string_db, "get_interactions", lambda g, required_score=400: [])
    monkeypatch.setattr(net.biogrid, "get_interactions",
                        lambda g, api_key=None: called.append("biogrid") or [])
    net.collect_ppi_partners("APP", sources=["biogrid"], biogrid_api_key="")
    assert called == []


def test_string_threshold_is_passed_through(monkeypatch):
    seen = {}
    monkeypatch.setattr(net.signor, "get_interactions", lambda g: [])
    monkeypatch.setattr(net.string_db, "get_interactions",
                        lambda g, required_score=400: seen.update(score=required_score) or [])
    net.collect_ppi_partners("APP", sources=["string"], string_required_score=900)
    assert seen["score"] == 900


def test_sources_are_merged(monkeypatch):
    patch_sources(monkeypatch,
                  signor=signor_rows(("APP", "S1"), score=0.5),
                  string=string_rows("APP", ("T1", 0.9)))
    result = net.collect_ppi_partners("APP", sources=["signor", "string"])
    assert set(result["partners"]) == {"S1", "T1"}


def test_one_failing_source_does_not_break_the_others(monkeypatch):
    def boom(g):
        raise RuntimeError("signor down")
    monkeypatch.setattr(net.signor, "get_interactions", boom)
    monkeypatch.setattr(net.string_db, "get_interactions",
                        lambda g, required_score=400: string_rows("APP", ("T1", 0.9)))
    monkeypatch.setattr(net.biogrid, "get_interactions", lambda g, api_key=None: [])
    assert net.collect_ppi_partners("APP", sources=["signor", "string"])["partners"] == ["T1"]


# --- ranking ---------------------------------------------------------------

def test_multi_source_partners_rank_above_single_source(monkeypatch):
    patch_sources(monkeypatch,
                  signor=signor_rows(("APP", "BOTH"), score=0.1),
                  string=string_rows("APP", ("BOTH", 0.1), ("STRING_ONLY", 0.99)))
    partners = net.collect_ppi_partners("APP", sources=["signor", "string"])["partners"]
    assert partners[0] == "BOTH"


def test_representative_score_prefers_signor_over_string(monkeypatch):
    # Both partners are single-source, so the score decides.
    patch_sources(monkeypatch,
                  signor=signor_rows(("APP", "SIG_HIGH"), score=0.9),
                  string=string_rows("APP", ("STR_LOW", 0.2)))
    partners = net.collect_ppi_partners("APP", sources=["signor", "string"])["partners"]
    assert partners == ["SIG_HIGH", "STR_LOW"]


def test_higher_score_ranks_first(monkeypatch):
    patch_sources(monkeypatch, string=string_rows("APP", ("LOW", 0.2), ("HIGH", 0.95), ("MID", 0.5)))
    assert net.collect_ppi_partners("APP", sources=["string"])["partners"] == ["HIGH", "MID", "LOW"]


def test_max_nodes_caps_the_result(monkeypatch):
    patch_sources(monkeypatch, string=string_rows("APP", *[(f"P{i:02d}", i / 100) for i in range(50)]))
    assert len(net.collect_ppi_partners("APP", sources=["string"], max_nodes=5)["partners"]) == 5


# --- filtering -------------------------------------------------------------

def test_min_score_drops_weak_edges(monkeypatch):
    patch_sources(monkeypatch, string=string_rows("APP", ("WEAK", 0.2), ("STRONG", 0.9)))
    partners = net.collect_ppi_partners("APP", sources=["string"], min_score=0.5)["partners"]
    assert partners == ["STRONG"]


def test_min_score_none_keeps_everything(monkeypatch):
    patch_sources(monkeypatch, string=string_rows("APP", ("WEAK", 0.01), ("STRONG", 0.9)))
    assert len(net.collect_ppi_partners("APP", sources=["string"], min_score=None)["partners"]) == 2


def test_non_gene_partners_are_excluded(monkeypatch):
    patch_sources(monkeypatch, signor=(
        signor_rows(("APP", "REALGENE"), score=0.5)
        + signor_rows(("APP", "SOMEDRUG"), score=0.5, partner_type="chemical")))
    assert net.collect_ppi_partners("APP", sources=["signor"])["partners"] == ["REALGENE"]


def test_the_gene_itself_is_never_a_partner(monkeypatch):
    patch_sources(monkeypatch, signor=signor_rows(("APP", "APP"), ("APP", "PSEN1"), score=0.5))
    assert net.collect_ppi_partners("APP", sources=["signor"])["partners"] == ["PSEN1"]


# --- hub exclusion ---------------------------------------------------------

@pytest.mark.parametrize("symbol", ["UBC", "ACTB", "TUBB", "GNAS", "HSP90AA1", "YWHAZ", "SUMO1"])
def test_known_hub_families_are_recognised(symbol):
    assert net.is_hub_family(symbol)


@pytest.mark.parametrize("symbol", ["APP", "PSEN1", "TP53", "MAPT"])
def test_ordinary_genes_are_not_hub_families(symbol):
    assert not net.is_hub_family(symbol)


def test_hub_family_partners_are_excluded_and_reported(monkeypatch):
    patch_sources(monkeypatch, signor=signor_rows(("APP", "UBC"), ("APP", "PSEN1"), score=0.5))
    result = net.collect_ppi_partners("APP", sources=["signor"])
    assert result["partners"] == ["PSEN1"]
    assert result["excluded_hubs"] == ["UBC"]


def test_high_degree_partners_are_excluded(monkeypatch):
    patch_sources(monkeypatch, signor=signor_rows(("APP", "PROMISCUOUS"), ("APP", "SPECIFIC"), score=0.5))
    monkeypatch.setattr(net, "global_interactor_count",
                        lambda g: 5000 if g == "PROMISCUOUS" else 10)
    result = net.collect_ppi_partners("APP", sources=["signor"], hub_threshold=1000)
    assert result["partners"] == ["SPECIFIC"]


def test_hub_threshold_is_configurable(monkeypatch):
    patch_sources(monkeypatch, signor=signor_rows(("APP", "MIDDLING"), score=0.5))
    monkeypatch.setattr(net, "global_interactor_count", lambda g: 500)
    assert net.collect_ppi_partners("APP", sources=["signor"], hub_threshold=100)["partners"] == []
    assert net.collect_ppi_partners("APP", sources=["signor"], hub_threshold=1000)["partners"] == ["MIDDLING"]


def test_hub_exclusion_can_be_turned_off(monkeypatch):
    patch_sources(monkeypatch, signor=signor_rows(("APP", "UBC"), ("APP", "PSEN1"), score=0.5))
    result = net.collect_ppi_partners("APP", sources=["signor"], exclude_hubs=False)
    assert set(result["partners"]) == {"UBC", "PSEN1"}


def test_unknown_interactor_count_is_not_hub_evidence(monkeypatch):
    patch_sources(monkeypatch, signor=signor_rows(("APP", "PSEN1"), score=0.5))
    monkeypatch.setattr(net, "global_interactor_count", lambda g: -1)  # lookup failed
    assert net.collect_ppi_partners("APP", sources=["signor"])["partners"] == ["PSEN1"]


# --- score fallback --------------------------------------------------------

def test_unscored_edges_fall_back_to_inverse_interactor_count(monkeypatch):
    patch_sources(monkeypatch, biogrid=[
        {"source": "APP", "target": "RARE", "partner": "RARE", "partner_type": "gene",
         "score": None, "db": "BioGRID"},
        {"source": "APP", "target": "COMMON", "partner": "COMMON", "partner_type": "gene",
         "score": None, "db": "BioGRID"},
    ])
    counts = {"RARE": 10, "COMMON": 900}
    monkeypatch.setattr(net, "global_interactor_count", lambda g: counts.get(g, -1))
    partners = net.collect_ppi_partners(
        "APP", sources=["biogrid"], biogrid_api_key="k", hub_threshold=1000)["partners"]
    # 1/10 > 1/900, so the more specific partner ranks first.
    assert partners == ["RARE", "COMMON"]


# --- reporting -------------------------------------------------------------

def test_source_counts_are_reported(monkeypatch):
    patch_sources(monkeypatch,
                  signor=signor_rows(("APP", "A"), ("APP", "B"), score=0.5),
                  string=string_rows("APP", ("A", 0.9)))
    counts = net.collect_ppi_partners("APP", sources=["signor", "string"])["source_counts"]
    assert counts == {"SIGNOR": 2, "STRING": 1}


def test_blank_gene_returns_nothing(monkeypatch):
    patch_sources(monkeypatch)
    assert net.collect_ppi_partners("  ", sources=["signor"])["partners"] == []


def test_no_partners_returns_empty(monkeypatch):
    patch_sources(monkeypatch)
    result = net.collect_ppi_partners("APP", sources=["signor", "string"])
    assert result["partners"] == []
    assert result["excluded_hubs"] == []
