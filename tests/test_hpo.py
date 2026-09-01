import pytest

from collectors import hpo
from tests.helpers import FakeResponse, FakeSession


@pytest.fixture(autouse=True)
def clean(monkeypatch):
    monkeypatch.setenv("PPI_CACHE_DISABLED", "1")
    hpo.reset_cache()
    yield
    hpo.reset_cache()


# The real files: phenotype.hpoa has a plain header, phenotype_to_genes.txt a
# "#Format:" one. Both spellings are exercised.
HPOA = """#description: HPO annotations
#version: 2025-01-15
database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration
OMIM:104300\tAlzheimer disease\t\tHP:0002354\tOMIM:104300\tTAS\t\t\t\t\tP\tHPO:x
OMIM:104300\tAlzheimer disease\t\tHP:0000726\tOMIM:104300\tTAS\t\t\t\t\tP\tHPO:x
OMIM:104300\tAlzheimer disease\t\tHP:0000006\tOMIM:104300\tTAS\t\t\t\t\tI\tHPO:x
ORPHA:1020\tAlzheimer disease\t\tHP:0002354\tORPHA:1020\tPCS\t\tHP:0040281\t\t\tP\tORPHA:x
ORPHA:1020\tAlzheimer disease\t\tHP:0001300\tORPHA:1020\tPCS\t\tHP:0040283\t\t\tP\tORPHA:x
ORPHA:1020\tAlzheimer disease\tNOT\tHP:0002960\tORPHA:1020\tPCS\t\t\t\t\tP\tORPHA:x
OMIM:999999\tOther disease\t\tHP:0009999\tOMIM:999999\tTAS\t\t\t\t\tP\tHPO:x
"""

GENES = """#Format: HPO-id<tab>HPO label<tab>entrez-gene-id<tab>entrez-gene-symbol<tab>Additional Info from G-D source<tab>G-D source<tab>disease-ID for link
HP:0002354\tMemory impairment\t351\tAPP\t-\tmim2gene\tOMIM:104300
HP:0002354\tMemory impairment\t4137\tMAPT\t-\tmim2gene\tOMIM:104300
HP:0002354\tMemory impairment\t351\tAPP\t-\torphadata\tORPHA:1020
HP:0000726\tDementia\t5663\tPSEN1\t-\tmim2gene\tOMIM:104300
HP:0001300\tParkinsonism\t6622\tSNCA\t-\torphadata\tORPHA:1020
"""
# The gene file's header actually uses literal tabs, not "<tab>".
GENES = GENES.replace("<tab>", "\t")


def session(hpoa=HPOA, genes=GENES):
    def handler(url, **kwargs):
        # Both files live under .../hpoa/, so match on the filename.
        if url.rstrip("/").endswith("phenotype.hpoa"):
            if hpoa is None:
                raise ConnectionError("hpoa down")
            return FakeResponse(body=hpoa)
        if url.rstrip("/").endswith("phenotype_to_genes.txt"):
            if genes is None:
                raise ConnectionError("genes down")
            return FakeResponse(body=genes)
        raise AssertionError(f"unexpected url {url}")
    return FakeSession(get=handler)


# --- disease -> phenotypes -------------------------------------------------

def test_phenotypes_are_read_for_an_omim_id(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    names = [p["hpo_id"] for p in hpo.get_disease_phenotypes(["OMIM:104300"])]
    assert set(names) == {"HP:0002354", "HP:0000726"}


def test_inheritance_terms_are_not_symptoms(monkeypatch):
    """Aspect I (inheritance) and C (course) are not phenotypic abnormalities."""
    monkeypatch.setattr(hpo, "SESSION", session())
    assert "HP:0000006" not in [p["hpo_id"] for p in hpo.get_disease_phenotypes(["OMIM:104300"])]


def test_orphanet_ids_are_read(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    got = {p["hpo_id"]: p for p in hpo.get_disease_phenotypes(["ORPHA:1020"])}
    # The NOT-qualified term is returned flagged, not dropped.
    assert set(got) == {"HP:0002354", "HP:0001300", "HP:0002960"}
    assert got["HP:0001300"]["resources"] == ["ORPHANET"]
    assert got["HP:0002960"]["excluded"] is True
    assert got["HP:0002354"]["excluded"] is False


@pytest.mark.parametrize("spelling", ["ORPHA:1020", "Orphanet:1020", "ORPHANET:1020",
                                      "orphacode:1020"])
def test_every_orphanet_spelling_resolves(monkeypatch, spelling):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert len(hpo.get_disease_phenotypes([spelling])) == 3


def test_omim_and_orphanet_annotations_are_merged(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    got = {p["hpo_id"]: p for p in hpo.get_disease_phenotypes(["OMIM:104300", "ORPHA:1020"])}
    assert set(got) == {"HP:0002354", "HP:0000726", "HP:0001300", "HP:0002960"}
    # Memory impairment is annotated by both.
    assert set(got["HP:0002354"]["resources"]) == {"OMIM", "ORPHANET"}


def test_frequency_comes_from_orphanet(monkeypatch):
    """OMIM rows usually carry no frequency; Orphanet is where it comes from."""
    monkeypatch.setattr(hpo, "SESSION", session())
    got = {p["hpo_id"]: p for p in hpo.get_disease_phenotypes(["OMIM:104300", "ORPHA:1020"])}
    assert got["HP:0002354"]["frequency"] == "HP:0040281"
    assert got["HP:0000726"]["frequency"] == ""      # OMIM-only, none available


def test_a_not_annotation_is_excluded(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    all_terms = {p["hpo_id"]: p for p in hpo.get_disease_phenotypes(["ORPHA:1020"])}
    assert all_terms["HP:0002960"]["excluded"] is True


def test_a_phenotype_only_some_sources_negate_is_kept(monkeypatch):
    hpoa = HPOA + "OMIM:104300\tAlzheimer disease\t\tHP:0002960\tOMIM:104300\tTAS\t\t\t\t\tP\tHPO:x\n"
    monkeypatch.setattr(hpo, "SESSION", session(hpoa=hpoa))
    got = {p["hpo_id"]: p for p in hpo.get_disease_phenotypes(["OMIM:104300", "ORPHA:1020"])}
    assert got["HP:0002960"]["excluded"] is False


def test_term_names_come_from_the_gene_file(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    got = {p["hpo_id"]: p for p in hpo.get_disease_phenotypes(["OMIM:104300"])}
    assert got["HP:0002354"]["name"] == "Memory impairment"


def test_an_unknown_disease_returns_nothing(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert hpo.get_disease_phenotypes(["OMIM:000000"]) == []


def test_the_limit_is_applied(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert len(hpo.get_disease_phenotypes(["OMIM:104300", "ORPHA:1020"], limit=2)) == 2
    assert len(hpo.get_disease_phenotypes(["OMIM:104300", "ORPHA:1020"])) == 4


def test_the_result_matches_the_open_targets_shape(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    p = hpo.get_disease_phenotypes(["ORPHA:1020"])[0]
    assert set(p) >= {"hpo_id", "ontology_id", "name", "description", "frequency",
                      "aspect", "resources", "resource", "excluded"}
    assert p["ontology_id"] == p["hpo_id"].replace(":", "_")


# --- phenotype -> genes ----------------------------------------------------

def test_genes_are_read_for_a_phenotype(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert [g["symbol"] for g in hpo.get_phenotype_genes("HP:0002354")] == ["APP", "MAPT"]


def test_duplicate_gene_rows_are_collapsed(monkeypatch):
    """APP is listed twice for HP:0002354, via OMIM and via Orphanet."""
    monkeypatch.setattr(hpo, "SESSION", session())
    symbols = [g["symbol"] for g in hpo.get_phenotype_genes("HP:0002354")]
    assert symbols.count("APP") == 1


def test_curated_links_are_unweighted(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert all(g["score"] == 1.0 for g in hpo.get_phenotype_genes("HP:0002354"))


def test_the_underscore_ontology_form_is_accepted(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert hpo.get_phenotype_genes("HP_0002354") == hpo.get_phenotype_genes("HP:0002354")


def test_a_phenotype_with_no_genes_returns_empty(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert hpo.get_phenotype_genes("HP:0009999") == []


def test_gene_limit_is_applied(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert len(hpo.get_phenotype_genes("HP:0002354", top_n=1)) == 1


# --- name lookup -----------------------------------------------------------

def test_disease_ids_can_be_found_by_exact_name(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert set(hpo.find_disease_ids_by_name("Alzheimer disease")) == {"OMIM:104300", "ORPHA:1020"}


def test_name_lookup_is_case_insensitive(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert hpo.find_disease_ids_by_name("alzheimer DISEASE")


def test_a_partial_name_does_not_match(monkeypatch):
    """A substring match would pick the wrong rare disease far too often."""
    monkeypatch.setattr(hpo, "SESSION", session())
    assert hpo.find_disease_ids_by_name("Alzheimer") == []


# --- disease registry search -----------------------------------------------

def test_the_registry_lists_every_annotated_disease(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    by_id = {d["id"]: d for d in hpo.list_diseases()}
    assert set(by_id) == {"OMIM:104300", "ORPHA:1020", "OMIM:999999"}
    assert by_id["ORPHA:1020"]["source"] == "ORPHANET"
    assert by_id["OMIM:104300"]["source"] == "OMIM"


def test_phenotype_counts_exclude_inheritance_and_negated_terms(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    by_id = {d["id"]: d for d in hpo.list_diseases()}
    # OMIM:104300: two P terms, plus an inheritance term that must not count.
    assert by_id["OMIM:104300"]["phenotype_count"] == 2
    # ORPHA:1020: two P terms, plus a NOT-qualified one that must not count.
    assert by_id["ORPHA:1020"]["phenotype_count"] == 2


def test_search_finds_a_disease_by_exact_name(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    results = hpo.search_diseases("Alzheimer disease")
    assert {r["id"] for r in results} == {"OMIM:104300", "ORPHA:1020"}
    assert all(r["exact"] for r in results)


def test_search_matches_a_partial_name(monkeypatch):
    """Unlike the automatic fallback, interactive search does partial matching."""
    monkeypatch.setattr(hpo, "SESSION", session())
    assert {r["id"] for r in hpo.search_diseases("alzheimer")} == {"OMIM:104300", "ORPHA:1020"}


def test_search_is_case_insensitive(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert hpo.search_diseases("ALZHEIMER DISEASE")


def test_exact_matches_rank_above_partial_ones(monkeypatch):
    hpoa = HPOA + "OMIM:111111\tAlzheimer disease, familial, type 3\t\tHP:0002354\tX\tTAS\t\t\t\t\tP\tHPO:x\n"
    monkeypatch.setattr(hpo, "SESSION", session(hpoa=hpoa))
    results = hpo.search_diseases("Alzheimer disease")
    assert results[0]["exact"] is True
    assert results[-1]["id"] == "OMIM:111111"


def test_better_annotated_diseases_rank_first_within_a_tier(monkeypatch):
    extra = "".join(
        f"ORPHA:1020\tAlzheimer disease\t\tHP:00{i:05d}\tX\tPCS\t\t\t\t\tP\tORPHA:x\n"
        for i in range(3)
    )
    monkeypatch.setattr(hpo, "SESSION", session(hpoa=HPOA + extra))
    results = hpo.search_diseases("Alzheimer disease")
    # Same name and tier, so the better-annotated registration comes first.
    assert [r["id"] for r in results] == ["ORPHA:1020", "OMIM:104300"]


def test_both_registrations_of_one_disease_are_offered(monkeypatch):
    """A disease is often registered in both OMIM and Orphanet."""
    monkeypatch.setattr(hpo, "SESSION", session())
    sources = {r["source"] for r in hpo.search_diseases("Alzheimer disease")}
    assert sources == {"OMIM", "ORPHANET"}


def test_search_respects_the_limit(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert len(hpo.search_diseases("disease", limit=1)) == 1


def test_search_of_a_blank_query_returns_nothing(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert hpo.search_diseases("   ") == []


def test_search_with_no_match_returns_nothing(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert hpo.search_diseases("no such disease anywhere") == []


# --- HP terms (HPO's own ids) ----------------------------------------------

@pytest.mark.parametrize("value,expected", [
    ("HP:0002354", "HP:0002354"), ("HP_0002354", "HP:0002354"),
    ("hp0002354", "HP:0002354"), ("  HP:0002354 ", "HP:0002354"),
    ("HP:234", ""), ("OMIM:104300", ""), ("APP", ""), ("", ""),
])
def test_hp_ids_are_normalised(value, expected):
    assert hpo.normalise_hpo_id(value) == expected


def test_only_terms_with_genes_are_listed(monkeypatch):
    """A term with no genes cannot become a matrix column."""
    monkeypatch.setattr(hpo, "SESSION", session())
    listed = {p["hpo_id"] for p in hpo.list_phenotypes()}
    assert listed == {"HP:0002354", "HP:0000726", "HP:0001300"}
    assert "HP:0009999" not in listed


def test_terms_carry_their_gene_count(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    by_id = {p["hpo_id"]: p for p in hpo.list_phenotypes()}
    assert by_id["HP:0002354"]["gene_count"] == 2      # APP, MAPT
    assert by_id["HP:0002354"]["name"] == "Memory impairment"


def test_a_term_can_be_looked_up_by_its_id(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    got = hpo.search_phenotypes("HP:0002354")
    assert [p["hpo_id"] for p in got] == ["HP:0002354"]
    assert got[0]["exact"] is True


def test_the_underscore_id_form_is_accepted_in_search(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert hpo.search_phenotypes("HP_0002354")[0]["hpo_id"] == "HP:0002354"


def test_an_unknown_id_returns_nothing(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert hpo.search_phenotypes("HP:9999999") == []


def test_terms_are_searched_by_name(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert [p["hpo_id"] for p in hpo.search_phenotypes("Memory impairment")] == ["HP:0002354"]


def test_term_search_matches_partially(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert [p["hpo_id"] for p in hpo.search_phenotypes("memory")] == ["HP:0002354"]


def test_term_search_respects_the_limit(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert len(hpo.search_phenotypes("i", limit=1)) == 1


def test_chosen_terms_come_back_in_the_standard_shape(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    got = hpo.get_phenotypes_by_id(["HP_0002354", "HP:0000726"])
    assert [p["hpo_id"] for p in got] == ["HP:0002354", "HP:0000726"]
    assert got[0]["name"] == "Memory impairment"
    assert got[0]["ontology_id"] == "HP_0002354"
    # Nothing is inferred: no disease annotation is involved.
    assert got[0]["excluded"] is False
    assert got[0]["frequency"] == "" and got[0]["resources"] == []


def test_chosen_terms_keep_the_given_order_and_deduplicate(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    got = hpo.get_phenotypes_by_id(["HP:0000726", "HP:0002354", "HP_0000726"])
    assert [p["hpo_id"] for p in got] == ["HP:0000726", "HP:0002354"]


def test_invalid_term_ids_are_dropped(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert hpo.get_phenotypes_by_id(["OMIM:104300", "", "nonsense", None]) == []


def test_an_unannotated_term_is_still_returned_with_zero_genes(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    got = hpo.get_phenotypes_by_id(["HP:0009999"])
    assert got[0]["hpo_id"] == "HP:0009999"
    assert got[0]["gene_count"] == 0
    assert got[0]["name"] == "HP:0009999"     # no label available, not invented


# --- degraded inputs -------------------------------------------------------

def test_files_are_downloaded_once(monkeypatch):
    s = session()
    monkeypatch.setattr(hpo, "SESSION", s)
    hpo.get_disease_phenotypes(["OMIM:104300"])
    hpo.get_disease_phenotypes(["ORPHA:1020"])
    hpo.get_phenotype_genes("HP:0002354")
    assert len(s.get_calls) == 2   # one hpoa, one gene file


def test_an_unreachable_hpoa_returns_nothing(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session(hpoa=None))
    assert hpo.get_disease_phenotypes(["OMIM:104300"]) == []


def test_an_unreachable_gene_file_still_yields_phenotypes(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session(genes=None))
    got = hpo.get_disease_phenotypes(["OMIM:104300"])
    assert [p["hpo_id"] for p in got]
    assert got[0]["name"] == got[0]["hpo_id"]   # no label available, not invented


def test_a_reordered_release_still_parses(monkeypatch):
    """Columns are matched by name, so reordering must not shift fields."""
    reordered = (
        "hpo_id\tdatabase_id\tqualifier\tdisease_name\taspect\tfrequency\n"
        "HP:0002354\tOMIM:104300\t\tAlzheimer disease\tP\tHP:0040282\n"
    )
    monkeypatch.setattr(hpo, "SESSION", session(hpoa=reordered))
    got = hpo.get_disease_phenotypes(["OMIM:104300"])
    assert [p["hpo_id"] for p in got] == ["HP:0002354"]
    assert got[0]["frequency"] == "HP:0040282"


def test_malformed_lines_are_skipped(monkeypatch):
    hpoa = HPOA + "garbage\n\t\t\nOMIM:104300\t\t\tNOTANHPOID\t\t\t\t\t\t\tP\t\n"
    monkeypatch.setattr(hpo, "SESSION", session(hpoa=hpoa))
    assert {p["hpo_id"] for p in hpo.get_disease_phenotypes(["OMIM:104300"])} == {
        "HP:0002354", "HP:0000726"}


def test_blank_ids_are_ignored(monkeypatch):
    monkeypatch.setattr(hpo, "SESSION", session())
    assert hpo.get_disease_phenotypes(["", None]) == []
