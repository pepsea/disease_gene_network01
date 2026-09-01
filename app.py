"""Gene-Disease Network Overlap Evaluator.

Give it a disease and a list of genes, and for each gene it builds a PPI
network from the selected sources, then scores how much of the disease's Open
Targets gene network that neighbourhood covers.
"""
from __future__ import annotations

import concurrent.futures
import logging
import os
from typing import Any, Optional

from flask import Flask, jsonify, render_template, request

from collectors import hpo
from collectors.gprofiler import enrich_gene_list
from collectors.hgnc import resolve_gene_symbols
from collectors.opentargets import (
    OpenTargetsError,
    get_disease_label,
    get_disease_phenotypes,
    get_disease_top_genes,
    get_disease_with_top_genes,
    get_disease_xrefs,
    is_ontology_id,
    normalise_ontology_id,
    resolve_disease_id,
    search_diseases,
)
from enrichment_overlap import calc_enrichment_overlap
from nw_overlap import MODERATE_THRESHOLD, STRONG_THRESHOLD, calc_network_overlap
from symptom_stats import DEFAULT_BACKGROUND, apply_fdr, score_symptom_breadth
from symptom_genes import (
    DEFAULT_GENES_PER_PHENOTYPE,
    DEFAULT_MAX_PHENOTYPES,
    build_symptom_gene_set,
)
from ppi_network import (
    DEFAULT_MAX_NODES,
    DEFAULT_SOURCES,
    DEFAULT_STRING_SCORE,
    HUB_DEGREE_THRESHOLD,
    collect_ppi_partners,
)

DEFAULT_MAX_GENES = 30
DEFAULT_DISEASE_TOP_N = 100
DEFAULT_MAX_WORKERS = 5
# How many of the disease's top genes define the pathway signature. Defaults to
# the whole disease network, so the two tables describe the same disease.
DEFAULT_ENRICH_GENE_N = 100
DEFAULT_ENRICH_MAX_TERMS = 50
# How many symptoms to list, before the excluded ones are dropped.
DEFAULT_SYMPTOM_LIST_N = 50

VALID_SOURCES = ("signor", "string", "biogrid")

log = logging.getLogger(__name__)


def _env_int(name: str, default: int) -> int:
    try:
        value = int(os.environ.get(name, default))
    except (TypeError, ValueError):
        return default
    return value if value > 0 else default


def parse_genes(raw: Any) -> list[str]:
    """Normalise the gene input into a de-duplicated, ordered symbol list.

    Accepts a list of strings or a single blob using newline, comma, semicolon,
    tab or whitespace separators. Duplicates (case-insensitively) are dropped,
    keeping the first spelling the user typed.
    """
    if raw is None:
        return []
    if isinstance(raw, str):
        items = [raw]
    elif isinstance(raw, (list, tuple)):
        items = [str(item) for item in raw]
    else:
        return []

    genes: list[str] = []
    seen: set[str] = set()
    for item in items:
        for token in item.replace(",", "\n").replace(";", "\n").split():
            symbol = token.strip()
            if not symbol or symbol.upper() in seen:
                continue
            seen.add(symbol.upper())
            genes.append(symbol)
    return genes


def parse_ppi_options(raw: Any, defaults: dict[str, Any]) -> dict[str, Any]:
    """Validate the PPI block of the request against the configured defaults.

    Mirrors the original project's ``ppi`` payload: ``sources``,
    ``string_score``, ``min_score``, ``hub_threshold``, ``max_nodes`` and
    ``exclude_hubs``.
    """
    raw = raw if isinstance(raw, dict) else {}

    sources = raw.get("sources")
    if isinstance(sources, (list, tuple)):
        selected = [str(s).strip().lower() for s in sources]
        selected = [s for s in VALID_SOURCES if s in selected]  # fixed order
    else:
        selected = list(defaults["sources"])
    if not selected:
        selected = list(defaults["sources"])

    def _num(key: str, cast, lo, hi, default):
        value = raw.get(key)
        if value in (None, "", "null"):
            return default
        try:
            value = cast(value)
        except (TypeError, ValueError):
            return default
        return max(lo, min(value, hi))

    min_score = raw.get("min_score")
    if min_score in (None, "", "null"):
        min_score = None
    else:
        try:
            min_score = max(0.0, min(float(min_score), 1.0))
        except (TypeError, ValueError):
            min_score = None

    return {
        "sources": selected,
        "string_score": _num("string_score", int, 0, 1000, defaults["string_score"]),
        "min_score": min_score,
        "hub_threshold": _num("hub_threshold", int, 1, 10_000_000, defaults["hub_threshold"]),
        "max_nodes": _num("max_nodes", int, 1, 5000, defaults["max_nodes"]),
        "exclude_hubs": bool(raw.get("exclude_hubs", defaults["exclude_hubs"])),
    }


def create_app(config: dict[str, Any] | None = None) -> Flask:
    app = Flask(__name__)
    app.config.update(
        BIOGRID_KEY=os.environ.get("BIOGRID_KEY") or os.environ.get("BIOGRID_API_KEY", ""),
        MAX_GENES=_env_int("NW_MAX_GENES", DEFAULT_MAX_GENES),
        DISEASE_TOP_N=_env_int("NW_DISEASE_TOP_N", DEFAULT_DISEASE_TOP_N),
        MAX_WORKERS=_env_int("NW_MAX_WORKERS", DEFAULT_MAX_WORKERS),
        ENRICHMENT=True,
        ENRICH_GENE_N=_env_int("NW_ENRICH_GENE_N", DEFAULT_ENRICH_GENE_N),
        ENRICH_MAX_TERMS=_env_int("NW_ENRICH_MAX_TERMS", DEFAULT_ENRICH_MAX_TERMS),
        SYMPTOMS=True,
        SYMPTOM_LIST_N=_env_int("NW_SYMPTOM_LIST_N", DEFAULT_SYMPTOM_LIST_N),
        SYMPTOM_MAX_PHENOTYPES=_env_int("NW_SYMPTOM_MAX", DEFAULT_MAX_PHENOTYPES),
        SYMPTOM_GENES_PER=_env_int("NW_SYMPTOM_GENES_PER", DEFAULT_GENES_PER_PHENOTYPE),
        SYMPTOM_SOURCE=os.environ.get("NW_SYMPTOM_SOURCE", "auto"),
        SYMPTOM_BACKGROUND=_env_int("NW_SYMPTOM_BACKGROUND", DEFAULT_BACKGROUND),
        PPI_DEFAULTS={
            "sources": list(DEFAULT_SOURCES),
            "string_score": _env_int("NW_STRING_SCORE", DEFAULT_STRING_SCORE),
            "hub_threshold": _env_int("NW_HUB_THRESHOLD", HUB_DEGREE_THRESHOLD),
            "max_nodes": _env_int("NW_MAX_NODES", DEFAULT_MAX_NODES),
            "exclude_hubs": True,
        },
        VALIDATE_SYMBOLS=True,
    )
    if config:
        app.config.update(config)

    @app.route("/")
    def index():
        return render_template(
            "index.html",
            max_genes=app.config["MAX_GENES"],
            strong_threshold=STRONG_THRESHOLD,
            moderate_threshold=MODERATE_THRESHOLD,
            ppi_defaults=app.config["PPI_DEFAULTS"],
            biogrid_enabled=bool(app.config["BIOGRID_KEY"]),
            enrichment_default=app.config["ENRICHMENT"],
            symptoms_default=app.config["SYMPTOMS"],
            symptom_source_default=app.config["SYMPTOM_SOURCE"],
        )

    @app.route("/healthz")
    def healthz():
        return jsonify(
            {
                "status": "ok",
                "max_genes": app.config["MAX_GENES"],
                "disease_top_n": app.config["DISEASE_TOP_N"],
                "biogrid_enabled": bool(app.config["BIOGRID_KEY"]),
                "ppi_defaults": app.config["PPI_DEFAULTS"],
                "enrichment": app.config["ENRICHMENT"],
                "symptoms": app.config["SYMPTOMS"],
                "symptom_source": app.config["SYMPTOM_SOURCE"],
            }
        )

    @app.route("/api/diseases")
    def diseases():
        """Search Open Targets for diseases matching a free-text query."""
        query = (request.args.get("q") or "").strip()
        if not query:
            return jsonify({"error": "q is required"}), 400
        try:
            limit = min(int(request.args.get("limit", 50)), 200)
        except (TypeError, ValueError):
            limit = 50

        try:
            if is_ontology_id(query):
                ontology_id = normalise_ontology_id(query)
                label = get_disease_label(ontology_id)
                if not label:
                    return jsonify({"query": query, "results": []})
                return jsonify(
                    {
                        "query": query,
                        "total": 1,
                        "results": [
                            {
                                "id": ontology_id,
                                "name": label,
                                "description": "",
                                "exact": True,
                            }
                        ],
                    }
                )
            results = search_diseases(query, limit=limit)
            return jsonify({"query": query, "total": len(results), "results": results})
        except OpenTargetsError as exc:
            return jsonify({"error": str(exc)}), 502

    @app.route("/api/hpo/diseases")
    def hpo_diseases():
        """Search HPO's own disease registry (OMIM / Orphanet / DECIPHER)."""
        query = (request.args.get("q") or "").strip()
        if not query:
            return jsonify({"error": "q is required"}), 400
        try:
            limit = min(int(request.args.get("limit", 50)), 500)
        except (TypeError, ValueError):
            limit = 50
        results = hpo.search_diseases(query, limit=limit)
        total = results[0].get("match_total", len(results)) if results else 0
        return jsonify({"query": query, "total": total, "results": results})

    @app.route("/api/hpo/phenotypes")
    def hpo_phenotypes():
        """Search HPO's own terms (HP ids), to name the symptoms directly."""
        query = (request.args.get("q") or "").strip()
        if not query:
            return jsonify({"error": "q is required"}), 400
        try:
            limit = min(int(request.args.get("limit", 50)), 500)
        except (TypeError, ValueError):
            limit = 50
        results = hpo.search_phenotypes(query, limit=limit)
        total = results[0].get("match_total", len(results)) if results else 0
        return jsonify({"query": query, "total": total, "results": results})

    @app.route("/api/genes/validate", methods=["POST"])
    def validate_genes():
        """Check submitted symbols against HGNC without running an analysis."""
        data = request.get_json(silent=True) or {}
        genes = parse_genes(data.get("genes"))
        if not genes:
            return jsonify({"error": "genes are required"}), 400
        max_genes = app.config["MAX_GENES"]
        if len(genes) > max_genes:
            return jsonify({"error": _too_many(len(genes), max_genes)}), 400

        resolved = resolve_gene_symbols(genes)
        return jsonify({"genes": resolved, "summary": _symbol_summary(resolved)})

    @app.route("/api/analyze", methods=["POST"])
    def analyze():
        data = request.get_json(silent=True) or {}
        disease_input = str(data.get("disease") or "").strip()
        disease_id_input = str(data.get("disease_id") or "").strip()
        genes = parse_genes(data.get("genes"))

        if not (disease_input or disease_id_input) or not genes:
            return jsonify({"error": "disease and genes are required"}), 400

        max_genes = app.config["MAX_GENES"]
        if len(genes) > max_genes:
            return jsonify({"error": _too_many(len(genes), max_genes)}), 400

        ppi_options = parse_ppi_options(data.get("ppi"), app.config["PPI_DEFAULTS"])
        biogrid_key = app.config.get("BIOGRID_KEY", "")
        if "biogrid" in ppi_options["sources"] and not biogrid_key:
            ppi_options["sources"] = [s for s in ppi_options["sources"] if s != "biogrid"]
            ppi_options["biogrid_skipped"] = True

        top_n = app.config["DISEASE_TOP_N"]
        if isinstance(data.get("top_n"), int) and data["top_n"] > 0:
            top_n = min(data["top_n"], 500)

        # ── Disease: an id picked from the search endpoint wins over free text.
        try:
            if disease_id_input:
                disease_id = (
                    normalise_ontology_id(disease_id_input)
                    if is_ontology_id(disease_id_input)
                    else disease_id_input
                )
                disease_label = disease_input or disease_id
            else:
                disease_id, disease_label = resolve_disease_id(disease_input)
                if not disease_id:
                    return jsonify({"error": f"Disease not found: {disease_input}"}), 404

            fetched_label, disease_genes = get_disease_with_top_genes(disease_id, top_n=top_n)
            disease_label = fetched_label or disease_label
        except OpenTargetsError as exc:
            return jsonify({"error": str(exc)}), 502

        if not disease_genes:
            return (
                jsonify({"error": f"No associated genes found for {disease_label} ({disease_id})."}),
                404,
            )

        # ── Disease pathway signature: enriched once, shared by every gene.
        want_enrichment = data.get("enrichment", app.config["ENRICHMENT"])
        want_enrichment = bool(want_enrichment) and app.config["ENRICHMENT"]
        disease_pathways: list[dict[str, Any]] = []
        if want_enrichment:
            enrich_n = min(app.config["ENRICH_GENE_N"], len(disease_genes))
            disease_pathways = enrich_gene_list(
                [g["symbol"] for g in disease_genes[:enrich_n]],
                max_results=app.config["ENRICH_MAX_TERMS"],
            )

        # ── Symptoms: the disease's HPO phenotypes, expanded into a gene set.
        want_symptoms = data.get("symptoms", app.config["SYMPTOMS"])
        want_symptoms = bool(want_symptoms) and app.config["SYMPTOMS"]
        phenotypes: list[dict[str, Any]] = []
        symptom_set: dict[str, Any] = {"genes": [], "expanded": [], "empty": [],
                                       "failed": [], "phenotype_count": 0}
        symptom_meta: dict[str, Any] = {"phenotype_source": "", "xrefs": []}
        symptom_source = str(data.get("symptom_source") or app.config["SYMPTOM_SOURCE"])
        if symptom_source not in VALID_SYMPTOM_SOURCES:
            symptom_source = app.config["SYMPTOM_SOURCE"]
        if want_symptoms:
            chosen_ids = data.get("hpo_disease_ids")
            chosen_ids = chosen_ids if isinstance(chosen_ids, list) else []
            chosen_terms = data.get("hpo_phenotype_ids")
            chosen_terms = chosen_terms if isinstance(chosen_terms, list) else []
            phenotypes, symptom_meta = collect_symptoms(
                disease_id, disease_label, symptom_source,
                app.config["SYMPTOM_LIST_N"],
                hpo_disease_ids=chosen_ids[:10],
                hpo_phenotype_ids=chosen_terms[: app.config["SYMPTOM_MAX_PHENOTYPES"]],
            )
            if phenotypes:
                symptom_set = build_symptom_gene_set(
                    phenotypes,
                    genes_per_phenotype=app.config["SYMPTOM_GENES_PER"],
                    max_phenotypes=app.config["SYMPTOM_MAX_PHENOTYPES"],
                    max_workers=app.config["MAX_WORKERS"],
                    gene_fetcher=symptom_meta["gene_fetcher"],
                )
        symptom_genes = symptom_set["genes"]
        symptom_per_phenotype = symptom_set.get("per_phenotype") or []

        # ── Genes: map to approved HGNC symbols before analysing.
        if app.config["VALIDATE_SYMBOLS"]:
            resolved = resolve_gene_symbols(genes)
        else:
            resolved = [
                {"input": g, "symbol": g, "status": "unchecked", "hgnc_id": "", "name": ""}
                for g in genes
            ]

        def process_gene(entry: dict[str, Any]) -> dict[str, Any]:
            symbol = entry["symbol"]
            ppi = collect_ppi_partners(
                symbol,
                sources=ppi_options["sources"],
                string_required_score=ppi_options["string_score"],
                min_score=ppi_options["min_score"],
                hub_threshold=ppi_options["hub_threshold"],
                exclude_hubs=ppi_options["exclude_hubs"],
                max_nodes=ppi_options["max_nodes"],
                biogrid_api_key=biogrid_key,
            )
            overlap = calc_network_overlap(symbol, ppi["partners"], disease_genes)

            symptom_overlap: dict[str, Any] = {}
            if symptom_genes:
                symptom_overlap = _prefix_symptom(
                    calc_network_overlap(symbol, ppi["partners"], symptom_genes)
                )
                # Score each symptom on its own, so the gene x symptom matrix
                # shows where the overlap actually sits. The overall figure is
                # the mean across symptoms.
                cells = []
                for phenotype in symptom_per_phenotype:
                    cell = calc_network_overlap(symbol, ppi["partners"], phenotype["genes"])
                    cells.append(
                        {
                            "hpo_id": phenotype.get("hpo_id", ""),
                            "name": phenotype.get("name", ""),
                            "percent": cell["weighted_percent"],
                            "matched_count": cell["matched_count"],
                            "gene_count": cell["disease_gene_count"],
                            "target_self": bool(cell["target_self"]),
                        }
                    )
                symptom_overlap["symptom_cells"] = cells
                symptom_overlap["symptom_mean_percent"] = (
                    round(sum(c["percent"] for c in cells) / len(cells), 1)
                    if cells
                    else 0.0
                )
                # Is this gene's spread across symptoms more than chance? A
                # hypergeometric test per symptom, combined with Fisher's
                # method, so consistent coverage of many symptoms outweighs one
                # strong hit.
                symptom_overlap.update(
                    score_symptom_breadth(
                        cells,
                        ppi_partner_count=overlap["ppi_partner_count"],
                        background=app.config["SYMPTOM_BACKGROUND"],
                    )
                )

            enrichment: dict[str, Any] = {}
            if want_enrichment and disease_pathways:
                gene_pathways = enrich_gene_list(
                    [symbol] + ppi["partners"],
                    max_results=app.config["ENRICH_MAX_TERMS"],
                )
                enrichment = calc_enrichment_overlap(symbol, gene_pathways, disease_pathways)

            return {
                "gene": symbol,
                "input_gene": entry["input"],
                "symbol_status": entry["status"],
                "hgnc_id": entry.get("hgnc_id", ""),
                "gene_name": entry.get("name", ""),
                "excluded_hub_count": len(ppi["excluded_hubs"]),
                "source_counts": ppi["source_counts"],
                **overlap,
                **enrichment,
                **symptom_overlap,
            }

        results: list[dict[str, Any]] = []
        workers = max(1, min(app.config["MAX_WORKERS"], len(resolved)))
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(process_gene, e): e for e in resolved}
            for future in concurrent.futures.as_completed(futures):
                entry = futures[future]
                try:
                    results.append(future.result())
                except Exception as exc:  # one bad gene must not sink the request
                    results.append(
                        {
                            "gene": entry["symbol"],
                            "input_gene": entry["input"],
                            "symbol_status": entry["status"],
                            "error": str(exc),
                        }
                    )

        # Every gene was tested against the same symptom set, so correct for
        # multiple testing across them.
        apply_fdr(results)

        results.sort(key=lambda r: r.get("weighted_score", 0), reverse=True)

        return jsonify(
            {
                "disease_id": disease_id,
                "disease_label": disease_label,
                "disease_gene_count": len(disease_genes),
                "ppi": ppi_options,
                "genes": {"resolved": resolved, "summary": _symbol_summary(resolved)},
                "enrichment": {
                    "enabled": want_enrichment,
                    "disease_pathway_count": len(disease_pathways),
                    "disease_gene_n": min(app.config["ENRICH_GENE_N"], len(disease_genes)),
                    "top_pathways": [
                        {"term_id": p.get("term_id", ""), "name": p.get("name", ""),
                         "source": p.get("source", ""), "p_value": p.get("p_value", 1.0)}
                        for p in disease_pathways[:10]
                    ],
                },
                "symptoms": {
                    "enabled": want_symptoms,
                    "requested_source": symptom_source,
                    "phenotype_source": symptom_meta.get("phenotype_source", ""),
                    "gene_sources": sorted(
                        {e.get("gene_source", "") for e in symptom_set["expanded"]
                         if e.get("gene_source")}
                    ),
                    "xrefs": symptom_meta.get("xrefs", []),
                    "xref_origin": symptom_meta.get("xref_origin", ""),
                    "background": app.config["SYMPTOM_BACKGROUND"],
                    "phenotype_count": len(phenotypes),
                    "expanded_count": len(symptom_set["expanded"]),
                    "gene_count": len(symptom_genes),
                    "excluded_count": sum(1 for p in phenotypes if p.get("excluded")),
                    "unindexed_count": len(symptom_set["failed"]) + len(symptom_set["empty"]),
                    "phenotypes": [
                        {"hpo_id": p["hpo_id"], "name": p["name"],
                         "frequency": p.get("frequency", ""),
                         "resources": p.get("resources") or [],
                         "excluded": p.get("excluded", False)}
                        for p in phenotypes[:30]
                    ],
                    "expanded": symptom_set["expanded"][:30],
                    "matrix_symptoms": [
                        {"hpo_id": e.get("hpo_id", ""), "name": e.get("name", ""),
                         "frequency": e.get("frequency", ""),
                         "resources": e.get("resources") or [],
                         "gene_source": e.get("gene_source", ""),
                         "gene_count": e.get("gene_count", 0)}
                        for e in symptom_per_phenotype
                    ],
                    "top_genes": [
                        {"symbol": g["symbol"], "score": g["score"],
                         "phenotype_count": g["phenotype_count"],
                         "phenotypes": g["phenotypes"]}
                        for g in symptom_genes[:10]
                    ],
                },
                "results": results,
            }
        )

    return app


VALID_SYMPTOM_SOURCES = ("auto", "opentargets", "hpo")


def _hpo_gene_fetcher(ontology_id: str, top_n: int) -> dict[str, Any]:
    return {"genes": hpo.get_phenotype_genes(ontology_id, top_n=top_n), "source": "hpo"}


def _hybrid_gene_fetcher(ontology_id: str, top_n: int) -> dict[str, Any]:
    """Open Targets first, HPO's curated links when it yields nothing.

    Open Targets only resolves a phenotype's genes when it indexes that HP term
    as a disease, which is not guaranteed; HPO's own phenotype-gene file always
    has it.
    """
    try:
        genes = get_disease_top_genes(ontology_id, top_n=top_n)
    except OpenTargetsError:
        genes = []
    if genes:
        return {"genes": genes, "source": "opentargets"}
    return _hpo_gene_fetcher(ontology_id, top_n)


def collect_symptoms(
    disease_id: str,
    disease_label: str,
    source: str,
    list_n: int,
    hpo_disease_ids: Optional[list[str]] = None,
    hpo_phenotype_ids: Optional[list[str]] = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Get the symptoms to analyse, falling back from Open Targets to HPO.

    Args:
        hpo_disease_ids: OMIM / ORPHA / DECIPHER ids from HPO's disease
            registry. HPO annotates diseases under those namespaces — it has no
            disease ids of its own — so this identifies the disease without
            going through Open Targets' cross references.
        hpo_phenotype_ids: HP term ids (``HP:0002354``). These are HPO's own
            ids, and giving them names the symptoms outright: no disease is
            involved, and the matrix columns are exactly what was asked for.

    Returns ``(phenotypes, meta)``. ``meta`` records which source supplied the
    symptom list, the ids used, and the gene fetcher to expand them with.
    """
    meta: dict[str, Any] = {"phenotype_source": "", "xrefs": [],
                            "xref_origin": "", "gene_fetcher": _hybrid_gene_fetcher}

    # Symptoms named directly outrank everything: nothing has to be inferred.
    terms = [str(i).strip() for i in (hpo_phenotype_ids or []) if str(i).strip()]
    if terms:
        phenotypes = hpo.get_phenotypes_by_id(terms)
        meta["gene_fetcher"] = _hpo_gene_fetcher
        meta["xrefs"] = [p["hpo_id"] for p in phenotypes]
        meta["xref_origin"] = "phenotypes"
        if phenotypes:
            meta["phenotype_source"] = "hpo"
        return phenotypes, meta

    # A disease picked directly in HPO wins over anything derived from EFO.
    chosen = [str(i).strip() for i in (hpo_disease_ids or []) if str(i).strip()]
    if chosen:
        meta["xrefs"] = chosen
        meta["xref_origin"] = "selected"
        meta["gene_fetcher"] = _hpo_gene_fetcher
        phenotypes = hpo.get_disease_phenotypes(chosen, limit=list_n)
        if phenotypes:
            meta["phenotype_source"] = "hpo"
        return phenotypes, meta

    if source in ("auto", "opentargets"):
        try:
            phenotypes = get_disease_phenotypes(disease_id, limit=list_n)
        except OpenTargetsError as exc:
            log.info("[symptoms] Open Targets phenotypes unavailable: %s", exc)
            phenotypes = []
        if phenotypes:
            meta["phenotype_source"] = "opentargets"
            return phenotypes, meta
        if source == "opentargets":
            return [], meta

    # HPO direct: needs the OMIM / Orphanet / DECIPHER ids HPO keys on.
    xrefs: list[str] = []
    try:
        xrefs = get_disease_xrefs(disease_id)
    except OpenTargetsError as exc:
        log.info("[symptoms] cross references unavailable: %s", exc)
    meta["xref_origin"] = "dbxrefs" if xrefs else ""
    if not xrefs and disease_label:
        # No xref: fall back to an exact disease-name match inside HPO.
        xrefs = hpo.find_disease_ids_by_name(disease_label)
        if xrefs:
            meta["xref_origin"] = "name"

    meta["xrefs"] = xrefs
    meta["gene_fetcher"] = _hpo_gene_fetcher
    if not xrefs:
        return [], meta

    phenotypes = hpo.get_disease_phenotypes(xrefs, limit=list_n)
    if phenotypes:
        meta["phenotype_source"] = "hpo"
    return phenotypes, meta


# Table 3 reuses the gene-level scorer, so its keys are namespaced to sit
# alongside table 1's in the same result object.
_SYMPTOM_KEYS = {
    "weighted_score": "symptom_weighted_score",
    "weighted_percent": "symptom_weighted_percent",
    "simple_ratio": "symptom_simple_ratio",
    "overlap_percent": "symptom_overlap_percent",
    "overlap_count": "symptom_overlap_count",
    "matched_count": "symptom_matched_count",
    "disease_gene_count": "symptom_gene_count",
    "target_self": "symptom_target_self",
    "target_self_score": "symptom_target_self_score",
    "overlapping_genes": "symptom_overlapping_genes",
    "interpretation": "symptom_interpretation",
}


def _prefix_symptom(overlap: dict[str, Any]) -> dict[str, Any]:
    """Namespace the shared overlap keys for the symptom table."""
    return {new: overlap[old] for old, new in _SYMPTOM_KEYS.items() if old in overlap}


def _too_many(count: int, limit: int) -> str:
    return f"Too many genes: {count} submitted, limit is {limit}."


def _symbol_summary(resolved: list[dict[str, Any]]) -> dict[str, int]:
    summary: dict[str, int] = {}
    for entry in resolved:
        status = entry.get("status", "unknown")
        summary[status] = summary.get(status, 0) + 1
    return summary


app = create_app()


if __name__ == "__main__":
    app.run(
        host=os.environ.get("HOST", "127.0.0.1"),
        port=_env_int("PORT", 5005),
        debug=os.environ.get("FLASK_DEBUG", "").lower() in {"1", "true", "yes"},
    )
