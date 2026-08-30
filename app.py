"""Gene-Disease Network Overlap Evaluator.

Give it a disease and a list of genes, and for each gene it builds a PPI
network from the selected sources, then scores how much of the disease's Open
Targets gene network that neighbourhood covers.
"""
from __future__ import annotations

import concurrent.futures
import os
from typing import Any

from flask import Flask, jsonify, render_template, request

from collectors.hgnc import resolve_gene_symbols
from collectors.opentargets import (
    OpenTargetsError,
    get_disease_label,
    get_disease_with_top_genes,
    is_ontology_id,
    normalise_ontology_id,
    resolve_disease_id,
    search_diseases,
)
from nw_overlap import MODERATE_THRESHOLD, STRONG_THRESHOLD, calc_network_overlap
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

VALID_SOURCES = ("signor", "string", "biogrid")


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
        "max_nodes": _num("max_nodes", int, 1, 500, defaults["max_nodes"]),
        "exclude_hubs": bool(raw.get("exclude_hubs", defaults["exclude_hubs"])),
    }


def create_app(config: dict[str, Any] | None = None) -> Flask:
    app = Flask(__name__)
    app.config.update(
        BIOGRID_KEY=os.environ.get("BIOGRID_KEY") or os.environ.get("BIOGRID_API_KEY", ""),
        MAX_GENES=_env_int("NW_MAX_GENES", DEFAULT_MAX_GENES),
        DISEASE_TOP_N=_env_int("NW_DISEASE_TOP_N", DEFAULT_DISEASE_TOP_N),
        MAX_WORKERS=_env_int("NW_MAX_WORKERS", DEFAULT_MAX_WORKERS),
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
            }
        )

    @app.route("/api/diseases")
    def diseases():
        """Search Open Targets for diseases matching a free-text query."""
        query = (request.args.get("q") or "").strip()
        if not query:
            return jsonify({"error": "q is required"}), 400
        try:
            limit = min(int(request.args.get("limit", 10)), 25)
        except (TypeError, ValueError):
            limit = 10

        try:
            if is_ontology_id(query):
                ontology_id = normalise_ontology_id(query)
                label = get_disease_label(ontology_id)
                if not label:
                    return jsonify({"query": query, "results": []})
                return jsonify(
                    {
                        "query": query,
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
            return jsonify({"query": query, "results": search_diseases(query, limit=limit)})
        except OpenTargetsError as exc:
            return jsonify({"error": str(exc)}), 502

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
            return {
                "gene": symbol,
                "input_gene": entry["input"],
                "symbol_status": entry["status"],
                "hgnc_id": entry.get("hgnc_id", ""),
                "gene_name": entry.get("name", ""),
                "excluded_hub_count": len(ppi["excluded_hubs"]),
                "source_counts": ppi["source_counts"],
                **overlap,
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

        results.sort(key=lambda r: r.get("weighted_score", 0), reverse=True)

        return jsonify(
            {
                "disease_id": disease_id,
                "disease_label": disease_label,
                "disease_gene_count": len(disease_genes),
                "ppi": ppi_options,
                "genes": {"resolved": resolved, "summary": _symbol_summary(resolved)},
                "results": results,
            }
        )

    return app


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
