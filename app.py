"""Gene-Disease Network Overlap Evaluator.

A small Flask app: give it a disease name and a list of genes, and for each
gene it fetches PPI partners, then scores how much of the disease's Open
Targets gene signal that neighbourhood covers.
"""
from __future__ import annotations

import concurrent.futures
import os
from typing import Any

from flask import Flask, jsonify, render_template, request

from collectors.opentargets import (
    OpenTargetsError,
    get_disease_top_genes,
    resolve_disease_id,
)
from collectors.ppi import get_ppi_partners
from nw_overlap import MODERATE_THRESHOLD, STRONG_THRESHOLD, calc_network_overlap

# Defaults; every one is overridable through the environment.
DEFAULT_MAX_GENES = 30
DEFAULT_DISEASE_TOP_N = 100
DEFAULT_PPI_TOP_N = 30
DEFAULT_MAX_WORKERS = 5


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


def create_app(config: dict[str, Any] | None = None) -> Flask:
    app = Flask(__name__)
    app.config.update(
        BIOGRID_KEY=os.environ.get("BIOGRID_KEY", ""),
        MAX_GENES=_env_int("NW_MAX_GENES", DEFAULT_MAX_GENES),
        DISEASE_TOP_N=_env_int("NW_DISEASE_TOP_N", DEFAULT_DISEASE_TOP_N),
        PPI_TOP_N=_env_int("NW_PPI_TOP_N", DEFAULT_PPI_TOP_N),
        MAX_WORKERS=_env_int("NW_MAX_WORKERS", DEFAULT_MAX_WORKERS),
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
        )

    @app.route("/healthz")
    def healthz():
        return jsonify(
            {
                "status": "ok",
                "max_genes": app.config["MAX_GENES"],
                "disease_top_n": app.config["DISEASE_TOP_N"],
                "ppi_top_n": app.config["PPI_TOP_N"],
                "biogrid_enabled": bool(app.config["BIOGRID_KEY"]),
            }
        )

    @app.route("/api/analyze", methods=["POST"])
    def analyze():
        data = request.get_json(silent=True) or {}
        disease_name = str(data.get("disease") or "").strip()
        genes = parse_genes(data.get("genes"))

        if not disease_name or not genes:
            return jsonify({"error": "disease and genes are required"}), 400

        max_genes = app.config["MAX_GENES"]
        if len(genes) > max_genes:
            return (
                jsonify(
                    {
                        "error": (
                            f"Too many genes: {len(genes)} submitted, "
                            f"limit is {max_genes}."
                        )
                    }
                ),
                400,
            )

        top_n = app.config["DISEASE_TOP_N"]
        if isinstance(data.get("top_n"), int) and data["top_n"] > 0:
            top_n = min(data["top_n"], 500)

        try:
            disease_id, disease_label = resolve_disease_id(disease_name)
            if not disease_id:
                return jsonify({"error": f"Disease not found: {disease_name}"}), 404
            disease_genes = get_disease_top_genes(disease_id, top_n=top_n)
        except OpenTargetsError as exc:
            return jsonify({"error": str(exc)}), 502

        if not disease_genes:
            return (
                jsonify(
                    {
                        "error": (
                            f"No associated genes found for {disease_label} ({disease_id})."
                        )
                    }
                ),
                404,
            )

        biogrid_key = app.config.get("BIOGRID_KEY", "")
        ppi_top_n = app.config["PPI_TOP_N"]

        def process_gene(gene: str) -> dict[str, Any]:
            partners = get_ppi_partners(gene, biogrid_key=biogrid_key, top_n=ppi_top_n)
            overlap = calc_network_overlap(gene, partners, disease_genes)
            return {"gene": gene, **overlap}

        results: list[dict[str, Any]] = []
        workers = max(1, min(app.config["MAX_WORKERS"], len(genes)))
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(process_gene, gene): gene for gene in genes}
            for future in concurrent.futures.as_completed(futures):
                try:
                    results.append(future.result())
                except Exception as exc:  # one bad gene must not sink the request
                    results.append({"gene": futures[future], "error": str(exc)})

        results.sort(key=lambda r: r.get("weighted_score", 0), reverse=True)

        return jsonify(
            {
                "disease_id": disease_id,
                "disease_label": disease_label,
                "disease_gene_count": len(disease_genes),
                "ppi_top_n": ppi_top_n,
                "results": results,
            }
        )

    return app


app = create_app()


if __name__ == "__main__":
    app.run(
        host=os.environ.get("HOST", "127.0.0.1"),
        port=_env_int("PORT", 5005),
        debug=os.environ.get("FLASK_DEBUG", "").lower() in {"1", "true", "yes"},
    )
