"""PPI network construction and partner ranking.

Ported from the original project's ``network.py``: interactions from the
selected sources are merged into one graph centred on the query gene, then
partners are ranked by how reproducible and how well-scored the interaction is.

Ranking order (descending):
  1. number of distinct databases supporting the edge (reproducibility)
  2. representative score — the first source-specific score found in
     ``DB_SCORE_PRIORITY`` (SIGNOR > BioGRID > STRING); when no source supplies
     one, the inverse of the partner's global interactor count is substituted
  3. edge weight (how many times the interaction was observed)
"""
from __future__ import annotations

import json
import logging
import os
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Optional

import networkx as nx

from collectors import biogrid, signor, string_db
from collectors._cache import cache_enabled, cache_path
from collectors._http import SESSION

log = logging.getLogger(__name__)

HUMAN_TAX_ID = 9606

# ── Hub detection ────────────────────────────────────────────────────────────
# Two objective rules, OR'd: a very high global interactor count, or membership
# of a family known to interact promiscuously (G proteins, ubiquitin, tubulin,
# ...) which the count alone does not always catch.
HUB_DEGREE_THRESHOLD = int(os.environ.get("HUB_DEGREE_THRESHOLD", "1000"))

HUB_FAMILY_PATTERNS = [
    r"^GNA[0-9A-Z]",     # G protein alpha subunits
    r"^GNB[0-9]",        # G protein beta subunits
    r"^GNG[0-9]",        # G protein gamma subunits
    r"^UB[ABC][0-9]?$",  # ubiquitin
    r"^RPS27A$",         # ubiquitin fusion
    r"^TUB[AB]",         # tubulin
    r"^ACT[BG][0-9]?$",  # actin
    r"^HSP(90|A)",       # major chaperones
    r"^YWHA[BEGHQZ]$",   # 14-3-3
    r"^SUMO[0-9]$",      # SUMO
]
_HUB_FAMILY_RE = re.compile("|".join(HUB_FAMILY_PATTERNS))

INTACT_COUNT_URL = os.environ.get(
    "INTACT_URL", "https://www.ebi.ac.uk/intact/ws/interaction/findInteractions/{gene}"
)
_HUB_CACHE_TTL = 30 * 24 * 3600  # 30 days; interactor counts change slowly
_hub_mem_cache: dict[str, int] = {}
_hub_lock = threading.Lock()

# Which source's score to prefer when several support the same partner.
DB_SCORE_PRIORITY = ["SIGNOR", "BioGRID", "STRING"]

DEFAULT_SOURCES = ("signor", "string")
DEFAULT_STRING_SCORE = 700
DEFAULT_MAX_NODES = 100


def is_hub_family(gene_symbol: str) -> bool:
    """True when the symbol belongs to a known promiscuous-hub family."""
    return bool(_HUB_FAMILY_RE.match((gene_symbol or "").upper()))


def reset_hub_cache() -> None:
    """Clear the in-memory interactor-count cache (used by tests)."""
    with _hub_lock:
        _hub_mem_cache.clear()


def global_interactor_count(gene_symbol: str) -> int:
    """Return the gene's global IntAct interactor count, or -1 if unknown.

    -1 means "could not determine" and is never treated as hub evidence.
    """
    key = (gene_symbol or "").upper()
    if not key:
        return -1

    with _hub_lock:
        if key in _hub_mem_cache:
            return _hub_mem_cache[key]

    if cache_enabled():
        path = cache_path("hub_degree", f"{key}.json")
        if path.exists() and time.time() - path.stat().st_mtime < _HUB_CACHE_TTL:
            try:
                count = int(json.loads(path.read_text())["count"])
                with _hub_lock:
                    _hub_mem_cache[key] = count
                return count
            except Exception:
                pass

    try:
        response = SESSION.get(
            INTACT_COUNT_URL.format(gene=gene_symbol),
            params={"page": 0, "pageSize": 1, "query": f"species:{HUMAN_TAX_ID}"},
            timeout=(5, 15),
        )
        response.raise_for_status()
        count = int(response.json().get("totalElements", -1))
    except Exception:
        count = -1

    if count >= 0 and cache_enabled():
        try:
            cache_path("hub_degree", f"{key}.json").write_text(json.dumps({"count": count}))
        except Exception:
            pass
    with _hub_lock:
        _hub_mem_cache[key] = count
    return count


# ── Network construction ─────────────────────────────────────────────────────

def build_ppi_network(
    gene_symbol: str,
    use_signor: bool = True,
    use_string: bool = False,
    use_biogrid: bool = False,
    biogrid_api_key: Optional[str] = None,
    string_required_score: int = DEFAULT_STRING_SCORE,
    min_score: Optional[float] = None,
) -> nx.Graph:
    """Build the interaction graph around ``gene_symbol`` from the chosen sources.

    Args:
        use_signor / use_string / use_biogrid: which sources to query.
        string_required_score: STRING's own 0-1000 confidence threshold.
        min_score: drop edges scoring below this after the sources are merged.

    Every source is best-effort: one failing leaves the others intact.
    """
    graph = nx.Graph()
    center = (gene_symbol or "").strip().upper()
    if not center:
        return graph
    graph.add_node(center, db="center", direct_partner=True, entity_type="gene")

    def _item_partners(item: dict) -> list[tuple[str, Optional[float], str]]:
        score = item.get("score", item.get("confidence", None))
        try:
            score = float(score) if score is not None else None
        except (TypeError, ValueError):
            score = None

        ptype = (item.get("partner_type") or "gene").strip().lower()
        out: list[tuple[str, Optional[float], str]] = []
        src = str(item.get("source") or "").strip().upper()
        tgt = str(item.get("target") or "").strip().upper()
        if src and tgt:
            partner = tgt if src == center else src
            if partner and partner != center:
                out.append((partner, score, ptype))
        for p in item.get("partners", []) or []:
            p = str(p).strip().upper()
            if p and p != center:
                out.append((p, score, ptype))
        return out

    def add_edges(interactions: list[dict], source_label: str) -> None:
        for item in interactions:
            for partner, score, ptype in _item_partners(item):
                if partner not in graph:
                    graph.add_node(
                        partner, db=source_label, direct_partner=True, entity_type=ptype
                    )
                else:
                    if source_label not in graph.nodes[partner].get("db", ""):
                        graph.nodes[partner]["db"] += f",{source_label}"
                    # "gene" wins if any source calls the partner a gene.
                    if ptype == "gene":
                        graph.nodes[partner]["entity_type"] = "gene"

                if graph.has_edge(center, partner):
                    edge = graph.edges[center, partner]
                    edge["weight"] = edge.get("weight", 1) + 1
                    edge.setdefault("dbs", set()).add(source_label)
                    if source_label not in edge.get("db", ""):
                        edge["db"] += f",{source_label}"
                    if score is not None:
                        prev = edge.get("score")
                        edge["score"] = score if prev is None else max(prev, score)
                        db_scores = edge.setdefault("db_scores", {})
                        prev_db = db_scores.get(source_label)
                        db_scores[source_label] = (
                            score if prev_db is None else max(prev_db, score)
                        )
                else:
                    graph.add_edge(
                        center,
                        partner,
                        weight=1,
                        effect=item.get("effect", ""),
                        mechanism=item.get("mechanism", ""),
                        db=source_label,
                        dbs={source_label},
                        score=score,
                        db_scores=({source_label: score} if score is not None else {}),
                    )

    if use_signor:
        try:
            add_edges(signor.get_interactions(gene_symbol), "SIGNOR")
        except Exception as exc:
            log.warning("[SIGNOR] %s: %s", gene_symbol, exc)

    if use_string:
        try:
            add_edges(
                string_db.get_interactions(
                    gene_symbol, required_score=string_required_score
                ),
                "STRING",
            )
        except Exception as exc:
            log.warning("[STRING] %s: %s", gene_symbol, exc)

    if use_biogrid:
        try:
            add_edges(biogrid.get_interactions(gene_symbol, api_key=biogrid_api_key), "BioGRID")
        except Exception as exc:
            log.warning("[BioGRID] %s: %s", gene_symbol, exc)

    # Sources such as BioGRID often leave SCORE empty. For those edges the
    # inverse of the partner's global interactor count stands in: a partner with
    # few interactions overall makes for a more specific, more informative edge.
    no_score_edges = [(u, v) for u, v, d in graph.edges(data=True) if d.get("score") is None]
    if no_score_edges:
        partners = [v if u == center else u for u, v in no_score_edges]
        with ThreadPoolExecutor(max_workers=10) as executor:
            counts = dict(zip(partners, executor.map(global_interactor_count, partners)))
        for u, v in no_score_edges:
            partner = v if u == center else u
            count = counts.get(partner, -1)
            if count and count > 0:
                graph.edges[u, v]["score"] = 1.0 / count
                graph.edges[u, v]["score_inferred"] = True

    if min_score is not None:
        drop = [
            (u, v)
            for u, v, d in graph.edges(data=True)
            if d.get("score") is not None and d["score"] < min_score
        ]
        graph.remove_edges_from(drop)
        isolated = [n for n in list(graph.nodes) if n != center and graph.degree(n) == 0]
        graph.remove_nodes_from(isolated)

    return graph


def _n_distinct_dbs(edge: dict) -> int:
    """How many distinct databases support this edge (a reproducibility proxy)."""
    dbs = edge.get("dbs")
    if not dbs:
        dbs = {d for d in (edge.get("db", "") or "").split(",") if d}
    return len(dbs)


def _representative_score(edge: dict) -> float:
    """Pick the edge's representative score by source priority."""
    db_scores = edge.get("db_scores") or {}
    for db in DB_SCORE_PRIORITY:
        score = db_scores.get(db)
        if score is not None:
            return score
    fallback = edge.get("score")
    return fallback if fallback is not None else -1.0


def rank_partners(
    graph: nx.Graph,
    center: str,
    exclude_hubs: bool = True,
    exclude_non_gene: bool = True,
    hub_threshold: Optional[int] = None,
) -> list[str]:
    """Return the centre's partners, most reproducible and best-scored first.

    ``exclude_hubs`` removes promiscuous hubs (known families, or a global
    interactor count above ``hub_threshold``) so that non-specific connectors do
    not inflate the overlap. ``exclude_non_gene`` always drops chemical and
    phenotype nodes.
    """
    center = (center or "").strip().upper()
    if center not in graph:
        return []

    def sort_key(node: str):
        edge = graph.edges[center, node]
        return (_n_distinct_dbs(edge), _representative_score(edge), edge.get("weight", 1))

    ranked = sorted(graph.neighbors(center), key=sort_key, reverse=True)
    candidates = [n for n in ranked if n != center]

    excluded: set[str] = set()
    if exclude_non_gene:
        excluded |= {
            n for n in candidates if graph.nodes[n].get("entity_type", "gene") != "gene"
        }

    if exclude_hubs:
        threshold = hub_threshold if hub_threshold is not None else HUB_DEGREE_THRESHOLD
        family_hubs = {n for n in candidates if is_hub_family(n)}
        to_check = [n for n in candidates if n not in family_hubs and n not in excluded]
        degree_hubs: set[str] = set()
        if to_check:
            with ThreadPoolExecutor(max_workers=10) as executor:
                counts = dict(zip(to_check, executor.map(global_interactor_count, to_check)))
            degree_hubs = {n for n, c in counts.items() if c > threshold}
        excluded |= family_hubs | degree_hubs

    return [n for n in ranked if n not in excluded]


def collect_ppi_partners(
    gene_symbol: str,
    sources: tuple[str, ...] | list[str] = DEFAULT_SOURCES,
    string_required_score: int = DEFAULT_STRING_SCORE,
    min_score: Optional[float] = None,
    hub_threshold: int = HUB_DEGREE_THRESHOLD,
    exclude_hubs: bool = True,
    max_nodes: int = DEFAULT_MAX_NODES,
    biogrid_api_key: str = "",
) -> dict[str, Any]:
    """Build the network for one gene and return its ranked partners.

    Returns ``{"partners", "excluded_hubs", "source_counts", "graph_size"}``.
    """
    wanted = {str(s).strip().lower() for s in (sources or ())}
    graph = build_ppi_network(
        gene_symbol,
        use_signor="signor" in wanted,
        use_string="string" in wanted,
        use_biogrid="biogrid" in wanted and bool(biogrid_api_key),
        biogrid_api_key=biogrid_api_key,
        string_required_score=string_required_score,
        min_score=min_score,
    )

    center = (gene_symbol or "").strip().upper()
    ranked = rank_partners(
        graph, center, exclude_hubs=exclude_hubs, hub_threshold=hub_threshold
    )
    all_partners = rank_partners(graph, center, exclude_hubs=False, exclude_non_gene=True)

    source_counts: dict[str, int] = {}
    for _, _, data in graph.edges(data=True):
        for db in data.get("dbs") or {d for d in (data.get("db") or "").split(",") if d}:
            source_counts[db] = source_counts.get(db, 0) + 1

    return {
        "partners": ranked[: max(1, int(max_nodes))],
        "excluded_hubs": [p for p in all_partners if p not in set(ranked)],
        "source_counts": source_counts,
        "graph_size": graph.number_of_nodes(),
    }
