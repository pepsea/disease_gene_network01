#!/usr/bin/env python3
"""Print the symptom (HPO phenotype) data Open Targets holds for a disease.

A probe, not part of the app: run it where there is network access to see what
`disease.phenotypes` actually returns before deciding how to use it.

    python3 scripts/fetch_symptoms.py
    python3 scripts/fetch_symptoms.py "cystic fibrosis" EFO_0000249
    python3 scripts/fetch_symptoms.py --limit 25 --json raw.json "Marfan syndrome"

The query is tried richest-first and falls back field-by-field, reporting which
variant the API accepted — so if the schema differs from what is assumed here,
the output says so instead of failing.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys

import requests

OT_API = os.environ.get(
    "OT_API_URL", "https://api.platform.opentargets.org/api/v4/graphql"
)
TIMEOUT = 40

# A spread of disease types: common/polygenic through rare/monogenic, where
# annotation density is expected to differ a lot.
DEFAULT_DISEASES = [
    "Alzheimer disease",
    "type 2 diabetes mellitus",
    "cystic fibrosis",
    "Marfan syndrome",
    "Duchenne muscular dystrophy",
]

ONTOLOGY_ID = re.compile(r"^(EFO|MONDO|HP|DOID|Orphanet|NCIT|GO|MP|OTAR)_\d+$", re.I)

# HPO frequency terms, so the raw HP: ids read as something.
FREQUENCY = {
    "HP:0040280": "Obligate (100%)",
    "HP:0040281": "Very frequent (80-99%)",
    "HP:0040282": "Frequent (30-79%)",
    "HP:0040283": "Occasional (5-29%)",
    "HP:0040284": "Very rare (1-4%)",
    "HP:0040285": "Excluded (0%)",
}
ASPECT = {
    "P": "表現型異常", "C": "臨床経過", "I": "遺伝形式",
    "M": "臨床的修飾因子", "H": "既往",
}

SEARCH = """
query S($q: String!) {
  search(queryString: $q, entityNames: ["disease"]) {
    hits { id name entity }
  }
}
"""

# Tried in order; the first the API accepts is used.
VARIANTS = [
    ("full", """
query D($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id name description
    therapeuticAreas { id name }
    synonyms { relation terms }
    phenotypes(page: {index: 0, size: $size}) {
      count
      rows {
        phenotypeHPO { id name description namespace }
        phenotypeEFO { id name }
        evidence { aspect frequency evidenceType qualifierNot resource }
      }
    }
  }
}
"""),
    ("no-evidence-detail", """
query D($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id name description
    therapeuticAreas { id name }
    phenotypes(page: {index: 0, size: $size}) {
      count
      rows { phenotypeHPO { id name description } }
    }
  }
}
"""),
    ("phenotypes-only", """
query D($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id name
    phenotypes(page: {index: 0, size: $size}) {
      count
      rows { phenotypeHPO { id name } }
    }
  }
}
"""),
    ("no-phenotypes", """
query D($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) { id name description therapeuticAreas { id name } }
}
"""),
]


def post(query: str, variables: dict) -> dict:
    r = requests.post(OT_API, json={"query": query, "variables": variables}, timeout=TIMEOUT)
    r.raise_for_status()
    return r.json()


def resolve(term: str) -> tuple[str | None, str]:
    if ONTOLOGY_ID.match(term):
        return term.upper().replace("ORPHANET", "Orphanet"), term
    data = post(SEARCH, {"q": term})
    if data.get("errors"):
        raise SystemExit(f"search failed: {data['errors']}")
    hits = [h for h in (data.get("data", {}).get("search") or {}).get("hits", [])
            if h.get("entity") == "disease"]
    if not hits:
        return None, term
    exact = [h for h in hits if (h.get("name") or "").casefold() == term.casefold()]
    best = exact[0] if exact else hits[0]
    return best["id"], best.get("name", term)


def fetch(disease_id: str, size: int) -> tuple[str, dict, list[str]]:
    """Return (variant_name, disease_payload, errors_seen)."""
    seen: list[str] = []
    for name, query in VARIANTS:
        data = post(query, {"efoId": disease_id, "size": size})
        if data.get("errors"):
            msgs = "; ".join(str(e.get("message", e)) for e in data["errors"])
            seen.append(f"[{name}] {msgs}")
            continue
        disease = (data.get("data") or {}).get("disease")
        if disease is None:
            seen.append(f"[{name}] disease not found")
            continue
        return name, disease, seen
    return "none", {}, seen


def show(term: str, size: int, raw_sink: list) -> None:
    print("=" * 78)
    disease_id, label = resolve(term)
    if not disease_id:
        print(f"{term}: 見つかりません")
        return

    variant, disease, errors = fetch(disease_id, size)
    raw_sink.append({"query": term, "id": disease_id, "variant": variant, "disease": disease})

    print(f"{disease.get('name', label)}  ({disease_id})   [query variant: {variant}]")
    for e in errors:
        print(f"  ! {e}")
    if not disease:
        return

    areas = ", ".join(a.get("name", "") for a in disease.get("therapeuticAreas") or [])
    if areas:
        print(f"  領域: {areas}")
    desc = (disease.get("description") or "").strip()
    if desc:
        print(f"  説明: {desc[:300]}{'…' if len(desc) > 300 else ''}")
    syn = disease.get("synonyms")
    if syn:
        terms = [t for group in syn for t in (group.get("terms") or [])]
        print(f"  同義語: {len(terms)} 件 — {', '.join(terms[:6])}")

    ph = disease.get("phenotypes")
    if ph is None:
        print("  症状: このクエリでは取得できませんでした")
        return
    rows = ph.get("rows") or []
    print(f"  症状(HPO): 全 {ph.get('count', len(rows))} 件 / 表示 {len(rows)} 件")
    if not rows:
        print("    （この疾患には表現型注釈がありません）")
        return

    for row in rows:
        hpo = row.get("phenotypeHPO") or {}
        bits = []
        for ev in row.get("evidence") or []:
            if ev.get("qualifierNot"):
                bits.append("NOT")
            f = ev.get("frequency")
            if f:
                bits.append(FREQUENCY.get(f, f))
            a = ev.get("aspect")
            if a:
                bits.append(ASPECT.get(a, a))
            r = ev.get("resource")
            if r:
                bits.append(str(r))
        extra = f"  [{' / '.join(dict.fromkeys(bits))}]" if bits else ""
        print(f"    - {hpo.get('id','?'):<12} {hpo.get('name','?')}{extra}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("diseases", nargs="*", default=None,
                    help="disease names or ontology ids (default: a sample spread)")
    ap.add_argument("--limit", type=int, default=15, help="phenotypes to show (default 15)")
    ap.add_argument("--json", metavar="PATH", help="also write the raw payloads here")
    args = ap.parse_args()

    targets = args.diseases or DEFAULT_DISEASES
    raw: list = []
    for term in targets:
        try:
            show(term, args.limit, raw)
        except Exception as exc:
            print(f"{term}: 取得失敗 — {type(exc).__name__}: {exc}", file=sys.stderr)
        print()

    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(raw, fh, ensure_ascii=False, indent=2)
        print(f"raw payloads -> {args.json}")


if __name__ == "__main__":
    main()
