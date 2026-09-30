#!/usr/bin/env python3
"""
Resolve any HGNC-known gene identifier to a gene of the 19,354-gene universe.

WHY THIS SCRIPT EXISTS
----------------------
Users and tools query genes by whatever name they have: a previous symbol from an
older paper (AARS, FAM58A), an alias (P53), an Entrez or UniProt ID, or the CURRENT
HGNC symbol of a gene the universe still lists under its former name (26 genes, e.g.
the universe's AATK is HGNC's LMTK1 since the universe was built). A plain symbol
lookup returns "not found" for all of these. This script precomputes, from the pinned
HGNC snapshot, one resolution record per identifier, so the server only looks keys up
and never applies mapping rules itself.

UNIVERSE -> HGNC ID
-------------------
Each universe gene is given its HGNC ID by approved symbol, otherwise by a previous
symbol that names exactly one HGNC entry (the 26 renamed genes). The script stops if
any gene is left unmapped or two genes share an HGNC ID.

KEYS AND PRECEDENCE
-------------------
Every identifier of every HGNC entry (all locus types, 45,140 entries in the
2026-09-30 snapshot) becomes an upper-cased key. When a key matches in several ways,
the first tier that matches decides:

  1. index_symbol     the symbol the universe uses for a gene (C9orf72, AATK)
  2. approved_symbol  current HGNC approved symbol (LMTK1 -> universe AATK)
  3. hgnc_id          HGNC:28337
  4. entrez_id        203228
  5. ensembl_gene_id  ENSG00000147894
  6. uniprot_id       Q96LT7
  7. previous_symbol  AARS -> AARS1
  8. alias_symbol     P53 -> TP53

STATUS
------
  resolved      the deciding tier names exactly one HGNC entry, which is in the universe
  not_in_index  the deciding tier names exactly one HGNC entry, outside the universe
                (e.g. a non-coding RNA or pseudogene symbol)
  ambiguous     the deciding tier names more than one HGNC entry (candidates listed)

Matches at LOWER tiers that point to other entries are kept in other_matches (e.g.
MDR1: previous symbol of ABCB1, but also an alias of TBC1D9).

auto_resolve (a lookup may answer with the gene record directly) = status resolved AND
  - match_type is a symbol or ID of tiers 1-6 (lower-tier matches cannot override an
    exact symbol or ID; they are only reported), or
  - match_type is previous_symbol AND other_matches is empty.
Aliases never auto-resolve. Every other key is answered with its resolution record so
the caller decides.

OUTPUT (--output, gzipped TSV, one row per key)
------
query_key, status, match_type, gene (universe symbol or empty), hgnc_id, hgnc_symbol,
locus_type, auto_resolve (0/1), candidates (JSON, ambiguous only), other_matches
(JSON), note. Candidates and other_matches are lists of
{"hgnc_id", "hgnc_symbol", "gene", "match_type"}; gene is empty outside the universe.
"""
import argparse
import csv
import gzip
import json
import sys
from collections import Counter, defaultdict

import pandas as pd

TIERS = ["index_symbol", "approved_symbol", "hgnc_id", "entrez_id",
         "ensembl_gene_id", "uniprot_id", "previous_symbol", "alias_symbol"]
MULTI_COLS = {"previous_symbol": "prev_symbol", "alias_symbol": "alias_symbol",
              "ensembl_gene_id": "ensembl_gene_id", "entrez_id": "entrez_id",
              "uniprot_id": "uniprot_ids"}
COLUMNS = ["query_key", "status", "match_type", "gene", "hgnc_id", "hgnc_symbol",
           "locus_type", "auto_resolve", "candidates", "other_matches", "note"]


def split_field(value):
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    return [x.strip() for x in str(value).strip('"').split("|") if x.strip()]


def map_universe(genes, hgnc):
    """Universe symbol -> HGNC ID (approved symbol, else unique previous symbol)."""
    by_symbol = dict(zip(hgnc["symbol"], hgnc["hgnc_id"]))
    by_prev = defaultdict(set)
    for r in hgnc.itertuples(index=False):
        for s in split_field(r.prev_symbol):
            by_prev[s.upper()].add(r.hgnc_id)
    mapping, via_prev, unmapped = {}, [], []
    for g in genes:
        if g in by_symbol:
            mapping[g] = by_symbol[g]
        elif len(by_prev.get(g.upper(), ())) == 1:
            mapping[g] = next(iter(by_prev[g.upper()]))
            via_prev.append(g)
        else:
            unmapped.append(g)
    if unmapped:
        sys.exit(f"ERROR: {len(unmapped)} universe genes have no HGNC ID: {unmapped[:20]}")
    dup = [h for h, n in Counter(mapping.values()).items() if n > 1]
    if dup:
        sys.exit(f"ERROR: HGNC IDs shared by several universe genes: {dup[:20]}")
    return mapping, via_prev


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--genes", required=True, help="universe gene list (data/gene_all_score.txt, first column)")
    ap.add_argument("--hgnc", required=True, help="HGNC complete set (data/raw/hgnc_complete_set_YYYYMMDD.txt.gz)")
    ap.add_argument("--output", required=True, help="output .tsv.gz")
    args = ap.parse_args()

    genes = pd.read_csv(args.genes, sep="\t", dtype=str).iloc[:, 0].tolist()
    hgnc = pd.read_csv(args.hgnc, sep="\t", dtype=str, low_memory=False)
    hgnc = hgnc[hgnc["status"] == "Approved"]
    gene2hgnc, via_prev = map_universe(genes, hgnc)
    hgnc2gene = {h: g for g, h in gene2hgnc.items()}
    info = {r.hgnc_id: (r.symbol, r.locus_type) for r in hgnc.itertuples(index=False)}

    # key -> tier -> set of HGNC IDs
    keys = defaultdict(lambda: defaultdict(set))
    for g, h in gene2hgnc.items():
        keys[g.upper()]["index_symbol"].add(h)
    for r in hgnc.itertuples(index=False):
        keys[r.symbol.upper()]["approved_symbol"].add(r.hgnc_id)
        keys[r.hgnc_id.upper()]["hgnc_id"].add(r.hgnc_id)
        for tier, col in MULTI_COLS.items():
            for x in split_field(getattr(r, col)):
                keys[x.upper()][tier].add(r.hgnc_id)

    def cand(h, tier):
        return {"hgnc_id": h, "hgnc_symbol": info[h][0], "gene": hgnc2gene.get(h, ""),
                "match_type": tier}

    rows, stats = [], Counter()
    for key in sorted(keys):
        tiers = keys[key]
        decide = next(t for t in TIERS if tiers.get(t))
        chosen = tiers[decide]
        other = []
        seen = set(chosen)
        for t in TIERS[TIERS.index(decide) + 1:]:
            for h in sorted(tiers.get(t, ())):
                if h not in seen:
                    other.append(cand(h, t))
                    seen.add(h)
        rec = dict.fromkeys(COLUMNS, "")
        rec.update(query_key=key, match_type=decide,
                   other_matches=json.dumps(other, separators=(",", ":")) if other else "")
        if len(chosen) > 1:
            rec["status"] = "ambiguous"
            cands = [cand(h, decide) for h in sorted(chosen)]
            rec["candidates"] = json.dumps(cands, separators=(",", ":"))
            n_in = sum(1 for c in cands if c["gene"])
            rec["note"] = f"{decide} of {len(cands)} HGNC genes ({n_in} in iNDDx)"
            rec["auto_resolve"] = 0
        else:
            h = next(iter(chosen))
            sym, locus = info[h]
            rec.update(hgnc_id=h, hgnc_symbol=sym, locus_type=locus)
            if h in hgnc2gene:
                g = hgnc2gene[h]
                rec.update(status="resolved", gene=g)
                rec["auto_resolve"] = int(decide in TIERS[:6] or
                                          (decide == "previous_symbol" and not other))
                if decide == "index_symbol" and sym != g:
                    rec["note"] = f"iNDDx lists this gene as {g}; the current HGNC symbol is {sym}"
                elif decide == "approved_symbol":
                    rec["note"] = f"current HGNC symbol; iNDDx lists this gene under its former symbol {g}"
                elif decide in ("previous_symbol", "alias_symbol"):
                    rec["note"] = f"{decide.replace('_', ' ')} of {sym}" + (
                        f" (iNDDx symbol {g})" if g != sym else "")
            else:
                rec.update(status="not_in_index", auto_resolve=0,
                           note=f"HGNC {locus}; not among the 19,354 protein-coding genes of iNDDx")
            if other and rec["status"] != "ambiguous":
                rec["note"] += ("; " if rec["note"] else "") + \
                    f"also matches as {other[0]['match_type'].replace('_', ' ')} of " + \
                    ", ".join(c["hgnc_symbol"] for c in other[:5]) + ("..." if len(other) > 5 else "")
        stats[(rec["status"], decide)] += 1
        rows.append(rec)

    # consistency checks
    idx_rows = [r for r in rows if r["match_type"] == "index_symbol"]
    assert len(idx_rows) == len(genes), (len(idx_rows), len(genes))
    assert all(r["status"] == "resolved" and r["auto_resolve"] == 1 and
               gene2hgnc[r["gene"]] == r["hgnc_id"] for r in idx_rows)

    with gzip.open(args.output, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print(f"universe: {len(genes)} genes -> HGNC IDs ({len(via_prev)} via previous symbol: "
          f"{', '.join(via_prev)})", file=sys.stderr)
    print(f"keys: {len(rows)}  auto_resolve: {sum(r['auto_resolve'] for r in rows)}", file=sys.stderr)
    for (st, mtype), n in sorted(stats.items()):
        print(f"  {st:13s} {mtype:16s} {n}", file=sys.stderr)
    print(f"wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
