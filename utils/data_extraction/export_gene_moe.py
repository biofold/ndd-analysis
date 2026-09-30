#!/usr/bin/env python3
"""
Write the per-gene MOE table served by iNDDx (results/main/gene_moe.tsv).

The MOE score counts the libraries (GO BP, GO CC, GO MF, KEGG, Reactome) in
which a gene belongs to at least one term significant in the high-confidence
set (adjusted p < 0.01, strict, on the unrounded values written by
scripts/1_enrichr_all.py). scripts/2_supercandidate.py computes it for the
candidate set from the candidate-set enrichment tables; a candidate appears in
a term's overlap exactly when it belongs to the term, so scoring by library
membership gives the same value for candidates and extends it to every gene,
which is what the iNDDx web resource shows. This script computes it for all
genes in data/gene_all_score.txt and FAILS unless every candidate's score and
library set equal results/main/supercandidate.tsv.

Symbols are matched to each library ignoring case (all libs/ are upper case by
convention; HGNC writes C9orf72); the output keeps the HGNC spelling of
gene_all_score.txt. Matched terms are listed per library in the order of the
high-confidence enrichment table (ascending adjusted p-value).

Columns: Gene, Class, MOE_score, MOE_GOBP, MOE_GOCC, MOE_GOMF, MOE_KEGG,
MOE_Reactome, MOE_matched_terms ("Library:term;term|Library:term").

Usage:
  python3 utils/data_extraction/export_gene_moe.py --ndd-path . \
      --output results/main/gene_moe.tsv
"""
import argparse
import csv
import os
import sys

LIBS = [("GO_Biological_Process_2026", "MOE_GOBP"),
        ("GO_Cellular_Component_2026", "MOE_GOCC"),
        ("GO_Molecular_Function_2026", "MOE_GOMF"),
        ("KEGG_2021_Human", "MOE_KEGG"),
        ("Reactome_Pathways_2024", "MOE_Reactome")]
THRESHOLD = 0.01


def read_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\r\n").split("\t")
            sets[f[0]] = {g.upper() for g in f[2:] if g}
    return sets


def significant_terms(path):
    """HC-significant terms in table order (gseapy writes ascending adjusted p)."""
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        return [row["Term"] for row in r if float(row["Adjusted P-value"]) < THRESHOLD]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ndd-path", default=".", help="ndd-analysis root")
    ap.add_argument("--results", default=None, help="results/main directory (default <ndd-path>/results/main)")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    res = a.results or os.path.join(a.ndd_path, "results", "main")

    genes = []  # (symbol as in gene_all_score.txt, class)
    with open(os.path.join(a.ndd_path, "data", "gene_all_score.txt")) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t")
            genes.append((f[0], f[1]))

    matched = {g.upper(): {} for g, _ in genes}   # key -> {lib: [terms]}
    for lib, _ in LIBS:
        gmt = read_gmt(os.path.join(a.ndd_path, "libs", f"{lib}.gmt"))
        for term in significant_terms(os.path.join(res, f"gene_set2_{lib}.tsv")):
            for key in gmt[term] & matched.keys():
                matched[key].setdefault(lib, []).append(term)

    # Candidates must reproduce scripts/2_supercandidate.py exactly.
    with open(os.path.join(res, "supercandidate.tsv")) as fh:
        sc = {row["Gene"].upper(): (int(float(row["Score"])), set(filter(None, (row["Terms"] or "").split("|"))))
              for row in csv.DictReader(fh, delimiter="\t")}
    bad = [k for k, (s, libs) in sc.items() if len(matched.get(k, {})) != s or set(matched.get(k, {})) != libs]
    if bad:
        sys.exit(f"Error: {len(bad)} candidates differ from supercandidate.tsv (e.g. {', '.join(bad[:5])})")
    n_cand = sum(1 for _, c in genes if c == "candidate")
    if len(sc) != n_cand:
        sys.exit(f"Error: supercandidate.tsv has {len(sc)} genes, gene_all_score.txt {n_cand} candidates")

    with open(a.output, "w", newline="") as fo:
        fo.write("\t".join(["Gene", "Class", "MOE_score"] + [c for _, c in LIBS] + ["MOE_matched_terms"]) + "\n")
        for g, cls in genes:
            m = matched[g.upper()]
            flags = [str(int(lib in m)) for lib, _ in LIBS]
            terms = "|".join(f"{lib}:" + ";".join(m[lib]) for lib, _ in LIBS if lib in m)
            fo.write("\t".join([g, cls, str(len(m))] + flags + [terms]) + "\n")
    print(f"{a.output}: {len(genes)} genes; candidates identical to supercandidate.tsv ({len(sc)})", file=sys.stderr)


if __name__ == "__main__":
    main()
