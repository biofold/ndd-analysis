#!/usr/bin/env python3
"""
Build the GO Slim gene-set libraries (libs/GOslim_{Biological_Process,
Cellular_Component,Molecular_Function}_2026.gmt) with the GO Consortium's
map2slim (OWLTools).

WHY THIS SCRIPT EXISTS
----------------------
Reviewer 1, comment 2 asks for a redundancy-reduced GO analysis (GO Slim) for
the main figures. Enrichr provides no GO Slim library (its GO libraries are the
full ontology), so the slim libraries are built here from GO's own files:

  1. map2slim. `owltools <go.obo> --gaf <goa_human.gaf> --map2slim --subset
     goslim_generic --write-gaf <out>` maps every human GO annotation to the
     generic GO Slim. map2slim follows is_a, part_of and regulates links and
     assigns each annotation to its MOST SPECIFIC slim term(s) only, so a gene
     mapped to "transferase activity" is not also given "catalytic activity".
  2. Propagation. An enrichment gene set for a slim term must contain every
     gene annotated to it at any depth (the true-path rule), so each gene is
     additionally assigned to every slim term above its mapped slim terms,
     following is_a and part_of in the full GO graph.
  3. Filtering. Annotations with a NOT qualifier and ND ("no biological data
     available") annotations are removed; all other evidence codes, including
     IEA, are kept.

Inputs are pinned in data/raw (see data/README.md): go-basic.obo.gz,
goslim_generic.obo.gz and goa_human.gaf.gz. The OWLTools executable and a Java
runtime are needed only to rebuild the libraries, not to run the pipeline.

Usage:
  python3 utils/enrichment/build_goslim_gmt.py \\
      --owltools /path/to/owltools \\
      --go-obo data/raw/go-basic.obo.gz --gaf data/raw/goa_human.gaf.gz \\
      --outdir libs --report libs/GOslim_2026_build_report.tsv
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

ASPECT = {"P": "Biological_Process", "C": "Cellular_Component", "F": "Molecular_Function"}
GAF_COLS = 17


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_obo(path):
    """Return id -> {name, namespace, parents (is_a + part_of), subsets}."""
    terms, cur = {}, None
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = {"parents": set(), "subsets": set()}
                continue
            if line.startswith("["):
                cur = None
                continue
            if cur is None or not line:
                continue
            if line.startswith("id: "):
                cur["id"] = line[4:]
                terms[cur["id"]] = cur
            elif line.startswith("name: "):
                cur["name"] = line[6:]
            elif line.startswith("namespace: "):
                cur["namespace"] = line[11:]
            elif line.startswith("subset: "):
                cur["subsets"].add(line[8:].split()[0])
            elif line.startswith("is_a: "):
                cur["parents"].add(line[6:].split()[0])
            elif line.startswith("relationship: part_of "):
                cur["parents"].add(line.split()[2])
            elif line.startswith("is_obsolete: true"):
                cur["obsolete"] = True
    return terms


def obo_version(path):
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("data-version:"):
                return line.split(":", 1)[1].strip()
            if line.startswith("[Term]"):
                break
    return "unknown"


def gaf_date(path):
    with open_text(path) as fh:
        for line in fh:
            if not line.startswith("!"):
                break
            if line.startswith("!date-generated:"):
                return line.split(":", 1)[1].strip()
    return "unknown"


def ancestors(terms, gid, cache):
    if gid in cache:
        return cache[gid]
    out = set()
    for p in terms.get(gid, {}).get("parents", ()):
        out.add(p)
        out |= ancestors(terms, p, cache)
    cache[gid] = out
    return out


def run_map2slim(owltools, go_obo, gaf, subset, java_home, memory):
    """Run map2slim on uncompressed copies of the inputs; return slim GAF path."""
    tmp = tempfile.mkdtemp(prefix="map2slim_")
    def plain(src):
        if not src.endswith(".gz"):
            return src
        dst = os.path.join(tmp, os.path.basename(src)[:-3])
        with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        return dst
    obo_p, gaf_p = plain(go_obo), plain(gaf)
    out = os.path.join(tmp, "goslim_mapped.gaf")
    env = dict(os.environ, OWLTOOLS_MEMORY=memory)
    if java_home:
        env["PATH"] = os.path.join(java_home, "bin") + os.pathsep + env["PATH"]
    cmd = [owltools, obo_p, "--gaf", gaf_p, "--map2slim", "--subset", subset,
           "--write-gaf", out]
    log = os.path.join(tmp, "map2slim.log")
    with open(log, "w") as lf:
        r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, env=env)
    if r.returncode != 0 or not os.path.exists(out):
        sys.exit(f"Error: map2slim failed (exit {r.returncode}); see {log}")
    n_unmapped = sum(1 for l in open(log) if l.startswith("UNMAPPED:"))
    return out, n_unmapped, tmp


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--owltools", required=True, help="OWLTools executable (github.com/owlcollab/owltools)")
    ap.add_argument("--go-obo", required=True, help="go-basic.obo[.gz]")
    ap.add_argument("--gaf", required=True, help="goa_human.gaf[.gz]")
    ap.add_argument("--subset", default="goslim_generic")
    ap.add_argument("--outdir", default="libs")
    ap.add_argument("--suffix", default="2026", help="library name suffix (default: %(default)s)")
    ap.add_argument("--report", help="write a per-library build report TSV here")
    ap.add_argument("--java-home", help="JAVA_HOME to use if java is not on PATH")
    ap.add_argument("--memory", default="8G", help="OWLTools heap (default: %(default)s)")
    args = ap.parse_args()

    terms = parse_obo(args.go_obo)
    slim_ids = {t for t, r in terms.items() if args.subset in r["subsets"] and not r.get("obsolete")}
    if not slim_ids:
        sys.exit(f"Error: no terms tagged with subset '{args.subset}' in {args.go_obo}")

    slim_gaf, n_unmapped, tmp = run_map2slim(args.owltools, args.go_obo, args.gaf,
                                             args.subset, args.java_home, args.memory)

    cache = {}
    direct = defaultdict(lambda: defaultdict(set))       # aspect -> slim id -> genes
    n_rows = n_not = n_nd = 0
    with open(slim_gaf) as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < GAF_COLS - 2:
                continue
            n_rows += 1
            symbol, qual, go_id, evid, aspect = f[2], f[3], f[4], f[6], f[8]
            if "NOT" in qual.split("|"):
                n_not += 1
                continue
            if evid == "ND":
                n_nd += 1
                continue
            if go_id in slim_ids and aspect in ASPECT:
                # Upper-case (Enrichr convention, as in the GO_*/KEGG/Reactome
                # libraries): the GAF uses HGNC case (C9orf72), and gseapy
                # matches case-sensitively. See normalize_gmt_case.py.
                direct[aspect][go_id].add(symbol.upper())

    report = []
    os.makedirs(args.outdir, exist_ok=True)
    for aspect, label in ASPECT.items():
        prop = defaultdict(set)
        for sid, genes in direct[aspect].items():
            prop[sid] |= genes
            for anc in ancestors(terms, sid, cache) & slim_ids:
                prop[anc] |= genes
        out_path = os.path.join(args.outdir, f"GOslim_{label}_{args.suffix}.gmt")
        with open(out_path, "w") as out:
            for sid in sorted(prop):
                out.write(f"{terms[sid]['name']} ({sid})\t\t" + "\t".join(sorted(prop[sid])) + "\n")
        genes = set().union(*prop.values()) if prop else set()
        report.append({
            "library": os.path.basename(out_path), "slim_terms": len(prop),
            "genes": len(genes),
            "gene_term_pairs_most_specific": sum(len(v) for v in direct[aspect].values()),
            "gene_term_pairs_propagated": sum(len(v) for v in prop.values()),
        })
        print(f"wrote {out_path}: {len(prop)} slim terms, {len(genes)} genes")

    shutil.rmtree(tmp, ignore_errors=True)
    meta = {
        "go_obo": f"{os.path.basename(args.go_obo)} ({obo_version(args.go_obo)})",
        "gaf": f"{os.path.basename(args.gaf)} (date-generated {gaf_date(args.gaf)})",
        "subset": f"{args.subset} ({len(slim_ids)} terms)",
        "slim_annotations": n_rows, "removed_NOT": n_not, "removed_ND": n_nd,
        "annotations_without_slim_term": n_unmapped,
    }
    for k, val in meta.items():
        print(f"{k}\t{val}")
    if args.report:
        with open(args.report, "w") as rf:
            for k, val in meta.items():
                rf.write(f"# {k}\t{val}\n")
            cols = list(report[0])
            rf.write("\t".join(cols) + "\n")
            for r in report:
                rf.write("\t".join(str(r[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
