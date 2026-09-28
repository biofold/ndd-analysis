"""
Build GOslim_{Biological_Process,Cellular_Component,Molecular_Function}_2026.gmt
by mapping the full GO annotations already used for the MOE enrichment pipeline
(GO_Biological_Process_2026.gmt, GO_Cellular_Component_2026.gmt,
GO_Molecular_Function_2026.gmt -- Enrichr-format libraries, one row per full GO
term with its associated gene symbols) down to the generic GO Slim term set.

Source files:
  - go-basic.obo         : full GO DAG (fetched from current.geneontology.org)
  - goslim_generic.obo   : generic GO Slim subset definition (same source)
  - GO_{BP,CC,MF}_2026.gmt : existing full-GO gene-set libraries (already built
                             for the 19,354-gene universe from the MOE pipeline)

For each full GO term T that a gene is annotated to, goatools' `mapslim` walks
T's ancestors in the full DAG and returns the nearest enclosing slim term(s).
Genes are then re-aggregated under each slim term.
"""
import sys
import re
from collections import defaultdict

from goatools.obo_parser import GODag
from goatools.mapslim import mapslim

NS_MAP = {
    "biological_process": "Biological_Process",
    "cellular_component": "Cellular_Component",
    "molecular_function": "Molecular_Function",
}

GMT_SOURCES = {
    "biological_process": "../ndd-analysis/libs/GO_Biological_Process_2026.gmt",
    "cellular_component": "../ndd-analysis/libs/GO_Cellular_Component_2026.gmt",
    "molecular_function": "../ndd-analysis/libs/GO_Molecular_Function_2026.gmt",
}

TERM_RE = re.compile(r"^(?P<name>.+) \((?P<go_id>GO:\d+)\)$")


def parse_gmt(path):
    """Return dict: go_id -> {"name": term_name, "genes": set(genes)}"""
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term_field = parts[0]
            genes = set(g for g in parts[2:] if g)
            m = TERM_RE.match(term_field)
            if not m:
                print(f"WARNING: could not parse term field: {term_field!r}", file=sys.stderr)
                continue
            go_id = m.group("go_id")
            name = m.group("name")
            out[go_id] = {"name": name, "genes": genes}
    return out


def main():
    print("Loading full GO DAG (go-basic.obo) ...")
    go_dag = GODag("go-basic.obo")

    print("Loading GO Slim subset (goslim_generic.obo) ...")
    goslim_dag = GODag("goslim_generic.obo")

    for ns_key, ns_label in NS_MAP.items():
        gmt_path = GMT_SOURCES[ns_key]
        print(f"\n=== {ns_label} ===")
        print(f"Parsing existing full-GO library: {gmt_path}")
        full_terms = parse_gmt(gmt_path)
        print(f"  {len(full_terms)} full GO terms, "
              f"{sum(len(v['genes']) for v in full_terms.values())} term-gene rows")

        # go_id -> set of gene symbols annotated (from the existing library)
        slim_to_genes = defaultdict(set)
        slim_id_to_name = {}
        unmapped_terms = 0
        unknown_go_ids = 0

        for go_id, rec in full_terms.items():
            if go_id not in go_dag:
                unknown_go_ids += 1
                continue
            try:
                direct_anc, all_anc = mapslim(go_id, go_dag, goslim_dag)
            except Exception as e:
                unmapped_terms += 1
                continue
            # direct_anc = nearest enclosing slim term(s); fall back to all_anc if empty
            targets = direct_anc if direct_anc else all_anc
            if not targets:
                unmapped_terms += 1
                continue
            for slim_id in targets:
                slim_to_genes[slim_id] |= rec["genes"]
                if slim_id in goslim_dag:
                    slim_id_to_name[slim_id] = goslim_dag[slim_id].name

        print(f"  {unknown_go_ids} full GO ids not found in go-basic.obo (skipped)")
        print(f"  {unmapped_terms} full GO ids with no slim mapping found (skipped)")
        print(f"  -> {len(slim_to_genes)} slim terms populated, "
              f"{sum(len(v) for v in slim_to_genes.values())} slim term-gene rows")

        out_path = f"GOslim_{ns_label}_2026.gmt"
        n_rows = 0
        with open(out_path, "w") as out:
            for slim_id, genes in sorted(slim_to_genes.items()):
                name = slim_id_to_name.get(slim_id, slim_id)
                genes_sorted = sorted(genes)
                out.write(f"{name} ({slim_id})\t\t" + "\t".join(genes_sorted) + "\n")
                n_rows += 1
        print(f"  wrote {out_path}: {n_rows} slim terms")


if __name__ == "__main__":
    main()
