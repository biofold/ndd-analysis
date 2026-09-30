#!/usr/bin/env python3
"""
Upper-case the gene symbols of .gmt libraries (Enrichr convention).

The Enrichr-distributed libraries in libs/ (GO_*, KEGG, Reactome, SynGO) write
every gene symbol in upper case, including open-reading-frame genes (C9ORF72),
while HGNC -- and therefore libraries built from HGNC-cased sources such as
MONDO or the GO annotation file used for GO Slim -- writes C9orf72. gseapy
matches symbols case-sensitively, so mixed conventions across libraries make
results depend on the case of the input list. This tool brings a library to the
upper-case convention: only gene columns (3rd onward) are changed; term names
and descriptions are left as they are. Duplicates created within a term by
upper-casing are collapsed (and reported); with --sort the genes of each term
are written in sorted order, as utils/enrichment/build_goslim_gmt.py writes them.

Usage:
  python3 utils/enrichment/normalize_gmt_case.py libs/MONDO_2026.gmt [...]
  python3 utils/enrichment/normalize_gmt_case.py --check libs/*.gmt   # exit 1 if any library is not upper case
"""
import argparse
import os
import sys


def normalize(path, sort=False):
    """Return (new_lines, n_symbols_changed, n_collapsed, examples)."""
    out, changed, collapsed, examples = [], set(), 0, []
    # newline="" keeps the file's own line endings; empty fields (e.g. a
    # trailing tab) are kept in place so only gene symbols change.
    with open(path, newline="") as fh:
        for line in fh:
            body = line.rstrip("\r\n")
            eol = line[len(body):]
            f = body.split("\t")
            seen, new = set(), []
            for g in f[2:]:
                if not g:
                    new.append(g)
                    continue
                u = g.upper()
                if u != g:
                    if g not in changed and len(examples) < 5:
                        examples.append(f"{g}->{u}")
                    changed.add(g)
                if u in seen:
                    collapsed += 1
                    continue
                seen.add(u)
                new.append(u)
            if sort:
                new = sorted(g for g in new if g) + [g for g in new if not g]
            out.append("\t".join(f[:2] + new) + eol)
    return out, len(changed), collapsed, examples


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gmt", nargs="+")
    ap.add_argument("--check", action="store_true", help="report only; exit 1 if any symbol is not upper case")
    ap.add_argument("--sort", action="store_true", help="write the genes of each term in sorted order")
    args = ap.parse_args()
    bad = 0
    for p in args.gmt:
        lines, n, c, ex = normalize(p, sort=args.sort)
        name = os.path.basename(p)
        if args.check:
            status = "OK" if n == 0 else f"{n} symbols not upper case (e.g. {', '.join(ex)})"
            print(f"{name}\t{status}")
            bad += n > 0
            continue
        with open(p + ".tmp", "w", newline="") as fo:
            fo.writelines(lines)
        os.replace(p + ".tmp", p)
        print(f"{name}\t{n} symbols upper-cased\t{c} within-term duplicates collapsed"
              + (f"\te.g. {', '.join(ex)}" if ex else ""))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
