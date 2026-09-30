#!/usr/bin/env python3
"""
GRCh38 genomic coordinates for the 19,354-gene universe, from a pinned GENCODE GTF.

ASSIGNMENT RULE (never by symbol)
---------------------------------
Each universe gene gets its HGNC ID as in build_gene_id_resolution.py (approved symbol,
else a unique previous symbol). Its GENCODE `gene` record is then found by:
  1. ensembl_id   HGNC's Ensembl gene ID == GENCODE gene_id (version stripped)
  2. hgnc_attr    otherwise, GENCODE's own hgnc_id gene attribute == the gene's HGNC ID
Genes with neither are left without coordinates (their cytoband and MANE Select from HGNC
are kept). When both rules name a record and disagree, rule 1 is used and the conflict
is logged. Pseudoautosomal (PAR) genes carry a chrX record and a chrY copy with its own
Ensembl ID and the same hgnc_id attribute: the chrX record is used, the copy is noted.

Every exception is written to --log (severity warning = no coordinates or conflicting
rules; info = assigned by rule 2 or PAR copy) and summarised on stderr.

OUTPUT (--output, TSV, one row per universe gene, universe order)
------
Gene, hgnc_id, ensembl_gene_id (versioned, as in the GTF), chrom, start, end (1-based,
inclusive, as in GTF), strand, cytoband (HGNC location), mane_select (HGNC; Ensembl and
RefSeq transcript), assembly, annotation, coord_match (ensembl_id | hgnc_attr | empty),
coord_note.
"""
import argparse
import csv
import gzip
import re
import sys
from collections import Counter, defaultdict

import pandas as pd

sys.path.insert(0, __import__("os").path.dirname(__file__))
from build_gene_id_resolution import map_universe  # noqa: E402

ATTR = re.compile(r'(\S+) "([^"]*)"')
COLUMNS = ["Gene", "hgnc_id", "ensembl_gene_id", "chrom", "start", "end", "strand",
           "cytoband", "mane_select", "assembly", "annotation", "coord_match", "coord_note"]


def read_gtf_genes(path):
    """GENCODE gene records: {unversioned ENSG: [record, ...]} and {HGNC ID: [ENSG, ...]}."""
    by_ensg, by_hgnc = defaultdict(list), defaultdict(set)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t", 8)
            if len(f) < 9 or f[2] != "gene":
                continue
            a = dict(ATTR.findall(f[8]))
            gid = a["gene_id"]
            base = gid.split(".")[0]
            rec = {"gene_id": gid, "chrom": f[0], "start": int(f[3]), "end": int(f[4]),
                   "strand": f[6], "gene_name": a.get("gene_name", ""),
                   "gene_type": a.get("gene_type", ""), "par": gid.endswith("_PAR_Y")}
            by_ensg[base].append(rec)
            if a.get("hgnc_id"):
                by_hgnc[a["hgnc_id"]].add(base)
    return by_ensg, by_hgnc


def pick(records):
    """One record per ENSG: chrX over the chrY PAR copy."""
    main = [r for r in records if not r["par"] and r["chrom"] != "chrY"] or \
           [r for r in records if not r["par"]] or records
    return main[0], len(records) > 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--genes", required=True)
    ap.add_argument("--hgnc", required=True)
    ap.add_argument("--gtf", required=True, help="GENCODE primary-assembly or CHR annotation GTF (.gz)")
    ap.add_argument("--annotation", required=True, help="label, e.g. 'GENCODE v50'")
    ap.add_argument("--output", required=True)
    ap.add_argument("--log", required=True)
    args = ap.parse_args()

    genes = pd.read_csv(args.genes, sep="\t", dtype=str).iloc[:, 0].tolist()
    hgnc = pd.read_csv(args.hgnc, sep="\t", dtype=str, low_memory=False)
    hgnc = hgnc[hgnc["status"] == "Approved"]
    gene2hgnc, _ = map_universe(genes, hgnc)
    h = hgnc.set_index("hgnc_id")
    by_ensg, by_hgnc = read_gtf_genes(args.gtf)

    rows, log, stats = [], [], Counter()
    for g in genes:
        hid = gene2hgnc[g]
        hr = h.loc[hid]
        ens = hr["ensembl_gene_id"] if isinstance(hr["ensembl_gene_id"], str) else ""
        mane = str(hr["mane_select"]).strip('"').replace("|", " / ") if isinstance(hr["mane_select"], str) else ""
        row = dict.fromkeys(COLUMNS, "")
        row.update(Gene=g, hgnc_id=hid, cytoband=hr["location"] if isinstance(hr["location"], str) else "",
                   mane_select=mane, assembly="GRCh38", annotation=args.annotation)
        attr = sorted(by_hgnc.get(hid, ()))
        chosen, how, notes = None, "", []
        if ens and ens in by_ensg:
            chosen, how = ens, "ensembl_id"
            extra = [e for e in attr if e != ens]
            par_y = [e for e in extra if all(r["chrom"] == "chrY" for r in by_ensg.get(e, [{"chrom": ""}]))]
            if extra and par_y == extra and ens in attr and \
                    all(r["chrom"] == "chrX" for r in by_ensg[ens]):
                notes.append(f"pseudoautosomal: chrY copy {','.join(extra)}")
                log.append((g, "info", "par_gene", f"chrX record {ens} used; chrY copy {','.join(extra)}"))
            elif extra or (attr and ens not in attr):
                notes.append(f"GENCODE hgnc_id attribute names {','.join(attr)}")
                log.append((g, "warning", "rules_disagree", f"HGNC Ensembl {ens}; GENCODE hgnc_id attribute on {','.join(attr)}"))
        elif len(attr) == 1:
            chosen, how = attr[0], "hgnc_attr"
            notes.append(f"Ensembl gene ID in GENCODE is {attr[0]}" +
                         (f" (HGNC lists {ens}, absent from this release)" if ens else " (HGNC lists none)"))
            log.append((g, "info", "assigned_by_hgnc_attr",
                        f"HGNC Ensembl ID {ens or 'none'} not in GTF; GENCODE hgnc_id attribute on {attr[0]}"))
        else:
            why = (f"GENCODE hgnc_id attribute on several genes: {','.join(attr)}" if attr else
                   f"HGNC Ensembl ID {ens or 'none'} not in GTF and no GENCODE record carries {hid}")
            log.append((g, "warning", "no_coordinates", why))
            notes.append("no GENCODE gene record on the GRCh38 primary chromosomes")
        if chosen:
            rec, has_par = pick(by_ensg[chosen])
            row.update(ensembl_gene_id=rec["gene_id"], chrom=rec["chrom"], start=rec["start"],
                       end=rec["end"], strand=rec["strand"], coord_match=how)
            if has_par:
                notes.append("pseudoautosomal: also on chrY")
                log.append((g, "info", "par_gene", "chrX record used; chrY PAR copy not listed"))
        row["coord_note"] = "; ".join(notes)
        stats[how or "none"] += 1
        rows.append(row)

    assert len(rows) == len(genes)
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(args.log, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "severity", "category", "detail"])
        w.writerows(log)
    print(f"{len(genes)} genes: " + ", ".join(f"{k} {v}" for k, v in sorted(stats.items())), file=sys.stderr)
    for (sev, cat), n in sorted(Counter((s, c) for _, s, c, _ in log).items()):
        ex = ", ".join(g for g, s, c, _ in log if (s, c) == (sev, cat))[:120]
        print(f"  [{sev}] {cat}: {n} ({ex})", file=sys.stderr)
    print(f"wrote {args.output} and {args.log}", file=sys.stderr)


if __name__ == "__main__":
    main()
