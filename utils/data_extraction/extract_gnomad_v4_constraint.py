#!/usr/bin/env python
"""
Extract per-gene gnomAD v4.1 constraint metrics (pLI, LOEUF, o/e LoF) from
gnomAD's own primary release, for serving ALONGSIDE -- not in place of -- the
dbNSFP-derived gnomAD v2.1.1 pLI already in ndd_master_table.tsv.

WHY ALONGSIDE, NOT INSTEAD
--------------------------
The master table's gnomAD_pLI column is gnomAD v2.1.1 (via dbNSFP 5.2) and is
the value every published analysis in this project used -- the Figure 5 violin
panels, the KS p-value matrix, the MOE-tier comparisons. v4.1 recalculated pLI
on a ~6x larger cohort: only 28.9% of genes match v2.1.1 exactly, though
Spearman rho = 0.834 (see verify_gnomad_pli_source.py and
data/gnomad_pli_v41_vs_v211_comparison.tsv). Replacing the column would
silently re-baseline published figures, so v4 is added as separate fields and
the v2.1.1 column is left untouched.

LOEUF IS THE HEADLINE METRIC, NOT pLI
-------------------------------------
gnomAD's v4 release notes recommend LOEUF (the upper bound of the 90% CI on the
observed/expected LoF ratio) over pLI for ranking constraint. LOEUF is
continuous and lower = more constrained, the opposite direction to pLI. It is
column `lof.oe_ci.upper`. v4 pLI (`lof.pLI`) is carried too, so the v2.1.1 vs
v4.1 comparison can be made on like-for-like terms.

TRANSCRIPT DEDUPLICATION
------------------------
The source file has one row per gene x transcript (211,523 rows), including
both RefSeq- and Ensembl-named rows for the same transcript. Restricting to
mane_select == True gives one transcript per gene; MANE Select was introduced
in v4.1 specifically to remove the canonical-transcript ambiguity that
v2.1.1's file had (all 25 otherwise-unexplained v2.1.1 mismatches traced to
multi-transcript genes).

Residual duplicate MANE rows for one gene are the same transcript under both
naming schemes. verify_gnomad_pli_source.py's docstring records that these were
checked once to agree on pLI and then drops them with keep="first". That
agreement is ENFORCED here rather than documented: --max-dup-spread asserts the
within-gene spread on every emitted metric, so a future gnomAD release whose
duplicate rows disagree fails loudly instead of silently keeping whichever row
sorted first.

COVERAGE: chrX AND chrY ARE NOT IN THE RELEASE
---------------------------------------------
The v4.1 constraint file contains chr1-22 only. Of 19,354 iNDDx genes, 17,461
have a value; 1,893 do not, and 884 of those are X (840) or Y (44) -- including
176 of the 205 uncovered HIGH-CONFIDENCE genes (MECP2, CDKL5, FMR1, ARX, DCX,
ATRX, KDM5C, DMD, SLC6A8 ...). NDD genes are enriched on X, so absence here is
class-biased. An absent gene means "gnomAD v4.1 did not score it", NOT
"unconstrained". The v2.1.1 pLI in ndd_master_table.tsv does cover X, which is
a further reason v4 is served alongside rather than instead.

JOIN KEY: ENSEMBL GENE ID, NOT SYMBOL
-------------------------------------
Symbols drift between releases. Matching iNDDx to this file by upper-cased
symbol covers 17,214 genes; by Ensembl gene ID it covers 17,461 -- 252 genes are
recovered and 5 are not. Output is therefore keyed on `Ensembl_ID` (gnomAD's
`gene_id`), with the symbol carried for readability only.

CONSTRAINT FLAGS
----------------
gnomAD sets `constraint_flags` on genes whose estimates it considers
unreliable (too few expected LoF variants, outliers, etc.). The flag is emitted
verbatim rather than used to drop rows, so a consumer can filter on it; a
constraint value on a flagged gene should not be read as comparable to an
unflagged one.

Usage
-----
python utils/data_extraction/extract_gnomad_v4_constraint.py \
    --input data/raw/gnomad.v4.1.constraint_metrics.tsv.gz \
    --output data/gene_gnomad_v4_constraint.tsv
"""
import argparse
import json
import sys

import pandas as pd

# Source columns. gnomAD's names are dotted; keep the mapping explicit so a
# renamed upstream column fails on a KeyError here rather than silently
# producing an all-NaN output field.
COL_GENE = "gene"
COL_GENE_ID = "gene_id"
COL_TRANSCRIPT = "transcript"
COL_MANE = "mane_select"
COL_FLAGS = "constraint_flags"
METRICS = {
    "lof.pLI": "gnomad_v4_pli",
    "lof.oe_ci.upper": "gnomad_v4_loeuf",
    "lof.oe": "gnomad_v4_oe_lof",
}


def normalize_mane(series):
    """mane_select arrives as bool or as a 'true'/'false'/'' string depending on
    how the release was exported; normalize both to a real boolean."""
    if series.dtype == bool:
        return series
    return (series.astype(str).str.strip().str.lower()
            .isin(["true", "t", "1", "yes"]))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True,
                    help="gnomAD v4.1 constraint_metrics TSV (.gz accepted)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--max-dup-spread", type=float, default=1e-6,
                    help="maximum tolerated within-gene spread across duplicate "
                         "MANE-Select rows, per metric (default 1e-6)")
    args = ap.parse_args()

    print(f"Reading {args.input} ...")
    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    print(f"  {len(df):,} rows x {len(df.columns)} columns")

    missing = [c for c in [COL_GENE, COL_TRANSCRIPT, COL_MANE, COL_FLAGS] + list(METRICS)
               if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: expected columns absent from input: {missing}")

    missing = [c for c in [COL_GENE_ID] if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: expected columns absent from input: {missing}")

    mane = df[normalize_mane(df[COL_MANE])].copy()
    print(f"  MANE-Select rows: {len(mane):,}")

    # gnomAD lists every MANE transcript twice: once under its Ensembl name
    # (ENST..., gene_id = ENSG...) and once under its RefSeq name (NM_..., gene_id
    # = NCBI Entrez number). Only the ENSG rows can be joined to iNDDx, so those are
    # the emitted rows. The RefSeq twins are not discarded on trust, though: first
    # verify they agree with their Ensembl twin, per gene symbol, on every metric.
    is_ens = mane[COL_GENE_ID].astype("string").str.startswith("ENSG").fillna(False)
    ens, ref = mane[is_ens].copy(), mane[~is_ens].copy()
    print(f"  Ensembl-named rows (ENSG): {len(ens):,} | RefSeq-named rows: {len(ref):,}")
    if ens[COL_GENE_ID].isna().any():
        sys.exit("ERROR: MANE row without a gene_id")

    both = mane[mane[COL_GENE].notna()]
    dup_syms = both[COL_GENE][both[COL_GENE].duplicated(keep=False)].unique()
    print(f"  gene symbols with >1 MANE-Select row: {len(dup_syms):,}")
    g = both[both[COL_GENE].isin(dup_syms)]
    offenders = {}
    for src in METRICS:
        spread = g.groupby(COL_GENE)[src].agg(lambda s: s.max() - s.min())
        bad = spread[spread > args.max_dup_spread].dropna()
        mixed = g.groupby(COL_GENE)[src].agg(lambda s: s.isna().any() and s.notna().any())
        if mixed.any():
            bad = pd.concat([bad, pd.Series(float("inf"), index=mixed[mixed].index)])
        if len(bad):
            offenders[src] = bad
    if offenders:
        for src, bad in offenders.items():
            print(f"\nERROR: {len(bad)} gene(s) have disagreeing MANE-Select rows "
                  f"for {src}:", file=sys.stderr)
            print(bad.head(10).to_string(), file=sys.stderr)
        sys.exit("Refusing to emit: RefSeq- and Ensembl-named MANE rows disagree, so "
                 "dropping one set would be an arbitrary choice. Inspect the release.")
    print(f"  RefSeq and Ensembl twins agree within {args.max_dup_spread:g} "
          f"(and on missingness) for all {len(METRICS)} metrics")

    # Now one row per Ensembl gene, which must already hold.
    if ens[COL_GENE_ID].duplicated().any():
        n = ens[COL_GENE_ID].duplicated().sum()
        print(f"  NOTE: {n} Ensembl gene IDs have >1 MANE row; keeping the first "
              f"after sorting by transcript")
    mane = ens

    out = (mane.sort_values([COL_GENE_ID, COL_TRANSCRIPT])
               .drop_duplicates(subset=COL_GENE_ID, keep="first")
               .rename(columns={COL_GENE_ID: "Ensembl_ID", COL_GENE: "Gene",
                                COL_TRANSCRIPT: "mane_transcript",
                                COL_FLAGS: "gnomad_v4_constraint_flags", **METRICS}))
    out["Gene"] = out["Gene"].astype("string").str.upper()

    # constraint_flags arrives as a JSON list string; "[]" means NO flags. Parse it
    # so "flagged" is countable and the served field is a real list downstream.
    def parse_flags(v):
        if pd.isna(v):
            return ""
        try:
            return ";".join(json.loads(v))
        except (TypeError, ValueError):
            sys.exit(f"ERROR: unparseable constraint_flags value: {v!r}")
    out["gnomad_v4_constraint_flags"] = out["gnomad_v4_constraint_flags"].map(parse_flags)

    out = out[["Ensembl_ID", "Gene", "mane_transcript", "gnomad_v4_pli",
               "gnomad_v4_loeuf", "gnomad_v4_oe_lof", "gnomad_v4_constraint_flags"]]
    out = out.sort_values("Ensembl_ID")

    assert out["Ensembl_ID"].is_unique, "Ensembl_ID not unique after deduplication"
    n_flagged = (out["gnomad_v4_constraint_flags"] != "").sum()
    print(f"\n  emitting {len(out):,} genes (one per Ensembl gene ID)")
    for c in ("gnomad_v4_pli", "gnomad_v4_loeuf", "gnomad_v4_oe_lof"):
        print(f"    {c:22} non-null {out[c].notna().sum():,}")
    print(f"    with >=1 constraint flag {n_flagged:,} (emitted, NOT dropped)")

    out.to_csv(args.output, sep="\t", index=False)
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
