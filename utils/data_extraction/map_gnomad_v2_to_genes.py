#!/usr/bin/env python3
"""
Attach gnomAD v2.1.1 constraint (pLI, LOEUF) to the 19,354-gene universe.

WHY THIS SCRIPT EXISTS
----------------------
The MOE validation analyses (R1.3, R2.10) used pLI from dbNSFP 5.2, which joined
gnomAD to genes by SYMBOL. Against gnomAD's own v2.1.1 file that left 1,385 scored
genes empty (almost all HGNC renames after gnomAD's 2018 naming, e.g. AARS1/AARS)
and gave 18 genes another gene's value (TUBB3, PI4K2A). dbNSFP also carries no
LOEUF. This script joins gnomAD's own table (data/gene_gnomad_v2_constraint.tsv,
from extract_gnomad_constraint.py) to the gene universe.

ASSIGNMENT RULE
---------------
A gnomAD constraint value is computed on ONE transcript, so it describes the gene
that owns that transcript. gnomAD v2.1.1 labels rows with 2018 (GENCODE v19) gene
IDs and symbols, some of which HGNC has since moved to other genes. Rows are
therefore assigned in two passes, one gnomAD row per gene and never two:

  1. Base: HGNC Ensembl gene ID == gnomAD gene_id; otherwise an unambiguous gnomAD
     symbol (98 v2.1.1 symbols name more than one gene_id and are never used);
     otherwise a PREVIOUS HGNC symbol of the gene, only when that previous symbol
     belongs to exactly one approved gene, is not the approved symbol of another
     gene, names exactly one gnomAD row, and that row is still unassigned
     (warning: assigned_by_previous_symbol).
  2. Transcript check. The current owner of each row's transcript is read from the
     gnomAD v4.1 constraint file (every GENCODE v39 transcript with its gene_id),
     which is pinned in data/raw, so no network lookup happens at run time.
       - owner is another gene of the universe WITHOUT a row of its own
         -> the row moves to the owner                 (warning: reassigned_by_transcript)
       - owner is another gene that already has its own row
         -> keep the base assignment                   (warning: transcript_owner_has_own_row)
       - transcript retired (absent from v4.1) and the gnomAD symbol names a different,
         still unmatched gene -> the owner cannot be decided, the value is dropped
                                                       (warning: unresolved_retired_transcript)
       - owner == holder, but the gnomAD symbol names a different unmatched gene
         -> keep, the transcript settles it            (info: symbol_conflict_resolved_by_transcript)

Every exception is written to --exceptions and summarised on stderr, and each gene
carries it in the gnomad_v2_warning column.

RENAMING CHECK
--------------
gnomAD v2.1.1 uses 2018 gene symbols. For every gene the symbol of its gnomAD row is
compared with the current HGNC symbol and, when they differ, classified with HGNC's
prev_symbol / alias_symbol fields (written to --renames, warned on stderr):
  renamed_since_gnomAD        gnomAD symbol is a previous HGNC symbol of the gene
  alias_in_gnomAD             gnomAD symbol is an HGNC alias of the gene
  symbol_differs_unverified   neither (e.g. clone names such as RP11-...)
  gnomAD_symbol_now_other_gene  gnomAD symbol is today the approved symbol of a
                              DIFFERENT gene -- the case most likely to mislead a
                              symbol-based lookup (e.g. the ADORA3 row held by TMIGD3)
  possible_missed_match       gene has no row, but one of its previous symbols names a
                              gnomAD row left unassigned

All exceptions and renaming checks are also written to ONE log (--log) with a
severity column (warning = value moved, dropped, or assigned by a weaker rule;
info = value unchanged, provenance only), so every run leaves a complete audit trail.

Output columns: Gene, Ensembl_ID, gnomad_v2_match (ensembl|symbol|previous_symbol|transcript),
gnomad_v2_transcript, gnomad_v2_warning, gnomAD_pLI, gnomAD_LOEUF,
gnomAD_LOEUF_decile, gnomAD_oe_lof, gnomAD_constraint_flag, gnomad_v2_warning_severity
(warning|info, the most severe of the gene's notes), gnomad_v2_note_text (plain language,
from NOTE_TEXT; served by iNDDx).

Usage:
  python3 utils/data_extraction/map_gnomad_v2_to_genes.py \\
      --genes data/gene_all_score.txt \\
      --hgnc data/raw/hgnc_complete_set_20260930.txt.gz \\
      --gnomad data/gene_gnomad_v2_constraint.tsv \\
      --transcript-owners data/raw/gnomad.v4.1.constraint_metrics.tsv.gz \\
      --output data/gene_gnomad_v2_mapped.tsv \\
      --exceptions data/gnomad_v2_mapping_exceptions.tsv \\
      --renames data/gnomad_v2_renamed_genes.tsv \\
      --log data/gnomad_v2_mapping_log.tsv
"""

import argparse
import sys
import warnings

import pandas as pd

# Severity of each exception category (warning = value moved, dropped, or assigned by a
# weaker rule; info = documented but not in doubt). Unknown categories count as warning.
SEVERITY = {"reassigned_by_transcript": "warning", "transcript_owner_has_own_row": "warning",
            "unresolved_retired_transcript": "warning", "assigned_by_previous_symbol": "warning",
            "symbol_conflict_resolved_by_transcript": "info",
            "gnomAD_symbol_now_other_gene": "warning", "possible_missed_match": "warning",
            "symbol_differs_unverified": "warning", "alias_in_gnomAD": "info",
            "renamed_since_gnomAD": "info"}


# Plain-language text for each category (shown by iNDDx); {x} is the detail after
# the colon, with the "from_", "vs_" and "row_kept_by_" prefixes removed. The
# definitions are the ones in the module docstring above.
NOTE_TEXT = {
    "renamed_since_gnomAD": "gnomAD v2.1.1 lists this gene under its previous HGNC symbol {x}.",
    "alias_in_gnomAD": "gnomAD v2.1.1 lists this gene under {x}, an HGNC alias of the gene.",
    "symbol_differs_unverified": "gnomAD v2.1.1 lists this gene as {x}, which is neither a previous "
                                 "symbol nor an alias in HGNC.",
    "gnomAD_symbol_now_other_gene": "gnomAD v2.1.1 lists this gene as {x}, today the approved symbol "
                                    "of a different gene.",
    "assigned_by_previous_symbol": "No gnomAD row matches this gene's Ensembl ID or symbol; the value is "
                                   "taken from the row of its previous HGNC symbol {x}.",
    "reassigned_by_transcript": "The value comes from the gnomAD row first assigned to {x}: that row's "
                                "transcript now belongs to this gene.",
    "lost_row_to_transcript_owner": "No gnomAD value: the row first assigned to this gene was computed on "
                                    "a transcript that now belongs to {x}, which receives it.",
    "transcript_owner_has_own_row": "The transcript of this gene's gnomAD row now belongs to {x}, which "
                                    "has its own row; the value is kept for this gene.",
    "unresolved_retired_transcript": "No gnomAD value: the row was computed on a retired transcript and "
                                     "is also claimed by {x}, so its owner cannot be decided.",
    "symbol_conflict_resolved_by_transcript": "No gnomAD value: the row labelled with this gene's symbol "
                                              "was computed on a transcript of {x}, which keeps it.",
    "possible_missed_match": "No gnomAD value, although a previous symbol names an unassigned gnomAD "
                             "row ({x}).",
}


def note_text(note):
    """Plain-language rendering of a gene's ';'-joined notes (unknown categories: the code)."""
    if not isinstance(note, str) or not note:
        return pd.NA
    out = []
    for c in (c for c in note.split(";") if c):
        cat, _, x = c.partition(":")
        for pre in ("from_", "vs_", "row_kept_by_"):
            if x.startswith(pre):
                x = x[len(pre):]
        t = NOTE_TEXT.get(cat)
        out.append(t.format(x=x) if t and x else (t.replace(" {x}", "").replace("({x})", "") if t else c))
    return " ".join(out)


def note_severity(note):
    """Most severe category among a gene's ';'-joined 'category:detail' notes."""
    if not isinstance(note, str) or not note:
        return pd.NA
    cats = [c.split(":", 1)[0] for c in note.split(";") if c]
    return "warning" if any(SEVERITY.get(c, "warning") == "warning" for c in cats) else "info"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", required=True, help="gene universe (first column = HGNC symbol)")
    ap.add_argument("--hgnc", required=True, help="HGNC hgnc_complete_set.txt[.gz]")
    ap.add_argument("--gnomad", required=True, help="data/gene_gnomad_v2_constraint.tsv")
    ap.add_argument("--transcript-owners", required=True,
                    help="gnomAD v4.1 constraint_metrics.tsv[.gz] (transcript -> current gene_id)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--exceptions", required=True, help="TSV listing every exception")
    ap.add_argument("--renames", required=True, help="TSV listing genes renamed since gnomAD v2.1.1")
    ap.add_argument("--log", required=True, help="single TSV log of every exception (with severity)")
    args = ap.parse_args()

    genes = pd.read_csv(args.genes, sep="\t", dtype=str)
    genes = genes.rename(columns={genes.columns[0]: "Gene"})[["Gene"]]
    genes["key"] = genes["Gene"].str.upper()

    hgnc = pd.read_csv(args.hgnc, sep="\t", dtype=str, low_memory=False)
    hgnc = hgnc[hgnc["status"] == "Approved"][["symbol", "ensembl_gene_id",
                                                  "prev_symbol", "alias_symbol"]]
    hgnc["key"] = hgnc["symbol"].str.upper()
    hgnc = hgnc.drop_duplicates("key")

    gn = pd.read_csv(args.gnomad, sep="\t", dtype=str)
    gn["key"] = gn["Gene"].str.upper()
    gn["tx"] = gn["transcript"].str.split(".").str[0]
    value_cols = ["gnomad_v2_pli", "gnomad_v2_loeuf", "gnomad_v2_loeuf_decile",
                  "gnomad_v2_oe_lof", "gnomad_v2_constraint_flag"]
    by_ens = gn.drop_duplicates("Ensembl_ID").set_index("Ensembl_ID")

    v4 = pd.read_csv(args.transcript_owners, sep="\t", dtype=str,
                     usecols=["gene", "gene_id", "transcript"])
    v4 = v4[v4["transcript"].str.startswith("ENST", na=False)].copy()
    v4["tx"] = v4["transcript"].str.split(".").str[0]
    if v4.groupby("tx")["gene_id"].nunique().gt(1).any():
        sys.exit("Error: a transcript has more than one owner in the transcript-owner file")
    tx_owner = v4.drop_duplicates("tx").set_index("tx")["gene_id"]

    # ---- pass 1: base assignment -------------------------------------------------
    g = genes.merge(hgnc[["key", "ensembl_gene_id"]], on="key", how="left")
    g["Ensembl_ID"] = g["ensembl_gene_id"]
    g["gnomad_gene_id"] = pd.NA
    g["gnomad_v2_match"] = pd.NA
    hit = g["Ensembl_ID"].isin(by_ens.index)
    g.loc[hit, "gnomad_gene_id"] = g.loc[hit, "Ensembl_ID"]
    g.loc[hit, "gnomad_v2_match"] = "ensembl"

    sym_counts = gn.groupby("key")["Ensembl_ID"].nunique()
    unambiguous = gn[gn["key"].map(sym_counts) == 1].drop_duplicates("key").set_index("key")
    claimed = set(g["gnomad_gene_id"].dropna())
    for i in g.index[g["gnomad_gene_id"].isna() & g["key"].isin(unambiguous.index)]:
        gid = unambiguous.loc[g.at[i, "key"], "Ensembl_ID"]
        if gid not in claimed:
            g.at[i, "gnomad_gene_id"] = gid
            g.at[i, "gnomad_v2_match"] = "symbol"
            claimed.add(gid)

    # previous HGNC symbol (renamed genes whose Ensembl ID also changed)
    split = lambda v: set() if pd.isna(v) else {x.strip().upper() for x in str(v).split("|")}
    prev = {k: split(v) for k, v in hgnc.set_index("key")["prev_symbol"].items()}
    alias = {k: split(v) for k, v in hgnc.set_index("key")["alias_symbol"].items()}
    approved = set(hgnc["key"])
    prev_owners = {}
    for k, ps in prev.items():
        for p_ in ps:
            prev_owners.setdefault(p_, set()).add(k)
    prev_exc = []
    for i in g.index[g["gnomad_gene_id"].isna()]:
        k = g.at[i, "key"]
        cands = [p_ for p_ in sorted(prev.get(k, set()))
                 if p_ in unambiguous.index and p_ not in approved
                 and len(prev_owners.get(p_, ())) == 1
                 and unambiguous.at[p_, "Ensembl_ID"] not in claimed]
        if len(cands) == 1:
            gid = unambiguous.at[cands[0], "Ensembl_ID"]
            g.at[i, "gnomad_gene_id"] = gid
            g.at[i, "gnomad_v2_match"] = "previous_symbol"
            claimed.add(gid)
            prev_exc.append((g.at[i, "Gene"], gid))

    # ---- pass 2: transcript check ------------------------------------------------
    g = g.set_index("Gene", drop=False)
    ens2gene = g.dropna(subset=["Ensembl_ID"]).drop_duplicates("Ensembl_ID")\
                .set_index("Ensembl_ID")["Gene"].to_dict()
    key2gene = g.set_index("key")["Gene"].to_dict()
    holder = {r: gene for gene, r in g["gnomad_gene_id"].dropna().items()}
    g["gnomad_v2_warning"] = pd.NA
    exc = []

    def note(gene, code):
        if gene is None:
            return
        prev = g.at[gene, "gnomad_v2_warning"]
        g.at[gene, "gnomad_v2_warning"] = code if pd.isna(prev) else f"{prev};{code}"

    for gene, rid in prev_exc:
        row = by_ens.loc[rid]
        note(gene, f"assigned_by_previous_symbol:{row['Gene']}")
        exc.append(dict(exception="assigned_by_previous_symbol", assigned_to=gene,
                        previous_holder=None, gnomad_row_gene_id=rid, gnomad_row_symbol=row["Gene"],
                        gnomad_row_transcript=row["transcript"],
                        transcript_owner_gene_id=tx_owner.get(row["tx"]),
                        transcript_owner=ens2gene.get(tx_owner.get(row["tx"])),
                        pLI=row["gnomad_v2_pli"], LOEUF=row["gnomad_v2_loeuf"]))

    for rid, row in by_ens.iterrows():
        owner_id = tx_owner.get(row["tx"])
        owner = ens2gene.get(owner_id) if owner_id is not None else None
        cur = holder.get(rid)
        sym_gene = key2gene.get(row["key"])
        base = dict(gnomad_row_gene_id=rid, gnomad_row_symbol=row["Gene"],
                    gnomad_row_transcript=row["transcript"], transcript_owner_gene_id=owner_id,
                    transcript_owner=owner, pLI=row["gnomad_v2_pli"], LOEUF=row["gnomad_v2_loeuf"])
        if owner is not None and owner != cur:
            if pd.isna(g.at[owner, "gnomad_gene_id"]):
                if cur is not None:
                    g.at[cur, "gnomad_gene_id"] = pd.NA
                    g.at[cur, "gnomad_v2_match"] = pd.NA
                    note(cur, f"lost_row_to_transcript_owner:{owner}")
                g.at[owner, "gnomad_gene_id"] = rid
                g.at[owner, "gnomad_v2_match"] = "transcript"
                holder[rid] = owner
                note(owner, "reassigned_by_transcript" + (f":from_{cur}" if cur else ""))
                exc.append(dict(exception="reassigned_by_transcript", assigned_to=owner,
                                previous_holder=cur, **base))
            elif cur is not None:
                note(cur, f"transcript_owner_has_own_row:{owner}")
                exc.append(dict(exception="transcript_owner_has_own_row", assigned_to=cur,
                                previous_holder=cur, **base))
        elif owner_id is None and cur is not None and sym_gene not in (None, cur) \
                and pd.isna(g.at[sym_gene, "gnomad_gene_id"]):
            g.at[cur, "gnomad_gene_id"] = pd.NA
            g.at[cur, "gnomad_v2_match"] = pd.NA
            holder.pop(rid, None)
            note(cur, f"unresolved_retired_transcript:vs_{sym_gene}")
            note(sym_gene, f"unresolved_retired_transcript:vs_{cur}")
            exc.append(dict(exception="unresolved_retired_transcript", assigned_to=None,
                            previous_holder=cur, symbol_gene=sym_gene, **base))
        elif owner is not None and owner == cur and sym_gene not in (None, cur) \
                and pd.isna(g.at[sym_gene, "gnomad_gene_id"]):
            note(sym_gene, f"symbol_conflict_resolved_by_transcript:row_kept_by_{cur}")
            exc.append(dict(exception="symbol_conflict_resolved_by_transcript", assigned_to=cur,
                            previous_holder=cur, symbol_gene=sym_gene, **base))

    # ---- renaming check ----------------------------------------------------------
    ren = []
    for gene, row in g.iterrows():
        rid = row["gnomad_gene_id"]
        k = row["key"]
        if pd.notna(rid):
            gsym = by_ens.at[rid, "Gene"]
            gk = gsym.upper()
            if gk == k:
                continue
            if gk in approved:
                code = "gnomAD_symbol_now_other_gene"
            elif gk in prev.get(k, set()):
                code = "renamed_since_gnomAD"
            elif gk in alias.get(k, set()):
                code = "alias_in_gnomAD"
            else:
                code = "symbol_differs_unverified"
            ren.append(dict(check=code, gene=gene, gnomad_symbol=gsym, gnomad_row_gene_id=rid,
                            match=row["gnomad_v2_match"]))
            note(gene, f"{code}:{gsym}")
        else:
            free = [p for p in prev.get(k, set()) if p in unambiguous.index
                    and unambiguous.at[p, "Ensembl_ID"] not in holder]
            for p in free:
                ren.append(dict(check="possible_missed_match", gene=gene,
                                gnomad_symbol=unambiguous.at[p, "Gene"],
                                gnomad_row_gene_id=unambiguous.at[p, "Ensembl_ID"], match=None))
                note(gene, f"possible_missed_match:{unambiguous.at[p, 'Gene']}")

    g = g.reset_index(drop=True)
    dup = g["gnomad_gene_id"].dropna().duplicated().sum()
    if dup:
        sys.exit(f"Error: {dup} gnomAD rows assigned to more than one gene")
    g = g.join(by_ens[value_cols + ["transcript"]], on="gnomad_gene_id")

    out = g.rename(columns={
        "gnomad_v2_pli": "gnomAD_pLI", "gnomad_v2_loeuf": "gnomAD_LOEUF",
        "gnomad_v2_loeuf_decile": "gnomAD_LOEUF_decile", "gnomad_v2_oe_lof": "gnomAD_oe_lof",
        "gnomad_v2_constraint_flag": "gnomAD_constraint_flag", "transcript": "gnomad_v2_transcript",
    })[["Gene", "Ensembl_ID", "gnomad_v2_match", "gnomad_v2_transcript", "gnomad_v2_warning",
        "gnomAD_pLI", "gnomAD_LOEUF", "gnomAD_LOEUF_decile", "gnomAD_oe_lof",
        "gnomAD_constraint_flag"]]
    # appended last so readers selecting columns by position are unaffected
    out = out.assign(gnomad_v2_warning_severity=out["gnomad_v2_warning"].map(note_severity),
                     gnomad_v2_note_text=out["gnomad_v2_warning"].map(note_text))
    out.to_csv(args.output, sep="\t", index=False)

    ex = pd.DataFrame(exc, columns=["exception", "assigned_to", "previous_holder", "symbol_gene",
                                    "gnomad_row_gene_id", "gnomad_row_symbol",
                                    "gnomad_row_transcript", "transcript_owner",
                                    "transcript_owner_gene_id", "pLI", "LOEUF"])
    ex.to_csv(args.exceptions, sep="\t", index=False)
    rn = pd.DataFrame(ren, columns=["check", "gene", "gnomad_symbol", "gnomad_row_gene_id", "match"])
    rn.to_csv(args.renames, sep="\t", index=False)

    severity = SEVERITY
    vals = by_ens[["gnomad_v2_pli", "gnomad_v2_loeuf", "transcript"]]
    log = pd.concat([
        ex.assign(category=ex["exception"], gene=ex["assigned_to"].fillna(ex["previous_holder"]),
                  other_gene=ex["symbol_gene"].fillna(ex["previous_holder"]).where(
                      ex["symbol_gene"].notna() | (ex["previous_holder"] != ex["assigned_to"])),
                  source="assignment"),
        rn.rename(columns={"gnomad_symbol": "gnomad_row_symbol"}).assign(
            category=rn["check"], other_gene=pd.NA, source="renaming")
          .join(vals, on="gnomad_row_gene_id")
          .rename(columns={"gnomad_v2_pli": "pLI", "gnomad_v2_loeuf": "LOEUF",
                           "transcript": "gnomad_row_transcript"}),
    ], ignore_index=True)
    log["severity"] = log["category"].map(severity).fillna("warning")
    log = log[["severity", "category", "source", "gene", "other_gene", "gnomad_row_symbol",
               "gnomad_row_gene_id", "gnomad_row_transcript", "transcript_owner", "pLI", "LOEUF"]]
    log = log.sort_values(["severity", "category", "gene"], ascending=[False, True, True])
    log.to_csv(args.log, sep="\t", index=False)

    m = out["gnomad_v2_match"]
    print(f"genes: {len(out):,}")
    print(f"matched to a gnomAD v2.1.1 row: {m.notna().sum():,} (ensembl {int((m == 'ensembl').sum()):,}, "
          f"symbol {int((m == 'symbol').sum()):,}, previous_symbol {int((m == 'previous_symbol').sum()):,}, "
          f"transcript {int((m == 'transcript').sum()):,})")
    print(f"with pLI and LOEUF: {int(out['gnomAD_pLI'].notna().sum()):,}")
    for code, grp in ex.groupby("exception"):
        genes_hit = sorted(set(grp["assigned_to"].dropna()) | set(grp["previous_holder"].dropna())
                           | set(grp["symbol_gene"].dropna()))
        warnings.warn(f"gnomAD v2.1.1 mapping exception '{code}': {len(grp)} row(s); genes: "
                      f"{', '.join(genes_hit)}", stacklevel=1)
    for code, grp in rn.groupby("check"):
        ex_list = ", ".join(f"{a}<-{b}" for a, b in zip(grp["gene"].head(8), grp["gnomad_symbol"].head(8)))
        warnings.warn(f"possible gene renaming '{code}': {len(grp)} gene(s), e.g. {ex_list}"
                      f"{' ...' if len(grp) > 8 else ''} (full list: {args.renames})", stacklevel=1)
    print(f"written: {args.output}\nexceptions ({len(ex)}): {args.exceptions}\n"
          f"renaming checks ({len(rn)}): {args.renames}\n"
          f"log ({len(log)} rows; {int((log['severity'] == 'warning').sum())} warnings, "
          f"{int((log['severity'] == 'info').sum())} info): {args.log}")


if __name__ == "__main__":
    warnings.simplefilter("always")
    warnings.formatwarning = lambda msg, cat, *a, **k: f"WARNING: {msg}\n"
    main()
