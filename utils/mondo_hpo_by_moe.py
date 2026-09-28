#!/usr/bin/env python3
"""
Count candidate genes carrying at least one MONDO disease term and at least
one HPO phenotype term, broken down by MOE score tier, and draw the grouped
bar chart used in the response to reviewer point R2.10.

WHY BOTH ONTOLOGIES
-------------------
The manuscript already reports MONDO disease-term coverage across MOE tiers
(scripts/3_count_supercandidate.py). Reviewer R2.10 asks whether HPO could
contribute to the MOE score or to its validation, so the natural first thing
to show is how HPO phenotype coverage behaves across the same tiers and how it
compares with MONDO. Plotting the two against the tier totals makes the
comparison direct: if HPO tracked MONDO exactly it would add nothing, and if
it covered a different slice of genes it is worth considering as a component.

Note the interpretive caveat that belongs with this figure: both MONDO and HPO
gene annotations derive largely from genes already linked to disease, so
coverage rising with MOE is consistent with the score working, but neither
ontology is independent enough of curated disease knowledge to serve as
validation on its own. That is why the independent validation in step 6 uses
gnomAD constraint and ClinVar variant data instead.

Restricted to the candidate class by default, matching every other MOE-tier
analysis in this pipeline (the score is only defined for candidates).

Usage:
  python3 utils/mondo_hpo_by_moe.py results/main/ndd_master_table.tsv \
      --output-prefix results/docs/figures/mondo_hpo_by_moe --png
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

# Fixed colour convention for this figure family: MONDO orange, HPO red/pink,
# tier totals blue.
COLOR_MONDO = "orange"
COLOR_HPO = "#F17C7C"
COLOR_TOTAL = "#7FA8D9"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("master_table", help="Path to ndd_master_table.tsv")
    ap.add_argument("--class-filter", default="candidate",
                    help="Restrict to this Class (default: candidate)")
    ap.add_argument("--all-classes", action="store_true",
                    help="Override --class-filter and use every gene with a MOE score")
    ap.add_argument("--output-prefix", default="mondo_hpo_by_moe",
                    help="Output path without extension")
    ap.add_argument("--png", action="store_true",
                    help="Also write a PNG alongside the PDF")
    args = ap.parse_args()

    if not os.path.exists(args.master_table):
        sys.exit(f"Error: file not found: {args.master_table}")

    df = pd.read_csv(args.master_table, sep="\t", dtype=str)
    for col in ("MOE_score", "MONDO_terms", "n_hpo_terms"):
        if col not in df.columns:
            sys.exit(f"Error: master table is missing required column: {col}")

    df["MOE_score"] = pd.to_numeric(df["MOE_score"], errors="coerce")
    df = df.dropna(subset=["MOE_score"])
    df["MOE_score"] = df["MOE_score"].astype(int)

    if not args.all_classes:
        if "Class" not in df.columns:
            sys.exit("Error: --class-filter given but master table has no 'Class' column")
        df = df[df["Class"] == args.class_filter]
        if df.empty:
            sys.exit(f"Error: no genes with Class == '{args.class_filter}'")

    has_mondo = df["MONDO_terms"].notna() & (df["MONDO_terms"].str.strip() != "")
    has_hpo = pd.to_numeric(df["n_hpo_terms"], errors="coerce").fillna(0) > 0

    rows = []
    for tier in sorted(df["MOE_score"].unique()):
        sel = df["MOE_score"] == tier
        total = int(sel.sum())
        n_mondo = int((sel & has_mondo).sum())
        n_hpo = int((sel & has_hpo).sum())
        rows.append({"MOE_score": tier, "genes_with_MONDO": n_mondo,
                     "genes_with_HPO": n_hpo, "total_genes": total,
                     "fraction_MONDO": n_mondo / total,
                     "fraction_HPO": n_hpo / total})
    counts = pd.DataFrame(rows)

    print(counts.to_string(index=False))
    out_tsv = f"{args.output_prefix}.tsv"
    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    counts.to_csv(out_tsv, sep="\t", index=False, float_format="%.4f")
    print(f"\nCounts saved to: {out_tsv}")

    x = np.arange(len(counts))
    width = 0.27
    fig, ax = plt.subplots(figsize=(10, 6), facecolor="white")
    ax.set_facecolor("white")
    ax.bar(x - width, counts["genes_with_MONDO"], width,
           color=COLOR_MONDO, edgecolor="#333333", linewidth=0.6)
    ax.bar(x, counts["genes_with_HPO"], width,
           color=COLOR_HPO, edgecolor="#333333", linewidth=0.6)
    ax.bar(x + width, counts["total_genes"], width,
           color=COLOR_TOTAL, edgecolor="#333333", linewidth=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels(counts["MOE_score"].astype(str), fontsize=13)
    ax.set_xlabel("MOE Score", fontsize=14, labelpad=10)
    ax.set_ylabel("Number of genes", fontsize=14, labelpad=10)
    label = "candidate genes" if not args.all_classes else "all scored genes"
    ax.set_title(f"MONDO / HPO term coverage by MOE tier ({label})", fontsize=15, pad=14)
    # Square colour swatches, top-left, naming the three gene sets.
    handles = [
        Patch(facecolor=COLOR_MONDO, edgecolor="#333333", linewidth=0.6, label="MONDO genes"),
        Patch(facecolor=COLOR_HPO, edgecolor="#333333", linewidth=0.6, label="HPO genes"),
        Patch(facecolor=COLOR_TOTAL, edgecolor="#333333", linewidth=0.6, label="MOE genes"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=12, frameon=False,
              handlelength=1.1, handleheight=1.1, labelspacing=0.6,
              borderaxespad=0.8)
    ax.tick_params(labelsize=12)

    # Full plot border, matching the other MOE-tier panels.
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#333333")
        spine.set_linewidth(0.8)

    fig.tight_layout()
    out_pdf = f"{args.output_prefix}.pdf"
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight", facecolor="white")
    print(f"PDF figure saved to: {out_pdf}")
    if args.png:
        out_png = f"{args.output_prefix}.png"
        fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"PNG figure saved to: {out_png}")


if __name__ == "__main__":
    main()
