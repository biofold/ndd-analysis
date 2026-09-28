#!/usr/bin/env python3
"""
Classification-rule check and per-source contribution (Reviewer 1 comment 1,
Reviewer 2 comment 3).

WHY THIS SCRIPT EXISTS
----------------------
R1.1 asks for a plain-language version of the Table 1 Boolean rules; the
response states that the plain-language rule reproduces the original class of
every gene. R2.3 asks what SFARI (and each resource) actually contributes, and
which five SFARI genes were left out of the candidate Venn diagram. This
script produces the evidence for all three statements from the master table's
seven source flags:

  1. rule_check.tsv  -- the Table 1 rule, re-implemented in evaluation order
     (high-confidence rules first, then "any flag" for candidate), cross-
     tabulated against the class assigned by the pipeline. Every off-diagonal
     cell must be zero.
  2. contribution.tsv -- leave-one-resource-out reclassification: for each of
     the four resources (Orphanet, SFARI Gene, GeneTrek, Sanchis-Juan), remove
     all of its flags, reclassify with the same rule, and count the genes that
     change class. This is each resource's unique contribution.
  3. sfari_only.tsv  -- candidate genes whose only source flag is SFARI (the
     genes that cannot be placed in a three-set Venn diagram of Orphanet /
     GeneTrek / Sanchis-Juan).

Rule (Table 1), with HC = GeneTrek high-confidence, LC = GeneTrek candidate:
  high-confidence  <- HC or SFARI score 1 ("sfari+") or (Orphanet neuro and
                      Orphanet develop) or (Sanchis-Juan and (neuro or develop
                      or LC or any SFARI))
  candidate        <- not high-confidence and at least one flag
  no reported evidence <- no flag

Usage:
  python3 utils/table_generation/source_contribution.py \\
      results/main/ndd_master_table.tsv --output-prefix results/docs/tables/source
"""

import argparse
import sys

import numpy as np
import pandas as pd

FLAGS = ["GeneTrek_HC", "GeneTrek_LC", "Orphanet_neuro", "Orphanet_develop",
         "SanchisJuan", "SFARI_plus", "SFARI_other"]
RESOURCES = {
    "Orphanet": ["Orphanet_neuro", "Orphanet_develop"],
    "SFARI Gene": ["SFARI_plus", "SFARI_other"],
    "GeneTrek": ["GeneTrek_HC", "GeneTrek_LC"],
    "Sanchis-Juan et al.": ["SanchisJuan"],
}
LABEL = {"curated": "high-confidence", "candidate": "candidate",
         "no_evidence": "no reported evidence"}


def classify(f):
    sfari_any = f["SFARI_plus"] | f["SFARI_other"]
    high = (f["GeneTrek_HC"] | f["SFARI_plus"]
            | (f["Orphanet_neuro"] & f["Orphanet_develop"])
            | (f["SanchisJuan"] & (f["Orphanet_neuro"] | f["Orphanet_develop"]
                                   | f["GeneTrek_LC"] | sfari_any)))
    return pd.Series(np.where(high, "curated",
                              np.where(f.any(axis=1), "candidate", "no_evidence")),
                     index=f.index)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("master_table")
    ap.add_argument("--output-prefix", required=True)
    args = ap.parse_args()

    mt = pd.read_csv(args.master_table, sep="\t", low_memory=False)
    missing = [c for c in FLAGS + ["Gene", "Class"] if c not in mt.columns]
    if missing:
        sys.exit(f"Error: master table lacks columns {missing}")
    f = mt[FLAGS].fillna(0).astype(int).astype(bool)

    # 1. rule check
    rule = classify(f)
    xt = pd.crosstab(mt["Class"].map(LABEL), rule.map(LABEL),
                     rownames=["pipeline class"], colnames=["Table 1 rule"])
    xt.to_csv(f"{args.output_prefix}_rule_check.tsv", sep="\t")
    n_agree = int((rule == mt["Class"]).sum())

    # 2. leave-one-resource-out
    rows = []
    for res, cols in RESOURCES.items():
        g = f.copy()
        g[cols] = False
        new = classify(g)
        changed = new != rule
        tr = pd.crosstab(rule[changed].map(LABEL), new[changed].map(LABEL))
        row = {"resource_removed": res,
               "genes_flagged_by_resource": int(f[cols].any(axis=1).sum()),
               "genes_changing_class": int(changed.sum())}
        for a in ["high-confidence", "candidate"]:
            for b in ["candidate", "no reported evidence"]:
                if a == b:
                    continue
                row[f"{a} -> {b}"] = int(tr.loc[a, b]) if a in tr.index and b in tr.columns else 0
        rows.append(row)
    contrib = pd.DataFrame(rows)
    contrib.to_csv(f"{args.output_prefix}_contribution.tsv", sep="\t", index=False)

    # 3. SFARI-only candidates
    others = [c for c in FLAGS if c not in RESOURCES["SFARI Gene"]]
    sf_only = mt[(mt["Class"] == "candidate")
                 & (f["SFARI_plus"] | f["SFARI_other"]) & ~f[others].any(axis=1)]
    keep = [c for c in ["Gene", "SFARI_score", "MOE_score"] if c in mt.columns]
    sf_only[keep].to_csv(f"{args.output_prefix}_sfari_only.tsv", sep="\t", index=False)

    print(f"rule_agreement\t{n_agree}/{len(mt)}")
    print(f"sfari_only_candidates\t{len(sf_only)}\t{','.join(sf_only['Gene'])}")
    print(contrib.to_string(index=False))


if __name__ == "__main__":
    main()
