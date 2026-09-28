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
  3. rule_counts.tsv -- for high-confidence genes, how many satisfy each of
     the four rules, how many satisfy only that rule, and the distribution of
     the number of rules and of source flags per gene in every class (the
     flags are not mutually exclusive; this quantifies the overlap).
  4. decision_tree.tsv / decision_tree.png -- the rule as a decision tree in
     the same evaluation order as the pseudocode, with the number of genes
     leaving at every branch (each gene is counted once, at the first rule it
     satisfies, so these differ from the per-rule counts in rule_counts.tsv).
  5. sfari_only.tsv  -- candidate genes whose only source flag is SFARI (the
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


def plot_decision_tree(steps, out_png):
    """steps: list of (question, n_yes, yes_label) in evaluation order, the
    last entry being the any-flag test; plus the final 'no' count."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    COL = {"high-confidence": "#4689a3", "candidate": "#e08a2c",
           "no reported NDD evidence": "#9a9a9a"}
    plt.rcParams.update({"font.size": 8})
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    ax.set_xlim(0, 10); ax.set_ylim(-0.7, 10.45); ax.axis("off")

    def box(x, y, w, h, text, fc, ec, tc="black", bold=False):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                    boxstyle="round,pad=0.02,rounding_size=0.12",
                                    fc=fc, ec=ec, lw=1))
        ax.text(x, y, text, ha="center", va="center", color=tc,
                fontweight="bold" if bold else "normal", fontsize=8, linespacing=1.3)

    def arrow(x0, y0, x1, y1, label=None, lx=0, ly=0, ha="center"):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color="#333333", lw=0.9,
                                    shrinkA=0, shrinkB=0))
        if label:
            ax.text((x0 + x1) / 2 + lx, (y0 + y1) / 2 + ly, label, ha=ha, va="center",
                    fontsize=7, color="#333333")

    xq, wq, hq = 3.2, 5.0, 0.9
    ys = [9.2, 7.6, 6.0, 4.4, 2.8]
    n_total = steps["n_total"]
    box(xq, 10.15, 4.6, 0.42, f"All HGNC protein-coding genes (n = {n_total:,})",
        "white", "#333333")
    arrow(xq, 9.93, xq, ys[0] + hq / 2)

    remaining = n_total
    xh = 8.3
    for i, (q, n_yes, lab) in enumerate(steps["questions"]):
        y = ys[i]
        box(xq, y, wq, hq, f"{lab}\n{q}", "#f4f4f4", "#333333")
        remaining -= n_yes
        if i < 4:
            arrow(xq + wq / 2, y, xh - 1.25, y, f"yes: {n_yes:,}", ly=0.16)
        if i < len(steps["questions"]) - 1:
            arrow(xq, y - hq / 2, xq, ys[i + 1] + hq / 2, f"no: {remaining:,}",
                  lx=0.12, ha="left")

    # high-confidence exit box spanning the four rule rows
    hc_n = sum(n for _, n, _ in steps["questions"][:4])
    ytop, ybot = ys[0] + 0.45, ys[3] - 0.45
    ax.add_patch(FancyBboxPatch((xh - 1.25, ybot), 2.5, ytop - ybot,
                                boxstyle="round,pad=0.02,rounding_size=0.15",
                                fc=COL["high-confidence"], ec="#333333", lw=1))
    ax.text(xh, (ytop + ybot) / 2, f"High-confidence\n(score 2)\n\nn = {hc_n:,}",
            ha="center", va="center", color="white", fontweight="bold", fontsize=9)

    # final split: candidate vs no evidence
    q5, n_cand, _ = steps["questions"][4]
    n_none = steps["n_none"]
    yb = 1.0
    arrow(xq + wq / 2, ys[4], xh - 1.25, ys[4], f"yes: {n_cand:,}", ly=0.16)
    box(xh, ys[4], 2.5, 0.9, f"Candidate\n(score 1)\nn = {n_cand:,}",
        COL["candidate"], "#333333", tc="white", bold=True)
    arrow(xq, ys[4] - hq / 2, xq, yb + 0.45, f"no: {n_none:,}", lx=0.12, ha="left")
    box(xq, yb, 3.4, 0.9, f"No reported NDD evidence\n(score 0)\nn = {n_none:,}",
        COL["no reported NDD evidence"], "#333333", tc="white", bold=True)

    ax.text(0.05, -0.65,
            "HC/LC: GeneTrek high-confidence / candidate; NEU/DEV: Orphanet neurological / "
            "developmental;\nSJ: Sanchis-Juan et al.; SF1: SFARI Gene score 1; SF: other "
            "SFARI Gene entry. Counts: genes leaving at each branch.",
            fontsize=6.5, color="#555555", va="bottom")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


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

    # 3. per-rule counts and flag multiplicity
    sfari_any = f["SFARI_plus"] | f["SFARI_other"]
    rules = pd.DataFrame({
        "R1: GeneTrek high-confidence": f["GeneTrek_HC"],
        "R2: SFARI Gene score 1": f["SFARI_plus"],
        "R3: Orphanet neurological AND developmental": f["Orphanet_neuro"] & f["Orphanet_develop"],
        "R4: Sanchis-Juan AND >=1 supporting source": f["SanchisJuan"] & (
            f["Orphanet_neuro"] | f["Orphanet_develop"] | f["GeneTrek_LC"] | sfari_any),
    })
    hc = rule == "curated"
    n_rules = rules[hc].sum(axis=1)
    rc = []
    for name in rules.columns:
        rc.append({"section": "rule", "item": name,
                   "high_conf_genes_satisfying": int(rules.loc[hc, name].sum()),
                   "satisfying_only_this_rule": int((rules.loc[hc, name] & (n_rules == 1)).sum())})
    for k, cnt in n_rules.value_counts().sort_index().items():
        rc.append({"section": "rules_per_high_conf_gene", "item": int(k),
                   "high_conf_genes_satisfying": int(cnt), "satisfying_only_this_rule": ""})
    nflag = f.sum(axis=1)
    for cls in ["curated", "candidate", "no_evidence"]:
        vc = nflag[rule == cls].value_counts().sort_index()
        for k, cnt in vc.items():
            rc.append({"section": f"flags_per_gene:{LABEL[cls]}", "item": int(k),
                       "high_conf_genes_satisfying": int(cnt), "satisfying_only_this_rule": ""})
    sanchis_alone = f["SanchisJuan"] & ~f.drop(columns="SanchisJuan").any(axis=1)
    rc.append({"section": "sanchis_juan_without_support", "item": "candidate",
               "high_conf_genes_satisfying": int(sanchis_alone.sum()), "satisfying_only_this_rule": ""})
    pd.DataFrame(rc).rename(columns={"high_conf_genes_satisfying": "n_genes",
                                     "satisfying_only_this_rule": "n_only_this_rule"}) \
        .to_csv(f"{args.output_prefix}_rule_counts.tsv", sep="\t", index=False)

    # 4. decision tree, evaluation order of the pseudocode
    left = pd.Series(True, index=f.index)
    qs = [
        ("GeneTrek high-confidence (HC)?", f["GeneTrek_HC"], "Rule R1"),
        ("SFARI Gene score 1 (SF1)?", f["SFARI_plus"], "Rule R2"),
        ("Orphanet neurological AND developmental\n(NEU and DEV)?",
         f["Orphanet_neuro"] & f["Orphanet_develop"], "Rule R3"),
        ("Sanchis-Juan (SJ) AND at least one of\nNEU, DEV, LC, SF1, SF?",
         f["SanchisJuan"] & (f["Orphanet_neuro"] | f["Orphanet_develop"]
                             | f["GeneTrek_LC"] | sfari_any), "Rule R4"),
        # HC genes have all left at R1, so HC is not re-tested here
        ("Any source flag\n(LC, NEU, DEV, SJ, SF1, SF)?",
         f.drop(columns="GeneTrek_HC").any(axis=1), "Remaining genes"),
    ]
    tree = []
    for q, cond, lab in qs:
        yes = left & cond
        tree.append((q, int(yes.sum()), lab))
        left = left & ~cond
    n_none = int(left.sum())
    assert sum(n for _, n, _ in tree[:4]) == int((rule == "curated").sum())
    assert tree[4][1] == int((rule == "candidate").sum()) and n_none == int((rule == "no_evidence").sum())
    pd.DataFrame([{"step": lab, "question": q.replace("\n", " "), "n_yes": n}
                  for q, n, lab in tree] + [{"step": "end", "question": "no flag",
                                              "n_yes": n_none}]) \
        .to_csv(f"{args.output_prefix}_decision_tree.tsv", sep="\t", index=False)
    plot_decision_tree({"questions": tree, "n_none": n_none, "n_total": len(mt)},
                       f"{args.output_prefix}_decision_tree.png")

    # 5. SFARI-only candidates
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
