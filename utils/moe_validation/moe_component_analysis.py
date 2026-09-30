#!/usr/bin/env python3
"""
Component-level analyses of the MOE score, candidate set only.

WHY THIS SCRIPT EXISTS
----------------------
Reviewer 2, comment 10 asks whether the small differences between non-zero
MOE tiers reflect redundancy among the five annotation sources (GO BP/MF,
KEGG/Reactome "have a large overlap"), and whether HPO should be added as a
sixth component. Reviewer 1, comment 3 asks for the MOE >= 4 threshold to be
justified. Reviewer 2, comment 9 asks whether MOE merely tracks annotation
density. This script answers those directly from files already in the repo:

  1. Redundancy between the five binary MOE indicators: phi coefficient and
     Jaccard index for every pair, and a PCA of the five indicators
     (explained-variance ratio per component).                 -> redundancy.tsv, pca.tsv
  2. What each component adds: ClinVar P/LP ROC-AUC of the full score, of
     each single component, and of every leave-one-component-out 4-point
     score, each compared with the full score by DeLong's test for correlated
     ROC curves.                                              -> score_variants.tsv
  3. HPO as a sixth component, built exactly like the other five: a Fisher/
     hypergeometric enrichment of the high-confidence set against the
     19,354-gene background over HPO terms (BH-adjusted p < 0.01), and one
     point to a candidate annotated to at least one enriched term. MOE6 =
     MOE + HPO point, compared with MOE by DeLong's test. HPO gene-phenotype
     annotation is derived from known disease genes, so it is circular with
     ClinVar; the comparison is therefore also run against gnomAD pLI >= 0.9
     and ClinGen haploinsufficiency, which are not derived from disease
     annotation.                                              -> score_variants.tsv
  4. Threshold choice: Youden's J and F1 at every MOE cutoff, for ClinVar
     P/LP and pLI >= 0.9.                                     -> thresholds.tsv
  5. ClinGen haploinsufficiency (score 3, "sufficient evidence") by MOE tier
     as a third outcome external to the score's construction, with a
     Cochran-Armitage trend test.                             -> clingen_hi.tsv
  6. Annotation density: Spearman correlation of MOE with the number of
     KEGG and Reactome pathways per gene (the third density measure named in
     the response draft), and MOE's ClinVar ROC-AUC WITHIN quartiles of
     PubMed record count, i.e. among genes of comparable study volume.
                                                              -> annotation_strata.tsv

All analyses use the candidate class only (n = 6,856): MOE is defined only
for candidates, and including high-confidence genes would be circular.

Usage:
  python3 utils/moe_validation/moe_component_analysis.py \\
      results/main/ndd_master_table.tsv \\
      --hpo data/gene_hpo_terms.tsv --libs libs \\
      --output-prefix results/docs/validation/moe_component
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

COMPONENTS = ["MOE_GOBP", "MOE_GOCC", "MOE_GOMF", "MOE_KEGG", "MOE_Reactome"]
COMPONENT_LABEL = {
    "MOE_GOBP": "GO BP", "MOE_GOCC": "GO CC", "MOE_GOMF": "GO MF",
    "MOE_KEGG": "KEGG", "MOE_Reactome": "Reactome", "HPO": "HPO",
}

# HPO terms that describe inheritance, onset, frequency or other clinical
# modifiers rather than a phenotypic abnormality. They sit outside the
# "Phenotypic abnormality" branch (HP:0000118) of the ontology; hp.obo is not
# a repo input, so they are excluded by name, and the excluded count is
# reported.
HPO_NON_PHENOTYPE_PATTERN = (
    r"inheritance|onset|^sporadic$|mosaicism|penetrance|heterogeneity|"
    r"^contiguous gene syndrome$|^somatic mutation$|^de novo$|anticipation"
)


# --------------------------------------------------------------------------
# DeLong test for two correlated ROC AUCs (Sun & Xu 2014 fast algorithm)
# --------------------------------------------------------------------------
def _midrank(x):
    order = np.argsort(x)
    xs = x[order]
    n = len(x)
    t = np.zeros(n)
    i = 0
    while i < n:
        j = i
        while j < n and xs[j] == xs[i]:
            j += 1
        t[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    out = np.empty(n)
    out[order] = t
    return out


def _delong_components(y, scores):
    """scores: (k, n) array. Returns aucs (k,), covariance (k, k)."""
    y = np.asarray(y).astype(bool)
    pos, neg = scores[:, y], scores[:, ~y]
    m, n = pos.shape[1], neg.shape[1]
    k = scores.shape[0]
    tx = np.array([_midrank(p) for p in pos])
    ty = np.array([_midrank(q) for q in neg])
    tz = np.array([_midrank(s) for s in np.hstack([pos, neg])])
    aucs = tz[:, :m].sum(axis=1) / (m * n) - (m + 1.0) / (2.0 * n)
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    s01 = np.atleast_2d(np.cov(v01))
    s10 = np.atleast_2d(np.cov(v10))
    return aucs, s01 / m + s10 / n


def delong_auc_ci(y, s):
    auc, cov = _delong_components(y, np.atleast_2d(np.asarray(s, float)))
    se = np.sqrt(cov[0, 0])
    return auc[0], auc[0] - 1.96 * se, auc[0] + 1.96 * se


def delong_test(y, s1, s2):
    """Two-sided DeLong p-value for AUC(s1) - AUC(s2)."""
    aucs, cov = _delong_components(y, np.vstack([s1, s2]).astype(float))
    var = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    if var <= 0:
        return aucs[0] - aucs[1], 1.0
    z = (aucs[0] - aucs[1]) / np.sqrt(var)
    return aucs[0] - aucs[1], 2 * stats.norm.sf(abs(z))


# --------------------------------------------------------------------------
def fmt_p(p):
    return "<2.2e-308" if p == 0 else f"{p:.3g}"


def phi_and_jaccard(a, b):
    a, b = a.astype(bool), b.astype(bool)
    n11 = int((a & b).sum()); n10 = int((a & ~b).sum())
    n01 = int((~a & b).sum()); n00 = int((~a & ~b).sum())
    den = np.sqrt(float(n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
    phi = (n11 * n00 - n10 * n01) / den if den > 0 else np.nan
    union = n11 + n10 + n01
    jac = n11 / union if union else np.nan
    return phi, jac, n11, n10, n01, n00


def read_gmt_counts(path):
    """Number of gene sets each (upper-cased) gene belongs to."""
    counts = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            for g in set(x.strip().upper() for x in parts[2:] if x.strip()):
                counts[g] = counts.get(g, 0) + 1
    return counts


def hpo_component(master, hpo_file, alpha=0.01):
    """Build the HPO MOE component exactly like the five existing ones."""
    hpo = pd.read_csv(hpo_file, sep="\t", dtype=str)
    hpo = hpo[hpo["Gene"].notna() & (hpo["Gene"].str.strip() != "-")]
    universe = set(master["Gene"].str.upper())
    term_genes, term_name = {}, {}
    n_excluded_terms = set()
    for _, r in hpo.iterrows():
        g = r["Gene"].strip().upper()
        if g not in universe or pd.isna(r["HPO_id"]):
            continue
        ids = r["HPO_id"].split(";")
        names = r["HPO_name"].split(";") if isinstance(r["HPO_name"], str) else [""] * len(ids)
        for tid, tname in zip(ids, names):
            tid, tname = tid.strip(), tname.strip()
            if not tid:
                continue
            if pd.Series([tname]).str.contains(HPO_NON_PHENOTYPE_PATTERN, case=False, regex=True).iloc[0]:
                n_excluded_terms.add(tid)
                continue
            term_genes.setdefault(tid, set()).add(g)
            term_name[tid] = tname

    hc = set(master.loc[master["Class"] == "curated", "Gene"].str.upper())
    N, L = len(universe), len(hc)
    rows = []
    for tid, genes in term_genes.items():
        k = len(genes & hc)
        if k == 0:
            continue
        n = len(genes)
        p = stats.hypergeom.sf(k - 1, N, n, L)  # P(X >= k), Enrichr's test
        rows.append((tid, term_name[tid], k, n, p))
    enr = pd.DataFrame(rows, columns=["HPO_id", "HPO_name", "k_high_conf", "n_term", "p"])
    # Benjamini-Hochberg
    enr = enr.sort_values("p").reset_index(drop=True)
    m = len(enr)
    q = enr["p"].values * m / np.arange(1, m + 1)
    enr["p_adj"] = np.minimum.accumulate(q[::-1])[::-1].clip(max=1)
    sig = set(enr.loc[enr["p_adj"] < alpha, "HPO_id"])

    point = {}
    for tid in sig:
        for g in term_genes[tid]:
            point[g] = 1
    comp = master["Gene"].str.upper().map(point).fillna(0).astype(int)
    info = {
        "hpo_terms_tested": m,
        "hpo_terms_significant": len(sig),
        "hpo_terms_excluded_non_phenotype": len(n_excluded_terms),
    }
    return comp, enr, info


def cochran_armitage(tiers, y):
    """Two-sided Cochran-Armitage trend test, scores = tier values."""
    t = np.asarray(tiers, float); y = np.asarray(y, float)
    N = len(y); pbar = y.mean()
    tbar = t.mean()
    num = np.sum(y * (t - tbar))
    var = pbar * (1 - pbar) * np.sum((t - tbar) ** 2)
    z = num / np.sqrt(var)
    return z, 2 * stats.norm.sf(abs(z))


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


def plot_figure(red, variants_df, out_png):
    """Panel a: phi matrix of the five MOE indicators. Panels b-d: ROC-AUC
    (DeLong 95% CI) of the full score, each leave-one-component-out score and
    the HPO-augmented score, for three outcomes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                         "xtick.labelsize": 7, "ytick.labelsize": 7,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "xtick.direction": "out", "ytick.direction": "out"})
    labels = [COMPONENT_LABEL[c] for c in COMPONENTS]
    M = np.eye(len(labels))
    for _, r in red.iterrows():
        i, j = labels.index(r["component_1"]), labels.index(r["component_2"])
        M[i, j] = M[j, i] = r["phi"]
    fig = plt.figure(figsize=(7.2, 2.9))
    gl = fig.add_gridspec(1, 1, left=0.10, right=0.275, bottom=0.22, top=0.86)
    gs = fig.add_gridspec(1, 4, wspace=0.12, left=0.415, right=0.99, bottom=0.22, top=0.86)
    ax = fig.add_subplot(gl[0])
    Mplot = np.where(np.eye(len(labels)) == 1, np.nan, M)
    im = ax.imshow(Mplot, cmap="Blues", vmin=0, vmax=1)
    for i in range(len(labels)):
        for j in range(len(labels)):
            if i != j:
                ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if M[i, j] > 0.6 else "black")
    ax.set_xticks(range(len(labels)), labels, rotation=45, ha="right")
    ax.set_yticks(range(len(labels)), labels)
    ax.set_title("Pairwise phi coefficient", loc="left")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.text(-0.36, 1.06, "a", transform=ax.transAxes, fontweight="bold", fontsize=10)

    order = ["MOE (5 components)"] + [f"MOE without {COMPONENT_LABEL[c]}" for c in COMPONENTS] \
            + ["MOE6 (MOE + HPO)"]
    ylab = ["MOE (all five)"] + [f"without {COMPONENT_LABEL[c]}" for c in COMPONENTS] + ["MOE + HPO"]
    y = np.arange(len(order))[::-1]
    outcomes = ["ClinVar P/LP", "gnomAD pLI >= 0.9", "gnomAD LOEUF decile 1", "ClinGen HI = 3"]
    titles = ["ClinVar P/LP", "pLI \u2265 0.9", "LOEUF top decile", "ClinGen HI = 3"]
    focal, grey, hpo_c = "#1f5f8b", "#8a8a8a", "#c0504d"
    first = None
    for k, (oname, title) in enumerate(zip(outcomes, titles)):
        axk = fig.add_subplot(gs[k], sharey=first)
        first = first or axk
        sub = variants_df[variants_df["outcome"] == oname].set_index("score").loc[order]
        ref = sub.loc["MOE (5 components)", "roc_auc"]
        axk.axvline(ref, color=focal, lw=0.8, ls=":", zorder=0)
        for yi, (name, r) in zip(y, sub.iterrows()):
            c = focal if name.startswith("MOE (5") else (hpo_c if "HPO" in name else grey)
            axk.errorbar(r["roc_auc"], yi, xerr=[[r["roc_auc"] - r["auc_ci_low"]],
                                                  [r["auc_ci_high"] - r["roc_auc"]]],
                         fmt="o", ms=4, color=c, ecolor=c, elinewidth=1, capsize=0)
        axk.set_title(title, loc="left")
        axk.set_xlabel("ROC-AUC (95% CI)")
        axk.margins(x=0.12)
        axk.xaxis.set_major_locator(plt.MaxNLocator(3))
        if k == 0:
            axk.set_yticks(y, ylab)
            axk.text(-0.62, 1.06, "b", transform=axk.transAxes, fontweight="bold", fontsize=10,
                     ha="right")
        else:
            plt.setp(axk.get_yticklabels(), visible=False)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("master_table")
    ap.add_argument("--hpo", required=True, help="data/gene_hpo_terms.tsv")
    ap.add_argument("--libs", required=True, help="libs/ directory with the GMT files")
    ap.add_argument("--kegg-gmt", default="KEGG_2021_Human.gmt",
                    help="KEGG library used by the MOE score (default: %(default)s)")
    ap.add_argument("--reactome-gmt", default="Reactome_Pathways_2024.gmt")
    ap.add_argument("--output-prefix", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.output_prefix)), exist_ok=True)
    full = pd.read_csv(args.master_table, sep="\t", low_memory=False)
    cand = full[full["Class"] == "candidate"].copy()
    if len(cand) == 0:
        sys.exit("Error: no candidate rows in the master table")
    for c in COMPONENTS + ["MOE_score"]:
        cand[c] = pd.to_numeric(cand[c], errors="coerce").fillna(0).astype(int)
    assert (cand[COMPONENTS].sum(axis=1) == cand["MOE_score"]).all(), \
        "MOE_score is not the sum of its five components"

    y_clin = (pd.to_numeric(cand["n_pathogenic_likely_pathogenic"], errors="coerce").fillna(0) >= 1).astype(int).values
    pli = pd.to_numeric(cand["gnomAD_pLI"], errors="coerce")
    # gnomAD's own LOEUF decile (0 = the 10% most constrained genes)
    loeuf_dec = pd.to_numeric(cand["gnomAD_LOEUF_decile"], errors="coerce")
    hi = pd.to_numeric(cand["ClinGen_HI_score"], errors="coerce")

    # ---- 1. redundancy ----------------------------------------------------
    red = []
    for i, a in enumerate(COMPONENTS):
        for b in COMPONENTS[i + 1:]:
            phi, jac, n11, n10, n01, n00 = phi_and_jaccard(cand[a].values, cand[b].values)
            red.append({"component_1": COMPONENT_LABEL[a], "component_2": COMPONENT_LABEL[b],
                        "phi": round(phi, 3), "jaccard": round(jac, 3),
                        "both": n11, "only_1": n10, "only_2": n01, "neither": n00})
    prev = {COMPONENT_LABEL[c]: round(cand[c].mean(), 3) for c in COMPONENTS}
    red = pd.DataFrame(red)
    red.to_csv(f"{args.output_prefix}_redundancy.tsv", sep="\t", index=False)

    X = cand[COMPONENTS].values.astype(float)
    Xs = (X - X.mean(0)) / X.std(0, ddof=0)
    eigval, eigvec = np.linalg.eigh(np.cov(Xs, rowvar=False))
    order = np.argsort(eigval)[::-1]
    eigval, eigvec = eigval[order], eigvec[:, order]
    evr = eigval / eigval.sum()
    pca = pd.DataFrame({
        "principal_component": [f"PC{i+1}" for i in range(len(evr))],
        "explained_variance_ratio": np.round(evr, 3),
        "cumulative": np.round(np.cumsum(evr), 3),
    })
    for j, c in enumerate(COMPONENTS):
        pca[f"loading_{COMPONENT_LABEL[c].replace(' ', '_')}"] = np.round(eigvec[j, :] * np.sign(eigvec[:, 0].sum()), 3)
    pca.to_csv(f"{args.output_prefix}_pca.tsv", sep="\t", index=False)

    # ---- 3. HPO component -------------------------------------------------
    hpo_point_all, hpo_enr, hpo_info = hpo_component(full, args.hpo)
    cand["HPO_point"] = hpo_point_all.loc[cand.index].values
    cand["MOE6"] = cand["MOE_score"] + cand["HPO_point"]
    hpo_enr.to_csv(f"{args.output_prefix}_hpo_enrichment.tsv", sep="\t", index=False)
    has_hpo = pd.to_numeric(cand["n_hpo_terms"], errors="coerce").fillna(0) > 0

    # ---- 2 + 3. score variants vs three outcomes ---------------------------
    outcomes = {
        "ClinVar P/LP": (np.ones(len(cand), bool), y_clin),
        "gnomAD pLI >= 0.9": (pli.notna().values, (pli >= 0.9).astype(int).values),
        "gnomAD LOEUF decile 1": (loeuf_dec.notna().values, (loeuf_dec == 0).astype(int).values),
        # genes not curated by ClinGen are counted as "not HI 3"
        "ClinGen HI = 3": (np.ones(len(cand), bool), (hi == 3).astype(int).values),
    }
    variants = {"MOE (5 components)": cand["MOE_score"].values}
    for c in COMPONENTS:
        variants[f"single: {COMPONENT_LABEL[c]}"] = cand[c].values
    for c in COMPONENTS:
        variants[f"MOE without {COMPONENT_LABEL[c]}"] = (cand["MOE_score"] - cand[c]).values
    variants["MOE6 (MOE + HPO)"] = cand["MOE6"].values
    variants["single: HPO"] = cand["HPO_point"].values

    rows = []
    for oname, (mask, y) in outcomes.items():
        ym = y[mask]
        ref = variants["MOE (5 components)"][mask]
        for vname, s in variants.items():
            sm = s[mask]
            auc, lo, hi_ = delong_auc_ci(ym, sm)
            if vname == "MOE (5 components)":
                d, p = 0.0, np.nan
            else:
                d, p = delong_test(ym, sm, ref)
            rows.append({"outcome": oname, "n": int(mask.sum()), "n_positive": int(ym.sum()),
                         "score": vname, "roc_auc": round(auc, 3),
                         "auc_ci_low": round(lo, 3), "auc_ci_high": round(hi_, 3),
                         "delta_auc_vs_MOE": round(d, 3),
                         "delong_p_vs_MOE": "" if np.isnan(p) else fmt_p(p)})
    variants_df = pd.DataFrame(rows)
    variants_df.to_csv(f"{args.output_prefix}_score_variants.tsv", sep="\t", index=False)
    plot_figure(red, variants_df, f"{args.output_prefix}_figure.png")

    # ---- 4. threshold choice ---------------------------------------------
    trows = []
    for oname in ["ClinVar P/LP", "gnomAD pLI >= 0.9", "gnomAD LOEUF decile 1"]:
        mask, y = outcomes[oname]
        s, ym = cand["MOE_score"].values[mask], y[mask]
        for cut in range(1, 6):
            pred = s >= cut
            tp = int((pred & (ym == 1)).sum()); fp = int((pred & (ym == 0)).sum())
            fn = int((~pred & (ym == 1)).sum()); tn = int((~pred & (ym == 0)).sum())
            sens = tp / (tp + fn); spec = tn / (tn + fp)
            prec = tp / (tp + fp) if tp + fp else np.nan
            f1 = 2 * prec * sens / (prec + sens) if prec + sens else np.nan
            trows.append({"outcome": oname, "moe_cutoff": f">={cut}", "n_selected": tp + fp,
                          "sensitivity": round(sens, 3), "specificity": round(spec, 3),
                          "precision": round(prec, 3), "youden_j": round(sens + spec - 1, 3),
                          "f1": round(f1, 3)})
    thr = pd.DataFrame(trows)
    for oname in thr["outcome"].unique():
        sub = thr[thr["outcome"] == oname]
        thr.loc[sub["youden_j"].idxmax(), "note"] = "max Youden J"
        thr.loc[sub["f1"].idxmax(), "note"] = (thr.loc[sub["f1"].idxmax(), "note"] + "; max F1"
                                               if isinstance(thr.loc[sub["f1"].idxmax(), "note"], str)
                                               else "max F1")
    thr["note"] = thr["note"].fillna("")
    thr.to_csv(f"{args.output_prefix}_thresholds.tsv", sep="\t", index=False)

    # ---- 5. ClinGen HI by tier -------------------------------------------
    y_hi = (hi == 3).astype(int).values
    hrows = []
    for t in range(6):
        m = cand["MOE_score"].values == t
        k, n = int(y_hi[m].sum()), int(m.sum())
        lo, up = wilson(k, n)
        hrows.append({"moe_score": t, "n_genes": n, "n_hi3": k, "fraction": round(k / n, 4),
                      "ci_low": round(lo, 4), "ci_high": round(up, 4)})
    z, p_trend = cochran_armitage(cand["MOE_score"].values, y_hi)
    pred = cand["MOE_score"].values >= 4
    tab = [[int((pred & (y_hi == 1)).sum()), int((pred & (y_hi == 0)).sum())],
           [int((~pred & (y_hi == 1)).sum()), int((~pred & (y_hi == 0)).sum())]]
    or4, p4 = stats.fisher_exact(tab, alternative="greater")
    a, b, c_, d = [x + 0.5 if 0 in sum(tab, []) else x for x in sum(tab, [])]
    se = np.sqrt(1 / a + 1 / b + 1 / c_ + 1 / d)
    hdf = pd.DataFrame(hrows)
    hdf.attrs = {}
    with open(f"{args.output_prefix}_clingen_hi.tsv", "w") as fh:
        hdf.to_csv(fh, sep="\t", index=False)
        fh.write(f"# Cochran-Armitage trend z = {z:.3f}, two-sided p = {fmt_p(p_trend)}\n")
        fh.write(f"# MOE>=4 vs <4: OR = {or4:.3f} (95% CI {np.exp(np.log(or4)-1.96*se):.3f}-"
                 f"{np.exp(np.log(or4)+1.96*se):.3f}), one-sided Fisher p = {fmt_p(p4)}\n")

    # ---- 6. annotation density -------------------------------------------
    kegg = read_gmt_counts(os.path.join(args.libs, args.kegg_gmt))
    reac = read_gmt_counts(os.path.join(args.libs, args.reactome_gmt))
    gu = cand["Gene"].str.upper()
    cand["n_kegg_pathways"] = gu.map(kegg).fillna(0).astype(int)
    cand["n_reactome_pathways"] = gu.map(reac).fillna(0).astype(int)
    arows = []
    for col, label in [("n_kegg_pathways", "KEGG pathways per gene"),
                       ("n_reactome_pathways", "Reactome pathways per gene")]:
        rho, p = stats.spearmanr(cand["MOE_score"], cand[col])
        n = len(cand)
        zf = np.arctanh(rho); se_z = 1 / np.sqrt(n - 3)
        arows.append({"analysis": "spearman", "measure": label, "stratum": "all candidates",
                      "n": n, "n_positive": "", "value": round(rho, 3),
                      "ci_low": round(np.tanh(zf - 1.96 * se_z), 3),
                      "ci_high": round(np.tanh(zf + 1.96 * se_z), 3), "p": fmt_p(p)})
    npub = pd.to_numeric(cand["n_pubmed"], errors="coerce")
    ok = npub.notna().values
    q = pd.qcut(npub[ok].rank(method="first"), 4, labels=["Q1 (least studied)", "Q2", "Q3", "Q4 (most studied)"])
    sub = cand[ok].assign(pub_q=q.values)
    ysub = y_clin[ok]
    for lab in q.cat.categories:
        m = (sub["pub_q"] == lab).values
        auc, lo, up = delong_auc_ci(ysub[m], sub["MOE_score"].values[m])
        rng = npub[ok][m]
        arows.append({"analysis": "MOE ROC-AUC for ClinVar P/LP within PubMed-count quartile",
                      "measure": f"PubMed records {int(rng.min())}-{int(rng.max())}",
                      "stratum": lab, "n": int(m.sum()), "n_positive": int(ysub[m].sum()),
                      "value": round(auc, 3), "ci_low": round(lo, 3), "ci_high": round(up, 3),
                      "p": ""})
    pd.DataFrame(arows).to_csv(f"{args.output_prefix}_annotation_strata.tsv", sep="\t", index=False)

    # ---- summary ----------------------------------------------------------
    summary = {
        "n_candidates": len(cand),
        "component_prevalence": prev,
        **hpo_info,
        "candidates_with_any_hpo_annotation": int(has_hpo.sum()),
        "candidates_receiving_hpo_point": int(cand["HPO_point"].sum()),
        "hpo_point_among_hpo_annotated": round(cand.loc[has_hpo, "HPO_point"].mean(), 3),
        "clingen_hi3_total": int(y_hi.sum()),
        "clingen_trend_z": round(z, 3), "clingen_trend_p": fmt_p(p_trend),
    }
    pd.Series(summary).to_csv(f"{args.output_prefix}_summary.tsv", sep="\t", header=False)
    for k, v in summary.items():
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
