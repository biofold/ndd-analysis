#!/usr/bin/env python3
"""
TADA-NDD association of the MOE score adjusted for gnomAD v2.1.1 constraint
(R1.3 robustness check; companion of moe_pli_validation.py for ClinVar).

Fits, on the candidate genes with a TADA-NDD result (moe_tada_input.tsv,
built by build_tada_input.py) and a non-missing constraint value:

    logit P(TADA_NDD = 1) = b0 + b1 * MOE_score              (unadjusted)
    logit P(TADA_NDD = 1) = b0 + b1 * MOE_score + b2 * Z      (Z = pLI or LOEUF)

TADA_NDD = 1 when the Fu et al. 2022 TADA-NDD FDR <= 0.1.

Outputs (in --output-dir), one row per term, same columns as
moe_<metric>_logistic_regression.tsv:
    moe_tada_unadjusted_logistic_regression.tsv
    moe_<metric>_tada_logistic_regression.tsv   (+ _summary.txt)
"""
import argparse
import os

import numpy as np
import pandas as pd
import statsmodels.api as sm


def fit(y, X):
    model = sm.Logit(y, sm.add_constant(X.astype(float))).fit(disp=False)
    ci = model.conf_int(alpha=0.05)
    tab = pd.DataFrame({
        "term": model.params.index,
        "coef": model.params.values,
        "std_err": model.bse.values,
        "z": model.tvalues.values,
        "p_value": model.pvalues.values,
    })
    tab["odds_ratio"] = np.exp(tab["coef"])
    tab["or_ci_low"] = np.exp(ci[0].values)
    tab["or_ci_high"] = np.exp(ci[1].values)
    tab["n_genes"] = int(len(y))
    tab["n_positive"] = int(y.sum())
    return model, tab


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("master_table", help="ndd_master_table.tsv (Gene, gnomAD_pLI, gnomAD_LOEUF)")
    ap.add_argument("--tada-input", required=True, help="moe_tada_input.tsv from build_tada_input.py")
    ap.add_argument("--outcome", default="TADA_NDD", help="Outcome column (TADA_NDD or TADA_NDD_fdr005)")
    ap.add_argument("--output-dir", default=".")
    args = ap.parse_args()

    tada = pd.read_csv(args.tada_input, sep="\t").dropna(subset=["FDR_TADA_NDD", "MOE_score"])
    master = pd.read_csv(args.master_table, sep="\t", dtype=str, usecols=["Gene", "gnomAD_pLI", "gnomAD_LOEUF"])
    for c in ("gnomAD_pLI", "gnomAD_LOEUF"):
        master[c] = pd.to_numeric(master[c], errors="coerce")
    df = tada.merge(master, on="Gene", how="left")
    y_all = df[args.outcome].astype(int)

    _, tab = fit(y_all, df[["MOE_score"]])
    out = os.path.join(args.output_dir, "moe_tada_unadjusted_logistic_regression.tsv")
    tab.to_csv(out, sep="\t", index=False)
    r = tab.set_index("term").loc["MOE_score"]
    print(f"TADA-NDD ~ MOE: n = {len(df):,}, OR = {r.odds_ratio:.4f} "
          f"({r.or_ci_low:.4f}-{r.or_ci_high:.4f}), p = {r.p_value:.3g}")

    for metric in ("pLI", "LOEUF"):
        col = f"gnomAD_{metric}"
        d = df.dropna(subset=[col])
        model, tab = fit(d[args.outcome].astype(int), d[["MOE_score", col]])
        out = os.path.join(args.output_dir, f"moe_{metric.lower()}_tada_logistic_regression.tsv")
        tab.to_csv(out, sep="\t", index=False)
        with open(out.replace(".tsv", "_summary.txt"), "w") as f:
            f.write(str(model.summary()))
        r = tab.set_index("term").loc["MOE_score"]
        print(f"TADA-NDD ~ MOE + {col}: n = {len(d):,}, OR = {r.odds_ratio:.4f} "
              f"({r.or_ci_low:.4f}-{r.or_ci_high:.4f}), p = {r.p_value:.3g}")


if __name__ == "__main__":
    main()
