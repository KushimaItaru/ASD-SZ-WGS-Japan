#!/usr/bin/env python3
# run_demo_burden.py
# - Minimal, self-contained DEMO of the primary burden model used in the manuscript.
# - Fits the B'-aligned logistic burden regression on the SYNTHETIC demo table:
#       case ~ burden_count + sex + PC1..PC10 + log1p(total_bases) + log1p(total_gene)
#   (same covariate structure as the primary case-control burden analyses).
# - Reports the odds ratio, 95% CI, Wald P value, degrees of freedom and N for the
#   burden term, plus total run time.
# - Input : demo_burden_input.tsv (created by make_demo_data.py; fully synthetic)
# - Depends only on declared dependencies: pandas, numpy, statsmodels.
# - Runs in < 5 s on a normal desktop computer; no controlled-access data required.

import time
import numpy as np
import pandas as pd
import statsmodels.api as sm

t0 = time.time()

df = pd.read_csv("demo_burden_input.tsv", sep="\t")

# B'-style covariates: log1p-transform span and gene-overlap burden
df["log1p_total_bases"] = np.log1p(df["total_bases"])
df["log1p_total_gene"]  = np.log1p(df["total_gene"])

predictors = ["burden_count", "sex"] + [f"PC{k}" for k in range(1, 11)] + \
             ["log1p_total_bases", "log1p_total_gene"]

X = sm.add_constant(df[predictors].astype(float))
y = df["case"].astype(int)

model = sm.Logit(y, X).fit(disp=0)

beta  = model.params["burden_count"]
se    = model.bse["burden_count"]
p     = model.pvalues["burden_count"]
or_   = np.exp(beta)
ci_lo = np.exp(beta - 1.96 * se)
ci_hi = np.exp(beta + 1.96 * se)

print("=" * 64)
print("DEMO: B'-aligned logistic burden regression (synthetic data)")
print("=" * 64)
print(f"N (controls / cases) : {(y==0).sum()} / {(y==1).sum()}")
print(f"Model                : case ~ burden_count + sex + PC1-10 + "
      f"log1p(total_bases) + log1p(total_gene)")
print(f"Burden term (Wald)   : df = 1")
print(f"  Odds ratio (OR)    : {or_:.3f}")
print(f"  95% CI             : {ci_lo:.3f} - {ci_hi:.3f}")
print(f"  P value            : {p:.3e}")
print(f"Total run time       : {time.time() - t0:.2f} s")
print("=" * 64)
