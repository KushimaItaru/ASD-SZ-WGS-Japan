#!/usr/bin/env python3
# make_demo_data.py
# - Generates a SMALL, fully SYNTHETIC per-individual burden table for the demo.
# - No real participant data are used; all values are simulated with a fixed seed.
# - Output columns mirror the B'-model burden inputs used in the manuscript:
#     sample_id, group (Control/Case), case (0/1), burden_count,
#     sex (0/1), PC1..PC10, total_bases, total_gene
# - A modest case-enriched burden (ground-truth OR ~ 1.5 per unit) is planted so
#   the demo regression returns a clear, positive, significant estimate.
# - Writes: demo_burden_input.tsv

import numpy as np
import pandas as pd

SEED = 42
rng = np.random.default_rng(SEED)

N_CTRL = 400   # synthetic controls
N_CASE = 200   # synthetic cases

def make_group(n, is_case):
    lam = 1.5 if is_case else 1.0          # planted case enrichment
    df = pd.DataFrame({
        "sample_id":  [f"{'CASE' if is_case else 'CTRL'}{i:04d}" for i in range(n)],
        "group":      ["Case" if is_case else "Control"] * n,
        "case":       np.full(n, 1 if is_case else 0, dtype=int),
        "burden_count": rng.poisson(lam, n).astype(int),
        "sex":        rng.integers(0, 2, n),                       # 0/1 synthetic sex
        "total_bases": rng.integers(2_000_000, 6_000_000, n),      # synthetic SV span
        "total_gene":  rng.integers(50, 400, n),                   # synthetic gene overlap
    })
    for k in range(1, 11):                                         # synthetic PC1..PC10
        df[f"PC{k}"] = rng.normal(0, 1, n).round(4)
    return df

demo = pd.concat([make_group(N_CTRL, False), make_group(N_CASE, True)], ignore_index=True)
cols = ["sample_id", "group", "case", "burden_count", "sex"] + \
       [f"PC{k}" for k in range(1, 11)] + ["total_bases", "total_gene"]
demo = demo[cols].sample(frac=1, random_state=SEED).reset_index(drop=True)
demo.to_csv("demo_burden_input.tsv", sep="\t", index=False)
print(f"Wrote demo_burden_input.tsv  ({len(demo)} synthetic individuals: "
      f"{(demo.case==0).sum()} controls, {(demo.case==1).sum()} cases)")
