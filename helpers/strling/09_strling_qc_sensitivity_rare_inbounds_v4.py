#!/usr/bin/env python3
# ============================================================================
# NOTE: This script contains hardcoded default paths specific to the original
# analysis environment (NIG supercomputer). To run in a different environment,
# update the paths below or set the corresponding environment variables
# (e.g. SAMPLE_INFO, PCA_EIGENVEC, CRAM_BASE_DIR1, CRAM_BASE_DIR2).
# ============================================================================
# 09_strling_qc_sensitivity_rare_inbounds_v4.py
# - Description:
#   - Load v9 per_sample.tsv (with rare_outlier_count/rare_any) and run QC + sensitivity analysis
#   - Per-group summary (Depth/PC/rare_outlier_count/rare_any, etc.)
#   - SMD (ASD vs Healthy / SZ vs Healthy)
#   - Models (if statsmodels available):
#       1) Logistic on rare_any (robust to zero-inflation and extreme values)
#       2) Poisson with offset on rare_outlier_count (exposure=tested_loci_total)
#       3) Poisson with q99 winsorization on rare_outlier_count (extreme value sensitivity)
#       4) Poisson with trimming (thr<10/20/30) on rare_outlier_count (extreme value sensitivity)
#   - Output top samples by rare_outlier_count
#   - Record execution time
#
# Usage (no arguments):
#   cd <repo_root>/helpers/strling
#   python 09_strling_qc_sensitivity_rare_inbounds_v4.py

from __future__ import annotations

import argparse
import math
import time
from datetime import datetime
from pathlib import Path
from typing import Iterable, List, Optional

import numpy as np
import pandas as pd

try:
    import statsmodels.api as sm
except Exception:
    sm = None


def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def detect_col(cols: Iterable[str], candidates: List[str]) -> Optional[str]:
    low = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    return None


def smd(x1: np.ndarray, x0: np.ndarray) -> float:
    x1 = x1.astype(float)
    x0 = x0.astype(float)
    m1 = np.nanmean(x1)
    m0 = np.nanmean(x0)
    v1 = np.nanvar(x1, ddof=1)
    v0 = np.nanvar(x0, ddof=1)
    s = math.sqrt((v1 + v0) / 2.0) if (v1 + v0) > 0 else np.nan
    return float((m1 - m0) / s) if s and not math.isnan(s) and s > 0 else np.nan


def group_summary(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    g = df.groupby("Group")[metric]
    return pd.DataFrame({
        "metric": metric,
        "Group": g.mean().index,
        "N": g.size().values,
        "mean": g.mean().values,
        "sd": g.std(ddof=1).values,
        "median": g.median().values,
        "q1": g.quantile(0.25).values,
        "q3": g.quantile(0.75).values,
        "min": g.min().values,
        "max": g.max().values,
    })


def fit_logit(df: pd.DataFrame, y_col: str, cov_cols: List[str]) -> pd.DataFrame:
    if sm is None:
        return pd.DataFrame()
    y = df[y_col].astype(int).values
    X = df[cov_cols].copy()
    X = sm.add_constant(X, has_constant="add")
    m = sm.Logit(y, X).fit(disp=False)
    out = []
    for term in m.params.index:
        beta = float(m.params[term])
        se = float(m.bse[term])
        OR = math.exp(beta)
        lo = math.exp(beta - 1.96 * se)
        hi = math.exp(beta + 1.96 * se)
        p = float(m.pvalues[term])
        out.append((term, beta, se, OR, lo, hi, p))
    return pd.DataFrame(out, columns=["term", "coef", "se", "OR", "OR_L", "OR_U", "p"])


def fit_poisson_offset(df: pd.DataFrame, y_col: str, exposure_col: str, cov_cols: List[str]) -> pd.DataFrame:
    if sm is None:
        return pd.DataFrame()
    y = df[y_col].astype(float).values
    X = df[cov_cols].copy()
    X = sm.add_constant(X, has_constant="add")
    offset = np.log(df[exposure_col].astype(float).values)
    m = sm.GLM(y, X, family=sm.families.Poisson(), offset=offset).fit()
    out = []
    for term in m.params.index:
        beta = float(m.params[term])
        se = float(m.bse[term])
        irr = math.exp(beta)
        lo = math.exp(beta - 1.96 * se)
        hi = math.exp(beta + 1.96 * se)
        p = float(m.pvalues[term])
        out.append((term, beta, se, irr, lo, hi, p))
    return pd.DataFrame(out, columns=["term", "coef", "se", "IRR", "IRR_L", "IRR_U", "p"])


def main() -> None:
    project_root = Path(__file__).resolve().parents[2]
    outdir = project_root / "analysis_results_strling"
    outdir.mkdir(parents=True, exist_ok=True)

    ap = argparse.ArgumentParser(add_help=True)
    # Changes: default set to v9
    ap.add_argument(
        "--per_sample_tsv",
        default=str(outdir / "strling_outlier_burden_rare_crossfit_inbounds_v9.per_sample.tsv"),
    )
    # Changes: output prefix also set to v9
    ap.add_argument("--out_prefix", default="qc_sensitivity_strling_inbounds_v9", help="Output prefix")
    ap.add_argument("--topn", type=int, default=50)
    ap.add_argument("--trim_thresholds", default="10,20,30")
    args = ap.parse_args()

    t0 = time.time()
    print(f"[{ts()}] [INFO] Start {Path(__file__).name}")
    print(f"[{ts()}] [INFO] per_sample_tsv={args.per_sample_tsv}")
    print(f"[{ts()}] [INFO] statsmodels_available={sm is not None}")

    per = Path(args.per_sample_tsv)
    if not per.exists():
        raise SystemExit(f"[ERROR] per_sample_tsv not found: {per}")

    df = pd.read_csv(per, sep="\t", dtype=str)

    # required columns
    req = ["SampleID", "Group", "Pedigree_No", "Depth", "Sex_M", "tested_loci_total", "rare_outlier_count", "rare_any"]
    for r in req:
        if r not in df.columns:
            raise SystemExit(f"[ERROR] Missing {r} in per_sample.tsv")
    for i in range(1, 11):
        if f"PC{i}" not in df.columns:
            raise SystemExit(f"[ERROR] Missing PC{i} in per_sample.tsv")

    # numeric
    num_cols = ["Depth", "Sex_M", "tested_loci_total", "rare_outlier_count", "rare_any"] + [f"PC{i}" for i in range(1, 11)]
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df[df["Group"].isin(["Healthy", "ASD", "SZ"])].copy()
    df = df.dropna(subset=num_cols).copy()
    df = df[df["tested_loci_total"] > 0].copy()

    print(f"[{ts()}] [INFO] N (complete-case) = {len(df)}")
    print(f"[{ts()}] [INFO] Group counts:\n{df['Group'].value_counts().to_string()}")

    # outputs
    report = outdir / f"{args.out_prefix}.report.txt"
    tables = outdir / f"{args.out_prefix}.tables.tsv"
    top = outdir / f"{args.out_prefix}.top{args.topn}_samples.tsv"
    models = outdir / f"{args.out_prefix}.models.tsv"

    # group summaries
    gs = []
    for m in ["tested_loci_total", "Depth", "Sex_M", "rare_outlier_count", "rare_any"]:
        gs.append(group_summary(df, m))
    gs_df = pd.concat(gs, ignore_index=True)

    # SMDs
    smd_rows = []
    h = df[df["Group"] == "Healthy"]
    for m in ["Depth", "Sex_M"] + [f"PC{i}" for i in range(1, 11)]:
        a = df[df["Group"] == "ASD"][m].to_numpy(dtype=float)
        s = df[df["Group"] == "SZ"][m].to_numpy(dtype=float)
        hh = h[m].to_numpy(dtype=float)
        smd_rows.append((m, smd(a, hh), smd(s, hh)))
    smd_df = pd.DataFrame(smd_rows, columns=["metric", "SMD_ASD_vs_Healthy", "SMD_SZ_vs_Healthy"])

    # top samples
    df_top = df.sort_values(["rare_outlier_count", "Depth"], ascending=[False, False]).head(args.topn).copy()
    df_top[["SampleID", "Group", "Pedigree_No", "Depth", "rare_outlier_count", "rare_any", "tested_loci_total"] + [f"PC{i}" for i in range(1, 11)]].to_csv(
        top, sep="\t", index=False
    )

    # models
    out_models = []

    # design covariates
    df["Group_ASD"] = (df["Group"] == "ASD").astype(int)
    df["Group_SZ"] = (df["Group"] == "SZ").astype(int)
    base_cov = ["Group_ASD", "Group_SZ", "Depth", "Sex_M"] + [f"PC{i}" for i in range(1, 11)]

    if sm is not None:
        # 1) logistic rare_any
        tab_logit = fit_logit(df, y_col="rare_any", cov_cols=base_cov)
        if not tab_logit.empty:
            tab_logit.insert(0, "_model", "logit_rare_any")
            out_models.append(tab_logit)

        # 2) poisson offset count
        tab_pois = fit_poisson_offset(df, y_col="rare_outlier_count", exposure_col="tested_loci_total", cov_cols=base_cov)
        if not tab_pois.empty:
            tab_pois.insert(0, "_model", "poisson_offset_count")
            out_models.append(tab_pois)

        # 3) winsorize q99
        q99 = float(df["rare_outlier_count"].quantile(0.99))
        d2 = df.copy()
        d2["count_clip_q99"] = np.minimum(d2["rare_outlier_count"], q99)
        tab_clip = fit_poisson_offset(d2, y_col="count_clip_q99", exposure_col="tested_loci_total", cov_cols=base_cov)
        if not tab_clip.empty:
            tab_clip.insert(0, "_model", f"poisson_offset_count_clip_q99_{q99:g}")
            out_models.append(tab_clip)

        # 4) trim thresholds
        thr_list = []
        for x in args.trim_thresholds.split(","):
            x = x.strip()
            if x:
                thr_list.append(int(x))
        for thr in thr_list:
            d3 = df[df["rare_outlier_count"] < thr].copy()
            if len(d3) < 1000:
                continue
            tab_trim = fit_poisson_offset(d3, y_col="rare_outlier_count", exposure_col="tested_loci_total", cov_cols=base_cov)
            if not tab_trim.empty:
                tab_trim.insert(0, "_model", f"poisson_offset_count_trim_lt_{thr}")
                out_models.append(tab_trim)

    if out_models:
        allm = pd.concat(out_models, ignore_index=True)
        allm.to_csv(models, sep="\t", index=False)
    else:
        pd.DataFrame({"_model": [], "term": [], "coef": [], "se": [], "p": []}).to_csv(models, sep="\t", index=False)

    # tables
    gs_df2 = gs_df.copy()
    gs_df2.insert(0, "_table", "group_summary")
    smd_df2 = smd_df.copy()
    smd_df2.insert(0, "_table", "smd")
    pd.concat([gs_df2, smd_df2], ignore_index=True).to_csv(tables, sep="\t", index=False)

    # report
    with report.open("w") as w:
        w.write(f"[{ts()}] QC/Sensitivity report\n")
        w.write(f"Input per_sample: {per}\n")
        w.write(f"N complete-case: {len(df)}\n")
        w.write(f"Group counts:\n{df['Group'].value_counts().to_string()}\n\n")

        w.write("== Group summaries ==\n")
        w.write(gs_df.to_string(index=False))
        w.write("\n\n== SMD (ASD vs Healthy / SZ vs Healthy) ==\n")
        w.write(smd_df.to_string(index=False))
        w.write("\n\n== Top samples by rare_outlier_count ==\n")
        w.write(df_top[["SampleID", "Group", "Pedigree_No", "Depth", "rare_outlier_count", "rare_any"]].to_string(index=False))
        w.write("\n\n== Models ==\n")
        if out_models:
            keep_terms = {"Group_ASD", "Group_SZ", "Depth", "Sex_M"}
            w.write(f"models.tsv: {models}\n")
            allm2 = allm.copy()
            for mname in sorted(set(allm2["_model"].tolist())):
                w.write(f"\n[{mname}]\n")
                sub = allm2[allm2["_model"] == mname].copy()
                sub = sub[sub["term"].isin(list(keep_terms))].copy()
                w.write(sub.to_string(index=False))
                w.write("\n")
        else:
            w.write("statsmodels not available; models skipped.\n")

        w.write(f"\nOutputs:\n  {report}\n  {tables}\n  {models}\n  {top}\n")
        w.write(f"\nElapsed_sec: {time.time() - t0:.1f}\n")

    print(f"[{ts()}] [DONE] Wrote report: {report}")
    print(f"[{ts()}] [DONE] Wrote tables: {tables}")
    print(f"[{ts()}] [DONE] Wrote models: {models}")
    print(f"[{ts()}] [DONE] Wrote top samples: {top}")
    print(f"[{ts()}] [DONE] Elapsed={time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
