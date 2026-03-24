#!/usr/bin/env python3
# ============================================================================
# NOTE: This script contains hardcoded default paths specific to the original
# analysis environment (NIG supercomputer). To run in a different environment,
# update the paths below or set the corresponding environment variables
# (e.g. SAMPLE_INFO, PCA_EIGENVEC, CRAM_BASE_DIR1, CRAM_BASE_DIR2).
# ============================================================================
# 20_tre_case_case_comparison_v3.py
#
# Description:
#   - Load per_sample.tsv from EHdn (v19) and STRling (v9) and run ASD vs SZ
#     direct comparison (case-case)
#   - Primary: logistic regression on rare_any (binary, >=1 carrier) -> OR (ASD vs SZ), 95%CI, P
#   - Secondary: Poisson GLM on rare_outlier_count (count) -> RR, 95%CI, P
#     Note (v2): EHdn uses observed_clusters_total, STRling uses tested_loci_total as offset
#     Note: samples with offset=0 are excluded from Poisson GLM (log(0) undefined)
#   - Auxiliary: Mann-Whitney U test (nonparametric)
#   - Heterogeneity test: evaluate homogeneity of two case-control ORs (ASD vs Ctrl, SZ vs Ctrl)
#     using Cochran's Q / Breslow-Day-like approach
#   - Covariates: Sex_M, Depth, PC1-PC10
#   - Record execution time
#
# v1->v2 changes:
#   - EHdn input changed from v18 to v19 per_sample (includes observed_clusters_total column)
#   - Poisson GLM offset dynamically detected from observed_clusters_total (EHdn) or
#     tested_loci_total (STRling)
#   - Added exclusion and logging of samples with offset=0
#   - Added exposure summary to output
#
# Output:
#   - case_case_comparison_v3.tsv: ASD vs SZ statistical results (one row per caller)
#   - case_case_comparison_v3.model_summary.txt: detailed report
#   - case_case_comparison_v3.heterogeneity.tsv: heterogeneity test results
#
# Usage:
#   Run via the top-level wrapper: crosscaller/04_tre_crosscaller_compare.sh

from __future__ import annotations

import argparse
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

try:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
except ImportError:
    sys.exit("[ERROR] statsmodels is required.")

try:
    from scipy.stats import mannwhitneyu, chi2
except ImportError:
    sys.exit("[ERROR] scipy is required.")


def ts():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def detect_col(cols, candidates):
    low = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    return None


def detect_offset_col(cols):
    """Detect offset column: observed_clusters_total (EHdn) or tested_loci_total (STRling)."""
    candidates = ["observed_clusters_total", "tested_loci_total"]
    for cand in candidates:
        if cand in cols:
            return cand
    return None


def load_per_sample(path, tool_name):
    df = pd.read_csv(path, sep="\t")
    print(f"[{ts()}] [{tool_name}] Loaded {len(df)} samples from {Path(path).name}")

    # Check for offset column
    offset_col = detect_offset_col(df.columns)
    if offset_col:
        print(f"[{ts()}] [{tool_name}] Offset column found: {offset_col}")
    else:
        print(f"[{ts()}] [{tool_name}] WARNING: No offset column found (observed_clusters_total / tested_loci_total)")

    df = df[df["Group"].isin(["ASD", "SZ"])].copy()
    print(f"[{ts()}] [{tool_name}] ASD+SZ samples: {len(df)}")
    covariates = ["Sex_M", "Depth"] + [f"PC{i}" for i in range(1, 11)]
    c_count = detect_col(df.columns, ["rare_outlier_count"])
    if c_count is None:
        raise ValueError(f"rare_outlier_count column not found in {path}")
    if "rare_any" not in df.columns:
        df["rare_any"] = (df[c_count].astype(float) >= 1).astype(int)
    required = [c_count] + covariates
    existing = [c for c in required if c in df.columns]
    len_before = len(df)
    df = df.dropna(subset=existing)
    if len(df) < len_before:
        print(f"[{ts()}] [{tool_name}] Dropped {len_before - len(df)} samples with missing covariates")
    for grp in ["ASD", "SZ"]:
        n = (df["Group"] == grp).sum()
        print(f"[{ts()}] [{tool_name}] {grp}: N={n}")

    # Report offset distribution for case-case subset
    if offset_col and offset_col in df.columns:
        for grp in ["ASD", "SZ"]:
            vals = df.loc[df["Group"] == grp, offset_col]
            n_zero = (vals == 0).sum()
            print(f"[{ts()}] [{tool_name}] {grp} {offset_col}: "
                  f"mean={vals.mean():.1f}, median={vals.median():.0f}, "
                  f"min={vals.min():.0f}, max={vals.max():.0f}, N_zero={n_zero}")

    return df


def load_per_sample_full(path, tool_name):
    """Load full dataset (Healthy + ASD + SZ) for heterogeneity test."""
    df = pd.read_csv(path, sep="\t")
    df = df[df["Group"].isin(["Healthy", "ASD", "SZ"])].copy()
    if "rare_any" not in df.columns:
        c_count = detect_col(df.columns, ["rare_outlier_count"])
        if c_count:
            df["rare_any"] = (df[c_count].astype(float) >= 1).astype(int)
    covariates = ["Sex_M", "Depth"] + [f"PC{i}" for i in range(1, 11)]
    df = df.dropna(subset=["rare_outlier_count"] + covariates)
    print(f"[{ts()}] [{tool_name}] Full dataset: {len(df)} samples "
          f"(Healthy={sum(df['Group']=='Healthy')}, "
          f"ASD={sum(df['Group']=='ASD')}, "
          f"SZ={sum(df['Group']=='SZ')})")
    return df


def run_case_case_tests(df, tool_name):
    cov_formula = "Sex_M + Depth + " + " + ".join([f"PC{i}" for i in range(1, 11)])
    df_asd = df[df["Group"] == "ASD"].copy()
    df_sz = df[df["Group"] == "SZ"].copy()
    df_sub = df.copy()
    df_sub["IsSZ"] = (df_sub["Group"] == "SZ").astype(int)

    row = {
        "Tool": tool_name,
        "Comparison": "ASD vs SZ",
        "N_ASD": len(df_asd),
        "N_SZ": len(df_sz),
        "rare_any_ASD": df_asd["rare_any"].mean(),
        "rare_any_SZ": df_sz["rare_any"].mean(),
        "mean_count_ASD": df_asd["rare_outlier_count"].mean(),
        "mean_count_SZ": df_sz["rare_outlier_count"].mean(),
    }

    # Mann-Whitney U (two-sided)
    try:
        stat, p_mwu = mannwhitneyu(
            df_asd["rare_outlier_count"].values,
            df_sz["rare_outlier_count"].values,
            alternative="two-sided")
        row["MWU_stat"] = stat
        row["MWU_P"] = p_mwu
    except Exception as e:
        print(f"[WARN] MWU failed for {tool_name}: {e}")
        row["MWU_stat"] = np.nan
        row["MWU_P"] = np.nan

    # Primary: Logistic Regression rare_any ~ IsSZ + covariates
    formula_logit = f"rare_any ~ IsSZ + {cov_formula}"
    try:
        model_logit = smf.logit(formula=formula_logit, data=df_sub).fit(
            disp=0, method="bfgs", maxiter=200)
        coef = model_logit.params["IsSZ"]
        ci = model_logit.conf_int().loc["IsSZ"]
        row["Logit_OR_SZvsASD"] = np.exp(coef)
        row["Logit_OR_CI_low"] = np.exp(ci[0])
        row["Logit_OR_CI_high"] = np.exp(ci[1])
        row["Logit_P"] = model_logit.pvalues["IsSZ"]
    except Exception as e:
        print(f"[WARN] Logistic regression failed for {tool_name}: {e}")
        row["Logit_OR_SZvsASD"] = np.nan
        row["Logit_OR_CI_low"] = np.nan
        row["Logit_OR_CI_high"] = np.nan
        row["Logit_P"] = np.nan

    # Secondary: Poisson GLM rare_outlier_count ~ IsSZ + covariates
    # v2: Use observed_clusters_total or tested_loci_total as offset
    formula_pois = f"rare_outlier_count ~ IsSZ + {cov_formula}"
    offset_col = detect_offset_col(df_sub.columns)

    try:
        if offset_col:
            # Filter out samples with offset=0 (log(0) is undefined)
            df_pois = df_sub[df_sub[offset_col] > 0].copy()
            n_excluded = len(df_sub) - len(df_pois)
            if n_excluded > 0:
                excluded_ids = df_sub[df_sub[offset_col] == 0]["SampleID"].tolist()
                print(f"[{ts()}] [{tool_name}] Poisson: excluded {n_excluded} samples "
                      f"with {offset_col}=0: {excluded_ids}")
            else:
                print(f"[{ts()}] [{tool_name}] Poisson: no samples excluded ({offset_col} all > 0)")

            offset_vals = np.log(df_pois[offset_col].astype(float).values)
            model_pois = smf.glm(
                formula=formula_pois, data=df_pois,
                family=sm.families.Poisson(),
                offset=offset_vals).fit()

            row["offset_col"] = offset_col
            row["N_Poisson"] = len(df_pois)
            row["N_excluded_zero"] = n_excluded
            row["exposure_mean"] = df_pois[offset_col].mean()
            row["exposure_median"] = df_pois[offset_col].median()
            print(f"[{ts()}] [{tool_name}] Poisson GLM with offset={offset_col} "
                  f"(N={len(df_pois)}, mean_exposure={df_pois[offset_col].mean():.1f})")
        else:
            print(f"[{ts()}] [{tool_name}] Poisson GLM WITHOUT offset (no offset column available)")
            model_pois = smf.glm(
                formula=formula_pois, data=df_sub,
                family=sm.families.Poisson()).fit()
            row["offset_col"] = "NONE"
            row["N_Poisson"] = len(df_sub)
            row["N_excluded_zero"] = 0
            row["exposure_mean"] = np.nan
            row["exposure_median"] = np.nan

        coef_p = model_pois.params["IsSZ"]
        ci_p = model_pois.conf_int().loc["IsSZ"]
        row["Poisson_RR_SZvsASD"] = np.exp(coef_p)
        row["Poisson_RR_CI_low"] = np.exp(ci_p[0])
        row["Poisson_RR_CI_high"] = np.exp(ci_p[1])
        row["Poisson_P"] = model_pois.pvalues["IsSZ"]
    except Exception as e:
        print(f"[WARN] Poisson GLM failed for {tool_name}: {e}")
        row["Poisson_RR_SZvsASD"] = np.nan
        row["Poisson_RR_CI_low"] = np.nan
        row["Poisson_RR_CI_high"] = np.nan
        row["Poisson_P"] = np.nan
        row["offset_col"] = offset_col if offset_col else "NONE"
        row["N_Poisson"] = np.nan
        row["N_excluded_zero"] = np.nan
        row["exposure_mean"] = np.nan
        row["exposure_median"] = np.nan

    return row


def run_heterogeneity_test(df_full, tool_name):
    cov_formula = "Sex_M + Depth + " + " + ".join([f"PC{i}" for i in range(1, 11)])
    df_ctrl = df_full[df_full["Group"] == "Healthy"].copy()
    results = {}
    log_ors = []
    se_log_ors = []

    for case_grp in ["ASD", "SZ"]:
        df_case = df_full[df_full["Group"] == case_grp].copy()
        df_sub = pd.concat([df_ctrl, df_case]).copy()
        df_sub["IsCase"] = (df_sub["Group"] == case_grp).astype(int)
        formula_logit = f"rare_any ~ IsCase + {cov_formula}"
        try:
            model = smf.logit(formula=formula_logit, data=df_sub).fit(
                disp=0, method="bfgs", maxiter=200)
            log_or = model.params["IsCase"]
            se = model.bse["IsCase"]
            log_ors.append(log_or)
            se_log_ors.append(se)
            results[f"log_OR_{case_grp}"] = log_or
            results[f"SE_{case_grp}"] = se
            results[f"OR_{case_grp}"] = np.exp(log_or)
        except Exception as e:
            print(f"[WARN] Heterogeneity: logistic failed for {case_grp}: {e}")
            return {"Tool": tool_name, "Q_stat": np.nan, "Q_P": np.nan, "I2_pct": np.nan}

    if len(log_ors) != 2:
        return {"Tool": tool_name, "Q_stat": np.nan, "Q_P": np.nan, "I2_pct": np.nan}

    log_ors = np.array(log_ors)
    se_log_ors = np.array(se_log_ors)
    weights = 1.0 / (se_log_ors ** 2)
    weighted_mean = np.sum(weights * log_ors) / np.sum(weights)
    Q = np.sum(weights * (log_ors - weighted_mean) ** 2)
    df_Q = len(log_ors) - 1
    p_Q = 1 - chi2.cdf(Q, df_Q)
    I2 = max(0, (Q - df_Q) / Q * 100) if Q > 0 else 0.0

    results["Tool"] = tool_name
    results["weighted_mean_log_OR"] = weighted_mean
    results["weighted_mean_OR"] = np.exp(weighted_mean)
    results["Q_stat"] = Q
    results["Q_df"] = df_Q
    results["Q_P"] = p_Q
    results["I2_pct"] = I2
    return results


def main():
    start_time = time.time()
    print(f"[{ts()}] === TRE Case-Case Comparison (ASD vs SZ) v3 ===")
    print(f"[{ts()}] v2→v3: EHdn Poisson GLM uses offset=log(observed_clusters_total)")

    repo_root = Path(__file__).resolve().parents[2]

    ap = argparse.ArgumentParser()
    # v2: EHdn uses v19 per_sample (with observed_clusters_total)
    ap.add_argument("--ehdn_per_sample",
                    default=str(repo_root / "analysis_results_novel" / "outlier_burden_rare_crossfit_v19.per_sample.tsv"))
    ap.add_argument("--strling_per_sample",
                    default=str(repo_root / "analysis_results_strling" / "strling_outlier_burden_rare_crossfit_inbounds_v9.per_sample.tsv"))
    ap.add_argument("--outdir", default=str(repo_root / "crosscaller_results"))
    ap.add_argument("--prefix", default="case_case_comparison_v3")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # === Part 1: ASD vs SZ Direct Comparison ===
    print(f"\n[{ts()}] === Part 1: ASD vs SZ Direct Comparison ===")
    all_results = []

    ehdn_path = Path(args.ehdn_per_sample)
    if ehdn_path.exists():
        df_ehdn = load_per_sample(ehdn_path, "EHdn")
        res_ehdn = run_case_case_tests(df_ehdn, "EHdn")
        all_results.append(res_ehdn)
    else:
        print(f"[WARN] EHdn per_sample not found: {ehdn_path}")

    strling_path = Path(args.strling_per_sample)
    if strling_path.exists():
        df_strling = load_per_sample(strling_path, "STRling")
        res_strling = run_case_case_tests(df_strling, "STRling")
        all_results.append(res_strling)
    else:
        print(f"[WARN] STRling per_sample not found: {strling_path}")

    if not all_results:
        sys.exit("[ERROR] No data loaded.")

    res_df = pd.DataFrame(all_results)
    out_tsv = outdir / f"{args.prefix}.tsv"
    res_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\n[{ts()}] Wrote {out_tsv}")

    # === Part 2: Heterogeneity ===
    print(f"\n[{ts()}] === Part 2: Heterogeneity of Case-Control ORs (ASD vs SZ) ===")
    het_results = []

    if ehdn_path.exists():
        df_ehdn_full = load_per_sample_full(ehdn_path, "EHdn")
        het_ehdn = run_heterogeneity_test(df_ehdn_full, "EHdn")
        het_results.append(het_ehdn)

    if strling_path.exists():
        df_strling_full = load_per_sample_full(strling_path, "STRling")
        het_strling = run_heterogeneity_test(df_strling_full, "STRling")
        het_results.append(het_strling)

    if het_results:
        het_df = pd.DataFrame(het_results)
        het_tsv = outdir / f"{args.prefix}.heterogeneity.tsv"
        het_df.to_csv(het_tsv, sep="\t", index=False)
        print(f"[{ts()}] Wrote {het_tsv}")

    # === Part 3: Report ===
    report_path = outdir / f"{args.prefix}.model_summary.txt"
    with report_path.open("w") as w:
        w.write(f"[{ts()}] TRE Case-Case Comparison: ASD vs SZ (v3)\n")
        w.write(f"[{ts()}] v2→v3: EHdn Poisson GLM with offset=log(observed_clusters_total)\n")
        w.write("=" * 70 + "\n\n")
        w.write("RATIONALE:\n")
        w.write("  The manuscript demonstrates significant TRE burden enrichment in both\n")
        w.write("  ASD and SZ vs controls. To formally support the 'shared burden' claim,\n")
        w.write("  this analysis directly tests whether ASD and SZ differ from each other.\n")
        w.write("  A non-significant case-case comparison indicates no evidence of\n")
        w.write("  heterogeneity, strengthening the cross-disorder interpretation.\n\n")
        w.write("=" * 70 + "\n")
        w.write("PART 1: DIRECT ASD vs SZ COMPARISON\n")
        w.write("  (OR/RR > 1 means SZ > ASD; OR/RR < 1 means ASD > SZ)\n")
        w.write("=" * 70 + "\n\n")
        for _, r in res_df.iterrows():
            tool = r["Tool"]
            w.write(f"--- {tool} ---\n")
            w.write(f"  N_ASD={r['N_ASD']}, N_SZ={r['N_SZ']}\n")
            w.write(f"  rare_any: ASD={r['rare_any_ASD']:.3f}, SZ={r['rare_any_SZ']:.3f}\n")
            w.write(f"  mean_count: ASD={r['mean_count_ASD']:.4f}, SZ={r['mean_count_SZ']:.4f}\n")
            w.write(f"  Poisson offset column: {r.get('offset_col', 'N/A')}\n")
            if not pd.isna(r.get('exposure_mean')):
                w.write(f"  Exposure mean: {r['exposure_mean']:.1f}, median: {r['exposure_median']:.0f}\n")
            w.write(f"  N in Poisson model: {r.get('N_Poisson', 'N/A')}, "
                    f"N excluded (zero): {r.get('N_excluded_zero', 'N/A')}\n\n")
            w.write(f"  Primary (Logistic, rare_any):\n")
            or_val = r.get("Logit_OR_SZvsASD", np.nan)
            ci_lo = r.get("Logit_OR_CI_low", np.nan)
            ci_hi = r.get("Logit_OR_CI_high", np.nan)
            p_val = r.get("Logit_P", np.nan)
            w.write(f"    OR(SZ vs ASD) = {or_val:.4f} (95%CI {ci_lo:.4f}-{ci_hi:.4f}), P = {p_val:.6e}\n")
            if not np.isnan(p_val) and p_val > 0.05:
                w.write(f"    -> NOT significant (P > 0.05): No evidence ASD != SZ\n")
            elif not np.isnan(p_val):
                w.write(f"    -> Significant (P < 0.05): Evidence of difference between ASD and SZ\n")
            w.write("\n")
            w.write(f"  Secondary (Poisson GLM, count):\n")
            rr_val = r.get("Poisson_RR_SZvsASD", np.nan)
            ci_lo_p = r.get("Poisson_RR_CI_low", np.nan)
            ci_hi_p = r.get("Poisson_RR_CI_high", np.nan)
            p_val_p = r.get("Poisson_P", np.nan)
            w.write(f"    RR(SZ vs ASD) = {rr_val:.4f} (95%CI {ci_lo_p:.4f}-{ci_hi_p:.4f}), P = {p_val_p:.6e}\n\n")
            w.write(f"  Mann-Whitney U (two-sided):\n")
            w.write(f"    P = {r.get('MWU_P', np.nan):.6e}\n\n")
        w.write("=" * 70 + "\n")
        w.write("PART 2: HETEROGENEITY OF CASE-CONTROL ORs\n")
        w.write("  Cochran's Q test: Are the ASD-vs-Ctrl and SZ-vs-Ctrl log(OR)s\n")
        w.write("  significantly different from each other?\n")
        w.write("=" * 70 + "\n\n")
        if het_results:
            for h in het_results:
                tool = h["Tool"]
                w.write(f"--- {tool} ---\n")
                w.write(f"  OR(ASD vs Ctrl) = {h.get('OR_ASD', np.nan):.4f}\n")
                w.write(f"  OR(SZ vs Ctrl)  = {h.get('OR_SZ', np.nan):.4f}\n")
                w.write(f"  Weighted mean OR = {h.get('weighted_mean_OR', np.nan):.4f}\n")
                w.write(f"  Cochran's Q = {h.get('Q_stat', np.nan):.4f}, df={h.get('Q_df', 1)}, "
                        f"P = {h.get('Q_P', np.nan):.6e}\n")
                w.write(f"  I2 = {h.get('I2_pct', np.nan):.1f}%\n")
                q_p = h.get("Q_P", np.nan)
                if not np.isnan(q_p) and q_p > 0.05:
                    w.write(f"  -> No significant heterogeneity (P > 0.05)\n")
                    w.write(f"    Supports the 'shared burden' interpretation.\n")
                elif not np.isnan(q_p):
                    w.write(f"  -> Significant heterogeneity (P < 0.05)\n")
                w.write("\n")
        w.write("=" * 70 + "\n")
        w.write("INTERPRETATION GUIDE:\n")
        w.write("  - If both Part 1 (direct) and Part 2 (heterogeneity) are\n")
        w.write("    non-significant -> strong support for 'shared TRE burden'\n")
        w.write("  - Non-significant does NOT mean 'equal'; it means 'no evidence of\n")
        w.write("    difference' given the sample size.\n")
        w.write("=" * 70 + "\n")
    print(f"[{ts()}] Wrote {report_path}")

    # Console summary
    print(f"\n{'=' * 60}")
    print("=== SUMMARY: ASD vs SZ Direct Comparison (v3) ===")
    print(f"{'=' * 60}")
    for _, r in res_df.iterrows():
        tool = r["Tool"]
        print(f"\n[{tool}] (offset={r.get('offset_col', 'N/A')})")
        or_val = r.get("Logit_OR_SZvsASD", np.nan)
        ci_lo = r.get("Logit_OR_CI_low", np.nan)
        ci_hi = r.get("Logit_OR_CI_high", np.nan)
        p_val = r.get("Logit_P", np.nan)
        print(f"  Logistic: OR(SZ vs ASD) = {or_val:.2f} "
              f"(95%CI {ci_lo:.2f}-{ci_hi:.2f}), P = {p_val:.2e}")
        rr_val = r.get("Poisson_RR_SZvsASD", np.nan)
        ci_lo_p = r.get("Poisson_RR_CI_low", np.nan)
        ci_hi_p = r.get("Poisson_RR_CI_high", np.nan)
        p_val_p = r.get("Poisson_P", np.nan)
        print(f"  Poisson:  RR(SZ vs ASD) = {rr_val:.2f} "
              f"(95%CI {ci_lo_p:.2f}-{ci_hi_p:.2f}), P = {p_val_p:.2e}")
        print(f"  MWU:      P = {r.get('MWU_P', np.nan):.2e}")

    if het_results:
        print(f"\n{'=' * 60}")
        print("=== SUMMARY: Heterogeneity of Case-Control ORs ===")
        print(f"{'=' * 60}")
        for h in het_results:
            tool = h["Tool"]
            print(f"\n[{tool}]")
            print(f"  OR(ASD vs Ctrl) = {h.get('OR_ASD', np.nan):.2f}, "
                  f"OR(SZ vs Ctrl) = {h.get('OR_SZ', np.nan):.2f}")
            print(f"  Cochran's Q = {h.get('Q_stat', np.nan):.3f}, "
                  f"P = {h.get('Q_P', np.nan):.2e}, "
                  f"I2 = {h.get('I2_pct', np.nan):.1f}%")

    elapsed = time.time() - start_time
    print(f"\n[{ts()}] [DONE] Total elapsed: {elapsed:.1f}s ({elapsed / 60:.1f}min)")


if __name__ == "__main__":
    main()
