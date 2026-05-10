#!/usr/bin/env python3
# compute_strict_loco_top_exclusions_v1.py
#
# Purpose:
#   R18 round Phase A+++ : aggregate Diff_any architecture-level primary
#   endpoint に対する 3 種類の追加 sensitivity analysis を 1 スクリプトに統合する。
#
#   B2. Strict LOCO sensitivity (10 fits):
#       LOCO_k strict = Diff_any union から membership_k=1 の bin を全削除した
#       bin set (他 class にも属していても削除)
#       「single L2 class dependence を否定」を強く言うための解析。
#       Standard LOCO は class-unique contributions のみ削除。Strict LOCO は
#       class-participating bins を全削除する。
#
#   B3. Top-N sample exclusion sensitivity (N=1, 2, 3, 5):
#       Per-individual Diff_any DEL exposure を計算し、降順 sort で上位 N 個の
#       sample を除外して ASD-DEL B' logistic を再走。
#       「single sample dependence を否定」のための解析。
#
#   B4. Top-N event exclusion sensitivity (N=1, 2, 3, 5):
#       Per-event "n_bins_disrupted in Diff_any union" を計算し、降順 sort で
#       上位 N 個の event を overlap table から除外、per-individual exposure を
#       再計算して ASD-DEL B' logistic を再走。
#       「single event dependence を否定」のための解析。
#
#   全解析共通:
#     - Step 03 v9 Global flag を継承 (subset-level 99% filter は適用しない)
#     - B' covariates: Sex_numeric + PC1-10 + log1p_total_del_bases + log1p_total_gene_DEL
#     - ASD vs Healthy (case=ASD=1, ctrl=Healthy=0)
#
# Inputs:
#   /home/kushima-pg/tad04212026/02_bin_l2_annotation/output_v2/
#     bin_l2_annotation_v2.tsv.gz  (4980-bin universe; bin_id + 10 membership_X cols)
#   /home/kushima-pg/tad04212026/04_wgs_sv_boundary_overlap/output_v10/
#     sample_boundary_event_overlap_v10.tsv.gz  (sample x event x bin overlap)
#   /home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v3/
#     sample_burden_L2_and_specificity_v3.tsv  (covariates + Diagnosis)
#
# Outputs:
#   /home/kushima-pg/tad04212026/16_aggregate_primary_sensitivity_v1/output/
#     strict_loco_v1.tsv                   (10 strict LOCO ORs + baseline)
#     top_sample_exclusion_v1.tsv          (baseline + top-1/2/3/5 sample drops)
#     top_event_exclusion_v1.tsv           (baseline + top-1/2/3/5 event drops)
#     strict_loco_top_exclusions_v1.log    (combined log)
#
# Run:
#   cd /home/kushima-pg/tad04212026/16_aggregate_primary_sensitivity_v1/
#   python3 compute_strict_loco_top_exclusions_v1.py
#
# Elapsed time logged.
import os
import sys
import time
import pandas as pd
import numpy as np
import statsmodels.api as sm

t0 = time.time()

BIN_ANNOT = "/home/kushima-pg/tad04212026/02_bin_l2_annotation/output_v2/bin_l2_annotation_v2.tsv.gz"
OVERLAP_TSV = "/home/kushima-pg/tad04212026/04_wgs_sv_boundary_overlap/output_v10/sample_boundary_event_overlap_v10.tsv.gz"
BURDEN_TSV = "/home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v3/sample_burden_L2_and_specificity_v3.tsv"

OUTDIR = "/home/kushima-pg/tad04212026/16_aggregate_primary_sensitivity_v1/output"
os.makedirs(OUTDIR, exist_ok=True)
OUT_STRICT = os.path.join(OUTDIR, "strict_loco_v1.tsv")
OUT_SAMPLE = os.path.join(OUTDIR, "top_sample_exclusion_v1.tsv")
OUT_EVENT = os.path.join(OUTDIR, "top_event_exclusion_v1.tsv")
OUT_LOG = os.path.join(OUTDIR, "strict_loco_top_exclusions_v1.log")

PRIMARY_BOUNDARY_CLASSES = [
    "HPC_Exc-CA", "HPC_Exc-DG", "HPC_Exc-ENT",
    "HPC_Inh-CGE", "HPC_Inh-MGE",
    "PFC_Astro", "PFC_Exc-DL", "PFC_Exc-UL",
    "PFC_Inh-CGE", "PFC_Inh-MGE",
]
COVS = ["Sex_numeric"] + ["PC%d" % i for i in range(1, 11)] + ["log1p_total_del_bases", "log1p_total_gene_DEL"]

TOP_N_LIST = [1, 2, 3, 5]


def log(msg, fh=None):
    line = "[%.1fs] %s" % (time.time() - t0, msg)
    print(line)
    if fh is not None:
        fh.write(line + "\n")
        fh.flush()


def fit_logistic_generic(exposure_series, sub_filt, label="exposure"):
    X = pd.DataFrame({label: exposure_series.values})
    for c in COVS:
        X[c] = sub_filt[c].values
    X = sm.add_constant(X)
    y = sub_filt["is_case"].values
    try:
        m = sm.Logit(y, X).fit(disp=0, maxiter=200)
        beta = float(m.params[label])
        se = float(m.bse[label])
        return {
            "beta": beta, "SE": se,
            "OR": float(np.exp(beta)),
            "ci_lo": float(np.exp(beta - 1.96 * se)),
            "ci_hi": float(np.exp(beta + 1.96 * se)),
            "wald_p": float(m.pvalues[label]),
            "status": "ok",
        }
    except Exception as e:
        return {"status": "error", "error": str(e)}


def compute_per_sample_burden(bin_set, ov_del, samples_index):
    """Per-sample count of unique bins in bin_set disrupted by DEL events."""
    ov_filt = ov_del[ov_del["bin_id"].isin(bin_set)]
    per_sample_bins = ov_filt.groupby("sample_id")["bin_id"].nunique()
    return per_sample_bins.reindex(samples_index, fill_value=0).astype(int)


with open(OUT_LOG, "w") as logfh:
    log("=" * 70, logfh)
    log("compute_strict_loco_top_exclusions_v1.py", logfh)
    log("(B2: strict LOCO + B3: top-N sample + B4: top-N event)", logfh)
    log("Step 03 v9 Global flag inherited; no subset-level 99% filter.", logfh)
    log("=" * 70, logfh)

    # --- Load bin_l2_annotation_v2 ---
    log("Loading bin_l2_annotation_v2: %s" % BIN_ANNOT, logfh)
    anno = pd.read_csv(BIN_ANNOT, sep="\t", compression="gzip", low_memory=False,
                       dtype={"bin_id": "string"})
    log("  annotation shape: %s" % str(anno.shape), logfh)

    membership_cols = ["membership_%s" % c for c in PRIMARY_BOUNDARY_CLASSES]
    missing = [c for c in membership_cols if c not in anno.columns]
    if missing:
        sys.exit("ERROR: missing membership columns: %s" % missing)

    diff_any_full = anno[anno[membership_cols].sum(axis=1) > 0].copy()
    full_bin_set = set(diff_any_full["bin_id"].tolist())
    log("  Full Diff_any union (10 classes): %d bins" % len(full_bin_set), logfh)
    assert len(full_bin_set) == 4980, "Expected 4,980 Diff_any bins; got %d" % len(full_bin_set)

    # Strict LOCO bin sets: drop ALL bins with membership_k=1
    strict_loco_sets = {}
    for k in PRIMARY_BOUNDARY_CLASSES:
        col_k = "membership_%s" % k
        # bins where (membership_k = 0) AND (any of other 9 classes = 1)
        other_cols = [c for c in membership_cols if c != col_k]
        mask = (anno[col_k] == 0) & (anno[other_cols].sum(axis=1) > 0)
        sl_bins = anno[mask]["bin_id"].tolist()
        strict_loco_sets[k] = set(sl_bins)
        n_sl = len(strict_loco_sets[k])
        n_dropped = len(full_bin_set) - n_sl
        log("  STRICT_LOCO_%s: %d bins (-%d from full)" % (k, n_sl, n_dropped), logfh)

    # --- Load Step 4 v10 overlap (DEL only) ---
    log("Loading Step 4 v10 overlap table: %s" % OVERLAP_TSV, logfh)
    ov = pd.read_csv(OVERLAP_TSV, sep="\t", compression="gzip", low_memory=False,
                     dtype={"sample_id": "string", "bin_id": "string", "event_id": "string"})
    log("  overlap shape: %s" % str(ov.shape), logfh)
    ov_del = ov[ov["sv_type_norm"] == "DEL"].copy()
    log("  DEL-filtered overlap: %d rows" % len(ov_del), logfh)

    # --- Load Step 5 v3 burden TSV ---
    log("Loading Step 5 v3 burden TSV: %s" % BURDEN_TSV, logfh)
    burden = pd.read_csv(BURDEN_TSV, sep="\t", low_memory=False)
    log("  burden shape: %s" % str(burden.shape), logfh)

    sub = burden[burden["Diagnosis"].isin(["ASD", "Healthy"])].copy()
    sub["is_case"] = (sub["Diagnosis"] == "ASD").astype(int)
    sub_filt = sub.dropna(subset=COVS).copy()
    log("  After covariates dropna [Step 03 v9 Global flag]: N=%d (ASD=%d, Healthy=%d)" % (
        len(sub_filt), int((sub_filt["is_case"] == 1).sum()),
        int((sub_filt["is_case"] == 0).sum())
    ), logfh)

    samples_idx = pd.Index(sub_filt["sample_id"])
    analytic_sample_set = set(sub_filt["sample_id"].tolist())

    # Restrict overlap to analytic samples (avoid contamination from non-analytic samples)
    ov_del_analytic = ov_del[ov_del["sample_id"].isin(analytic_sample_set)].copy()
    log("  ov_del restricted to analytic samples: %d rows (%d unique events, %d unique samples)" % (
        len(ov_del_analytic),
        ov_del_analytic["event_id"].nunique(),
        ov_del_analytic["sample_id"].nunique()
    ), logfh)

    # --- BASELINE: full Diff_any union ASD-DEL B' logistic ---
    log("=" * 70, logfh)
    log("BASELINE: full Diff_any union ASD-DEL B' logistic", logfh)
    full_burden = compute_per_sample_burden(full_bin_set, ov_del_analytic, samples_idx)
    sub_filt["Diff_any_full"] = full_burden.values
    n_exp_case = int(((sub_filt["is_case"] == 1) & (sub_filt["Diff_any_full"] >= 1)).sum())
    n_exp_ctrl = int(((sub_filt["is_case"] == 0) & (sub_filt["Diff_any_full"] >= 1)).sum())
    log("  Baseline: exp_case=%d, exp_ctrl=%d" % (n_exp_case, n_exp_ctrl), logfh)
    baseline_res = fit_logistic_generic(sub_filt["Diff_any_full"], sub_filt)
    log("  Baseline OR=%.3f [%.3f-%.3f], Wald P=%.3e" % (
        baseline_res["OR"], baseline_res["ci_lo"], baseline_res["ci_hi"], baseline_res["wald_p"]
    ), logfh)

    # ============================================================
    # B2. STRICT LOCO (10 fits)
    # ============================================================
    log("=" * 70, logfh)
    log("B2. STRICT LOCO: dropping ALL bins with membership_k=1 (10 fits)", logfh)
    log("=" * 70, logfh)
    strict_results = [{
        "loco_dropped_class": "(none; full Diff_any baseline)",
        "n_bins": len(full_bin_set),
        "n_dropped": 0,
        "n_exp_case": n_exp_case, "n_exp_ctrl": n_exp_ctrl,
        **baseline_res,
    }]
    for k in PRIMARY_BOUNDARY_CLASSES:
        log("STRICT_LOCO drop %s ..." % k, logfh)
        sl_burden = compute_per_sample_burden(strict_loco_sets[k], ov_del_analytic, samples_idx)
        sub_filt["Diff_any_strict_loco"] = sl_burden.values
        n_exp_case_k = int(((sub_filt["is_case"] == 1) & (sub_filt["Diff_any_strict_loco"] >= 1)).sum())
        n_exp_ctrl_k = int(((sub_filt["is_case"] == 0) & (sub_filt["Diff_any_strict_loco"] >= 1)).sum())
        sl_res = fit_logistic_generic(sub_filt["Diff_any_strict_loco"], sub_filt)
        log("  STRICT_LOCO %s: %d bins (-%d), exp_case=%d, exp_ctrl=%d, OR=%.3f [%.3f-%.3f], Wald P=%.3e" % (
            k, len(strict_loco_sets[k]), len(full_bin_set) - len(strict_loco_sets[k]),
            n_exp_case_k, n_exp_ctrl_k,
            sl_res["OR"], sl_res["ci_lo"], sl_res["ci_hi"], sl_res["wald_p"]
        ), logfh)
        strict_results.append({
            "loco_dropped_class": k,
            "n_bins": len(strict_loco_sets[k]),
            "n_dropped": len(full_bin_set) - len(strict_loco_sets[k]),
            "n_exp_case": n_exp_case_k, "n_exp_ctrl": n_exp_ctrl_k,
            **sl_res,
        })

    strict_df = pd.DataFrame(strict_results)
    strict_df.to_csv(OUT_STRICT, sep="\t", index=False, na_rep="NA")
    log("Saved strict LOCO: %s" % OUT_STRICT, logfh)
    sl_only = strict_df[strict_df["loco_dropped_class"] != "(none; full Diff_any baseline)"]
    if not sl_only.empty:
        sl_ors = sl_only["OR"].dropna()
        log("  10 STRICT_LOCO ORs: min=%.3f, max=%.3f, mean=%.3f, range=%.3f" % (
            sl_ors.min(), sl_ors.max(), sl_ors.mean(), sl_ors.max() - sl_ors.min()
        ), logfh)
        log("  ±%.1f%% of baseline OR=%.3f" % (
            100 * max(abs(sl_ors.max() - baseline_res["OR"]), abs(sl_ors.min() - baseline_res["OR"])) / baseline_res["OR"],
            baseline_res["OR"]
        ), logfh)

    # ============================================================
    # B3. TOP-N SAMPLE EXCLUSION (N = 1, 2, 3, 5)
    # ============================================================
    log("=" * 70, logfh)
    log("B3. TOP-N SAMPLE EXCLUSION (N=1, 2, 3, 5)", logfh)
    log("=" * 70, logfh)

    # Identify top-burden samples by Diff_any_full
    sub_filt_sorted = sub_filt.sort_values("Diff_any_full", ascending=False).reset_index(drop=True)
    top_samples_ranked = sub_filt_sorted[["sample_id", "Diagnosis", "Diff_any_full"]].head(10)
    log("Top 10 high-burden samples by Diff_any DEL count:", logfh)
    for i, row in top_samples_ranked.iterrows():
        log("  rank %d: %s (%s, Diff_any=%d)" % (
            i + 1, row["sample_id"], row["Diagnosis"], int(row["Diff_any_full"])
        ), logfh)

    sample_results = [{
        "scenario": "baseline (no exclusion)",
        "n_excluded": 0, "excluded_samples": "",
        "n_remaining": len(sub_filt),
        "n_exp_case": n_exp_case, "n_exp_ctrl": n_exp_ctrl,
        **baseline_res,
    }]
    for n_drop in TOP_N_LIST:
        excluded = top_samples_ranked["sample_id"].head(n_drop).tolist()
        sub_drop = sub_filt[~sub_filt["sample_id"].isin(excluded)].copy()
        # Re-fit with remaining samples
        X = pd.DataFrame({"exposure": sub_drop["Diff_any_full"].values})
        for c in COVS:
            X[c] = sub_drop[c].values
        X = sm.add_constant(X)
        y = sub_drop["is_case"].values
        n_e_case = int(((sub_drop["is_case"] == 1) & (sub_drop["Diff_any_full"] >= 1)).sum())
        n_e_ctrl = int(((sub_drop["is_case"] == 0) & (sub_drop["Diff_any_full"] >= 1)).sum())
        try:
            m = sm.Logit(y, X).fit(disp=0, maxiter=200)
            beta = float(m.params["exposure"])
            se = float(m.bse["exposure"])
            res = {
                "beta": beta, "SE": se,
                "OR": float(np.exp(beta)),
                "ci_lo": float(np.exp(beta - 1.96 * se)),
                "ci_hi": float(np.exp(beta + 1.96 * se)),
                "wald_p": float(m.pvalues["exposure"]),
                "status": "ok",
            }
        except Exception as e:
            res = {"status": "error", "error": str(e)}
        log("Top-%d samples excluded: %s" % (n_drop, ", ".join(excluded)), logfh)
        if res["status"] == "ok":
            log("  N=%d, exp_case=%d, exp_ctrl=%d, OR=%.3f [%.3f-%.3f], Wald P=%.3e" % (
                len(sub_drop), n_e_case, n_e_ctrl,
                res["OR"], res["ci_lo"], res["ci_hi"], res["wald_p"]
            ), logfh)
        sample_results.append({
            "scenario": "top-%d sample exclusion" % n_drop,
            "n_excluded": n_drop,
            "excluded_samples": ";".join(excluded),
            "n_remaining": len(sub_drop),
            "n_exp_case": n_e_case, "n_exp_ctrl": n_e_ctrl,
            **res,
        })

    sample_df = pd.DataFrame(sample_results)
    sample_df.to_csv(OUT_SAMPLE, sep="\t", index=False, na_rep="NA")
    log("Saved top-sample exclusion: %s" % OUT_SAMPLE, logfh)

    # ============================================================
    # B4. TOP-N EVENT EXCLUSION (N = 1, 2, 3, 5)
    # ============================================================
    log("=" * 70, logfh)
    log("B4. TOP-N EVENT EXCLUSION (N=1, 2, 3, 5)", logfh)
    log("=" * 70, logfh)

    # Identify events disrupting Diff_any bins, count unique bins per event
    ov_diff = ov_del_analytic[ov_del_analytic["bin_id"].isin(full_bin_set)].copy()
    event_bin_counts = ov_diff.groupby("event_id").agg(
        n_bins_disrupted=("bin_id", "nunique"),
        sample_id=("sample_id", "first"),
    ).sort_values("n_bins_disrupted", ascending=False).reset_index()
    log("Top 10 events by # Diff_any bins disrupted:", logfh)
    for i, row in event_bin_counts.head(10).iterrows():
        log("  rank %d: event=%s sample=%s n_bins=%d" % (
            i + 1, row["event_id"], row["sample_id"], int(row["n_bins_disrupted"])
        ), logfh)

    event_results = [{
        "scenario": "baseline (no exclusion)",
        "n_excluded": 0, "excluded_events": "",
        "n_exp_case": n_exp_case, "n_exp_ctrl": n_exp_ctrl,
        **baseline_res,
    }]
    for n_drop in TOP_N_LIST:
        excluded_events = event_bin_counts["event_id"].head(n_drop).tolist()
        ov_drop = ov_del_analytic[~ov_del_analytic["event_id"].isin(excluded_events)].copy()
        new_burden = compute_per_sample_burden(full_bin_set, ov_drop, samples_idx)
        sub_filt["Diff_any_event_drop"] = new_burden.values
        n_e_case = int(((sub_filt["is_case"] == 1) & (sub_filt["Diff_any_event_drop"] >= 1)).sum())
        n_e_ctrl = int(((sub_filt["is_case"] == 0) & (sub_filt["Diff_any_event_drop"] >= 1)).sum())
        ev_res = fit_logistic_generic(sub_filt["Diff_any_event_drop"], sub_filt)
        log("Top-%d events excluded: %s" % (n_drop, ";".join(excluded_events)), logfh)
        if ev_res["status"] == "ok":
            log("  exp_case=%d, exp_ctrl=%d, OR=%.3f [%.3f-%.3f], Wald P=%.3e" % (
                n_e_case, n_e_ctrl,
                ev_res["OR"], ev_res["ci_lo"], ev_res["ci_hi"], ev_res["wald_p"]
            ), logfh)
        event_results.append({
            "scenario": "top-%d event exclusion" % n_drop,
            "n_excluded": n_drop,
            "excluded_events": ";".join(excluded_events),
            "n_exp_case": n_e_case, "n_exp_ctrl": n_e_ctrl,
            **ev_res,
        })

    event_df = pd.DataFrame(event_results)
    event_df.to_csv(OUT_EVENT, sep="\t", index=False, na_rep="NA")
    log("Saved top-event exclusion: %s" % OUT_EVENT, logfh)

    # ============================================================
    # FINAL SUMMARY
    # ============================================================
    log("=" * 70, logfh)
    log("FINAL SUMMARY (all 3 sensitivity analyses)", logfh)
    log("=" * 70, logfh)
    log("Baseline ASD-DEL: OR=%.3f [%.3f-%.3f], Wald P=%.3e (full Diff_any union, N=%d)" % (
        baseline_res["OR"], baseline_res["ci_lo"], baseline_res["ci_hi"], baseline_res["wald_p"],
        len(sub_filt)
    ), logfh)
    log("", logfh)
    if not sl_only.empty:
        log("B2. STRICT LOCO (10 fits): OR range %.3f-%.3f (mean %.3f, ±%.1f%% of baseline)" % (
            sl_ors.min(), sl_ors.max(), sl_ors.mean(),
            100 * max(abs(sl_ors.max() - baseline_res["OR"]), abs(sl_ors.min() - baseline_res["OR"])) / baseline_res["OR"]
        ), logfh)
    log("", logfh)
    sample_only = sample_df[sample_df["scenario"] != "baseline (no exclusion)"]
    if not sample_only.empty:
        s_ors = sample_only["OR"].dropna()
        log("B3. TOP-N SAMPLE EXCLUSION: OR range %.3f-%.3f (across N=1,2,3,5)" % (
            s_ors.min(), s_ors.max()
        ), logfh)
        log("   max deviation from baseline: ±%.1f%%" % (
            100 * max(abs(s_ors.max() - baseline_res["OR"]), abs(s_ors.min() - baseline_res["OR"])) / baseline_res["OR"]
        ), logfh)
    log("", logfh)
    event_only = event_df[event_df["scenario"] != "baseline (no exclusion)"]
    if not event_only.empty:
        e_ors = event_only["OR"].dropna()
        log("B4. TOP-N EVENT EXCLUSION: OR range %.3f-%.3f (across N=1,2,3,5)" % (
            e_ors.min(), e_ors.max()
        ), logfh)
        log("   max deviation from baseline: ±%.1f%%" % (
            100 * max(abs(e_ors.max() - baseline_res["OR"]), abs(e_ors.min() - baseline_res["OR"])) / baseline_res["OR"]
        ), logfh)
    log("[Total elapsed: %.1fs]" % (time.time() - t0), logfh)
