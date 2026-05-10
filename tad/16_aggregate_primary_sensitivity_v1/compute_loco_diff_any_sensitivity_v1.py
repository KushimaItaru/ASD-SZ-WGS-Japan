#!/usr/bin/env python3
# compute_loco_diff_any_sensitivity_v1.py
#
# Purpose:
#   Leave-One-Class-Out (LOCO) Diff_any union sensitivity for R18 aggregate primary.
#   Test whether aggregate Diff_any signal is driven by a single L2 class.
#
#   For each LOCO_k (k=1..10):
#     LOCO_k bin set = bins where any of (10 PRIMARY_BOUNDARY_CLASSES \ {k}) has membership=1
#     Per-sample exposure = LOCO_k bins disrupted by rare DEL events
#     Outcome: ASD vs Healthy, B' logistic regression (Sex_numeric + PC1-10 +
#              log1p_total_del_bases + log1p_total_gene_DEL)
#     Report: OR + 95% CI + Wald P
#
#   Baseline: full Diff_any union (all 10 classes), same pipeline.
#   Expectation: 10 LOCO ORs within tight range of baseline -> single-class
#                dependence ruled out.
#
# Inputs:
#   /home/kushima-pg/tad04212026/02_bin_l2_annotation/output_v2/
#     bin_l2_annotation_v2.tsv.gz  (4980-bin universe; bin_id + 10 membership_X cols)
#   /home/kushima-pg/tad04212026/04_wgs_sv_boundary_overlap/output_v10/
#     sample_boundary_event_overlap_v10.tsv.gz  (sample x event x bin overlap;
#                                                 Step 04 v10 = top-1% CNV QC applied)
#   /home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v3/
#     sample_burden_L2_and_specificity_v3.tsv  (covariates + Diagnosis +
#                                                total_del_count for 99% filter)
#
# Outputs:
#   /home/kushima-pg/tad04212026/16_aggregate_primary_sensitivity_v1/
#     output/loco_diff_any_sensitivity_v1.tsv
#     output/loco_diff_any_sensitivity_v1.log
#
# Run:
#   cd /home/kushima-pg/tad04212026/16_aggregate_primary_sensitivity_v1/
#   python3 compute_loco_diff_any_sensitivity_v1.py
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
OUT_TSV = os.path.join(OUTDIR, "loco_diff_any_sensitivity_v1.tsv")
OUT_LOG = os.path.join(OUTDIR, "loco_diff_any_sensitivity_v1.log")

PRIMARY_BOUNDARY_CLASSES = [
    "HPC_Exc-CA", "HPC_Exc-DG", "HPC_Exc-ENT",
    "HPC_Inh-CGE", "HPC_Inh-MGE",
    "PFC_Astro", "PFC_Exc-DL", "PFC_Exc-UL",
    "PFC_Inh-CGE", "PFC_Inh-MGE",
]
COVS = ["Sex_numeric"] + ["PC%d" % i for i in range(1, 11)] + ["log1p_total_del_bases", "log1p_total_gene_DEL"]


def log(msg, fh=None):
    line = "[%.1fs] %s" % (time.time() - t0, msg)
    print(line)
    if fh is not None:
        fh.write(line + "\n")
        fh.flush()


with open(OUT_LOG, "w") as logfh:
    log("=" * 70, logfh)
    log("compute_loco_diff_any_sensitivity_v1.py", logfh)
    log("=" * 70, logfh)

    # --- Step 1: Load bin_l2_annotation_v2 ---
    log("Loading bin_l2_annotation_v2: %s" % BIN_ANNOT, logfh)
    anno = pd.read_csv(BIN_ANNOT, sep="\t", compression="gzip", low_memory=False,
                       dtype={"bin_id": "string"})
    log("  annotation shape: %s" % str(anno.shape), logfh)

    membership_cols = ["membership_%s" % c for c in PRIMARY_BOUNDARY_CLASSES]
    missing = [c for c in membership_cols if c not in anno.columns]
    if missing:
        sys.exit("ERROR: missing membership columns: %s" % missing)

    # Diff_any (any of 10 PRIMARY classes) = 4,980 bins
    diff_any_full = anno[anno[membership_cols].sum(axis=1) > 0].copy()
    full_bin_set = set(diff_any_full["bin_id"].tolist())
    log("  Full Diff_any union (10 classes): %d bins" % len(full_bin_set), logfh)
    assert len(full_bin_set) == 4980, "Expected 4,980 Diff_any bins; got %d" % len(full_bin_set)

    # LOCO_k bin sets
    loco_sets = {}
    for k_idx, k in enumerate(PRIMARY_BOUNDARY_CLASSES):
        # Drop class k: bins where any of the OTHER 9 classes have membership=1
        other_cols = [c for c in membership_cols if c != "membership_%s" % k]
        loco_bins = anno[anno[other_cols].sum(axis=1) > 0]["bin_id"].tolist()
        loco_sets[k] = set(loco_bins)
        n_loco = len(loco_sets[k])
        n_dropped = len(full_bin_set) - n_loco
        log("  LOCO_%s (drop class): %d bins (-%d from full)" % (k, n_loco, n_dropped), logfh)

    # --- Step 2: Load Step 4 v10 overlap table (DEL only) ---
    log("Loading Step 4 v10 overlap table: %s" % OVERLAP_TSV, logfh)
    ov = pd.read_csv(OVERLAP_TSV, sep="\t", compression="gzip", low_memory=False,
                     dtype={"sample_id": "string", "bin_id": "string", "event_id": "string"})
    log("  overlap shape: %s" % str(ov.shape), logfh)
    # Filter to DEL events
    ov_del = ov[ov["sv_type_norm"] == "DEL"].copy()
    log("  DEL-filtered overlap: %d rows" % len(ov_del), logfh)
    log("  unique DEL events: %d, unique bins: %d, unique samples: %d" % (
        ov_del["event_id"].nunique(), ov_del["bin_id"].nunique(),
        ov_del["sample_id"].nunique()
    ), logfh)

    # --- Step 3: Load Step 5 v3 burden TSV for covariates ---
    log("Loading Step 5 v3 burden TSV: %s" % BURDEN_TSV, logfh)
    burden = pd.read_csv(BURDEN_TSV, sep="\t", low_memory=False)
    log("  burden shape: %s" % str(burden.shape), logfh)
    log("  Diagnosis distribution: %s" % burden["Diagnosis"].value_counts().to_dict(), logfh)

    # ASD vs Healthy filter
    sub = burden[burden["Diagnosis"].isin(["ASD", "Healthy"])].copy()
    sub["is_case"] = (sub["Diagnosis"] == "ASD").astype(int)

    # 99-percentile filter on total_del_count
    sub_c = sub.dropna(subset=COVS + ["total_del_count"]).copy()
    thr = float(sub_c["total_del_count"].quantile(0.99))
    sub_filt = sub_c[sub_c["total_del_count"] <= thr].copy()
    log("  After 99%% filter (total_del_count <= %.0f): N=%d (ASD=%d, Healthy=%d)" % (
        thr, len(sub_filt), int((sub_filt["is_case"] == 1).sum()),
        int((sub_filt["is_case"] == 0).sum())
    ), logfh)

    analytic_samples = set(sub_filt["sample_id"].tolist())

    # --- Step 4: For each LOCO bin set + full Diff_any, compute per-sample disrupted bin count ---
    def compute_burden(bin_set, ov_del, samples_index):
        """Return per-sample count of unique bins in bin_set disrupted by DEL events."""
        ov_filt = ov_del[ov_del["bin_id"].isin(bin_set)]
        # unique (sample, bin) pairs (deduplicate multiple events hitting same bin)
        per_sample_bins = ov_filt.groupby("sample_id")["bin_id"].nunique()
        burden_series = per_sample_bins.reindex(samples_index, fill_value=0).astype(int)
        return burden_series

    samples_idx = pd.Index(sub_filt["sample_id"])

    # Full Diff_any union baseline
    log("=" * 70, logfh)
    log("Computing Full Diff_any union (baseline)", logfh)
    full_burden = compute_burden(full_bin_set, ov_del, samples_idx)
    sub_filt["Diff_any_full"] = full_burden.values

    n_exp_case_full = int(((sub_filt["is_case"] == 1) & (sub_filt["Diff_any_full"] >= 1)).sum())
    n_exp_ctrl_full = int(((sub_filt["is_case"] == 0) & (sub_filt["Diff_any_full"] >= 1)).sum())
    log("  Full Diff_any union: exp_case=%d, exp_ctrl=%d" % (n_exp_case_full, n_exp_ctrl_full), logfh)

    def fit_logistic(exposure_series):
        X = pd.DataFrame({"exposure": exposure_series.values})
        for c in COVS:
            X[c] = sub_filt[c].values
        X = sm.add_constant(X)
        y = sub_filt["is_case"].values
        try:
            m = sm.Logit(y, X).fit(disp=0, maxiter=200)
            beta = float(m.params["exposure"])
            se = float(m.bse["exposure"])
            return {
                "beta": beta, "SE": se,
                "OR": float(np.exp(beta)),
                "ci_lo": float(np.exp(beta - 1.96 * se)),
                "ci_hi": float(np.exp(beta + 1.96 * se)),
                "wald_p": float(m.pvalues["exposure"]),
                "status": "ok",
            }
        except Exception as e:
            return {"status": "error", "error": str(e)}

    full_res = fit_logistic(sub_filt["Diff_any_full"])
    log("  Full Diff_any baseline: OR=%.3f [%.3f-%.3f], Wald P=%.3e" % (
        full_res["OR"], full_res["ci_lo"], full_res["ci_hi"], full_res["wald_p"]
    ), logfh)

    # --- Step 5: LOCO computations ---
    log("=" * 70, logfh)
    log("Computing 10 LOCO Diff_any union ORs", logfh)
    results = [{
        "loco_dropped_class": "(none; full Diff_any baseline)",
        "n_bins": len(full_bin_set),
        "n_dropped": 0,
        "n_exp_case": n_exp_case_full,
        "n_exp_ctrl": n_exp_ctrl_full,
        **full_res
    }]
    for k in PRIMARY_BOUNDARY_CLASSES:
        log("LOCO drop %s ..." % k, logfh)
        loco_burden = compute_burden(loco_sets[k], ov_del, samples_idx)
        sub_filt["Diff_any_loco"] = loco_burden.values
        n_exp_case_k = int(((sub_filt["is_case"] == 1) & (sub_filt["Diff_any_loco"] >= 1)).sum())
        n_exp_ctrl_k = int(((sub_filt["is_case"] == 0) & (sub_filt["Diff_any_loco"] >= 1)).sum())
        loco_res = fit_logistic(sub_filt["Diff_any_loco"])
        log("  LOCO drop %s: %d bins (-%d), exp_case=%d, exp_ctrl=%d, OR=%.3f [%.3f-%.3f], Wald P=%.3e" % (
            k, len(loco_sets[k]), len(full_bin_set) - len(loco_sets[k]),
            n_exp_case_k, n_exp_ctrl_k,
            loco_res["OR"], loco_res["ci_lo"], loco_res["ci_hi"], loco_res["wald_p"]
        ), logfh)
        results.append({
            "loco_dropped_class": k,
            "n_bins": len(loco_sets[k]),
            "n_dropped": len(full_bin_set) - len(loco_sets[k]),
            "n_exp_case": n_exp_case_k,
            "n_exp_ctrl": n_exp_ctrl_k,
            **loco_res,
        })

    # --- Step 6: Save TSV ---
    out_df = pd.DataFrame(results)
    out_df.to_csv(OUT_TSV, sep="\t", index=False, na_rep="NA")
    log("=" * 70, logfh)
    log("Saved: %s" % OUT_TSV, logfh)
    log("=" * 70, logfh)
    log("Summary (LOCO range):", logfh)
    loco_only = out_df[out_df["loco_dropped_class"] != "(none; full Diff_any baseline)"]
    if not loco_only.empty:
        ors = loco_only["OR"].dropna()
        log("  10 LOCO ORs: min=%.3f, max=%.3f, mean=%.3f, range=%.3f" % (
            ors.min(), ors.max(), ors.mean(), ors.max() - ors.min()
        ), logfh)
        log("  Full Diff_any baseline OR=%.3f" % full_res["OR"], logfh)
        log("  All LOCO ORs ranging within +/-%.1f%% of baseline: %s" % (
            100 * max(abs(ors.max() - full_res["OR"]), abs(ors.min() - full_res["OR"])) / full_res["OR"],
            "single-class dependence ruled out" if ors.min() > 1.0 else "check classes with OR<1"
        ), logfh)
    log("[Total elapsed: %.1fs]" % (time.time() - t0), logfh)
