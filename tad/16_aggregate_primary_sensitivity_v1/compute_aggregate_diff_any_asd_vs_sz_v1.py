#!/usr/bin/env python3
# compute_aggregate_diff_any_asd_vs_sz_v1.py
#
# 処理内容:
#   R18 改訂提案の Phase A++ 補強解析。
#   aggregate Diff_any union (4,980-bin universe) における
#   ASD vs SZ heterogeneity を architecture-level で 2 tests:
#     1. ASD vs SZ DEL  (architecture-level ASD-preferential 主張の primary 検証)
#     2. ASD vs SZ DUP  (control comparator; 差なしを期待)
#
#   Aggregate Diff_any union exposure を Step 5 v3 burden TSV から取得し、
#   B' logistic regression (case=ASD=1, ctrl=SZ=0) で
#   OR + 95% CI + Wald P を出力。
#   case-control subset: 99-percentile high-burden CNV outlier 除外。
#
#   この解析の目的:
#   - 現状の per-class Stouffer P_global (6.1e-3) は v1 構造（per-class primary）の遺物。
#   - architecture-level primary に再整理した場合、ASD-vs-SZ も同じ aggregate endpoint で
#     直接定量化することで logical consistency を確保。
#   - Reviewer が「なぜ primary endpoint が ASD-vs-Healthy では aggregate で
#     ASD-vs-SZ では per-class なのか」と問う risk を排除。
#
# 入力:
#   /home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v3/
#     sample_burden_L2_and_specificity_v3.tsv
#       (4,980-universe per-class burden + group_primary aggregate columns;
#        Step 04 v10 sample-level top-1% CNV count QC 適用後)
#
# 出力:
#   /home/kushima-pg/tad04212026/16_aggregate_primary_sensitivity_v1/
#     output/aggregate_diff_any_asd_vs_sz_v1.tsv
#     output/aggregate_diff_any_asd_vs_sz_v1.log
#
# Aggregate Diff_any 構築:
#   Diff_any = Diff_specific_n1_DEL + Diff_shared_n2plus_DEL
#       (これは bin_l2_annotation_v2 の group_primary annotation 由来; 4,980-universe)
#
# B' logistic regression covariates (SV-type-specific):
#   DEL: Sex_numeric + PC1-10 + log1p_total_del_bases + log1p_total_gene_DEL
#   DUP: Sex_numeric + PC1-10 + log1p_total_dup_bases + log1p_total_gene_DUP
#
# 99-percentile filter:
#   DEL: total_del_count <= q99(total_del_count) on ASD+SZ subset
#   DUP: total_dup_count <= q99(total_dup_count) on ASD+SZ subset
#
# 実行方法:
#   cd /home/kushima-pg/tad04212026/16_aggregate_primary_sensitivity_v1/
#   python3 compute_aggregate_diff_any_asd_vs_sz_v1.py
#
# 実行時間ログ出力。
import os
import sys
import time
import pandas as pd
import numpy as np
import statsmodels.api as sm

t0 = time.time()

INPUT_BURDEN = "/home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v3/sample_burden_L2_and_specificity_v3.tsv"
OUTDIR = "/home/kushima-pg/tad04212026/16_aggregate_primary_sensitivity_v1/output"
os.makedirs(OUTDIR, exist_ok=True)
OUT_TSV = os.path.join(OUTDIR, "aggregate_diff_any_asd_vs_sz_v1.tsv")
OUT_LOG = os.path.join(OUTDIR, "aggregate_diff_any_asd_vs_sz_v1.log")


def log(msg, fh=None):
    line = "[%.1fs] %s" % (time.time() - t0, msg)
    print(line)
    if fh is not None:
        fh.write(line + "\n")
        fh.flush()


with open(OUT_LOG, "w") as logfh:
    log("=" * 70, logfh)
    log("compute_aggregate_diff_any_asd_vs_sz_v1.py", logfh)
    log("=" * 70, logfh)
    log("Input burden TSV: %s" % INPUT_BURDEN, logfh)

    df = pd.read_csv(INPUT_BURDEN, sep="\t", low_memory=False)
    log("Loaded burden TSV: %d rows x %d cols" % (df.shape[0], df.shape[1]), logfh)
    log("Diagnosis distribution: %s" % df["Diagnosis"].value_counts().to_dict(), logfh)

    COVS_BASE = ["Sex_numeric"] + ["PC%d" % i for i in range(1, 11)]
    SV_COVS_MAP = {
        "DEL": (COVS_BASE + ["log1p_total_del_bases", "log1p_total_gene_DEL"], "total_del_count"),
        "DUP": (COVS_BASE + ["log1p_total_dup_bases", "log1p_total_gene_DUP"], "total_dup_count"),
    }

    results = []
    for sv in ["DEL", "DUP"]:
        log("-" * 70, logfh)
        log("Test: ASD vs SZ %s (architecture-level cross-disorder heterogeneity)" % sv, logfh)

        # ASD vs SZ subset (no Healthy)
        sub = df[df["Diagnosis"].isin(["ASD", "SZ"])].copy()
        sub["is_case"] = (sub["Diagnosis"] == "ASD").astype(int)  # ASD=1, SZ=0

        diff_spec = "n_boundary_group_primary__Diff_specific_n1_%s" % sv
        diff_shr = "n_boundary_group_primary__Diff_shared_n2plus_%s" % sv
        if diff_spec not in sub.columns or diff_shr not in sub.columns:
            log("  MISSING columns: %s or %s" % (diff_spec, diff_shr), logfh)
            results.append({
                "comparison": "ASD_vs_SZ", "sv_type": sv,
                "status": "missing_columns",
            })
            continue

        sub["Diff_any"] = sub[diff_spec] + sub[diff_shr]

        sv_covs, thr_col = SV_COVS_MAP[sv]
        sub_c = sub.dropna(subset=sv_covs + ["Diff_any"]).copy()

        if thr_col in sub_c.columns:
            thr = float(sub_c[thr_col].quantile(0.99))
            sub_filt = sub_c[sub_c[thr_col] <= thr].copy()
            log("  99%% filter on %s (within ASD+SZ subset): cutoff=%.0f, n_before=%d, n_after=%d" %
                (thr_col, thr, len(sub_c), len(sub_filt)), logfh)
        else:
            log("  WARNING: %s column not found; skipping 99%% filter" % thr_col, logfh)
            sub_filt = sub_c.copy()

        n = sub_filt.shape[0]
        n_case = int((sub_filt["is_case"] == 1).sum())  # ASD
        n_ctrl = int((sub_filt["is_case"] == 0).sum())  # SZ
        n_exp_case = int(((sub_filt["is_case"] == 1) & (sub_filt["Diff_any"] >= 1)).sum())
        n_exp_ctrl = int(((sub_filt["is_case"] == 0) & (sub_filt["Diff_any"] >= 1)).sum())
        log("  N=%d (ASD=%d, SZ=%d, exp_ASD=%d, exp_SZ=%d)" %
            (n, n_case, n_ctrl, n_exp_case, n_exp_ctrl), logfh)

        if n_case < 10 or n_ctrl < 10 or n_exp_case + n_exp_ctrl < 5:
            log("  insufficient cells; skipping", logfh)
            results.append({
                "comparison": "ASD_vs_SZ", "sv_type": sv,
                "n": n, "n_case": n_case, "n_ctrl": n_ctrl,
                "n_exp_case": n_exp_case, "n_exp_ctrl": n_exp_ctrl,
                "status": "insufficient_cells",
            })
            continue

        X = sub_filt[["Diff_any"] + sv_covs].copy()
        X.columns = ["exposure"] + sv_covs
        X = sm.add_constant(X)
        y = sub_filt["is_case"].values

        try:
            model = sm.Logit(y, X).fit(disp=0, maxiter=200)
            beta = float(model.params["exposure"])
            se = float(model.bse["exposure"])
            OR = float(np.exp(beta))
            ci_lo = float(np.exp(beta - 1.96 * se))
            ci_hi = float(np.exp(beta + 1.96 * se))
            pval = float(model.pvalues["exposure"])
            log("  beta=%.4f, SE=%.4f, OR=%.3f [%.3f-%.3f], Wald P=%.3e" %
                (beta, se, OR, ci_lo, ci_hi, pval), logfh)
            results.append({
                "comparison": "ASD_vs_SZ", "sv_type": sv,
                "n": n, "n_case": n_case, "n_ctrl": n_ctrl,
                "n_exp_case": n_exp_case, "n_exp_ctrl": n_exp_ctrl,
                "beta": beta, "SE": se, "OR": OR,
                "ci_lo": ci_lo, "ci_hi": ci_hi, "wald_p": pval,
                "status": "ok",
            })
        except Exception as e:
            log("  ERROR: %s" % str(e), logfh)
            results.append({
                "comparison": "ASD_vs_SZ", "sv_type": sv,
                "n": n, "n_case": n_case, "n_ctrl": n_ctrl,
                "n_exp_case": n_exp_case, "n_exp_ctrl": n_exp_ctrl,
                "status": "fit_error", "error": str(e),
            })

    out_df = pd.DataFrame(results)
    out_df.to_csv(OUT_TSV, sep="\t", index=False, na_rep="NA")
    log("=" * 70, logfh)
    log("Saved: %s" % OUT_TSV, logfh)
    log("=" * 70, logfh)
    log("Summary table (ASD as case, SZ as control):", logfh)
    for _, row in out_df.iterrows():
        if row["status"] == "ok":
            log("  ASD_vs_SZ %s: OR=%.3f [%.3f-%.3f], Wald P=%.3e (n_ASD=%d, n_SZ=%d, exp_ASD=%d, exp_SZ=%d)" % (
                row["sv_type"],
                row["OR"], row["ci_lo"], row["ci_hi"],
                row["wald_p"],
                int(row["n_case"]), int(row["n_ctrl"]),
                int(row["n_exp_case"]), int(row["n_exp_ctrl"])
            ), logfh)
    log("=" * 70, logfh)
    log("Interpretation guide:", logfh)
    log("  ASD_vs_SZ DEL OR > 1, P < 0.05: ASD-preferential at architecture level (primary)", logfh)
    log("  ASD_vs_SZ DUP OR ~ 1.0, P > 0.05: control comparator, differential-DEL specificity preserved", logfh)
    log("  Compare with per-class Stouffer P_global = 6.1e-3 for consistency check.", logfh)
    log("[Total elapsed: %.1fs]" % (time.time() - t0), logfh)
