#!/usr/bin/env python3
# compute_aggregate_diff_any_asd_vs_sz_v2.py
#
# 処理内容:
#   R18 改訂提案の Phase A++ 補強解析。
#   aggregate Diff_any union (4,980-bin universe) における
#   ASD vs SZ heterogeneity を architecture-level で 2 tests:
#     1. ASD vs SZ DEL  (architecture-level ASD-preferential 主張の main support)
#     2. ASD vs SZ DUP  (control comparator; 差なしを期待)
#
#   Aggregate Diff_any union exposure を Step 5 v3 burden TSV から取得し、
#   B' logistic regression (case=ASD=1, ctrl=SZ=0) で
#   OR + 95% CI + Wald P を出力。
#
# v1 -> v2 変更点 (重要):
#   v1 では "subset 内 99-percentile high-burden CNV outlier filter" を
#   ASD+SZ subset 内で計算していた。これは ChatGPT review で
#   "case-case comparison 専用 filter は post hoc に見える" と却下された。
#   v2 では Step 03 v9 (34a_extract_rare_sv_events_v9.py) の Global flag に依存し、
#   subset-level 99-percentile filter を削除した。
#   Step 03 v9 で全コホート (ASD+SZ+Healthy) の cnv_count から q99 を計算し、
#   高 burden サンプルを全 event tableから永続的に除外しているため、
#   Step 05 v3 burden TSV は既に Global flag 通過後の値である。
#
# 入力:
#   /home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v3/
#     sample_burden_L2_and_specificity_v3.tsv
#       (4,980-universe per-class burden + group_primary aggregate columns;
#        Step 03 v9 + Step 04 v10 sample-level top-1% CNV count QC 適用後)
#
# 出力:
#   /home/kushima-pg/tad04212026/16_aggregate_primary_sensitivity_v1/
#     output/aggregate_diff_any_asd_vs_sz_v2.tsv
#     output/aggregate_diff_any_asd_vs_sz_v2.log
#
# B' logistic regression covariates (SV-type-specific):
#   DEL: Sex_numeric + PC1-10 + log1p_total_del_bases + log1p_total_gene_DEL
#   DUP: Sex_numeric + PC1-10 + log1p_total_dup_bases + log1p_total_gene_DUP
#
# 実行方法:
#   cd /home/kushima-pg/tad04212026/16_aggregate_primary_sensitivity_v1/
#   python3 compute_aggregate_diff_any_asd_vs_sz_v2.py
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
OUT_TSV = os.path.join(OUTDIR, "aggregate_diff_any_asd_vs_sz_v2.tsv")
OUT_LOG = os.path.join(OUTDIR, "aggregate_diff_any_asd_vs_sz_v2.log")


def log(msg, fh=None):
    line = "[%.1fs] %s" % (time.time() - t0, msg)
    print(line)
    if fh is not None:
        fh.write(line + "\n")
        fh.flush()


with open(OUT_LOG, "w") as logfh:
    log("=" * 70, logfh)
    log("compute_aggregate_diff_any_asd_vs_sz_v2.py", logfh)
    log("(v1 -> v2: removed subset-level 99% filter; rely on Step 03 v9 Global flag)", logfh)
    log("=" * 70, logfh)
    log("Input burden TSV: %s" % INPUT_BURDEN, logfh)

    df = pd.read_csv(INPUT_BURDEN, sep="\t", low_memory=False)
    log("Loaded burden TSV: %d rows x %d cols" % (df.shape[0], df.shape[1]), logfh)
    log("Diagnosis distribution: %s" % df["Diagnosis"].value_counts().to_dict(), logfh)

    COVS_BASE = ["Sex_numeric"] + ["PC%d" % i for i in range(1, 11)]
    SV_COVS_MAP = {
        "DEL": COVS_BASE + ["log1p_total_del_bases", "log1p_total_gene_DEL"],
        "DUP": COVS_BASE + ["log1p_total_dup_bases", "log1p_total_gene_DUP"],
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

        sv_covs = SV_COVS_MAP[sv]
        # v2: 99-percentile filter は削除。Step 03 v9 Global flag に依存。
        sub_filt = sub.dropna(subset=sv_covs + ["Diff_any"]).copy()

        n = sub_filt.shape[0]
        n_case = int((sub_filt["is_case"] == 1).sum())  # ASD
        n_ctrl = int((sub_filt["is_case"] == 0).sum())  # SZ
        n_exp_case = int(((sub_filt["is_case"] == 1) & (sub_filt["Diff_any"] >= 1)).sum())
        n_exp_ctrl = int(((sub_filt["is_case"] == 0) & (sub_filt["Diff_any"] >= 1)).sum())
        log("  N=%d (ASD=%d, SZ=%d, exp_ASD=%d, exp_SZ=%d) [Step 03 v9 Global flag only]" %
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

    # BH-FDR within 2-test ASD-vs-SZ family
    out_df = pd.DataFrame(results)
    if (out_df["status"] == "ok").any():
        ok = out_df[out_df["status"] == "ok"].copy()
        ok = ok.sort_values("wald_p").reset_index(drop=True)
        m = len(ok)
        ok["bh_rank"] = ok.index + 1
        ok["bh_q"] = (ok["wald_p"] * m / ok["bh_rank"]).cummax().clip(upper=1.0)
        out_df = out_df.merge(ok[["comparison", "sv_type", "bh_rank", "bh_q"]],
                              on=["comparison", "sv_type"], how="left")

    out_df.to_csv(OUT_TSV, sep="\t", index=False, na_rep="NA")
    log("=" * 70, logfh)
    log("Saved: %s" % OUT_TSV, logfh)
    log("=" * 70, logfh)
    log("Summary table (ASD as case, SZ as control; BH-FDR within 2-test family):", logfh)
    for _, row in out_df.iterrows():
        if row["status"] == "ok":
            log("  ASD_vs_SZ %s: OR=%.3f [%.3f-%.3f], Wald P=%.3e, BH q=%.3e (n_ASD=%d, n_SZ=%d, exp_ASD=%d, exp_SZ=%d)" % (
                row["sv_type"],
                row["OR"], row["ci_lo"], row["ci_hi"],
                row["wald_p"], row["bh_q"],
                int(row["n_case"]), int(row["n_ctrl"]),
                int(row["n_exp_case"]), int(row["n_exp_ctrl"])
            ), logfh)
    log("=" * 70, logfh)
    log("Interpretation guide:", logfh)
    log("  ASD_vs_SZ DEL OR > 1, P < 0.05: ASD-preferential at architecture level (main support)", logfh)
    log("  ASD_vs_SZ DUP OR ~ 1.0, P > 0.05: control comparator, differential-DEL specificity preserved", logfh)
    log("  Compare with per-class Stouffer P_global = 6.1e-3 for consistency check.", logfh)
    log("[Total elapsed: %.1fs]" % (time.time() - t0), logfh)
