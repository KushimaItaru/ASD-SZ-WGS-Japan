#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
40_run_replication_full_pipeline_v5.py
=====================================================
ファイル名: 40_run_replication_full_pipeline_v5.py
処理内容:
  - v4 -> v5 変更点 (WGS discovery sample-level top-1% CNV count QC 伝播のための
    path update のみ。arrayCGH v22 / MSSNG 入力は従来通り):
      * WGS_PRIMARY_OUT_DIR:
          tad04292026/06_wgs_primary_L2/output_v3
          -> tad04292026/06_wgs_primary_L2/output_v4
      * WGS_L2_RESULTS:
          B_prime_L2_classes_results_v3.tsv -> _v4.tsv
      * WGS_SPEC_RESULTS:
          B_prime_specificity_groups_results_v3.tsv -> _v4.tsv
      * OUT_DIR:
          tad04292026/10_replication_2way_meta/output_v4
          -> tad04292026/10_replication_2way_meta/output_v5
      * 全出力 TSV の suffix: _v4 -> _v5
      * ACGH/MSSNG/BIN_ANNOTATION は v4 と同一 (read-only 据え置き)
      * 解析ロジック (Step A/B/C/D) は v4 と完全一致

  - Replication (arrayCGH + MSSNG) の burden 計算 + B' logistic/GEE fit
    + IVW fixed-effect meta を 1 本にまとめた統合パイプライン。

  - Replication は DEL-only (arrayCGH は sv_type=="DEL" & pattern=="A" で明示 filter、
    MSSNG は元 overlap が既に DEL-only で抽出済み)。

  - Step A (burden):
      arrayCGH: sample_event_bin_overlap_v22.tsv.gz (tad04292026 v22) を
                bin_l2_annotation_v2 と join
      MSSNG   : mssng_sample_event_bin_overlap_v1.tsv.gz を bin_l2_annotation_v2 と join

  - Step B (fit):
      arrayCGH: B' pooled logistic (statsmodels.Logit)
      MSSNG   : GEE (Binomial + Independence, groups=FAMILYID)

  - Step C (2-way meta):
      arrayCGH ASD_vs_CONT beta/se と MSSNG ASD_vs_unaffSib beta/se を IVW.

  - Step D (3-way meta):
      WGS primary v4 (tad04292026 06_wgs_primary_L2 output_v4, ASD_vs_HC, DEL)
        + arrayCGH v22 (ASD_vs_CONT, DEL)
        + MSSNG (ASD_vs_unaffSib, DEL)
      IVW fixed-effect. N-cohort Cochran's Q / I^2 付き.

  - 実行時間記録あり

使い方:
  python 40_run_replication_full_pipeline_v5.py

出力 (OUT_DIR = tad04292026/10_replication_2way_meta/output_v5 に配置):
  [burden]
    arraycgh_burden_L2_and_specificity_v5.tsv
    mssng_burden_L2_and_specificity_v5.tsv
    replication_burden_summary_v5.tsv
  [fit cohort-wise]
    arraycgh_L2_classes_results_v5.tsv
    arraycgh_specificity_groups_results_v5.tsv
    mssng_L2_classes_results_v5.tsv
    mssng_specificity_groups_results_v5.tsv
  [fit 2-way meta]
    meta_L2_classes_results_v5.tsv
    meta_specificity_groups_results_v5.tsv
  [fit 3-way meta (WGS v4 + arrayCGH v22 + MSSNG)]
    meta3way_L2_classes_results_v5.tsv
    meta3way_specificity_groups_results_v5.tsv
  [summary]
    replication_fit_summary_v5.tsv
"""
from __future__ import annotations
import sys
import time
import warnings
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
from scipy import stats as sp_stats
import statsmodels.api as sm
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.genmod.families import Binomial
from statsmodels.genmod.cov_struct import Independence

warnings.filterwarnings("ignore")


# =========================================================
# PATHS  (固定。必要に応じて編集してください)
# =========================================================
# arrayCGH 入力 (v4 と同一; WGS top-1% QC は WGS discovery のみ適用、arrayCGH は不変)
ACGH_OVERLAP = Path(
    "/lustre12/home/kushima-pg/tad04292026/08_arraycgh_sample_burden/"
    "output_v22/sample_event_bin_overlap_v22.tsv.gz"
)
ACGH_COVARIATES = Path(
    "/lustre12/home/kushima-pg/tad04292026/08_arraycgh_sample_burden/"
    "output_v22/sample_covariates_v22.tsv"
)
# MSSNG 入力 (v4 と同一)
MSSNG_OVERLAP = Path(
    "/lustre12/home/kushima-pg/noncoding_tad_specificity_04022026/"
    "output_mssng_factorial_overlap_v2/mssng_sample_event_bin_overlap_v1.tsv.gz"
)
MSSNG_COVARIATES = Path(
    "/lustre12/home/kushima-pg/noncoding_tad_specificity_04022026/"
    "output_mssng_factorial_overlap_v2/mssng_sample_covariates_v1.tsv"
)
# bin_annotation (v4 と同一)
BIN_ANNOTATION = Path(
    "/lustre12/home/kushima-pg/tad04292026/02_bin_l2_annotation/"
    "output_v2/bin_l2_annotation_v2.tsv.gz"
)
# v5 CHANGE: WGS primary を tad04292026 Step 6 v4 出力へ切り替え
# (v4 = sample-level top-1% CNV count QC 適用済みの WGS discovery fit 結果)
# 39_fit_B_prime_L2_and_specificity_v4.R の出力を参照する。
# comparison=ASD_vs_HC, sv_type=DEL のフィルタ仕様は v3 と同一。
WGS_PRIMARY_OUT_DIR = Path(
    "/lustre12/home/kushima-pg/tad04292026/06_wgs_primary_L2/output_v4"
)
WGS_L2_RESULTS = WGS_PRIMARY_OUT_DIR / "B_prime_L2_classes_results_v4.tsv"
WGS_SPEC_RESULTS = WGS_PRIMARY_OUT_DIR / "B_prime_specificity_groups_results_v4.tsv"

# v5 CHANGE: OUT_DIR を tad04292026/10_replication_2way_meta/output_v5 へ変更
OUT_DIR = Path(
    "/lustre12/home/kushima-pg/tad04292026/10_replication_2way_meta/output_v5"
)

MIN_CELL_COUNT = 5


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


# =========================================================
# 定数
# =========================================================
L2_CLASSES = [
    "HPC_Exc-CA", "HPC_Exc-DG", "HPC_Exc-ENT",
    "HPC_Inh-CGE", "HPC_Inh-MGE",
    "PFC_Astro", "PFC_Exc-DL", "PFC_Exc-UL",
    "PFC_Inh-CGE", "PFC_Inh-MGE",
]
GROUP_COLUMNS = ["group_primary", "group_s2", "group_s3", "group_s4", "group_s5"]
EXPOSURES = ["n_boundary", "n_events", "carrier_boundary"]
SV_TAG = "DEL"  # replication は DEL-only

ACGH_COMPARISONS = [
    ("ASD_vs_CONT", "ASD", "CONT"),
    ("ASD_vs_SZ",   "ASD", "SCZ"),
]
MSSNG_COMPARISONS = [
    ("ASD_vs_unaffSib", "ASD", "unaffected_sibling"),
]
# 3-way メタで用いる WGS primary 結果のフィルタ条件
WGS_PRIMARY_COMPARISON = "ASD_vs_HC"
WGS_PRIMARY_SV = "DEL"

ACGH_COV_COLS = [
    "sample_id", "diagnosis", "sex", "platform_nimblegen",
    "log1p_total_del_bases_A", "log1p_total_gene_DEL_A",
]
MSSNG_COV_COLS = [
    "sample_id", "Status", "Sex_numeric", "FAMILYID",
    "log1p_total_del_bases", "log1p_total_gene_DEL",
    "anc_OTH", "anc_SAS", "anc_EAS", "anc_AMR", "anc_AFR",
    "plat_NovaSeq", "plat_HiSeq", "plat_HiSeq2000", "plat_HiSeq2500", "plat_CG",
]
MSSNG_COVS_BASE = ["Sex_numeric",
                   "log1p_total_del_bases", "log1p_total_gene_DEL"]
MSSNG_ANCESTRY = ["anc_OTH", "anc_SAS", "anc_EAS", "anc_AMR", "anc_AFR"]
MSSNG_PLATFORM = ["plat_NovaSeq", "plat_HiSeq", "plat_HiSeq2000",
                  "plat_HiSeq2500", "plat_CG"]


# =========================================================
# Step A shared: bin annotation loader
# =========================================================
def load_bin_annotation(path: Path) -> pd.DataFrame:
    log(f"Loading bin annotation: {path}")
    anno = pd.read_csv(path, sep="\t", compression="gzip",
                       dtype={"bin_id": "string",
                              "chrom": "string",
                              "best_matching_level": "string"})
    expected = (["bin_id", "overlaps_diffbound_any", "n_L2_diff_support"]
                + [f"membership_{c}" for c in L2_CLASSES]
                + GROUP_COLUMNS)
    missing = [c for c in expected if c not in anno.columns]
    if missing:
        raise RuntimeError(f"bin_annotation missing columns: {missing}")
    log(f"  annotation shape: {anno.shape}")
    needed = (["bin_id"]
              + [f"membership_{c}" for c in L2_CLASSES]
              + GROUP_COLUMNS)
    return anno[needed].copy().set_index("bin_id")


# =========================================================
# Step A: burden computation (DEL only)
# =========================================================
def compute_burden(
    long_df: pd.DataFrame,
    all_samples: pd.Index,
    sv_tag: str,
) -> pd.DataFrame:
    out = pd.DataFrame(index=all_samples)

    # 10 L2 classes
    for c in L2_CLASSES:
        col_mem = f"membership_{c}"
        sub_c = long_df[long_df[col_mem] == 1]
        if len(sub_c) == 0:
            out[f"n_boundary_{c}_{sv_tag}"] = 0
            out[f"n_events_{c}_{sv_tag}"] = 0
            out[f"carrier_boundary_{c}_{sv_tag}"] = 0
            continue
        nbin = sub_c.groupby("sample_id")["bin_id"].nunique()
        nevt = sub_c.groupby("sample_id")["event_id"].nunique()
        out[f"n_boundary_{c}_{sv_tag}"] = nbin.reindex(all_samples, fill_value=0).astype(int)
        out[f"n_events_{c}_{sv_tag}"] = nevt.reindex(all_samples, fill_value=0).astype(int)
        out[f"carrier_boundary_{c}_{sv_tag}"] = (out[f"n_boundary_{c}_{sv_tag}"] >= 1).astype(int)

    # 5 group definitions
    for gcol in GROUP_COLUMNS:
        labels = long_df[gcol].dropna().unique().tolist()
        for lbl in labels:
            sub_l = long_df[long_df[gcol] == lbl]
            col_b = f"n_boundary_{gcol}__{lbl}_{sv_tag}"
            col_e = f"n_events_{gcol}__{lbl}_{sv_tag}"
            col_c = f"carrier_boundary_{gcol}__{lbl}_{sv_tag}"
            if len(sub_l) == 0:
                out[col_b] = 0
                out[col_e] = 0
                out[col_c] = 0
                continue
            nbin = sub_l.groupby("sample_id")["bin_id"].nunique()
            nevt = sub_l.groupby("sample_id")["event_id"].nunique()
            out[col_b] = nbin.reindex(all_samples, fill_value=0).astype(int)
            out[col_e] = nevt.reindex(all_samples, fill_value=0).astype(int)
            out[col_c] = (out[col_b] >= 1).astype(int)

    return out.reset_index().rename(columns={"index": "sample_id"})


def process_arraycgh(overlap_path: Path, cov_path: Path, anno: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    log("=" * 50)
    log("Processing arrayCGH ...")
    log(f"  overlap: {overlap_path}")
    ov = pd.read_csv(overlap_path, sep="\t", compression="gzip", low_memory=False,
                     dtype={"sample_id": "string", "bin_id": "string",
                            "event_id": "string", "sv_type": "string",
                            "pattern": "string"})
    log(f"  raw overlap rows: {len(ov):,}")
    before = len(ov)
    if "pattern" in ov.columns:
        ov = ov.loc[(ov["sv_type"] == "DEL") & (ov["pattern"] == "A")].copy()
    else:
        ov = ov.loc[ov["sv_type"] == "DEL"].copy()
    log(f"  after DEL+PatternA filter: {len(ov):,} (dropped {before - len(ov):,})")

    cov = pd.read_csv(cov_path, sep="\t", low_memory=False,
                      dtype={"sample_id": "string"})
    missing = [c for c in ACGH_COV_COLS if c not in cov.columns]
    if missing:
        log(f"  [WARN] arrayCGH covariates missing: {missing}")
    keep = [c for c in ACGH_COV_COLS if c in cov.columns]
    cov = cov[keep].drop_duplicates("sample_id").copy()
    log(f"  covariates shape: {cov.shape}")
    if "diagnosis" in cov.columns:
        log(f"  diagnosis distribution:")
        for k, v in cov["diagnosis"].value_counts(dropna=False).items():
            log(f"    {k}: {v}")

    ov_slim = ov[["sample_id", "bin_id", "event_id"]].drop_duplicates()
    merged = ov_slim.join(anno, on="bin_id", how="inner")
    log(f"  after inner join: {len(merged):,} "
        f"(lost {len(ov_slim) - len(merged):,} rows, unique bins = {merged['bin_id'].nunique():,})")

    samples = pd.Index(cov["sample_id"].unique(), name="sample_id")
    log(f"  unique samples: {len(samples):,}")
    burden = compute_burden(merged, samples, sv_tag=SV_TAG)

    final = cov.merge(burden, how="left", on="sample_id")
    num_cols = [c for c in final.columns
                if c.startswith(("n_boundary_", "n_events_", "carrier_boundary_"))]
    final[num_cols] = final[num_cols].fillna(0).astype(int)
    log(f"  final arrayCGH burden shape: {final.shape}")

    qc = {
        "dataset": "arraycgh",
        "n_samples": int(len(final)),
        "n_ASD": int((final.get("diagnosis", pd.Series(dtype=str)).astype(str) == "ASD").sum()),
        "n_SCZ": int((final.get("diagnosis", pd.Series(dtype=str)).astype(str) == "SCZ").sum()),
        "n_SZ": int((final.get("diagnosis", pd.Series(dtype=str)).astype(str) == "SZ").sum()),
        "n_CONT": int((final.get("diagnosis", pd.Series(dtype=str)).astype(str) == "CONT").sum()),
    }
    for lbl in ["Diff_specific_n1", "Diff_shared_n2plus", "Static"]:
        col = f"carrier_boundary_group_primary__{lbl}_{SV_TAG}"
        if col in final.columns and "diagnosis" in final.columns:
            by_dx = final.groupby("diagnosis")[col].sum()
            qc[f"carrier_{lbl}"] = {str(k): int(v) for k, v in by_dx.items()}
    return final, qc


def process_mssng(overlap_path: Path, cov_path: Path, anno: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    log("=" * 50)
    log("Processing MSSNG ...")
    log(f"  overlap: {overlap_path}")
    ov = pd.read_csv(overlap_path, sep="\t", compression="gzip", low_memory=False,
                     dtype={"sample_id": "string", "bin_id": "string",
                            "event_id": "string"})
    log(f"  overlap rows: {len(ov):,} (DEL-only, pre-extracted)")

    cov = pd.read_csv(cov_path, sep="\t", low_memory=False,
                      dtype={"sample_id": "string", "FAMILYID": "string"})
    missing = [c for c in MSSNG_COV_COLS if c not in cov.columns]
    if missing:
        log(f"  [WARN] MSSNG covariates missing: {missing}")
    keep = [c for c in MSSNG_COV_COLS if c in cov.columns]
    cov = cov[keep].drop_duplicates("sample_id").copy()
    log(f"  covariates shape: {cov.shape}")
    if "Status" in cov.columns:
        log(f"  Status distribution:")
        for k, v in cov["Status"].value_counts(dropna=False).items():
            log(f"    {k}: {v}")

    before = len(ov)
    ov = ov[ov["sample_id"].isin(cov["sample_id"])].copy()
    log(f"  overlap filtered to covariate samples: {len(ov):,} "
        f"(dropped {before - len(ov):,})")

    ov_slim = ov[["sample_id", "bin_id", "event_id"]].drop_duplicates()
    merged = ov_slim.join(anno, on="bin_id", how="inner")
    log(f"  after inner join: {len(merged):,} "
        f"(lost {len(ov_slim) - len(merged):,} rows, unique bins = {merged['bin_id'].nunique():,})")

    samples = pd.Index(cov["sample_id"].unique(), name="sample_id")
    log(f"  unique samples: {len(samples):,}")
    burden = compute_burden(merged, samples, sv_tag=SV_TAG)

    final = cov.merge(burden, how="left", on="sample_id")
    num_cols = [c for c in final.columns
                if c.startswith(("n_boundary_", "n_events_", "carrier_boundary_"))]
    final[num_cols] = final[num_cols].fillna(0).astype(int)
    log(f"  final MSSNG burden shape: {final.shape}")

    qc = {
        "dataset": "mssng",
        "n_samples": int(len(final)),
        "n_families": int(final["FAMILYID"].nunique()) if "FAMILYID" in final.columns else None,
        "n_ASD": int((final.get("Status", pd.Series(dtype=str)).astype(str) == "ASD").sum()),
        "n_unaffected_sibling": int((final.get("Status", pd.Series(dtype=str)).astype(str)
                                    == "unaffected_sibling").sum()),
    }
    for lbl in ["Diff_specific_n1", "Diff_shared_n2plus", "Static"]:
        col = f"carrier_boundary_group_primary__{lbl}_{SV_TAG}"
        if col in final.columns and "Status" in final.columns:
            by_st = final.groupby("Status")[col].sum()
            qc[f"carrier_{lbl}"] = {str(k): int(v) for k, v in by_st.items()}
    return final, qc


# =========================================================
# Step B-1: BH FDR
# =========================================================
def bh_fdr(p: pd.Series) -> pd.Series:
    p = p.copy()
    mask = ~p.isna()
    if mask.sum() == 0:
        return p
    adj = pd.Series(np.nan, index=p.index)
    vals = p[mask].values.astype(float)
    order = np.argsort(vals)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(vals) + 1)
    n = len(vals)
    bh = vals * n / ranks
    bh_sorted = bh[order]
    for i in range(len(bh_sorted) - 2, -1, -1):
        bh_sorted[i] = min(bh_sorted[i], bh_sorted[i + 1])
    bh[order] = np.clip(bh_sorted, 0, 1)
    adj.loc[mask[mask].index] = bh
    return adj


# =========================================================
# Step B-2: fit helpers
# =========================================================
def _empty_res(status, n_case, n_ctrl, cc, cn):
    return {
        "n_case": int(n_case), "n_ctrl": int(n_ctrl),
        "carrier_case": int(cc), "carrier_ctrl": int(cn),
        "beta": np.nan, "se": np.nan,
        "or": np.nan, "or_lo95": np.nan, "or_hi95": np.nan,
        "p_value": np.nan, "fit_status": status,
    }


def _empty_res_gee(status, n_case, n_ctrl, cc, cn, n_fam):
    d = _empty_res(status, n_case, n_ctrl, cc, cn)
    d["n_families"] = int(n_fam)
    return d


def fit_bprime_logit_acgh(
    df: pd.DataFrame,
    exposure_col: str,
    case_label: str,
    control_label: str,
    min_cell: int = 5,
) -> dict:
    sub = df[df["diagnosis"].astype(str).isin([case_label, control_label])].copy()
    sub["is_case"] = (sub["diagnosis"].astype(str) == case_label).astype(int)
    if exposure_col not in sub.columns:
        return _empty_res("missing_exposure_column", 0, 0, 0, 0)
    needed = [exposure_col, "sex", "log1p_total_del_bases_A",
              "log1p_total_gene_DEL_A", "platform_nimblegen", "is_case"]
    missing_cov = [c for c in needed if c not in sub.columns]
    if missing_cov:
        return _empty_res(f"missing_cov:{','.join(missing_cov)}", 0, 0, 0, 0)
    sub = sub[needed].replace([np.inf, -np.inf], np.nan).dropna().copy()
    sub[exposure_col] = pd.to_numeric(sub[exposure_col], errors="coerce")
    sub = sub.dropna(subset=[exposure_col])
    n_case = int(sub["is_case"].sum())
    n_ctrl = int(len(sub) - n_case)
    carrier_case = int(((sub["is_case"] == 1) & (sub[exposure_col] >= 1)).sum())
    carrier_ctrl = int(((sub["is_case"] == 0) & (sub[exposure_col] >= 1)).sum())

    if sub[exposure_col].nunique() < 2 or sub["is_case"].nunique() < 2:
        return _empty_res("no_variance", n_case, n_ctrl, carrier_case, carrier_ctrl)
    if (carrier_case + carrier_ctrl) < min_cell:
        return _empty_res("insufficient_carriers", n_case, n_ctrl, carrier_case, carrier_ctrl)

    X_cols = [exposure_col, "sex", "log1p_total_del_bases_A",
              "log1p_total_gene_DEL_A", "platform_nimblegen"]
    X = sub[X_cols].astype(float)
    X = sm.add_constant(X, has_constant="add")
    y = sub["is_case"].astype(int).to_numpy()

    last_err = ""
    for method in ["newton", "lbfgs", "bfgs"]:
        try:
            fit = sm.Logit(y, X).fit(disp=0, maxiter=200, method=method)
            if np.isfinite(fit.params[exposure_col]):
                beta = float(fit.params[exposure_col])
                se = float(fit.bse[exposure_col])
                p = float(fit.pvalues[exposure_col])
                status = "ok" if fit.mle_retvals.get("converged", False) else "not_converged"
                return {
                    "n_case": n_case, "n_ctrl": n_ctrl,
                    "carrier_case": carrier_case, "carrier_ctrl": carrier_ctrl,
                    "beta": beta, "se": se,
                    "or": float(np.exp(beta)),
                    "or_lo95": float(np.exp(beta - 1.959964 * se)),
                    "or_hi95": float(np.exp(beta + 1.959964 * se)),
                    "p_value": p,
                    "fit_status": status,
                }
        except Exception as e:
            last_err = str(e)[:120]
    return _empty_res(f"glm_error:{last_err}", n_case, n_ctrl, carrier_case, carrier_ctrl)


def fit_gee_mssng(
    df: pd.DataFrame,
    exposure_col: str,
    case_label: str,
    control_label: str,
    min_cell: int = 5,
) -> dict:
    sub = df[df["Status"].astype(str).isin([case_label, control_label])].copy()
    sub["is_case"] = (sub["Status"].astype(str) == case_label).astype(int)
    anc_cols = [c for c in MSSNG_ANCESTRY if c in sub.columns]
    plat_cols = [c for c in MSSNG_PLATFORM if c in sub.columns]
    needed = [exposure_col, "is_case", "FAMILYID"] + MSSNG_COVS_BASE + anc_cols + plat_cols
    missing_cov = [c for c in needed if c not in sub.columns]
    if missing_cov:
        return _empty_res_gee(f"missing_cov:{','.join(missing_cov)}", 0, 0, 0, 0, 0)
    sub = sub[needed].replace([np.inf, -np.inf], np.nan).dropna().copy()
    sub[exposure_col] = pd.to_numeric(sub[exposure_col], errors="coerce")
    sub = sub.dropna(subset=[exposure_col])
    n_case = int(sub["is_case"].sum())
    n_ctrl = int(len(sub) - n_case)
    carrier_case = int(((sub["is_case"] == 1) & (sub[exposure_col] >= 1)).sum())
    carrier_ctrl = int(((sub["is_case"] == 0) & (sub[exposure_col] >= 1)).sum())
    n_fam = int(sub["FAMILYID"].nunique())

    if sub[exposure_col].nunique() < 2 or sub["is_case"].nunique() < 2:
        return _empty_res_gee("no_variance", n_case, n_ctrl, carrier_case, carrier_ctrl, n_fam)
    if (carrier_case + carrier_ctrl) < min_cell:
        return _empty_res_gee("insufficient_carriers",
                              n_case, n_ctrl, carrier_case, carrier_ctrl, n_fam)

    X_cols = [exposure_col] + MSSNG_COVS_BASE + anc_cols + plat_cols
    X = sub[X_cols].astype(float).to_numpy()
    X = np.column_stack([np.ones(len(X)), X])
    y = sub["is_case"].astype(int).to_numpy()

    try:
        model = GEE(
            endog=y, exog=X, groups=sub["FAMILYID"].values,
            family=Binomial(), cov_struct=Independence(),
        )
        fit = model.fit(maxiter=200)
        if fit.converged:
            beta = float(fit.params[1])
            se = float(fit.bse[1])
            p = float(fit.pvalues[1])
            return {
                "n_case": n_case, "n_ctrl": n_ctrl,
                "carrier_case": carrier_case, "carrier_ctrl": carrier_ctrl,
                "n_families": n_fam,
                "beta": beta, "se": se,
                "or": float(np.exp(beta)),
                "or_lo95": float(np.exp(beta - 1.959964 * se)),
                "or_hi95": float(np.exp(beta + 1.959964 * se)),
                "p_value": p,
                "fit_status": "ok",
            }
        else:
            return _empty_res_gee("not_converged",
                                  n_case, n_ctrl, carrier_case, carrier_ctrl, n_fam)
    except Exception as e:
        return _empty_res_gee(f"gee_error:{str(e)[:120]}",
                              n_case, n_ctrl, carrier_case, carrier_ctrl, n_fam)


# =========================================================
# Step C: IVW fixed-effect meta (2-way)
# =========================================================
def ivw_meta(beta_a, se_a, beta_b, se_b) -> dict:
    if any(pd.isna([beta_a, se_a, beta_b, se_b])):
        return {"beta_meta": np.nan, "se_meta": np.nan,
                "or_meta": np.nan, "or_meta_lo95": np.nan, "or_meta_hi95": np.nan,
                "p_meta": np.nan, "q_het": np.nan, "p_het": np.nan, "i2_het": np.nan,
                "meta_status": "missing_inputs"}
    w_a = 1.0 / (se_a ** 2)
    w_b = 1.0 / (se_b ** 2)
    beta_m = (w_a * beta_a + w_b * beta_b) / (w_a + w_b)
    se_m = np.sqrt(1.0 / (w_a + w_b))
    z = beta_m / se_m
    p_m = 2.0 * sp_stats.norm.sf(abs(z))
    q = w_a * (beta_a - beta_m) ** 2 + w_b * (beta_b - beta_m) ** 2
    p_q = 1.0 - sp_stats.chi2.cdf(q, df=1)
    i2 = max(0.0, (q - 1) / q) if q > 0 else 0.0
    return {
        "beta_meta": float(beta_m), "se_meta": float(se_m),
        "or_meta": float(np.exp(beta_m)),
        "or_meta_lo95": float(np.exp(beta_m - 1.959964 * se_m)),
        "or_meta_hi95": float(np.exp(beta_m + 1.959964 * se_m)),
        "p_meta": float(p_m),
        "q_het": float(q), "p_het": float(p_q), "i2_het": float(i2),
        "meta_status": "ok",
    }


# =========================================================
# Step D support: IVW fixed-effect meta (N-way, generalized)
# =========================================================
def ivw_meta_nway(betas: List[float], ses: List[float]) -> dict:
    """
    Inverse-variance-weighted fixed-effect meta-analysis for N >= 2 cohorts.
    Inputs with any NaN are dropped; if fewer than 2 valid cohorts remain,
    returns a degraded / missing-inputs result.
    Returns Cochran's Q, df, p_het, I^2 in addition to meta beta/se/OR/p.
    """
    pairs = [(float(b), float(s)) for b, s in zip(betas, ses)
             if not (pd.isna(b) or pd.isna(s))]
    k = len(pairs)
    if k < 2:
        return {"beta_meta": np.nan, "se_meta": np.nan,
                "or_meta": np.nan, "or_meta_lo95": np.nan, "or_meta_hi95": np.nan,
                "p_meta": np.nan,
                "q_het": np.nan, "df_het": np.nan,
                "p_het": np.nan, "i2_het": np.nan,
                "k_included": int(k), "meta_status": "missing_inputs"}
    bs = np.array([p[0] for p in pairs], dtype=float)
    ss = np.array([p[1] for p in pairs], dtype=float)
    ws = 1.0 / (ss ** 2)
    w_sum = float(ws.sum())
    beta_m = float((ws * bs).sum() / w_sum)
    se_m = float(np.sqrt(1.0 / w_sum))
    z = beta_m / se_m
    p_m = 2.0 * sp_stats.norm.sf(abs(z))
    q = float((ws * (bs - beta_m) ** 2).sum())
    df = k - 1
    p_q = 1.0 - sp_stats.chi2.cdf(q, df=df) if df > 0 else np.nan
    i2 = max(0.0, (q - df) / q) if q > 0 else 0.0
    return {
        "beta_meta": beta_m, "se_meta": se_m,
        "or_meta": float(np.exp(beta_m)),
        "or_meta_lo95": float(np.exp(beta_m - 1.959964 * se_m)),
        "or_meta_hi95": float(np.exp(beta_m + 1.959964 * se_m)),
        "p_meta": float(p_m),
        "q_het": q, "df_het": int(df),
        "p_het": float(p_q) if df > 0 else np.nan,
        "i2_het": float(i2),
        "k_included": int(k), "meta_status": "ok",
    }


# =========================================================
# Step D support: WGS primary loader
# =========================================================
def load_wgs_primary_results(l2_path: Path, spec_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load tad04292026 06_wgs_primary_L2 output_v4 (R script
    39_fit_B_prime_L2_and_specificity_v4) outputs and filter to
    comparison = ASD_vs_HC, sv_type = DEL.
    """
    log(f"Loading WGS primary L2 results:   {l2_path}")
    wgs_l2 = pd.read_csv(l2_path, sep="\t", low_memory=False)
    wgs_l2 = wgs_l2.loc[
        (wgs_l2["comparison"] == WGS_PRIMARY_COMPARISON)
        & (wgs_l2["sv_type"] == WGS_PRIMARY_SV)
    ].copy()
    log(f"  WGS L2 filtered rows (ASD_vs_HC, DEL): {len(wgs_l2):,}")

    log(f"Loading WGS primary spec results: {spec_path}")
    wgs_sp = pd.read_csv(spec_path, sep="\t", low_memory=False)
    wgs_sp = wgs_sp.loc[
        (wgs_sp["comparison"] == WGS_PRIMARY_COMPARISON)
        & (wgs_sp["sv_type"] == WGS_PRIMARY_SV)
    ].copy()
    log(f"  WGS spec filtered rows (ASD_vs_HC, DEL): {len(wgs_sp):,}")
    return wgs_l2, wgs_sp


# =========================================================
# Step B drivers
# =========================================================
def run_acgh_l2(acgh: pd.DataFrame, min_cell: int) -> pd.DataFrame:
    rows = []
    for l2 in L2_CLASSES:
        for exp_kind in EXPOSURES:
            col = f"{exp_kind}_{l2}_{SV_TAG}"
            for cmp_name, case_lbl, ctrl_lbl in ACGH_COMPARISONS:
                res = fit_bprime_logit_acgh(acgh, col, case_lbl, ctrl_lbl, min_cell)
                res.update({
                    "dataset": "arraycgh", "analysis": "L2_class",
                    "L2_class": l2, "sv_type": SV_TAG,
                    "exposure": exp_kind, "exposure_col": col,
                    "comparison": cmp_name, "case_label": case_lbl,
                    "control_label": ctrl_lbl,
                })
                rows.append(res)
    return pd.DataFrame(rows)


def run_mssng_l2(mssng: pd.DataFrame, min_cell: int) -> pd.DataFrame:
    rows = []
    for l2 in L2_CLASSES:
        for exp_kind in EXPOSURES:
            col = f"{exp_kind}_{l2}_{SV_TAG}"
            for cmp_name, case_lbl, ctrl_lbl in MSSNG_COMPARISONS:
                res = fit_gee_mssng(mssng, col, case_lbl, ctrl_lbl, min_cell)
                res.update({
                    "dataset": "mssng", "analysis": "L2_class",
                    "L2_class": l2, "sv_type": SV_TAG,
                    "exposure": exp_kind, "exposure_col": col,
                    "comparison": cmp_name, "case_label": case_lbl,
                    "control_label": ctrl_lbl,
                })
                rows.append(res)
    return pd.DataFrame(rows)


def _collect_labels(df: pd.DataFrame, gcol: str) -> List[str]:
    prefix = f"carrier_boundary_{gcol}__"
    labs = []
    suffix = f"_{SV_TAG}"
    for c in df.columns:
        if c.startswith(prefix) and c.endswith(suffix):
            lb = c[len(prefix):-len(suffix)]
            labs.append(lb)
    return sorted(set(labs))


def run_acgh_groups(acgh: pd.DataFrame, min_cell: int) -> pd.DataFrame:
    rows = []
    for gcol in GROUP_COLUMNS:
        labels = _collect_labels(acgh, gcol)
        for lbl in labels:
            for exp_kind in EXPOSURES:
                col = f"{exp_kind}_{gcol}__{lbl}_{SV_TAG}"
                if col not in acgh.columns:
                    continue
                for cmp_name, case_lbl, ctrl_lbl in ACGH_COMPARISONS:
                    res = fit_bprime_logit_acgh(acgh, col, case_lbl, ctrl_lbl, min_cell)
                    res.update({
                        "dataset": "arraycgh", "analysis": "specificity_group",
                        "group_scheme": gcol, "group_label": lbl,
                        "sv_type": SV_TAG,
                        "exposure": exp_kind, "exposure_col": col,
                        "comparison": cmp_name, "case_label": case_lbl,
                        "control_label": ctrl_lbl,
                    })
                    rows.append(res)
    return pd.DataFrame(rows)


def run_mssng_groups(mssng: pd.DataFrame, min_cell: int) -> pd.DataFrame:
    rows = []
    for gcol in GROUP_COLUMNS:
        labels = _collect_labels(mssng, gcol)
        for lbl in labels:
            for exp_kind in EXPOSURES:
                col = f"{exp_kind}_{gcol}__{lbl}_{SV_TAG}"
                if col not in mssng.columns:
                    continue
                for cmp_name, case_lbl, ctrl_lbl in MSSNG_COMPARISONS:
                    res = fit_gee_mssng(mssng, col, case_lbl, ctrl_lbl, min_cell)
                    res.update({
                        "dataset": "mssng", "analysis": "specificity_group",
                        "group_scheme": gcol, "group_label": lbl,
                        "sv_type": SV_TAG,
                        "exposure": exp_kind, "exposure_col": col,
                        "comparison": cmp_name, "case_label": case_lbl,
                        "control_label": ctrl_lbl,
                    })
                    rows.append(res)
    return pd.DataFrame(rows)


# =========================================================
# Step C driver: 2-way meta
# =========================================================
def build_meta(
    acgh_res: pd.DataFrame,
    mssng_res: pd.DataFrame,
    key_cols: List[str],
) -> pd.DataFrame:
    a = acgh_res[acgh_res["comparison"] == "ASD_vs_CONT"].copy()
    m = mssng_res[mssng_res["comparison"] == "ASD_vs_unaffSib"].copy()
    a = a[key_cols + ["beta", "se", "n_case", "n_ctrl", "carrier_case",
                      "carrier_ctrl", "or", "or_lo95", "or_hi95",
                      "p_value", "fit_status"]].rename(
        columns={
            "beta": "beta_acgh", "se": "se_acgh",
            "n_case": "n_case_acgh", "n_ctrl": "n_ctrl_acgh",
            "carrier_case": "carrier_case_acgh",
            "carrier_ctrl": "carrier_ctrl_acgh",
            "or": "or_acgh", "or_lo95": "or_lo95_acgh", "or_hi95": "or_hi95_acgh",
            "p_value": "p_acgh", "fit_status": "fit_status_acgh",
        })
    m = m[key_cols + ["beta", "se", "n_case", "n_ctrl", "carrier_case",
                      "carrier_ctrl", "or", "or_lo95", "or_hi95",
                      "p_value", "fit_status"]].rename(
        columns={
            "beta": "beta_mssng", "se": "se_mssng",
            "n_case": "n_case_mssng", "n_ctrl": "n_ctrl_mssng",
            "carrier_case": "carrier_case_mssng",
            "carrier_ctrl": "carrier_ctrl_mssng",
            "or": "or_mssng", "or_lo95": "or_lo95_mssng", "or_hi95": "or_hi95_mssng",
            "p_value": "p_mssng", "fit_status": "fit_status_mssng",
        })
    merged = a.merge(m, how="outer", on=key_cols)
    meta_records = []
    for _, r in merged.iterrows():
        meta = ivw_meta(r["beta_acgh"], r["se_acgh"],
                        r["beta_mssng"], r["se_mssng"])
        meta_records.append(meta)
    meta_df = pd.DataFrame(meta_records)
    return pd.concat([merged.reset_index(drop=True),
                      meta_df.reset_index(drop=True)], axis=1)


# =========================================================
# Step D driver: 3-way meta (WGS + arrayCGH + MSSNG)
# =========================================================
def build_meta3way(
    wgs_res: pd.DataFrame,
    acgh_res: pd.DataFrame,
    mssng_res: pd.DataFrame,
    key_cols: List[str],
) -> pd.DataFrame:
    """
    Merge 3 cohort results by key_cols and compute per-key IVW meta.
    WGS is pre-filtered (ASD_vs_HC, DEL) upstream. arrayCGH filter:
    ASD_vs_CONT ; MSSNG filter: ASD_vs_unaffSib.
    """
    w = wgs_res[key_cols + ["beta", "se", "n_case", "n_ctrl", "carrier_case",
                            "carrier_ctrl", "or", "or_lo95", "or_hi95",
                            "p_value", "fit_status"]].rename(
        columns={
            "beta": "beta_wgs", "se": "se_wgs",
            "n_case": "n_case_wgs", "n_ctrl": "n_ctrl_wgs",
            "carrier_case": "carrier_case_wgs",
            "carrier_ctrl": "carrier_ctrl_wgs",
            "or": "or_wgs", "or_lo95": "or_lo95_wgs", "or_hi95": "or_hi95_wgs",
            "p_value": "p_wgs", "fit_status": "fit_status_wgs",
        })
    a = acgh_res[acgh_res["comparison"] == "ASD_vs_CONT"].copy()
    a = a[key_cols + ["beta", "se", "n_case", "n_ctrl", "carrier_case",
                      "carrier_ctrl", "or", "or_lo95", "or_hi95",
                      "p_value", "fit_status"]].rename(
        columns={
            "beta": "beta_acgh", "se": "se_acgh",
            "n_case": "n_case_acgh", "n_ctrl": "n_ctrl_acgh",
            "carrier_case": "carrier_case_acgh",
            "carrier_ctrl": "carrier_ctrl_acgh",
            "or": "or_acgh", "or_lo95": "or_lo95_acgh", "or_hi95": "or_hi95_acgh",
            "p_value": "p_acgh", "fit_status": "fit_status_acgh",
        })
    m = mssng_res[mssng_res["comparison"] == "ASD_vs_unaffSib"].copy()
    m = m[key_cols + ["beta", "se", "n_case", "n_ctrl", "carrier_case",
                      "carrier_ctrl", "or", "or_lo95", "or_hi95",
                      "p_value", "fit_status"]].rename(
        columns={
            "beta": "beta_mssng", "se": "se_mssng",
            "n_case": "n_case_mssng", "n_ctrl": "n_ctrl_mssng",
            "carrier_case": "carrier_case_mssng",
            "carrier_ctrl": "carrier_ctrl_mssng",
            "or": "or_mssng", "or_lo95": "or_lo95_mssng", "or_hi95": "or_hi95_mssng",
            "p_value": "p_mssng", "fit_status": "fit_status_mssng",
        })
    merged = w.merge(a, how="outer", on=key_cols).merge(m, how="outer", on=key_cols)
    meta_records = []
    for _, r in merged.iterrows():
        meta = ivw_meta_nway(
            [r.get("beta_wgs", np.nan),
             r.get("beta_acgh", np.nan),
             r.get("beta_mssng", np.nan)],
            [r.get("se_wgs", np.nan),
             r.get("se_acgh", np.nan),
             r.get("se_mssng", np.nan)],
        )
        meta_records.append(meta)
    meta_df = pd.DataFrame(meta_records)
    return pd.concat([merged.reset_index(drop=True),
                      meta_df.reset_index(drop=True)], axis=1)


# =========================================================
# main
# =========================================================
def main() -> None:
    t0 = time.time()
    log("=" * 60)
    log("Start 40_run_replication_full_pipeline_v5.py")
    log(f"  acgh_overlap:     {ACGH_OVERLAP}")
    log(f"  acgh_covariates:  {ACGH_COVARIATES}")
    log(f"  mssng_overlap:    {MSSNG_OVERLAP}")
    log(f"  mssng_covariates: {MSSNG_COVARIATES}")
    log(f"  bin_annotation:   {BIN_ANNOTATION}")
    log(f"  wgs_L2_results:   {WGS_L2_RESULTS}   [v5: from tad04292026 v4 (top-1% QC applied)]")
    log(f"  wgs_spec_results: {WGS_SPEC_RESULTS} [v5: from tad04292026 v4 (top-1% QC applied)]")
    log(f"  outdir:           {OUT_DIR}         [v5: tad04292026/10_replication_2way_meta/output_v5]")
    log(f"  min_cell_count:   {MIN_CELL_COUNT}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ================================================================
    # Step A: burden
    # ================================================================
    log("#" * 60)
    log("Step A: burden computation (arrayCGH + MSSNG)")
    log("#" * 60)

    anno = load_bin_annotation(BIN_ANNOTATION)

    acgh_burden, acgh_qc = process_arraycgh(ACGH_OVERLAP, ACGH_COVARIATES, anno)
    out_acgh = OUT_DIR / "arraycgh_burden_L2_and_specificity_v5.tsv"
    log(f"Writing: {out_acgh}")
    acgh_burden.to_csv(out_acgh, sep="\t", index=False)

    mssng_burden, mssng_qc = process_mssng(MSSNG_OVERLAP, MSSNG_COVARIATES, anno)
    out_mssng = OUT_DIR / "mssng_burden_L2_and_specificity_v5.tsv"
    log(f"Writing: {out_mssng}")
    mssng_burden.to_csv(out_mssng, sep="\t", index=False)

    # burden QC summary
    rows = []
    for qc in (acgh_qc, mssng_qc):
        for k, v in qc.items():
            if isinstance(v, dict):
                for kk, vv in v.items():
                    rows.append({"dataset": qc["dataset"],
                                 "metric": f"{k}__{kk}", "value": vv})
            else:
                rows.append({"dataset": qc["dataset"],
                             "metric": k, "value": v})
    summary_path = OUT_DIR / "replication_burden_summary_v5.tsv"
    log(f"Writing: {summary_path}")
    pd.DataFrame(rows).to_csv(summary_path, sep="\t", index=False)

    # free memory: anno no longer needed beyond this point
    del anno

    # ================================================================
    # Step B: fit cohort-wise
    # ================================================================
    log("#" * 60)
    log("Step B: cohort-wise fit (arrayCGH B' logit + MSSNG GEE)")
    log("#" * 60)

    log("---- Analysis (i): 10 L2 classes [arrayCGH] ----")
    acgh_l2 = run_acgh_l2(acgh_burden, MIN_CELL_COUNT)
    log(f"  arraycgh L2 rows: {len(acgh_l2)}")

    log("---- Analysis (i): 10 L2 classes [MSSNG] ----")
    mssng_l2 = run_mssng_l2(mssng_burden, MIN_CELL_COUNT)
    log(f"  mssng L2 rows: {len(mssng_l2)}")

    acgh_l2["p_fdr_within_stratum"] = np.nan
    for _, idx in acgh_l2.groupby(["comparison", "exposure", "sv_type"]).groups.items():
        acgh_l2.loc[idx, "p_fdr_within_stratum"] = bh_fdr(acgh_l2.loc[idx, "p_value"])
    mssng_l2["p_fdr_within_stratum"] = np.nan
    for _, idx in mssng_l2.groupby(["comparison", "exposure", "sv_type"]).groups.items():
        mssng_l2.loc[idx, "p_fdr_within_stratum"] = bh_fdr(mssng_l2.loc[idx, "p_value"])

    acgh_l2.to_csv(OUT_DIR / "arraycgh_L2_classes_results_v5.tsv",
                   sep="\t", index=False)
    mssng_l2.to_csv(OUT_DIR / "mssng_L2_classes_results_v5.tsv",
                    sep="\t", index=False)

    log("---- Analysis (ii): specificity groups [arrayCGH] ----")
    acgh_grp = run_acgh_groups(acgh_burden, MIN_CELL_COUNT)
    log(f"  arraycgh group rows: {len(acgh_grp)}")

    log("---- Analysis (ii): specificity groups [MSSNG] ----")
    mssng_grp = run_mssng_groups(mssng_burden, MIN_CELL_COUNT)
    log(f"  mssng group rows: {len(mssng_grp)}")

    acgh_grp["p_fdr_within_stratum"] = np.nan
    for _, idx in acgh_grp.groupby(
        ["group_scheme", "comparison", "exposure", "sv_type"]
    ).groups.items():
        acgh_grp.loc[idx, "p_fdr_within_stratum"] = bh_fdr(acgh_grp.loc[idx, "p_value"])
    mssng_grp["p_fdr_within_stratum"] = np.nan
    for _, idx in mssng_grp.groupby(
        ["group_scheme", "comparison", "exposure", "sv_type"]
    ).groups.items():
        mssng_grp.loc[idx, "p_fdr_within_stratum"] = bh_fdr(mssng_grp.loc[idx, "p_value"])

    acgh_grp.to_csv(OUT_DIR / "arraycgh_specificity_groups_results_v5.tsv",
                    sep="\t", index=False)
    mssng_grp.to_csv(OUT_DIR / "mssng_specificity_groups_results_v5.tsv",
                     sep="\t", index=False)

    # ================================================================
    # Step C: IVW meta (arrayCGH + MSSNG, 2-way)
    # ================================================================
    log("#" * 60)
    log("Step C: IVW fixed-effect meta (arrayCGH ASD_vs_CONT + MSSNG ASD_vs_unaffSib)")
    log("#" * 60)

    meta_l2 = build_meta(acgh_l2, mssng_l2,
                         key_cols=["L2_class", "exposure", "sv_type"])
    meta_l2["p_fdr_within_stratum"] = np.nan
    for _, idx in meta_l2.groupby(["exposure", "sv_type"]).groups.items():
        meta_l2.loc[idx, "p_fdr_within_stratum"] = bh_fdr(meta_l2.loc[idx, "p_meta"])
    meta_l2.to_csv(OUT_DIR / "meta_L2_classes_results_v5.tsv",
                   sep="\t", index=False)

    meta_grp = build_meta(acgh_grp, mssng_grp,
                          key_cols=["group_scheme", "group_label",
                                    "exposure", "sv_type"])
    meta_grp["p_fdr_within_stratum"] = np.nan
    for _, idx in meta_grp.groupby(["group_scheme", "exposure", "sv_type"]).groups.items():
        meta_grp.loc[idx, "p_fdr_within_stratum"] = bh_fdr(meta_grp.loc[idx, "p_meta"])
    meta_grp.to_csv(OUT_DIR / "meta_specificity_groups_results_v5.tsv",
                    sep="\t", index=False)

    # ================================================================
    # Step D: IVW meta (WGS + arrayCGH + MSSNG, 3-way, ASD-focused DEL)
    # ================================================================
    log("#" * 60)
    log("Step D: 3-way IVW meta (WGS v4 ASD_vs_HC + arrayCGH v22 ASD_vs_CONT + MSSNG ASD_vs_unaffSib, DEL)")
    log("#" * 60)

    try:
        wgs_l2, wgs_sp = load_wgs_primary_results(WGS_L2_RESULTS, WGS_SPEC_RESULTS)
    except FileNotFoundError as e:
        log(f"[WARN] WGS primary result TSV not found: {e}")
        log("[WARN] Step D (3-way meta) skipped.")
        wgs_l2, wgs_sp = pd.DataFrame(), pd.DataFrame()

    if not wgs_l2.empty:
        meta3_l2 = build_meta3way(
            wgs_l2, acgh_l2, mssng_l2,
            key_cols=["L2_class", "exposure", "sv_type"],
        )
        meta3_l2["p_fdr_within_stratum"] = np.nan
        for _, idx in meta3_l2.groupby(["exposure", "sv_type"]).groups.items():
            meta3_l2.loc[idx, "p_fdr_within_stratum"] = bh_fdr(meta3_l2.loc[idx, "p_meta"])
        meta3_l2.to_csv(OUT_DIR / "meta3way_L2_classes_results_v5.tsv",
                        sep="\t", index=False)
        log(f"  meta3way L2 rows: {len(meta3_l2)}")
    else:
        meta3_l2 = pd.DataFrame()

    if not wgs_sp.empty:
        meta3_grp = build_meta3way(
            wgs_sp, acgh_grp, mssng_grp,
            key_cols=["group_scheme", "group_label", "exposure", "sv_type"],
        )
        meta3_grp["p_fdr_within_stratum"] = np.nan
        for _, idx in meta3_grp.groupby(["group_scheme", "exposure", "sv_type"]).groups.items():
            meta3_grp.loc[idx, "p_fdr_within_stratum"] = bh_fdr(meta3_grp.loc[idx, "p_meta"])
        meta3_grp.to_csv(OUT_DIR / "meta3way_specificity_groups_results_v5.tsv",
                         sep="\t", index=False)
        log(f"  meta3way spec rows: {len(meta3_grp)}")
    else:
        meta3_grp = pd.DataFrame()

    # ---- fit summary ----
    all_fit = pd.concat([
        acgh_l2[["dataset", "fit_status"]].assign(block="L2_class"),
        mssng_l2[["dataset", "fit_status"]].assign(block="L2_class"),
        acgh_grp[["dataset", "fit_status"]].assign(block="specificity_group"),
        mssng_grp[["dataset", "fit_status"]].assign(block="specificity_group"),
    ])
    summary_fit = (all_fit
                   .groupby(["dataset", "block", "fit_status"])
                   .size().reset_index(name="n"))
    summary_fit.to_csv(OUT_DIR / "replication_fit_summary_v5.tsv",
                       sep="\t", index=False)

    # ---- previews ----
    log("=" * 50)
    log("Preview: arrayCGH group_primary carrier ASD_vs_CONT")
    cond = ((acgh_grp["group_scheme"] == "group_primary")
            & (acgh_grp["exposure"] == "carrier_boundary")
            & (acgh_grp["comparison"] == "ASD_vs_CONT"))
    print(acgh_grp.loc[cond, ["group_label", "n_case", "n_ctrl",
                              "carrier_case", "carrier_ctrl",
                              "or", "or_lo95", "or_hi95", "p_value",
                              "p_fdr_within_stratum", "fit_status"]]
          .to_string(index=False))

    log("Preview: MSSNG group_primary carrier ASD_vs_unaffSib")
    cond = ((mssng_grp["group_scheme"] == "group_primary")
            & (mssng_grp["exposure"] == "carrier_boundary")
            & (mssng_grp["comparison"] == "ASD_vs_unaffSib"))
    print(mssng_grp.loc[cond, ["group_label", "n_case", "n_ctrl",
                               "carrier_case", "carrier_ctrl",
                               "or", "or_lo95", "or_hi95", "p_value",
                               "p_fdr_within_stratum", "fit_status"]]
          .to_string(index=False))

    log("Preview: META (2-way) group_primary carrier (arrayCGH + MSSNG IVW)")
    cond = ((meta_grp["group_scheme"] == "group_primary")
            & (meta_grp["exposure"] == "carrier_boundary"))
    print(meta_grp.loc[cond, ["group_label",
                              "or_acgh", "p_acgh",
                              "or_mssng", "p_mssng",
                              "or_meta", "or_meta_lo95", "or_meta_hi95",
                              "p_meta", "p_het", "i2_het",
                              "p_fdr_within_stratum", "meta_status"]]
          .to_string(index=False))

    if not meta3_grp.empty:
        log("Preview: META3way group_primary carrier (WGS v4 + arrayCGH v22 + MSSNG IVW)")
        cond = ((meta3_grp["group_scheme"] == "group_primary")
                & (meta3_grp["exposure"] == "carrier_boundary"))
        print(meta3_grp.loc[cond, ["group_label",
                                   "or_wgs", "p_wgs",
                                   "or_acgh", "p_acgh",
                                   "or_mssng", "p_mssng",
                                   "or_meta", "or_meta_lo95", "or_meta_hi95",
                                   "p_meta", "q_het", "df_het", "p_het",
                                   "i2_het", "k_included",
                                   "p_fdr_within_stratum", "meta_status"]]
              .to_string(index=False))

    elapsed = time.time() - t0
    log(f"Done. Elapsed: {elapsed:.2f}s = {elapsed/60:.2f}min")
    log("=" * 60)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        import traceback
        traceback.print_exc()
        sys.exit(1)
