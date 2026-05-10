#!/usr/bin/env python3
# ============================================================================
# ファイル名: tad_replication_mssng_v18.py
#
# 処理概要:
#   MSSNGデータ（ASD vs unaffected_sibling）を用いた
#   GRIFIN TAD boundary disruption解析の再現検証パイプライン
#   *** Full-family版: one-per-family選択を廃止し、全eligible ASD + sibling を使用 ***
#
#   - Phase 0 : (条件付き) Heffel boundary master を h5ad + diffbound から構築 (static_all用)
#   - Phase 1 : NAHR GDロカス読み込み
#   - Phase 2 : CNV抽出(cnv_merged.csv) + フィルタリング (full-family, one-per-family廃止済み)
#   - Phase 3 : L2 diffbound BED files から 10 lineage class + static_all のbin setsを構築
#   - Phase 4 : unique disrupted bin count per sample (Pattern A + Pattern C)
#   - Phase 5 : サンプルデータフレーム構築 (FAMILYID for GEE)
#   - Phase 6 : GEE logistic回帰 (Independence, family-clustered)
#   - Phase 7 : 結果出力 (両パターン, FDR補正なし)
#
# v17 から v18 への変更点 (2026-04-21 tad04212026 パイプライン移植):
#   1. _BASEDIR を tad04212026/09_mssng_sample_burden に変更。
#      MSSNG 原データ (cnv_merged.csv, subject_table.csv, Sample.tsv, 各種 BED) と
#      Heffel boundary / L2 diffbound は noncoding_tad_mssng_03132026 / heffel_deep_analysis_03242026
#      のまま read-only 参照する。
#   2. _OUTDIR を {_BASEDIR}/output_v18 に変更 (スケルトンと整合)。
#      MIN_SV_LEN != 25000 の場合は {_BASEDIR}/output_v18_{MIN_SV_LEN}bp。
#   3. 出力 TSV ファイル名サフィックスを _v17 → _v18 に更新。
#   4. script metadata (version key)、log banner、SBATCH 例を v18 に更新。
#   5. 解析ロジック (Phase 0-7, GEE Independence, 共変量, Pattern A/C, negative control) は
#      v17 から一切変更なし。
#
# v13からの主要変更点:
#   1. PRIMARY_BOUNDARY_CLASSES: 4旧クラス → 10 L2 lineage classes
#      - hpc_exc_ca, hpc_exc_dg, hpc_exc_ent, hpc_inh_cge, hpc_inh_mge,
#      - pfc_astro, pfc_exc_dl, pfc_exc_ul, pfc_inh_cge, pfc_inh_mge
#   2. DISCOVERY_SIG_CLASSES: Discovery段階で有意だった8クラス
#      - hpc_exc_dg, hpc_exc_ent, hpc_inh_cge, hpc_inh_mge,
#      - pfc_exc_dl, pfc_exc_ul, pfc_inh_cge, pfc_inh_mge
#   3. Phase 3 完全改写: L2 diffbound BED files から直接読み込み
#   4. Phase 4 パターン追加: Pattern A (NAHR excluded) + Pattern C (NAHR excluded + abs_len ≤ 1MB)
#      - Pattern C 専用共変量 (sample_total_base_c, sample_total_gene_c)
#   5. Phase 6 パターン分岐: PRIMARY (Pattern A) + SECONDARY (Pattern C) + negative control
#      - Discovery有意8クラス: 片側P, 非有意2クラス: 両側P
#   6. Phase 7 出力: pattern列を追加
#
# PRIMARY replication (10検定):
#   - ASD_vs_unaffected_sibling × DEL × 10 boundary classes × bin_count
#   - Pattern A (NAHR excluded)
#   - GEE logistic回帰 (family-clustered, Independence)
#   - Discovery有意クラス: one-sided P (beta > 0), 非有意クラス: two-sided P
#   - FDR補正なし
#
# SECONDARY Pattern C (10検定):
#   - Same as PRIMARY but using Pattern C (abs_len ≤ 1MB)
#   - analysis_type = "SECONDARY_pattern_c"
#
# Negative control / specificity (同様に両パターン):
#   - ASD × DUP × 10 boundary classes × bin_count (DEL特異性)
#   - ASD × DEL × static_all × bin_count (発達期TAD特異性)
#   - two-sided P で報告 (非有意を期待)
#
# SBATCH example:
#   sbatch --wrap="python tad_replication_mssng_v18.py" \
#     -p ncbn-cpu --account=ncbn-cpu --cpus-per-task=8 --mem=64G \
#     --job-name=tad_mssng_v18 --output=logs/tad_mssng_v18_%j.out \
#     --error=logs/tad_mssng_v18_%j.err --time=24:00:00
#
#   感度解析 (10kb):
#   sbatch --wrap="MIN_SV_LEN=10000 python tad_replication_mssng_v18.py" \
#     -p ncbn-cpu --account=ncbn-cpu --cpus-per-task=8 --mem=64G \
#     --job-name=tad_mssng_v18_10k --output=logs/tad_mssng_v18_10k_%j.out \
#     --error=logs/tad_mssng_v18_10k_%j.err --time=24:00:00
# ============================================================================

import os
import sys
import re
import time
import gzip
import json
import csv as csv_mod
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from collections import defaultdict, Counter
from typing import List, Tuple, Dict, Optional

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.genmod.families import Binomial
from statsmodels.genmod.cov_struct import Independence
from scipy.stats import chi2
import warnings
warnings.filterwarnings("ignore")

# ============================================================================
# CONFIGURATION
# ============================================================================
_BASEDIR = "/lustre12/home/kushima-pg/tad04212026/09_mssng_sample_burden"
# 補足: MSSNG 原データ (cnv_merged.csv 等) は /lustre12/home/kushima-pg/resource/mssng/ と
# /lustre12/home/kushima-pg/noncoding_tad_mssng_03132026 を read-only で参照する
# (下記 _CNV_MERGED_FILE 等の絶対パス指定はそのまま)。
_HEFFEL_DIR = "/lustre12/home/kushima-pg/heffel_deep_analysis_03242026"

# --- MSSNG固有データ ---
_CNV_MERGED_FILE = "/lustre12/home/kushima-pg/resource/mssng/cnv_merged.csv"
_SUBJECT_FILE = "/lustre12/home/kushima-pg/resource/mssng/subject_table.csv"
_ANCESTRY_FILE = "/lustre12/home/kushima-pg/resource/mssng/Sample.tsv"
_EXCLUSION_BED = "/lustre12/home/kushima-pg/resource/hg38_cnv_exclusion_regions.bed"
_SEGDUP_BED = "/lustre12/home/kushima-pg/annotationInfo/segdup_hg38_sorted_merged.bed"

# --- 共有リソース ---
_BOUNDARY_MASTER = f"{_HEFFEL_DIR}/heffel_master_outputs/heffel_boundary_master_v7.tsv.gz"
_CURATED_GD_FILE = "/lustre12/home/kushima-pg/cnv_01012026/curated_genomic_disorder_cnv_loci_v3.txt"

# --- h5ad / diffbound (master生成用) ---
_RAW_H5AD = f"{_HEFFEL_DIR}/heffel_inspect/domain_boundaries/domain_boundaries/BrainDev_raw.boundary.h5ad"
_IMPUTE_H5AD = f"{_HEFFEL_DIR}/heffel_inspect/domain_boundaries/domain_boundaries/BrainDev_impute.boundary.h5ad"
_DIFF_DIR = f"{_HEFFEL_DIR}/heffel_inspect/diff_boundaries/L2_diff_domain_boundaries"

# --- L2 diffbound BED directory ---
_L2_DIFFBOUND_DIR = f"{_HEFFEL_DIR}/heffel_inspect/diff_boundaries/L2_diff_domain_boundaries" 

# --- 出力 ---
MIN_SV_LEN = int(os.environ.get("MIN_SV_LEN", 25000))
if MIN_SV_LEN == 25000:
    _OUTDIR = f"{_BASEDIR}/output_v18"
else:
    _OUTDIR = f"{_BASEDIR}/output_v18_{MIN_SV_LEN}bp"

# --- パラメータ ---
CASE_LABEL = "ASD"
CTRL_LABEL = "unaffected_sibling"
SV_TYPES = ["DEL", "DUP"]
# 10 L2 lineage classes (primary boundary classes)
PRIMARY_BOUNDARY_CLASSES = [
    "hpc_exc_ca", "hpc_exc_dg", "hpc_exc_ent",
    "hpc_inh_cge", "hpc_inh_mge",
    "pfc_astro", "pfc_exc_dl", "pfc_exc_ul",
    "pfc_inh_cge", "pfc_inh_mge",
]

# 8 classes that were FDR<0.05 in Discovery (Pattern A, ASD×DEL)
# Non-significant in Discovery: hpc_exc_ca (FDR=0.054), pfc_astro (FDR=0.065)
DISCOVERY_SIG_CLASSES = {
    "hpc_exc_dg", "hpc_exc_ent", "hpc_inh_cge", "hpc_inh_mge",
    "pfc_exc_dl", "pfc_exc_ul", "pfc_inh_cge", "pfc_inh_mge",
}

# 2 classes that were NOT significant in Discovery (Pattern A, ASD×DEL)
DISCOVERY_NONSIG_CLASSES = {"hpc_exc_ca", "pfc_astro"}

MAX_SV_LEN_PATTERN_C = 1000000  # Pattern C: abs_len <= 1MB

UNCLEAN_MAX_PCT = 50.0
COMMON_FREQ_THR = 0.1
EXCLUSION_OVERLAP_THR = 0.50
HIGH_BURDEN_PERCENTILE = 0.99
NAHR_RO_THR = 0.50
SEGDUP_MAX_PCT = 50.0
BIN_OVERLAP_THR = 0.10

FREQ_COLS = [
    "CNVnator_percent_freq_MSSNG_parents_HiSeqX",
    "ERDS_percent_freq_MSSNG_parents_HiSeq2000",
    "CNVnator_percent_freq_MSSNG_parents_HiSeq2000",
    "ERDS_percent_freq_MSSNG_parents_HiSeqX",
    "ERDS_percent_freq_1000g",
    "CNVnator_percent_freq_1000g",
    "CG_percent_freq_MSSNG_parents",
]
ANCESTRY_DUMMIES = [f"anc_{a}" for a in ["OTH", "SAS", "EAS", "AMR", "AFR"]]
PLATFORM_DUMMIES = ["plat_NovaSeq", "plat_HiSeq", "plat_HiSeq2000", "plat_HiSeq2500", "plat_CG"]
PLATFORM_NORM = {
    "Illumina HiSeqX": "HiSeqX",
    "Illumina NovaSeq": "NovaSeq",
    "Illumina HiSeq": "HiSeq",
    "Illumina HiSeq2000": "HiSeq2000",
    "Illumina HiSeq2500": "HiSeq2500",
    "Complete Genomics": "CG",
}

# --- Discovery結果: 期待効果方向 ---
DISCOVERY_EXPECTED_DIRECTION = "positive"  # beta > 0

# ============================================================================
# LOGGING
# ============================================================================
_T0 = time.time()

def log(msg):
    elapsed = time.time() - _T0
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================
def norm_chrom(c):
    s = str(c).replace("chr", "")
    if s.isdigit():
        return f"chr{int(s)}"
    if s.upper() in ("X", "Y", "M", "MT"):
        return f"chr{s.upper()}"
    return f"chr{s}"

def calc_reciprocal_overlap(s1, e1, s2, e2):
    ov = max(0, min(e1, e2) - max(s1, s2))
    len1, len2 = e1 - s1, e2 - s2
    if len1 <= 0 or len2 <= 0:
        return 0.0
    return min(ov / len1, ov / len2)

def calc_target_coverage(s1, e1, s2, e2):
    ov = max(0, min(e1, e2) - max(s1, s2))
    len2 = e2 - s2
    return ov / len2 if len2 > 0 else 0.0

def count_genes_from_symbol(gene_str):
    if pd.isna(gene_str) or str(gene_str).strip() == "":
        return 0
    return len(set(g.strip() for g in str(gene_str).split("|") if g.strip()))

def assert_unique_keys(df, key_cols, label, n_show=10):
    dup_mask = df.duplicated(subset=key_cols, keep=False)
    if dup_mask.any():
        n_dup = int(dup_mask.sum())
        examples = df.loc[dup_mask, key_cols].head(n_show).to_string(index=False)
        raise ValueError(f"[{label}] {n_dup} duplicate rows found.\n{examples}")

# ============================================================================
# BOUNDARY MASTER GENERATION (unified_heffel_tad_pipeline_v11.py Step 1 相当)
# ============================================================================
def normalize_group_key(group_key):
    x = str(group_key)
    x = x.replace("Inh-eCGE", "Inh-CGE")
    x = x.replace("Inh-eMGE", "Inh-MGE")
    return x

def stage_rank(colname):
    c = str(colname)
    if "2T" in c: return 1
    if "3T" in c: return 2
    if "infant" in c: return 3
    if "adult" in c: return 4
    return 99

def infer_trajectory_class_from_filename(basename):
    b = str(basename)
    if "Exc-" in b: return "neuronal_exc"
    if "Astro" in b: return "glial_astro"
    if "Inh-" in b: return "neuronal_inh"
    return "other"

def infer_lineage_key_from_filename(basename):
    b = str(basename)
    b = re.sub(r"_diffbound\.bed\.gz$", "", b)
    return b

def infer_region_from_lineage_key(lineage_key):
    lk = str(lineage_key)
    if lk.startswith("PFC_"): return "PFC"
    if lk.startswith("HPC_"): return "HPC"
    return "OTHER"

def extract_nonzero_sparse_long(adata, matrix_name, value_filter=None):
    from scipy import sparse as sp
    X = adata.X.tocoo() if sp.issparse(adata.X) else sp.coo_matrix(adata.X)
    df = pd.DataFrame({"row_idx": X.row, "col_idx": X.col, "value": X.data})
    if value_filter is not None:
        df = df.loc[df["value"] == value_filter].copy()
    obs_names = np.array(adata.obs_names)
    var_df = adata.var.reset_index(drop=True).copy()
    required_var_cols = ["chrom", "start", "end"]
    missing = [c for c in required_var_cols if c not in var_df.columns]
    if missing:
        raise ValueError(f"{matrix_name}: var に必要列がありません: {missing}")
    df["group_key"] = obs_names[df["row_idx"].to_numpy()]
    df["group_key_norm"] = df["group_key"].map(normalize_group_key)
    df["chrom"] = var_df.loc[df["col_idx"], "chrom"].to_numpy()
    df["start0"] = var_df.loc[df["col_idx"], "start"].to_numpy()
    df["end"] = var_df.loc[df["col_idx"], "end"].to_numpy()
    df["bin_id"] = (
        df["chrom"].astype(str) + ":" + df["start0"].astype(str) + "-" + df["end"].astype(str)
    )
    df["matrix_name"] = matrix_name
    return df[["matrix_name", "group_key", "group_key_norm", "chrom", "start0", "end",
               "bin_id", "value"]].reset_index(drop=True)

def read_diffbound_file(path):
    df = pd.read_csv(path, sep="\t", compression="gzip")
    expected = {"chrom", "start", "end", "chi2_sc"}
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(f"{path}: 必須列がありません: {missing}")
    return df

def summarize_diffbound_row_structure(df, basename):
    meta_cols = {"chrom", "start", "end", "chi2_sc"}
    value_cols = [c for c in df.columns if c not in meta_cols]
    rg_cols = [c for c in value_cols if "RG" in c]
    if len(rg_cols) > 0:
        early_cols = sorted(rg_cols, key=stage_rank)
        late_cols = sorted([c for c in value_cols if c not in rg_cols], key=stage_rank)
        anchor_type = "RG_anchor"
    else:
        ordered = sorted(value_cols, key=stage_rank)
        if len(ordered) >= 2:
            early_cols = [ordered[0]]
            late_cols = ordered[1:]
        else:
            early_cols = ordered[:]
            late_cols = ordered[:]
        anchor_type = "temporal_noRG"
    if len(early_cols) == 0 or len(late_cols) == 0:
        ordered = sorted(value_cols, key=stage_rank)
        early_cols = ordered[:1]
        late_cols = ordered[:]
        anchor_type = "fallback_single"
    return early_cols, late_cols, anchor_type

def build_diffbound_long(path):
    basename = path.name
    lineage_key = infer_lineage_key_from_filename(basename)
    trajectory_class = infer_trajectory_class_from_filename(basename)
    region = infer_region_from_lineage_key(lineage_key)
    df = read_diffbound_file(path)
    early_cols, late_cols, anchor_type = summarize_diffbound_row_structure(df, basename)
    df["early_mean_prob"] = df[early_cols].mean(axis=1)
    df["late_mean_prob"] = df[late_cols].mean(axis=1)
    df["late_max_prob"] = df[late_cols].max(axis=1)
    df["dev_delta_prob"] = df["late_mean_prob"] - df["early_mean_prob"]

    def direction(x):
        if pd.isna(x): return "NA"
        if x > 0: return "gain"
        if x < 0: return "loss"
        return "flat"

    df["dev_direction"] = df["dev_delta_prob"].map(direction)
    meta_cols_list = [
        "chrom", "start", "end", "chi2_sc",
        "early_mean_prob", "late_mean_prob", "late_max_prob",
        "dev_delta_prob", "dev_direction"
    ]
    exclude_cols = {
        "early_mean_prob", "late_mean_prob", "late_max_prob",
        "dev_delta_prob", "dev_direction"
    }
    value_cols = [c for c in df.columns if c not in set(meta_cols_list)]
    group_cols = [c for c in value_cols if c not in exclude_cols and c != "chi2_sc"]
    long_list = []
    for group_col in group_cols:
        sub = df[meta_cols_list].copy()
        sub["group_key_from_diff"] = group_col
        sub["group_key_norm"] = sub["group_key_from_diff"].map(normalize_group_key)
        sub["group_prob_in_diff"] = df[group_col].to_numpy()
        sub["lineage_key"] = lineage_key
        sub["trajectory_class"] = trajectory_class
        sub["anchor_type"] = anchor_type
        sub["region"] = region
        sub["diff_file"] = basename
        sub["bin_id"] = (
            sub["chrom"].astype(str) + ":" + sub["start"].astype(str) + "-" + sub["end"].astype(str)
        )
        sub = sub.rename(columns={"start": "start0"})
        long_list.append(sub)
    out = pd.concat(long_list, axis=0, ignore_index=True)
    return out[[
        "group_key_from_diff", "group_key_norm", "chrom", "start0", "end", "bin_id",
        "group_prob_in_diff", "chi2_sc",
        "early_mean_prob", "late_mean_prob", "late_max_prob",
        "dev_delta_prob", "dev_direction",
        "trajectory_class", "anchor_type", "region", "lineage_key", "diff_file"
    ]]

def deduplicate_diff_annotations(df):
    if df.empty:
        return df.copy()
    tmp = df.copy()
    tmp["abs_dev_delta_prob"] = tmp["dev_delta_prob"].abs()
    tmp = tmp.sort_values(
        ["group_key_norm", "bin_id", "abs_dev_delta_prob", "chi2_sc"],
        ascending=[True, True, False, False]
    )
    tmp = tmp.drop_duplicates(subset=["group_key_norm", "bin_id"], keep="first").copy()
    tmp = tmp.drop(columns=["abs_dev_delta_prob"])
    return tmp

def merge_diff_annotation_two_pass(master, diff_ann):
    diff_cols = [
        "group_key_from_diff", "group_prob_in_diff", "chi2_sc",
        "early_mean_prob", "late_mean_prob", "late_max_prob",
        "dev_delta_prob", "dev_direction",
        "trajectory_class", "anchor_type", "region", "lineage_key", "diff_file",
    ]
    diff_for_exact = diff_ann.copy()
    diff_for_exact = diff_for_exact.rename(columns={"group_key_from_diff": "_diff_gk"})
    pass1 = master[["group_key", "bin_id"]].merge(
        diff_for_exact[["_diff_gk", "bin_id"] + [c for c in diff_cols if c != "group_key_from_diff"] + ["group_key_norm"]],
        left_on=["group_key", "bin_id"],
        right_on=["_diff_gk", "bin_id"],
        how="inner"
    )
    pass1 = pass1.rename(columns={"_diff_gk": "group_key_from_diff"})
    pass1_keys = set(zip(pass1["group_key"], pass1["bin_id"]))
    log(f"  diff merge Pass 1 (exact): {len(pass1_keys)} master rows matched")
    master_unmatched_mask = ~pd.Series(
        list(zip(master["group_key"], master["bin_id"]))
    ).isin(pass1_keys).values
    master_unmatched = master.loc[master_unmatched_mask, ["group_key", "group_key_norm", "bin_id"]].copy()
    pass2 = master_unmatched.merge(
        diff_ann[["group_key_norm", "bin_id"] + diff_cols],
        on=["group_key_norm", "bin_id"],
        how="inner"
    )
    pass2_keys = set(zip(pass2["group_key"], pass2["bin_id"]))
    log(f"  diff merge Pass 2 (normalized): {len(pass2_keys)} additional master rows matched")
    matched = pd.concat([pass1, pass2], axis=0, ignore_index=True)
    assert_unique_keys(matched, ["group_key", "bin_id"], "diff_merge_matched")
    matched_for_join = matched[["group_key", "bin_id"] + diff_cols].copy()
    result = master.merge(matched_for_join, on=["group_key", "bin_id"], how="left", validate="1:1")
    n_matched = result["group_prob_in_diff"].notna().sum()
    n_total = result.shape[0]
    log(f"  diff merge total: {n_matched}/{n_total} master rows annotated ({n_matched / n_total * 100:.1f}%)")
    return result

def build_boundary_master(master_path):
    """Phase 0: boundary masterが存在しなければh5adから生成する。"""
    master_out = Path(master_path)
    if master_out.exists():
        log(f"Boundary master already exists: {master_out}")
        log("  Skipping master generation.")
        return str(master_out)

    t0_master = time.time()
    log("=" * 70)
    log("PHASE 0: Building Heffel boundary master (h5ad → master)")
    log("  Master file not found — generating from scratch.")
    log("=" * 70)

    import anndata as ad

    raw_h5ad = Path(_RAW_H5AD)
    impute_h5ad = Path(_IMPUTE_H5AD)
    diff_dir = Path(_DIFF_DIR)
    for p in [raw_h5ad, impute_h5ad, diff_dir]:
        if not p.exists():
            raise FileNotFoundError(f"Not found: {p}")

    outdir = master_out.parent
    outdir.mkdir(parents=True, exist_ok=True)

    log(f"Loading raw h5ad: {raw_h5ad}")
    adata_raw = ad.read_h5ad(raw_h5ad)
    log(f"Loading impute h5ad: {impute_h5ad}")
    adata_impute = ad.read_h5ad(impute_h5ad)

    log("Extracting raw support=2 bins")
    raw_support2 = extract_nonzero_sparse_long(adata_raw, matrix_name="raw", value_filter=2.0)
    raw_support2 = raw_support2.rename(columns={"value": "raw_value"})
    raw_support2["is_raw_support2"] = 1
    assert_unique_keys(raw_support2, ["group_key", "bin_id"], "raw_support2")
    log(f"raw_support2: {raw_support2.shape[0]} rows")

    log("Extracting impute nonzero bins")
    impute_nonzero = extract_nonzero_sparse_long(adata_impute, matrix_name="impute", value_filter=None)
    impute_nonzero = impute_nonzero.rename(columns={"value": "impute_value"})
    impute_nonzero = impute_nonzero[["group_key", "bin_id", "impute_value"]].copy()
    assert_unique_keys(impute_nonzero, ["group_key", "bin_id"], "impute_nonzero")

    master = raw_support2.merge(impute_nonzero, on=["group_key", "bin_id"], how="left", validate="1:1")
    master["impute_value"] = master["impute_value"].fillna(0.0)
    master["is_impute_support2"] = (master["impute_value"] == 2.0).astype(int)
    master["is_consensus_support2"] = ((master["raw_value"] == 2.0) & (master["impute_value"] == 2.0)).astype(int)
    log(f"master after raw×impute merge: {master.shape[0]} rows")

    log("Parsing diffbound files")
    diff_files = sorted(Path(_DIFF_DIR).glob("*_diffbound.bed.gz"))
    if len(diff_files) == 0:
        raise FileNotFoundError(f"No diffbound files found in: {_DIFF_DIR}")

    diff_long_list = []
    for f in diff_files:
        log(f"Reading diffbound: {f.name}")
        diff_long_list.append(build_diffbound_long(f))
    diff_long_all = pd.concat(diff_long_list, axis=0, ignore_index=True)
    diff_long_all = deduplicate_diff_annotations(diff_long_all)

    diff_ann = diff_long_all[[
        "group_key_from_diff", "group_key_norm", "bin_id",
        "group_prob_in_diff", "chi2_sc",
        "early_mean_prob", "late_mean_prob", "late_max_prob",
        "dev_delta_prob", "dev_direction",
        "trajectory_class", "anchor_type", "region", "lineage_key", "diff_file",
    ]].copy()
    assert_unique_keys(diff_ann, ["group_key_norm", "bin_id"], "diff_ann")

    log("Merging developmental annotations (two-pass)")
    master = merge_diff_annotation_two_pass(master, diff_ann)
    master["overlaps_diffbound"] = master["group_prob_in_diff"].notna().astype(int)
    assert_unique_keys(master, ["group_key", "bin_id"], "master_final")

    master["is_rg_anchored"] = (master["anchor_type"] == "RG_anchor").astype(int)
    master["is_temporal_norg"] = (master["anchor_type"] == "temporal_noRG").astype(int)
    master["is_dev_gain"] = (master["dev_direction"] == "gain").astype(int)
    master["is_dev_loss"] = (master["dev_direction"] == "loss").astype(int)
    master["is_dev_exc"] = (master["trajectory_class"] == "neuronal_exc").astype(int)
    master["is_dev_astro"] = (master["trajectory_class"] == "glial_astro").astype(int)
    master["is_dev_inh"] = (master["trajectory_class"] == "neuronal_inh").astype(int)
    master["is_region_pfc"] = (master["region"] == "PFC").astype(int)
    master["is_region_hpc"] = (master["region"] == "HPC").astype(int)
    master["is_lineage_inh_cge"] = master["lineage_key"].astype(str).str.contains("Inh-CGE", regex=False, na=False).astype(int)
    master["is_lineage_inh_mge"] = master["lineage_key"].astype(str).str.contains("Inh-MGE", regex=False, na=False).astype(int)

    master["is_exc_rg_gain"] = (
        (master["is_dev_exc"] == 1) & (master["is_rg_anchored"] == 1) & (master["is_dev_gain"] == 1)
    ).astype(int)
    master["is_astro_rg_loss"] = (
        (master["is_dev_astro"] == 1) & (master["is_rg_anchored"] == 1) & (master["is_dev_loss"] == 1)
    ).astype(int)
    master["is_inh_cge_gain"] = (
        (master["is_dev_inh"] == 1) & (master["is_lineage_inh_cge"] == 1) & (master["is_dev_gain"] == 1)
    ).astype(int)
    master["is_inh_mge_gain"] = (
        (master["is_dev_inh"] == 1) & (master["is_lineage_inh_mge"] == 1) & (master["is_dev_gain"] == 1)
    ).astype(int)

    log(f"Writing boundary master: {master_out}")
    master.to_csv(master_out, sep="\t", index=False, compression="gzip")

    elapsed_master = time.time() - t0_master
    log(f"Phase 0 done. {master.shape[0]} master rows. elapsed={elapsed_master:.1f}s")
    return str(master_out)

# ============================================================================
# PHASE 1: Load NAHR GD loci
# ============================================================================
def run_phase1_load_nahr():
    log("=" * 70)
    log("PHASE 1: Loading NAHR GD loci...")
    log("=" * 70)

    col_map_key = lambda cols: {c.lower(): c for c in cols}
    gd_df = pd.read_csv(_CURATED_GD_FILE, sep="\t", dtype=str)
    cm = col_map_key(gd_df.columns)
    gd_df["nahr_bool"] = gd_df[cm["nahr"]].str.strip().str.lower().isin(["true", "1", "yes"])
    gd_nahr = gd_df[gd_df["nahr_bool"]].copy()
    gd_nahr["chrom_norm"] = gd_nahr[cm["chr"]].apply(norm_chrom)
    gd_nahr["start_int"] = gd_nahr[cm["start"]].astype(int)
    gd_nahr["end_int"] = gd_nahr[cm["end"]].astype(int)
    gd_nahr["cnv_type"] = gd_nahr[cm["cnv"]].str.upper().str.strip()
    gd_nahr["gd_id"] = gd_nahr[cm["gd_id"]].str.strip()
    log(f"  Loaded {len(gd_nahr)} NAHR GD loci")
    return gd_nahr

def is_nahr_gd_cnv(chrom, start, end, svtype, gd_nahr):
    sub = gd_nahr[(gd_nahr["chrom_norm"] == chrom) & (gd_nahr["cnv_type"] == svtype)]
    for _, gd in sub.iterrows():
        gid = gd["gd_id"]
        gs, ge = gd["start_int"], gd["end_int"]
        if start < ge and end > gs:
            if "15q11.2" in gid and "15q13" not in gid:
                if calc_target_coverage(start, end, gs, ge) >= NAHR_RO_THR:
                    return True
            else:
                if calc_reciprocal_overlap(start, end, gs, ge) >= NAHR_RO_THR:
                    return True
    return False

# ============================================================================
# PHASE 2: Extract + filter MSSNG CNV data + one-per-family
# ============================================================================
def run_phase2_load_cnv(gd_nahr):
    """Phase 2: cnv_merged.csvからCNV抽出 + フィルタリング + one-per-family"""
    log("=" * 70)
    log("PHASE 2: Extracting and filtering MSSNG CNV data")
    log(f"  Source: {_CNV_MERGED_FILE}")
    log("=" * 70)

    # --- subject_table読み込み ---
    subjects = {}
    with open(_SUBJECT_FILE) as f:
        for row in csv_mod.DictReader(f):
            subjects[row['INDEXID']] = row

    def get_sex(sample_id):
        s = subjects.get(sample_id) or subjects.get(sample_id.rstrip('AB'))
        if s is None:
            return np.nan
        return 1.0 if s['SEX'] == 'M' else 0.0

    # 親IDセット
    parent_ids = set()
    for s in subjects.values():
        if s['FATHERID'] not in ('0', ''):
            parent_ids.add(s['FATHERID'])
        if s['MOTHERID'] not in ('0', ''):
            parent_ids.add(s['MOTHERID'])

    # ASD / unaffected_sibling 判定
    sample_status = {}
    for sid, s in subjects.items():
        if s['AFFECTION'] == '2':
            sample_status[sid] = 'ASD'
        elif s['AFFECTION'] == '1':
            has_parents = s['FATHERID'] not in ('0', '') or s['MOTHERID'] not in ('0', '')
            is_parent = sid in parent_ids
            if has_parents and not is_parent:
                sample_status[sid] = 'unaffected_sibling'

    def get_status(sample_id):
        st = sample_status.get(sample_id)
        if st:
            return st
        return sample_status.get(sample_id.rstrip('AB'))

    # --- Ancestry / Platform ---
    ancestry_map = {}
    platform_map = {}
    with open(_ANCESTRY_FILE, encoding="utf-8", errors="replace") as f:
        for row in csv_mod.DictReader(f, delimiter="\t"):
            sid = row["SUBMITTEDID"]
            anc = row.get("PREDICTED_ANCESTRY", "")
            plat = row.get("PLATFORM", "")
            if anc:
                ancestry_map[sid] = anc
                idx = row.get("INDEXID", "")
                if idx and idx != sid:
                    ancestry_map[idx] = anc
            if plat:
                platform_map[sid] = plat
                idx = row.get("INDEXID", "")
                if idx and idx != sid:
                    platform_map[idx] = plat

    def get_ancestry(sample_id):
        a = ancestry_map.get(sample_id) or ancestry_map.get(sample_id.rstrip("AB"))
        return a if a else "OTH"

    def get_platform_norm(sample_id):
        raw = platform_map.get(sample_id) or platform_map.get(sample_id.rstrip("AB"))
        raw = raw if raw else ""
        return PLATFORM_NORM.get(raw, "HiSeqX")

    log(f"  Loaded {len(subjects)} subjects, ancestry for {len(ancestry_map)}, platform for {len(platform_map)}")

    # --- CNV抽出 (extract_cnv_mssng_v2.py 統合) ---
    log("  Extracting CNVs from cnv_merged.csv (QC=ok, high_quality_rare, ASD/sibling)...")
    sv_rows = []
    n_skipped = 0
    with open(_CNV_MERGED_FILE) as f:
        reader = csv_mod.DictReader(f)
        for row in reader:
            if row['sample_QC'] != 'ok':
                n_skipped += 1
                continue
            if row['high_quality_rare'] != 'high quality rare':
                n_skipped += 1
                continue
            status = get_status(row['sample'])
            if status is None:
                n_skipped += 1
                continue
            row['status'] = status
            sv_rows.append(row)

    sv = pd.DataFrame(sv_rows)
    log(f"  Extracted: {len(sv)} CNVs (skipped: {n_skipped})")

    # --- 基本フィルタ ---
    sv["SV_chrom"] = sv["chr"].apply(norm_chrom)
    sv["abs_len"] = sv["size"].astype(float).abs()
    sv["sv_start"] = np.minimum(sv["start_position"].astype(int), sv["end_position"].astype(int))
    sv["sv_end"] = np.maximum(sv["start_position"].astype(int), sv["end_position"].astype(int))

    autosome_chroms = {f"chr{i}" for i in range(1, 23)}
    sv = sv[sv["SV_chrom"].isin(autosome_chroms)].copy()
    sv = sv[sv["abs_len"] >= MIN_SV_LEN].copy()
    sv = sv[sv["CNV_type"].isin(SV_TYPES)].copy()
    log(f"  After autosome/size({MIN_SV_LEN / 1000:.0f}kb)/type: {len(sv)}")

    # 頻度フィルタ
    for fc in FREQ_COLS:
        sv[fc + "_num"] = pd.to_numeric(sv[fc], errors="coerce").fillna(0)
    freq_max = sv[[fc + "_num" for fc in FREQ_COLS]].max(axis=1)
    sv = sv[freq_max < COMMON_FREQ_THR].copy()

    sv["Sex_numeric"] = sv["sample"].apply(get_sex)
    sv = sv.dropna(subset=["Sex_numeric"]).copy()

    sv["unclean_pct"] = pd.to_numeric(sv["unclean_genome_percent_overlap"], errors="coerce").fillna(0)
    sv = sv[sv["unclean_pct"] < UNCLEAN_MAX_PCT].copy()
    log(f"  After freq/sex/unclean: {len(sv)}")

    # --- Exclusion BED ---
    excl_raw = defaultdict(list)
    with open(_EXCLUSION_BED) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            c = norm_chrom(parts[0])
            try:
                excl_raw[c].append((int(parts[1]), int(parts[2])))
            except Exception:
                continue

    def check_excl_overlap(chrom, start, end, excl_dict, thr):
        sv_len = end - start
        if sv_len <= 0 or chrom not in excl_dict:
            return False
        cov = 0
        for es, ee in excl_dict[chrom]:
            ov = max(0, min(end, ee) - max(start, es))
            cov += ov
        return (cov / sv_len) >= thr

    excl_mask = np.array([check_excl_overlap(
        sv["SV_chrom"].values[i], sv["sv_start"].values[i], sv["sv_end"].values[i],
        excl_raw, EXCLUSION_OVERLAP_THR) for i in range(len(sv))])
    sv = sv[~excl_mask].reset_index(drop=True)
    log(f"  After exclusion BED: {len(sv)}")
    del excl_raw

    # --- Segdup filter ---
    log("  Loading segdup BED...")
    segdup_raw = defaultdict(list)
    with open(_SEGDUP_BED) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            c = norm_chrom(parts[0])
            try:
                segdup_raw[c].append((int(parts[1]), int(parts[2])))
            except Exception:
                continue

    segdup_merged = {}
    for chrom, intervals in segdup_raw.items():
        intervals.sort()
        merged = []
        for s, e in intervals:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        segdup_merged[chrom] = merged
    del segdup_raw

    def check_segdup_overlap(chrom, start, end, segdup_dict, thr_pct):
        sv_len = end - start
        if sv_len <= 0 or chrom not in segdup_dict:
            return False
        cov = 0
        for ss, se in segdup_dict[chrom]:
            if ss >= end:
                break
            ov = max(0, min(end, se) - max(start, ss))
            cov += ov
        return (cov / sv_len * 100.0) >= thr_pct

    n_before_segdup = len(sv)
    segdup_mask = np.array([check_segdup_overlap(
        sv["SV_chrom"].values[i], sv["sv_start"].values[i], sv["sv_end"].values[i],
        segdup_merged, SEGDUP_MAX_PCT) for i in range(len(sv))])
    sv = sv[~segdup_mask].reset_index(drop=True)
    log(f"  After segdup filter (<{SEGDUP_MAX_PCT}%): {n_before_segdup} -> {len(sv)} ({int(segdup_mask.sum())} removed)")
    del segdup_merged

    # --- NAHR flag ---
    nahr_flags = np.array([is_nahr_gd_cnv(
        sv["SV_chrom"].values[i], sv["sv_start"].values[i], sv["sv_end"].values[i],
        sv["CNV_type"].values[i], gd_nahr) for i in range(len(sv))])
    sv["is_nahr_gd"] = nahr_flags.astype(int)
    log(f"  NAHR GD flagged: {int(sv['is_nahr_gd'].sum())}/{len(sv)}")

    sv["gene_count"] = sv["gene_symbol"].apply(count_genes_from_symbol)

    # --- One-per-family (全サンプルベース: subject_table から構築) ---
    # v12バグ修正: CNV保有サンプルだけでなく、subject_tableの全ASD/siblingから選択
    log("  One-per-family selection (subject_table-based)...")

    # subject_tableから全ASD/sibling + FAMILYID を取得
    all_sample_status = {}  # sid -> "ASD" or "unaffected_sibling"
    all_fid_map = {}  # sid -> FAMILYID
    for sid, s in subjects.items():
        status = sample_status.get(sid)
        if status is None:
            continue
        if np.isnan(get_sex(sid)):
            continue
        all_sample_status[sid] = status
        fid = s.get("FAMILYID", "")
        all_fid_map[sid] = fid if fid else f"_solo_{sid}"

    log(f"  Total eligible samples (with sex): ASD={sum(1 for v in all_sample_status.values() if v == CASE_LABEL)}, "
        f"sib={sum(1 for v in all_sample_status.values() if v == CTRL_LABEL)}")

    # High-burden除外: CNVデータに基づく
    hi_sample_set = set()
    if len(sv) > 0:
        sample_cnv_counts = sv.groupby("sample").size()
        thr_99 = sample_cnv_counts.quantile(HIGH_BURDEN_PERCENTILE)
        hi_sample_set = set(sample_cnv_counts[sample_cnv_counts >= thr_99].index)
        log(f"  High-burden samples (>={thr_99:.0f} CNVs): {len(hi_sample_set)}")
        # Remove from all_sample_status and CNV data
        for sid in hi_sample_set:
            all_sample_status.pop(sid, None)
        sv = sv[~sv["sample"].isin(hi_sample_set)].copy()
        log(f"  After high-burden excl: {len(sv)} CNVs")

    family_asd, family_sib = defaultdict(list), defaultdict(list)
    for sid, status in all_sample_status.items():
        fid = all_fid_map.get(sid, f"_solo_{sid}")
        if status == CASE_LABEL:
            family_asd[fid].append(sid)
        elif status == CTRL_LABEL:
            family_sib[fid].append(sid)

    # v15: full-family - 全eligible ASD/siblingを保持 (GEEがFAMILYIDで相関処理)
    selected_samples = set(all_sample_status.keys())
    n_asd_all = sum(1 for v in all_sample_status.values() if v == CASE_LABEL)
    n_sib_all = sum(1 for v in all_sample_status.values() if v == CTRL_LABEL)
    log(f"  v18 full-family: ASD={n_asd_all}, sib={n_sib_all}, total={len(selected_samples)}")

    # CNVを選択サンプルに絞る
    sv = sv[sv["sample"].isin(selected_samples)].copy()
    log(f"  CNVs after eligible filter: {len(sv)}")

    # 最終sample_status_dictとfid_mapを選択サンプルに絞る
    final_sample_status = dict(all_sample_status)
    final_fid_map = dict(all_fid_map)

    all_fids = set(family_asd.keys()) | set(family_sib.keys())
    n_multi_asd_fam = sum(1 for sids in family_asd.values() if len(sids) > 1)
    n_multi_sib_fam = sum(1 for sids in family_sib.values() if len(sids) > 1)
    log(f"  Families: {len(all_fids)} total, {n_multi_asd_fam} w/ 2+ ASD, {n_multi_sib_fam} w/ 2+ sib")

    # Platform分布
    log("  Platform distribution:")
    asd_set = {sid for sid, s in all_sample_status.items() if s == CASE_LABEL}
    sib_set = {sid for sid, s in all_sample_status.items() if s == CTRL_LABEL}
    for label, sample_set in [("ASD", asd_set), ("Sibling", sib_set)]:
        pcounts = Counter()
        for sid in sample_set:
            pcounts[get_platform_norm(sid)] += 1
        log(f"    {label} (n={len(sample_set)}):")
        for p, n in sorted(pcounts.items(), key=lambda x: -x[1]):
            log(f"      {p}: {n} ({100 * n / len(sample_set):.1f}%)")

    return sv, final_sample_status, final_fid_map, get_sex, get_ancestry, get_platform_norm

# ============================================================================
# PHASE 3: Load boundary bins
# ============================================================================
def run_phase3_load_bins():
    """Phase 3: L2 diffbound BEDから10 lineage class + static_all のbin setsを構築"""
    log("=" * 70)
    log("PHASE 3: Loading boundary bins from L2 diffbound BED files")
    log("=" * 70)

    import anndata as ad

    # Read all diffbound BED files
    diff_dir = Path(_L2_DIFFBOUND_DIR)
    diff_files = sorted(diff_dir.glob("*_diffbound.bed.gz"))
    if len(diff_files) == 0:
        raise FileNotFoundError(f"No diffbound files found in: {diff_dir}")

    # Build bin_sets per lineage class
    all_diffbound_bins = set()
    raw_bin_sets = {}
    for f in diff_files:
        lineage_key = f.name.replace("_diffbound.bed.gz", "")
        safe_key = lineage_key.lower().replace("-", "_")
        log(f"  Reading {f.name}")
        df = pd.read_csv(f, sep="\t", compression="gzip", usecols=["chrom", "start", "end"])
        df["chrom"] = df["chrom"].apply(norm_chrom)
        df["bin_id"] = df["chrom"].astype(str) + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
        bins_df = df[["chrom", "start", "end", "bin_id"]].drop_duplicates("bin_id").copy()
        bins_df = bins_df.rename(columns={"start": "start0"})
        raw_bin_sets[safe_key] = bins_df
        all_diffbound_bins.update(bins_df["bin_id"].tolist())
        log(f"    {safe_key}: {bins_df['bin_id'].nunique()} unique bins")

    # Filter by min_bin_threshold
    min_bin_threshold = 500
    bin_sets = {}
    excluded_classes = []
    for key, df in raw_bin_sets.items():
        n = df["bin_id"].nunique()
        if n >= min_bin_threshold:
            bin_sets[key] = df
            log(f"  {key}: {n} bins")
        else:
            excluded_classes.append((key, n))
            log(f"  {key}: {n} bins -> EXCLUDED (< {min_bin_threshold})")

    # Verify PRIMARY_BOUNDARY_CLASSES are all present
    for bc in PRIMARY_BOUNDARY_CLASSES:
        if bc not in bin_sets:
            raise ValueError(f"Expected boundary class {bc} not found in L2 diffbound data. "
                             f"Excluded classes: {excluded_classes}")

    # Static bins from raw h5ad (support=2, exclude diffbound bins)
    log("  Reading static bins from raw h5ad (support=2, excluding diffbound)")
    adata_raw = ad.read_h5ad(_RAW_H5AD)
    raw_support2 = extract_nonzero_sparse_long(adata_raw, "raw", value_filter=2.0)
    raw_support2["chrom"] = raw_support2["chrom"].apply(norm_chrom)
    static_bins = raw_support2.loc[~raw_support2["bin_id"].isin(all_diffbound_bins),
                                    ["chrom", "start0", "end", "bin_id"]].drop_duplicates("bin_id")
    bin_sets["static_all"] = static_bins
    log(f"  static_all: {len(static_bins)} bins (after excluding diffbound)")

    # Build bin_index (sorted arrays per chromosome for fast lookup)
    bin_index = {}
    for bclass, df in bin_sets.items():
        by_chrom = {}
        for chrom, grp in df.groupby("chrom"):
            starts = grp["start0"].values.astype(np.int64)
            ends = grp["end"].values.astype(np.int64)
            bin_ids = grp["bin_id"].values
            order = np.argsort(starts)
            by_chrom[chrom] = (starts[order], ends[order], bin_ids[order])
        bin_index[bclass] = by_chrom

    return bin_index

# ============================================================================
# PHASE 4: Compute unique disrupted bin counts per sample
# ============================================================================
def run_phase4_bin_counts(sv, sample_status_dict, bin_index, fid_map=None,
                            get_sex=None, get_platform_norm=None, get_ancestry=None):
    """Phase 4: サンプルごとのunique disrupted bin count (Pattern A + Pattern C)
    v19 (patch v3): 末尾で constraint enrichment 用 dual dump を出力する。
    fid_map / get_sex / get_platform_norm / get_ancestry は dump で使うのみ、v18 の既存ロジックには影響しない。
    """
    log("=" * 70)
    log("PHASE 4: Computing unique disrupted bin counts (Pattern A + Pattern C)")
    log("=" * 70)

    ALL_BCLASSES = PRIMARY_BOUNDARY_CLASSES + ["static_all"]
    all_samples = sorted(sample_status_dict.keys())
    sample_to_idx = {s: i for i, s in enumerate(all_samples)}
    N_samples = len(all_samples)

    # Pattern A (NAHR excluded)
    disrupted_bins_a = {}
    for bclass in ALL_BCLASSES:
        disrupted_bins_a[bclass] = {svt: [set() for _ in range(N_samples)] for svt in SV_TYPES}

    # Pattern C (NAHR excluded + abs_len <= 1MB)
    disrupted_bins_c = {}
    for bclass in ALL_BCLASSES:
        disrupted_bins_c[bclass] = {svt: [set() for _ in range(N_samples)] for svt in SV_TYPES}

    sample_total_base = {svt: np.zeros(N_samples) for svt in SV_TYPES}
    sample_total_gene = {svt: np.zeros(N_samples) for svt in SV_TYPES}

    # Pattern C 専用共変量 (NAHR除外 + abs_len <= 1MB のCNVのみ)
    sample_total_base_c = {svt: np.zeros(N_samples) for svt in SV_TYPES}
    sample_total_gene_c = {svt: np.zeros(N_samples) for svt in SV_TYPES}

    # v19 patch v4: SV-interval-aware records for exact direct-overlap exclusion
    # (constraint enrichment v9 で使用; per-(sample, bin) で SV interval を保持)
    _sv_records_a = []

    for _, row in sv.iterrows():
        si = sample_to_idx.get(row["sample"])
        if si is None:
            continue
        svt = row["CNV_type"]
        if svt not in SV_TYPES:
            continue
        chrom = row["SV_chrom"]
        sv_s, sv_e = int(row["sv_start"]), int(row["sv_end"])
        abs_len = float(row["abs_len"])
        is_nahr = bool(row["is_nahr_gd"])

        # Pattern A only (NAHR excluded)
        if is_nahr:
            continue
        sample_total_base[svt][si] += abs_len
        sample_total_gene[svt][si] += int(row["gene_count"])

        # Pattern C 専用共変量: <= 1MB のCNVのみ加算
        if abs_len <= MAX_SV_LEN_PATTERN_C:
            sample_total_base_c[svt][si] += abs_len
            sample_total_gene_c[svt][si] += int(row["gene_count"])

        # Pattern A: all qualifying SVs
        # v19 patch v4: collect bin set across all bclasses for this SV (deduped),
        # then emit one SV-aware record per (sample, bin) tuple.
        _bins_this_sv_a = set()
        for bclass in ALL_BCLASSES:
            if chrom not in bin_index[bclass]:
                continue
            starts, ends, bid_arr = bin_index[bclass][chrom]
            lo = np.searchsorted(ends, sv_s, side="right")
            hi = np.searchsorted(starts, sv_e, side="left")
            for j in range(lo, hi):
                if sv_s < ends[j] and sv_e > starts[j]:
                    ov = min(sv_e, ends[j]) - max(sv_s, starts[j])
                    bin_len = ends[j] - starts[j]
                    if ov / bin_len >= BIN_OVERLAP_THR:
                        disrupted_bins_a[bclass][svt][si].add(bid_arr[j])
                        _bins_this_sv_a.add(bid_arr[j])
        # v19 patch v4: emit SV-aware records (sample, bin, sv interval, sv_type)
        for _bid in _bins_this_sv_a:
            _sv_records_a.append((row["sample"], _bid, chrom, sv_s, sv_e, svt))

        # Pattern C: NAHR excluded + abs_len <= 1MB
        if abs_len <= MAX_SV_LEN_PATTERN_C:
            for bclass in ALL_BCLASSES:
                if chrom not in bin_index[bclass]:
                    continue
                starts, ends, bid_arr = bin_index[bclass][chrom]
                lo = np.searchsorted(ends, sv_s, side="right")
                hi = np.searchsorted(starts, sv_e, side="left")
                for j in range(lo, hi):
                    if sv_s < ends[j] and sv_e > starts[j]:
                        ov = min(sv_e, ends[j]) - max(sv_s, starts[j])
                        bin_len = ends[j] - starts[j]
                        if ov / bin_len >= BIN_OVERLAP_THR:
                            disrupted_bins_c[bclass][svt][si].add(bid_arr[j])

    # Compute bin counts for both patterns
    bin_counts_a = {}
    for bclass in ALL_BCLASSES:
        bin_counts_a[bclass] = {}
        for svt in SV_TYPES:
            counts = np.array([len(disrupted_bins_a[bclass][svt][si]) for si in range(N_samples)])
            bin_counts_a[bclass][svt] = counts

    bin_counts_c = {}
    for bclass in ALL_BCLASSES:
        bin_counts_c[bclass] = {}
        for svt in SV_TYPES:
            counts = np.array([len(disrupted_bins_c[bclass][svt][si]) for si in range(N_samples)])
            bin_counts_c[bclass][svt] = counts

    # ========================================================================
    # === V19 ADDITION (patch v3): Dual dump for constraint enrichment analysis
    # === Output:
    # ===  1. mssng_event_bins_dumped_v1.tsv.gz (sample x bin x sv_type x pattern x Diagnosis)
    # ===  2. mssng_sample_covariates_dumped_v1.tsv.gz
    # === Robust Sex mapping: M/F, 1/0, 1/2, Male/Female, true/false 対応
    # ========================================================================
    log("=" * 70)
    log("V19 ADDITION (patch v3): Dumping event-bin + sample covariate")
    log("=" * 70)
    import os as _os
    _OUT_DIR = "/lustre12/home/kushima-pg/tad04212026/14_constraint_enrichment_v1/output_v1"
    _os.makedirs(_OUT_DIR, exist_ok=True)

    def _sex_to_numeric(x):
        """v19 patch v3 (bug fix): Map various Sex encodings to 1=Male, 0=Female, NaN=unknown.
        Robust to:
          - numeric (int/float/np.nan): 1.0/0.0/2.0 as returned by get_sex()
          - string: 'M'/'F', 'Male'/'Female', '1'/'0'/'2', 'TRUE'/'FALSE'
        v18 の get_sex() は 1.0/0.0/np.nan を返すため numeric pathway が primary.
        """
        if x is None:
            return float("nan")
        # numeric pathway (handles 1.0/0.0/np.nan from v18 get_sex)
        try:
            v = float(x)
            if np.isnan(v):
                return float("nan")
            vi = int(round(v))
            if vi == 1:
                return 1
            if vi in (0, 2):
                return 0
            return float("nan")
        except (TypeError, ValueError):
            pass
        # string pathway (raw subject_table 'M'/'F'/'Male'/'Female')
        s = str(x).strip().upper()
        if s.startswith("M"):
            return 1
        if s.startswith("F"):
            return 0
        if s == "TRUE":
            return 1
        if s == "FALSE":
            return 0
        return float("nan")

    def _diag_label(_status):
        # v19 patch v3 (bug fix): v18 sample_status は 'unaffected_sibling' を使うため
        # それを Sibling にマップする (line 611 of v18.py).
        if _status in ("Affected", "ASD", "Proband"):
            return "ASD"
        elif _status in ("Sibling", "Unaffected", "unaffected_sibling"):
            return "Sibling"
        return _status

    # ---- Dump A: event-bin records (sample, bin_id, sv_type, pattern, Diagnosis) ----
    # v19 patch v4: 既存 dump file (mssng_event_bins_dumped_v1.tsv.gz) に
    # sv_chr / sv_start / sv_end 列を追加して上書き出力。downstream で exact
    # direct-overlap exclusion (gene_start < sv_end AND gene_end > sv_start)
    # を実装するための情報。同じ (sample, bin) でも複数 SV record があり得る。
    _evbin_records = []
    for _sample, _bid, _chrom, _sv_s, _sv_e, _svt in _sv_records_a:
        _diag = _diag_label(sample_status_dict.get(_sample, ""))
        _evbin_records.append({
            "sample_id": _sample,
            "bin_id": _bid,
            "sv_type": _svt,
            "pattern": "A",
            "Diagnosis": _diag,
            "sv_chr": _chrom,
            "sv_start": _sv_s,
            "sv_end": _sv_e,
        })

    _evbin_path = _os.path.join(_OUT_DIR, "mssng_event_bins_dumped_v1.tsv.gz")
    pd.DataFrame(_evbin_records).to_csv(_evbin_path, sep="\t", index=False, compression="gzip")
    log(f"  Saved {len(_evbin_records)} event-bin records (with SV interval): {_evbin_path}")

    # ---- Dump B: sample-level covariate file (REQUIRED by 14_fit v5) ----
    _cov_records = []
    for _si, _sample in enumerate(all_samples):
        _diag = _diag_label(sample_status_dict.get(_sample, ""))
        _rec = {
            "sample_id": _sample,
            "Diagnosis": _diag,
            "log1p_total_del_bases": float(np.log1p(sample_total_base["DEL"][_si])),
            "log1p_total_gene_DEL": float(np.log1p(sample_total_gene["DEL"][_si])),
        }
        # FAMILYID
        if fid_map is not None:
            _rec["FAMILYID"] = fid_map.get(_sample, f"_solo_{_sample}")
        else:
            _rec["FAMILYID"] = f"_solo_{_sample}"
        # Sex (v3: robust mapping)
        if get_sex is not None:
            _rec["Sex"] = _sex_to_numeric(get_sex(_sample))
        else:
            _rec["Sex"] = float("nan")
        # Platform (string; downstream sanitizes + dummies)
        if get_platform_norm is not None:
            _rec["Platform"] = get_platform_norm(_sample) or ""
        else:
            _rec["Platform"] = ""
        # Ancestry (string; downstream sanitizes + dummies)
        if get_ancestry is not None:
            _rec["ancestry"] = get_ancestry(_sample) or ""
        else:
            _rec["ancestry"] = ""
        _cov_records.append(_rec)

    _cov_path = _os.path.join(_OUT_DIR, "mssng_sample_covariates_dumped_v1.tsv.gz")
    pd.DataFrame(_cov_records).to_csv(_cov_path, sep="\t", index=False, compression="gzip")
    log(f"  Saved {len(_cov_records)} sample covariate records: {_cov_path}")

    # ---- v3 Sanity: Sex NaN count ----
    _df_chk = pd.DataFrame(_cov_records)
    _n_sex_na = _df_chk["Sex"].isna().sum()
    log(f"  Sex NaN count: {_n_sex_na} / {len(_df_chk)}")
    if _n_sex_na > 0:
        log(f"  WARNING: {_n_sex_na} samples have Sex=NaN. "
            f"Verify subject_table.csv encoding before downstream pipeline.")

    del _evbin_records, _cov_records, _df_chk, _sv_records_a
    # === END V19 ADDITION (patch v3 + v4) ===

    del disrupted_bins_a, disrupted_bins_c

    log("--- Pattern A (NAHR excluded) ---")
    for bclass in ALL_BCLASSES:
        for svt in SV_TYPES:
            c = bin_counts_a[bclass][svt]
            n_exp = int(np.sum(c > 0))
            mean_c = float(c[c > 0].mean()) if n_exp > 0 else 0
            max_c = int(c.max())
            log(f"  {bclass}/{svt}: n_exposed={n_exp}, mean_bins={mean_c:.2f}, max={max_c}")

    log("--- Pattern C (NAHR excluded + abs_len <= 1MB) ---")
    for bclass in ALL_BCLASSES:
        for svt in SV_TYPES:
            c = bin_counts_c[bclass][svt]
            n_exp = int(np.sum(c > 0))
            mean_c = float(c[c > 0].mean()) if n_exp > 0 else 0
            max_c = int(c.max())
            log(f"  {bclass}/{svt}: n_exposed={n_exp}, mean_bins={mean_c:.2f}, max={max_c}")

    return (all_samples, sample_to_idx, bin_counts_a, bin_counts_c,
            sample_total_base, sample_total_gene,
            sample_total_base_c, sample_total_gene_c)

# ============================================================================
# PHASE 5: Build sample dataframe with FAMILYID for GEE
# ============================================================================
def run_phase5_sample_df(all_samples, sample_status_dict, fid_map,
                         get_sex, get_ancestry, get_platform_norm):
    log("=" * 70)
    log("PHASE 5: Building sample dataframe with FAMILYID")
    log("=" * 70)

    rows = []
    for si, sid in enumerate(all_samples):
        sex = get_sex(sid)
        if np.isnan(sex):
            continue
        anc = get_ancestry(sid)
        anc_d = [1.0 if anc == a else 0.0 for a in ["OTH", "SAS", "EAS", "AMR", "AFR"]]
        plat_norm = get_platform_norm(sid)
        plat_d = [1.0 if plat_norm == p.replace("plat_", "") else 0.0 for p in PLATFORM_DUMMIES]
        fid = fid_map.get(sid)
        fid = fid if pd.notna(fid) and fid else f"_solo_{sid}"
        rows.append([si, sid, sample_status_dict[sid], sex, fid] + anc_d + plat_d)

    cols = ["sample_idx", "SampleID", "Status", "Sex", "FAMILYID"] + ANCESTRY_DUMMIES + PLATFORM_DUMMIES
    sample_df = pd.DataFrame(rows, columns=cols)
    case_mask = (sample_df["Status"] == CASE_LABEL).values
    ctrl_mask = (sample_df["Status"] == CTRL_LABEL).values
    keep = case_mask | ctrl_mask
    y = case_mask[keep].astype(int)
    n_case, n_ctrl = int(y.sum()), int((1 - y).sum())
    log(f"  n_case={n_case}, n_ctrl={n_ctrl}, total={n_case + n_ctrl}")
    valid_si = sample_df.loc[keep, "sample_idx"].values
    family_ids = sample_df.loc[keep, "FAMILYID"].values
    n_clusters = len(np.unique(family_ids))
    log(f"  n_clusters (families): {n_clusters}")

    return sample_df, y, keep, valid_si, family_ids, n_case, n_ctrl, n_clusters

# ============================================================================
# PHASE 6: GEE logistic regression (PRIMARY replication + negative control)
# ============================================================================
def run_phase6_regression(sample_df, y, keep, valid_si, family_ids,
                          bin_counts_a, bin_counts_c,
                          sample_total_base, sample_total_gene,
                          sample_total_base_c, sample_total_gene_c,
                          n_case, n_ctrl, n_clusters):
    """Phase 6: GEE logistic回帰 (Independence, family-clustered)
    PRIMARY replication (Pattern A): DEL x 10 dev boundary classes
      - Discovery sig8: one-sided P (beta > 0), nonsig2: two-sided P
    SECONDARY Pattern C: Same but using Pattern C bin counts + covariates
    Negative control (both patterns): DUP x 10 dev classes + DEL x static_all
      - always two-sided P, P_onesided = NaN
    """
    log("=" * 70)
    log("PHASE 6: GEE regression analyses (Pattern A + Pattern C)")
    log("  PRIMARY replication (Pattern A): DEL x 10 dev boundary classes")
    log("  SECONDARY Pattern C: DEL x 10 dev boundary classes (Pattern C bins)")
    log("  Negative control: DUP x 10 dev classes + DEL x static_all (two-sided P)")
    log("  Inference: GEE (Binomial, Independence, family-clustered)")
    log("=" * 70)

    ALL_BCLASSES = PRIMARY_BOUNDARY_CLASSES + ["static_all"]
    all_results = []

    anc_vals = sample_df.loc[keep, ANCESTRY_DUMMIES].values
    sex_vals = sample_df.loc[keep, "Sex"].values
    plat_vals = sample_df.loc[keep, PLATFORM_DUMMIES].values

    for pattern, bin_counts in [("A", bin_counts_a), ("C", bin_counts_c)]:
        log(f"\n{'=' * 70}")
        log(f"Pattern {pattern}")
        log(f"{'=' * 70}")

        # Pattern-specific covariates
        if pattern == "C":
            _total_base = sample_total_base_c
            _total_gene = sample_total_gene_c
            log("  Using Pattern C-specific covariates (NAHR excluded + abs_len <= 1MB)")
        else:
            _total_base = sample_total_base
            _total_gene = sample_total_gene

        for svt in SV_TYPES:
            log(f"\n--- SV_type: {svt} ---")
            log1p_base = np.log1p(_total_base[svt][valid_si])
            log1p_gene = np.log1p(_total_gene[svt][valid_si])

            for bclass in ALL_BCLASSES:
                # Classification
                is_primary = (pattern == "A" and svt == "DEL" and bclass in PRIMARY_BOUNDARY_CLASSES)
                is_secondary_c = (pattern == "C" and svt == "DEL" and bclass in PRIMARY_BOUNDARY_CLASSES)
                is_negative_ctrl = (svt == "DUP" and bclass in PRIMARY_BOUNDARY_CLASSES) or \
                                   (svt == "DEL" and bclass == "static_all")

                if not is_primary and not is_secondary_c and not is_negative_ctrl:
                    # DUP x static_all is not needed
                    continue

                bsv_col = bin_counts[bclass][svt][valid_si].astype(float)
                n_exposed_case = int(np.sum(bsv_col[y == 1] > 0))
                n_exposed_ctrl = int(np.sum(bsv_col[y == 0] > 0))

                # X_full: intercept, bsv, sex, log1p_base, log1p_gene, ancestry, platform
                X_full = np.column_stack([np.ones(int(keep.sum())), bsv_col, sex_vals,
                                          log1p_base, log1p_gene, anc_vals, plat_vals])

                obs_beta = obs_se = obs_or = obs_p_twosided = np.nan

                try:
                    model = GEE(
                        endog=y,
                        exog=X_full,
                        groups=family_ids,
                        family=Binomial(),
                        cov_struct=Independence()
                    )
                    result = model.fit(maxiter=200)

                    if result.converged:
                        obs_beta = float(result.params[1])
                        obs_se = float(result.bse[1])
                        obs_or = float(np.exp(obs_beta))
                        obs_p_twosided = float(result.pvalues[1])
                    else:
                        log(f"  WARNING: GEE Independence did not converge for {bclass}/{svt} Pattern {pattern}")
                except Exception as e:
                    log(f"  WARNING: GEE Independence failed for {bclass}/{svt} Pattern {pattern}: {str(e)[:100]}")

                # One-sided / two-sided P based on Discovery significance
                if is_negative_ctrl:
                    # Negative control: always two-sided, P_onesided = NaN
                    obs_p_onesided = np.nan
                elif np.isfinite(obs_p_twosided) and np.isfinite(obs_beta):
                    if bclass in DISCOVERY_SIG_CLASSES:
                        # Discovery significant -> one-sided (expected direction: beta > 0)
                        if obs_beta > 0:
                            obs_p_onesided = obs_p_twosided / 2.0
                        else:
                            obs_p_onesided = 1.0 - obs_p_twosided / 2.0
                    else:
                        # Discovery non-significant -> two-sided
                        obs_p_onesided = obs_p_twosided
                else:
                    obs_p_onesided = np.nan

                ci_lo = float(np.exp(obs_beta - 1.96 * obs_se)) if np.isfinite(obs_se) else np.nan
                ci_hi = float(np.exp(obs_beta + 1.96 * obs_se)) if np.isfinite(obs_se) else np.nan

                # Replication significance
                if is_primary or is_secondary_c:
                    if bclass in DISCOVERY_SIG_CLASSES:
                        repl_sig = (np.isfinite(obs_p_onesided) and obs_p_onesided < 0.05 and obs_beta > 0)
                    else:
                        repl_sig = (np.isfinite(obs_p_onesided) and obs_p_onesided < 0.05)
                else:
                    repl_sig = False

                # Analysis type and log display
                if is_primary:
                    analysis_type = "PRIMARY_replication"
                    tag = "**" if repl_sig else "  "
                elif is_secondary_c:
                    analysis_type = "SECONDARY_pattern_c"
                    tag = "**" if repl_sig else "  "
                else:
                    analysis_type = "negative_control"
                    tag = "NC" if np.isfinite(obs_p_twosided) and obs_p_twosided < 0.05 else "  "

                log(f"  {tag} {bclass:20s}/{svt}: "
                    f"OR={obs_or:.3f} beta={obs_beta:+.4f} "
                    f"P_2s={obs_p_twosided:.2e} P_1s={obs_p_onesided:.2e} "
                    f"n_exp_case={n_exposed_case} n_exp_ctrl={n_exposed_ctrl} "
                    f"n_clusters={n_clusters} [{analysis_type}]")

                result_row = {
                    "Pattern": pattern,
                    "Boundary_class": bclass, "SV_type": svt,
                    "N_case": n_case, "N_control": n_ctrl, "N_complete": n_case + n_ctrl,
                    "N_exposed_case": n_exposed_case, "N_exposed_control": n_exposed_ctrl,
                    "N_clusters": n_clusters,
                    "Beta": obs_beta, "SE": obs_se, "OR": obs_or,
                    "CI_lo": ci_lo, "CI_hi": ci_hi,
                    "P_twosided": obs_p_twosided,
                    "P_onesided": obs_p_onesided,
                    "direction_match": bool(obs_beta > 0) if np.isfinite(obs_beta) else False,
                    "analysis_type": analysis_type,
                    "is_primary": is_primary,
                    "is_secondary_c": is_secondary_c,
                    "replicated": repl_sig,
                    "formal_test": bool(is_primary and bclass in DISCOVERY_SIG_CLASSES),
                    "exploratory_test": bool(is_primary and bclass in DISCOVERY_NONSIG_CLASSES),
                    "formal_replicated": bool(is_primary and bclass in DISCOVERY_SIG_CLASSES and repl_sig),
                    "exploratory_signal": bool(is_primary and bclass in DISCOVERY_NONSIG_CLASSES and repl_sig),
                }
                all_results.append(result_row)

    return all_results

# ============================================================================
# PHASE 7: Output (no FDR correction for replication)
# ============================================================================
def run_phase7_output(all_results, n_case, n_ctrl, n_clusters, gd_nahr_count):
    """Phase 7: 結果出力 (両パターン, FDR補正なし)"""
    log(f"\n{'=' * 70}")
    log("PHASE 7: Output (both patterns, no FDR)")
    log("=" * 70)

    outdir = Path(_OUTDIR)
    outdir.mkdir(parents=True, exist_ok=True)
    results_df = pd.DataFrame(all_results)

    # --- 全結果を1つのファイルに出力 ---
    out_df = pd.DataFrame({
        "cohort": "MSSNG_Replication",
        "comparison": "ASD_vs_unaffected_sibling",
        "pattern": results_df["Pattern"].values,
        "sv_type": results_df["SV_type"].values,
        "boundary_class": results_df["Boundary_class"].values,
        "model": "GEE_Independence_Bprime",
        "exposure_type": "bin_count",
        "analysis_type": results_df["analysis_type"].values,
        "n_case": results_df["N_case"].values,
        "n_control": results_df["N_control"].values,
        "n_complete": results_df["N_complete"].values,
        "n_exposed_case": results_df["N_exposed_case"].values,
        "n_exposed_control": results_df["N_exposed_control"].values,
        "n_clusters": results_df["N_clusters"].values,
        "beta": results_df["Beta"].values,
        "SE": results_df["SE"].values,
        "OR": results_df["OR"].values,
        "CI_lower": results_df["CI_lo"].values,
        "CI_upper": results_df["CI_hi"].values,
        "P_twosided": results_df["P_twosided"].values,
        "P_onesided": results_df["P_onesided"].values,
        "direction_match": results_df["direction_match"].values,
        "discovery_expected_direction": "positive",
        "replicated": results_df["replicated"].values,
        "formal_test": results_df["formal_test"].values,
        "exploratory_test": results_df["exploratory_test"].values,
        "formal_replicated": results_df["formal_replicated"].values,
        "exploratory_signal": results_df["exploratory_signal"].values,
    })
    all_out_path = str(outdir / "tad_replication_mssng_v18.tsv")
    out_df.to_csv(all_out_path, sep="\t", index=False, float_format="%.6g")
    log(f"  All results output: {all_out_path}")

    # --- Summary: PRIMARY replication (Pattern A) ---
    primary_df = out_df[out_df["analysis_type"] == "PRIMARY_replication"].copy()
    log("\n=== PRIMARY REPLICATION (Pattern A, DEL, bin_count, GEE Independence) ===")

    # Formal replication: sig8
    formal_df = primary_df[primary_df["formal_test"] == True].copy()
    n_formal_repl = int(formal_df["formal_replicated"].sum())
    log(f"  Formal replication (sig8): {n_formal_repl} / {len(formal_df)}")
    for _, r in formal_df.sort_values("P_onesided").iterrows():
        sig = "**" if r["formal_replicated"] else "  "
        dir_s = "+" if r["direction_match"] else "-"
        log(f"  {sig} {r['boundary_class']:20s} OR={r['OR']:.3f} [{r['CI_lower']:.3f}-{r['CI_upper']:.3f}] "
            f"P_1s={r['P_onesided']:.2e} P_2s={r['P_twosided']:.2e} dir={dir_s} "
            f"n_exp_case={int(r['n_exposed_case'])} n_exp_ctrl={int(r['n_exposed_control'])} "
            f"n_clusters={int(r['n_clusters'])}")

    # Exploratory: nonsig2
    explor_df = primary_df[primary_df["exploratory_test"] == True].copy()
    n_explor_sig = int(explor_df["exploratory_signal"].sum())
    log(f"  Exploratory (nonsig2): {n_explor_sig} / {len(explor_df)}")
    for _, r in explor_df.sort_values("P_twosided").iterrows():
        sig = "**" if r["exploratory_signal"] else "  "
        dir_s = "+" if r["direction_match"] else "-"
        log(f"  {sig} {r['boundary_class']:20s} OR={r['OR']:.3f} [{r['CI_lower']:.3f}-{r['CI_upper']:.3f}] "
            f"P_2s={r['P_twosided']:.2e} dir={dir_s} "
            f"n_exp_case={int(r['n_exposed_case'])} n_exp_ctrl={int(r['n_exposed_control'])} "
            f"n_clusters={int(r['n_clusters'])}")

    # --- Summary: SECONDARY Pattern C ---
    secondary_df = out_df[out_df["analysis_type"] == "SECONDARY_pattern_c"].copy()
    log("\n=== SECONDARY REPLICATION (Pattern C, DEL, abs_len <= 1MB, bin_count, GEE Independence) ===")

    sec_sig8 = secondary_df[secondary_df["boundary_class"].isin(DISCOVERY_SIG_CLASSES)].copy()
    n_sec_sig8_repl = int(sec_sig8["replicated"].sum())
    log(f"  Formal secondary (sig8): {n_sec_sig8_repl} / {len(sec_sig8)}")
    for _, r in sec_sig8.sort_values("P_onesided").iterrows():
        sig = "**" if r["replicated"] else "  "
        dir_s = "+" if r["direction_match"] else "-"
        log(f"  {sig} {r['boundary_class']:20s} OR={r['OR']:.3f} [{r['CI_lower']:.3f}-{r['CI_upper']:.3f}] "
            f"P_1s={r['P_onesided']:.2e} P_2s={r['P_twosided']:.2e} dir={dir_s} "
            f"n_exp_case={int(r['n_exposed_case'])} n_exp_ctrl={int(r['n_exposed_control'])} "
            f"n_clusters={int(r['n_clusters'])}")

    sec_nonsig2 = secondary_df[~secondary_df["boundary_class"].isin(DISCOVERY_SIG_CLASSES)].copy()
    n_sec_nonsig2_repl = int(sec_nonsig2["replicated"].sum())
    log(f"  Exploratory secondary (nonsig2): {n_sec_nonsig2_repl} / {len(sec_nonsig2)}")
    for _, r in sec_nonsig2.sort_values("P_twosided").iterrows():
        sig = "**" if r["replicated"] else "  "
        dir_s = "+" if r["direction_match"] else "-"
        log(f"  {sig} {r['boundary_class']:20s} OR={r['OR']:.3f} [{r['CI_lower']:.3f}-{r['CI_upper']:.3f}] "
            f"P_2s={r['P_twosided']:.2e} dir={dir_s} "
            f"n_exp_case={int(r['n_exposed_case'])} n_exp_ctrl={int(r['n_exposed_control'])} "
            f"n_clusters={int(r['n_clusters'])}")

    # --- Summary: Negative control ---
    negctrl_df = out_df[out_df["analysis_type"] == "negative_control"].copy()
    log("\n=== NEGATIVE CONTROL / SPECIFICITY (both patterns) ===")
    for pattern in sorted(negctrl_df["pattern"].unique()):
        pattern_data = negctrl_df[negctrl_df["pattern"] == pattern]
        log(f"\nPattern {pattern}:")
        for _, r in pattern_data.iterrows():
            tag = "NC" if np.isfinite(r["P_twosided"]) and r["P_twosided"] < 0.05 else "  "
            ctrl_type = "DUP_specificity" if r["sv_type"] == "DUP" else "static_specificity"
            log(f"  {tag} {r['boundary_class']:20s}/{r['sv_type']}: "
                f"OR={r['OR']:.3f} P_2s={r['P_twosided']:.2e} [{ctrl_type}]")

    # --- Config ---
    with open(str(outdir / "run_config.json"), "w") as f:
        json.dump({
            "script": "tad_replication_mssng_v18.py",
            "version": "v18",
            "changes_from_v17": [
                "tad04212026 pipeline: _BASEDIR relocated to 09_mssng_sample_burden",
                "_OUTDIR set to {_BASEDIR}/output_v18 (or output_v18_{MIN_SV_LEN}bp for sensitivity)",
                "Output TSV renamed: tad_replication_mssng_v17.tsv -> tad_replication_mssng_v18.tsv",
                "Preflight file renamed: preflight_exposed_counts_v17.tsv -> preflight_exposed_counts_v18.tsv",
                "All Phase 0-7 analysis logic (GEE Independence, covariates, Pattern A/C, negctrl) UNCHANGED from v17",
            ],
            "changes_from_v14_v15": [
                "One-per-family selection removed: all eligible ASD + sibling retained",
                "GEE handles family correlation via FAMILYID clustering",
                "15q11.2 NAHR logic fix: added and 15q13 not in gid (matching arrayCGH)",
                "DRY_RUN mode added (env DRY_RUN=1 -> Phase 4 only, output exposed counts)",
                "(v15 added Exchangeable primary + Independence sensitivity; v16 reversed; v17 removed Exchangeable)",
            ],
            "changes_from_v15": [
                "GEE Independence as primary (v15 Exchangeable failed to converge for all classes)",
                "fit(maxiter=200) explicit",
            ],
            "changes_from_v16": [
                "Exchangeable sensitivity removed from runtime (non-convergent in v15/v16)",
                "sens_exch_* columns removed from TSV",
                "formal_test / exploratory_test / formal_replicated / exploratory_signal columns added",
                "formal_replication / exploratory_replication columns removed",
                "Negative control P_onesided = NaN (always two-sided)",
                "model label changed to GEE_Independence_Bprime",
                "Summary log split: formal (sig8) / exploratory (nonsig2)",
            ],
            "dataset": "MSSNG",
            "comparison": "ASD vs unaffected_sibling",
            "primary_boundary_classes": PRIMARY_BOUNDARY_CLASSES,
            "discovery_sig_classes": sorted(list(DISCOVERY_SIG_CLASSES)),
            "discovery_nonsig_classes": sorted(list(DISCOVERY_NONSIG_CLASSES)),
            "discovery_nonsig_classes": ["hpc_exc_ca", "pfc_astro"],
            "primary_replication": {
                "pattern": "A",
                "sv_type": "DEL",
                "boundary_classes": PRIMARY_BOUNDARY_CLASSES,
                "exposure_type": "bin_count",
                "model": "GEE_Independence_Bprime",
                "p_value_type": "one-sided (beta > 0) for SIG classes, two-sided for non-SIG",
                "significance_threshold": 0.05,
                "fdr_correction": "none",
            },
            "secondary_pattern_c": {
                "pattern": "C",
                "sv_type": "DEL",
                "abs_len_max": MAX_SV_LEN_PATTERN_C,
                "boundary_classes": PRIMARY_BOUNDARY_CLASSES,
                "exposure_type": "bin_count",
                "model": "GEE_Independence_Bprime",
                "p_value_type": "one-sided (beta > 0) for SIG classes, two-sided for non-SIG",
                "significance_threshold": 0.05,
                "fdr_correction": "none",
            },
            "negative_control": {
                "dup_specificity": {
                    "sv_type": "DUP",
                    "boundary_classes": PRIMARY_BOUNDARY_CLASSES,
                    "purpose": "Confirm DEL-specificity of TAD disruption effect",
                },
                "static_specificity": {
                    "sv_type": "DEL",
                    "boundary_classes": ["static_all"],
                    "purpose": "Confirm developmental TAD boundary specificity",
                },
                "p_value_type": "two-sided",
            },
            "discovery_reference": {
                "cohort": "GRIFIN_Discovery (srWGS)",
                "significant_results": "8 boundary classes FDR < 0.05 (ASD_vs_Healthy, Pattern A, DEL, bin_count, beta > 0)",
            },
            "inference": "GEE logistic regression (family-clustered, Independence working correlation, full-family)",
            "gee_family": "Binomial",
            "gee_cov_struct_primary": "Independence",
            "gee_cov_struct_sensitivity": "none",
            "exchangeable_removed_from_runtime": True,
            "reason_for_removal": "Exchangeable failed to converge across all tested classes in v15/v16",
            "formal_test_classes": sorted(list(DISCOVERY_SIG_CLASSES)),
            "exploratory_classes": ["hpc_exc_ca", "pfc_astro"],
            "family_variable": "FAMILYID",
            "family_selection": "all eligible (no one-per-family filtering)",
            "covariates": "Sex + Ancestry(5 dummies, ref=EUR) + Platform(5 dummies, ref=HiSeqX) + log1p(total_bases) + log1p(total_gene) [Pattern A: all non-NAHR CNVs, Pattern C: non-NAHR + <=1MB CNVs]",
            "cnv_source": _CNV_MERGED_FILE,
            "min_sv_len": MIN_SV_LEN,
            "max_sv_len_pattern_c": MAX_SV_LEN_PATTERN_C,
            "bin_overlap_thr": BIN_OVERLAP_THR,
            "n_samples": {"ASD": n_case, "sibling": n_ctrl, "total": n_case + n_ctrl},
            "n_families": n_clusters,
            "n_nahr_loci": gd_nahr_count,
            "l2_diffbound_dir": _L2_DIFFBOUND_DIR,
            "raw_h5ad": _RAW_H5AD,
            "timestamp": datetime.now().isoformat(),
        }, f, indent=2)

    return all_out_path


# ============================================================================
# MAIN
# ============================================================================
def main():
    total_start = time.time()
    dry_run = os.environ.get("DRY_RUN", "0") == "1"

    log("=" * 70)
    log("tad_replication_mssng_v18.py")
    log("  MSSNG replication pipeline (ASD vs unaffected_sibling)")
    log("  *** Full-family version (no one-per-family filtering) ***")
    log("  Design: Discovery-Replication framework")
    log("  PRIMARY (Pattern A): DEL x 10 L2 lineage classes")
    log("  SECONDARY (Pattern C): DEL x 10 L2 lineage classes (abs_len <= 1MB)")
    log("  Negative control: DUP x 10 classes + DEL x static_all")
    log("  Inference: GEE Independence (family-clustered, no Exchangeable)")
    log(f"  MIN_SV_LEN={MIN_SV_LEN}")
    if dry_run:
        log("  *** DRY_RUN MODE: will stop after Phase 4 (exposed counts only) ***")
    log("=" * 70)

    # Phase 0: Boundary master (conditional - needed for static_all via raw h5ad)
    master_path = build_boundary_master(_BOUNDARY_MASTER)

    # Phase 1: Load NAHR GD loci
    gd_nahr = run_phase1_load_nahr()

    # Phase 2: Extract + filter CNVs (full-family, no one-per-family)
    (sv, sample_status_dict, fid_map,
     get_sex, get_ancestry, get_platform_norm) = run_phase2_load_cnv(gd_nahr)

    # Phase 3: Load boundary bins from L2 diffbound BED files
    bin_index = run_phase3_load_bins()

    # Phase 4: Compute disrupted bin counts (Pattern A + Pattern C)
    # v19 (patch v3): pass fid_map + getters for constraint-enrichment dual dump
    (all_samples, sample_to_idx, bin_counts_a, bin_counts_c,
     sample_total_base, sample_total_gene,
     sample_total_base_c, sample_total_gene_c) = run_phase4_bin_counts(
        sv, sample_status_dict, bin_index,
        fid_map=fid_map,
        get_sex=get_sex,
        get_platform_norm=get_platform_norm,
        get_ancestry=get_ancestry)

    # --- DRY_RUN: stop here, output exposed counts for preflight check ---
    if dry_run:
        log("\n" + "=" * 70)
        log("DRY_RUN: Phase 4 complete. Outputting exposed counts.")
        log("=" * 70)

        asd_samples = [s for s in all_samples if sample_status_dict.get(s) == CASE_LABEL]
        sib_samples = [s for s in all_samples if sample_status_dict.get(s) == CTRL_LABEL]
        log(f"  Total eligible ASD: {len(asd_samples)}")
        log(f"  Total eligible sibling: {len(sib_samples)}")

        fam_ids_set = set()
        for s in all_samples:
            fid = fid_map.get(s)
            if fid:
                fam_ids_set.add(fid)
        log(f"  Total unique families: {len(fam_ids_set)}")

        DRY_BCLASSES = PRIMARY_BOUNDARY_CLASSES + ["static_all"]
        outdir = Path(_OUTDIR)
        outdir.mkdir(parents=True, exist_ok=True)
        preflight_path = str(outdir / "preflight_exposed_counts_v18.tsv")
        rows = []
        for pattern, bin_counts in [("A", bin_counts_a), ("C", bin_counts_c)]:
            for svt in SV_TYPES:
                for bclass in DRY_BCLASSES:
                    is_relevant = (svt == "DEL" and bclass in PRIMARY_BOUNDARY_CLASSES) or \
                                  (svt == "DUP" and bclass in PRIMARY_BOUNDARY_CLASSES) or \
                                  (svt == "DEL" and bclass == "static_all")
                    if not is_relevant:
                        continue
                    bc = bin_counts[bclass][svt]
                    asd_idx = [sample_to_idx[s] for s in asd_samples if s in sample_to_idx]
                    sib_idx = [sample_to_idx[s] for s in sib_samples if s in sample_to_idx]
                    n_exp_asd = int(np.sum(bc[asd_idx] > 0)) if asd_idx else 0
                    n_exp_sib = int(np.sum(bc[sib_idx] > 0)) if sib_idx else 0
                    rows.append({
                        "pattern": pattern, "sv_type": svt, "boundary_class": bclass,
                        "n_asd": len(asd_samples), "n_sib": len(sib_samples),
                        "n_exposed_asd": n_exp_asd, "n_exposed_sib": n_exp_sib,
                    })
        preflight_df = pd.DataFrame(rows)
        preflight_df.to_csv(preflight_path, sep="\t", index=False)
        log(f"  Preflight exposed counts: {preflight_path}")

        log("\n  PRIMARY (Pattern A, DEL) exposed counts:")
        for _, r in preflight_df[(preflight_df["pattern"] == "A") &
                                  (preflight_df["sv_type"] == "DEL") &
                                  (preflight_df["boundary_class"].isin(PRIMARY_BOUNDARY_CLASSES))].iterrows():
            sig_mark = "*" if r["boundary_class"] in DISCOVERY_SIG_CLASSES else " "
            log(f"    {sig_mark} {r['boundary_class']:20s}: ASD={r['n_exposed_asd']}/{r['n_asd']} "
                f"sib={r['n_exposed_sib']}/{r['n_sib']}")

        total_elapsed = time.time() - total_start
        log(f"\nDRY_RUN complete. Elapsed: {total_elapsed:.1f}s ({total_elapsed / 60:.1f}min)")
        return

    # Phase 5: Build sample dataframe
    (sample_df, y, keep, valid_si, family_ids,
     n_case, n_ctrl, n_clusters) = run_phase5_sample_df(
        all_samples, sample_status_dict, fid_map,
        get_sex, get_ancestry, get_platform_norm)

    # Phase 6: GEE regression (PRIMARY + SECONDARY Pattern C + negative control)
    all_results = run_phase6_regression(
        sample_df, y, keep, valid_si, family_ids,
        bin_counts_a, bin_counts_c,
        sample_total_base, sample_total_gene,
        sample_total_base_c, sample_total_gene_c,
        n_case, n_ctrl, n_clusters)

    # Phase 7: Output (both patterns, no FDR)
    out_path = run_phase7_output(all_results, n_case, n_ctrl, n_clusters, len(gd_nahr))

    total_elapsed = time.time() - total_start
    log(f"\nPipeline complete. Total elapsed: {total_elapsed:.1f}s ({total_elapsed / 60:.1f}min)")
    log(f"Results output: {out_path}")

if __name__ == "__main__":
    main()
