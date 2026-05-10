#!/usr/bin/env python3
# ============================================================================
# ファイル名: tad_replication_arraycgh_v22.py
#
# 処理概要:
#   アレイCGHデータ（ASD, SZ, CONT）を用いた
#   GRIFIN TAD boundary disruption解析の再現検証 + ASD-SZ heterogeneity検証パイプライン
#
#   - Phase 0 : Validation of L2 diffbound BED files
#   - Phase 1 : CNVデータ読み込み + common CNV除外(platform+DEL/DUP別) + liftOver (hg18→hg38)
#   - Phase 2 : サンプルQC + discovery overlap除外（overlap_srWGS列ベース, ASD + SCZ 両方）
#   - Phase 3 : CNVフィルタ + segdup/exclusion BED + NAHR GD flagging + 遺伝子カウント
#   - Phase 4 : L2 diffbound BED files から 10 lineage class + static_all のbin setsを構築
#   - Phase 5 : unique disrupted bin count per sample (Pattern A + Pattern C)
#   - Phase 5B: per-sample disrupted bin IDs + sample covariatesをTSV出力
#   - Phase 6 : ASD vs CONT ロジスティック回帰 (PRIMARY replication + negative control)
#   - Phase 6B: SZ vs CONT ロジスティック回帰 (external SZ burden assessment)
#   - Phase 6C: ASD-SZ heterogeneity test (MNLogit + Stouffer + sign test)
#   - Phase 7 : 結果出力 (全比較, FDR補正なし)
#
# v18からの変更点 (v19):
#   1. Phase 2: discovery overlap除外方式を変更
#      - 旧方式: GRIFIN_srWGS_SampleInfo_11242025.txt の SampleID を CNV解析ID とマッチング
#      - 新方式: sample_data_01102023.xlsx の overlap_srWGS 列が "Yes" のサンプルを除外
#      - _SRWGS_SAMPLE_FILE 設定を削除
#   2. 出力ファイル名を _v19 に更新
#
# *** Phase 0-1, 3-7 の解析ロジックは v18 から変更なし ***
#
# v21 からの変更点 (v22, 2026-04-21 tad04212026 パイプライン移植):
#   1. _BASEDIR を tad04212026/08_arraycgh_sample_burden に変更。
#      arrayCGH 入力 (xlsx, chain, commonCNV) の _SCRIPT_DIR は
#      noncoding_tad_arraycgh_03302026 のままで参照（read-only）。
#   2. _OUTDIR を {_BASEDIR}/output_v22 に変更。
#   3. 出力 TSV ファイル名サフィックスを _v21 → _v22 に更新。
#   4. script metadata (version key)、log banner、SBATCH 例を v22 に更新。
#   5. 解析ロジック (Phase 0-7) は v21 から一切変更なし。
#
# SBATCH example:
#   sbatch --wrap="python tad_replication_arraycgh_v22.py" \
#     -p ncbn-cpu --account=ncbn-cpu --cpus-per-task=8 --mem=64G \
#     --job-name=tad_acgh_v22 --output=logs/tad_acgh_v22_%j.out \
#     --error=logs/tad_acgh_v22_%j.err --time=24:00:00
# ============================================================================

import os
import sys
import re
import time
import gzip
import json
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from bisect import bisect_left, bisect_right
from collections import defaultdict, Counter
from typing import List, Tuple, Dict, Optional

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2, norm
from scipy import stats as scipy_stats
import warnings
warnings.filterwarnings("ignore")

# ============================================================================
# CONFIGURATION
# ============================================================================
_BASEDIR = "/lustre12/home/kushima-pg/tad04212026/08_arraycgh_sample_burden"
_HEFFEL_DIR = "/lustre12/home/kushima-pg/heffel_deep_analysis_03242026"
_SCRIPT_DIR = "/lustre12/home/kushima-pg/noncoding_tad_arraycgh_03302026"

_CNV_FILE = f"{_SCRIPT_DIR}/cnv.annotation_gencode.v35_01042023_11h39m05.xlsx"
_SAMPLE_FILE = f"{_SCRIPT_DIR}/sample_data_01102023.xlsx"
_CHAIN_FILE = "/lustre12/home/kushima-pg/arrayCGH/hg18ToHg38.over.chain.gz"
_LIFTOVER_BIN = "/home/kushima-pg/bin/liftOver"
_COMMON_CNV_AGILENT = f"{_SCRIPT_DIR}/commonCNV/agilent_0.1percent_commonCNV.txt"
_COMMON_CNV_NIMBLEGEN = f"{_SCRIPT_DIR}/commonCNV/nimblegen_0.1percent_commonCNV.txt"

# F2 fix: unused dead code — commented out (MSSNG v17 still uses this variable)
# _BOUNDARY_MASTER = f"{_HEFFEL_DIR}/heffel_master_outputs/heffel_boundary_master_v7.tsv.gz"
_CURATED_GD_FILE = "/lustre12/home/kushima-pg/cnv_01012026/curated_genomic_disorder_cnv_loci_v3.txt"
_GTF_FILE = "/lustre12/home/kushima-pg/annotationInfo/gencode.v46.annotation.gtf.gz"
# v19: _SRWGS_SAMPLE_FILE は廃止。overlap_srWGS列（sample_data_01102023.xlsx内）で除外する

_RAW_H5AD = f"{_HEFFEL_DIR}/heffel_inspect/domain_boundaries/domain_boundaries/BrainDev_raw.boundary.h5ad"
_IMPUTE_H5AD = f"{_HEFFEL_DIR}/heffel_inspect/domain_boundaries/domain_boundaries/BrainDev_impute.boundary.h5ad"
_DIFF_DIR = f"{_HEFFEL_DIR}/heffel_inspect/diff_boundaries/L2_diff_domain_boundaries"
_L2_DIFFBOUND_DIR = f"{_HEFFEL_DIR}/heffel_inspect/diff_boundaries/L2_diff_domain_boundaries"

_OUTDIR = f"{_BASEDIR}/output_v22"

DISCOVERY_EXPECTED_DIRECTION = "positive"

CASE_LABEL = "ASD"
CTRL_LABEL = "CONT"
SCZ_LABEL = "SCZ"       # v16 新規
SV_TYPES = ["DEL", "DUP"]
STATE_MAP = {1: "DEL", 2: "DUP"}

PRIMARY_BOUNDARY_CLASSES = [
    "hpc_exc_ca", "hpc_exc_dg", "hpc_exc_ent",
    "hpc_inh_cge", "hpc_inh_mge",
    "pfc_astro", "pfc_exc_dl", "pfc_exc_ul",
    "pfc_inh_cge", "pfc_inh_mge",
]

DISCOVERY_SIG_CLASSES = {
    "hpc_exc_dg", "hpc_exc_ent", "hpc_inh_cge", "hpc_inh_mge",
    "pfc_exc_dl", "pfc_exc_ul", "pfc_inh_cge", "pfc_inh_mge",
}

DISCOVERY_NONSIG_CLASSES = {"hpc_exc_ca", "pfc_astro"}

NAHR_RO_THR = 0.50
COMMON_CNV_OVERLAP_THR = 0.30       # v15互換: one-direction overlap 30%
MIN_SV_LEN = 10000
MAX_SV_LEN_PATTERN_C = 1000000
BIN_OVERLAP_THR = 0.10

# v18追加: segdup / exclusion BED (WGS v18と同一)
_SEGDUP_BED = "/lustre12/home/kushima-pg/annotationInfo/segdup_hg38_sorted_merged.bed"
_EXCLUSION_BED = "/lustre12/home/kushima-pg/resource/hg38_cnv_exclusion_regions.bed"
SEGDUP_MAX_PCT = 50.0           # segdup coverage < 50% (WGS v18と同一)
EXCLUSION_OVERLAP_THR = 0.50    # exclusion BED coverage < 50% (WGS v18と同一)


# ============================================================================
# UTILITIES (v15と同一 + v18追加)
# ============================================================================
def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def norm_chrom(c):
    s = str(c).strip()
    if s.startswith("chr"):
        return s
    if s.isdigit() or s in ("X", "Y", "M", "MT"):
        return f"chr{s}"
    return s


def calc_reciprocal_overlap(s1, e1, s2, e2):
    ov = max(0, min(e1, e2) - max(s1, s2))
    l1, l2 = e1 - s1, e2 - s2
    if l1 <= 0 or l2 <= 0:
        return 0.0
    return min(ov / l1, ov / l2)


def calc_target_coverage(s1, e1, s2, e2):
    ov = max(0, min(e1, e2) - max(s1, s2))
    tlen = e2 - s2
    if tlen <= 0:
        return 0.0
    return ov / tlen


def assert_unique_keys(df: pd.DataFrame, key_cols: List[str], label: str, n_show: int = 10) -> None:
    dup = df[df.duplicated(subset=key_cols, keep=False)]
    if len(dup) > 0:
        log(f"FATAL: {label} has {len(dup)} duplicate rows on {key_cols}")
        log(dup.head(n_show).to_string())
        raise ValueError(f"{label}: {len(dup)} duplicate rows on {key_cols}. "
                         f"Pipeline halted to prevent data contamination.")


def parse_region_str(region_str):
    region_str = str(region_str).replace(",", "")
    m = re.match(r"(chr[\dXYMT]+):(\d+)-(\d+)", region_str)
    if not m:
        return None, None, None
    return m.group(1), int(m.group(2)), int(m.group(3))


def one_direction_overlap(cnv_start, cnv_end, common_start, common_end):
    """Fraction of CNV covered by common CNV region (v15互換)."""
    ov = max(0, min(cnv_end, common_end) - max(cnv_start, common_start))
    cnv_len = cnv_end - cnv_start
    return ov / cnv_len if cnv_len > 0 else 0.0


# --- v18追加: segdup / exclusion BED ユーティリティ (WGS v18 から移植) ---

def load_bed3(path: Path) -> Dict[str, list]:
    """BED3ファイルを読み込み、chrom -> [(start, end), ...] の辞書を返す。"""
    out = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = norm_chrom(parts[0])
            try:
                s = int(parts[1])
                e = int(parts[2])
            except ValueError:
                continue
            if e > s:
                out[chrom].append((s, e))
    return dict(out)


def _merge_intervals_sorted(starts, ends):
    if len(starts) == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    ms, me = [starts[0]], [ends[0]]
    for s, e in zip(starts[1:], ends[1:]):
        if s <= me[-1]:
            me[-1] = max(me[-1], e)
        else:
            ms.append(s)
            me.append(e)
    return np.array(ms, dtype=np.int64), np.array(me, dtype=np.int64)


def build_merged_interval_index(intervals_by_chrom: Dict) -> Dict:
    """chrom -> {start: np.array, end: np.array} のマージ済みインデックスを構築。"""
    index = {}
    for chrom, lst in intervals_by_chrom.items():
        arr = np.asarray(lst, dtype=np.int64)
        if arr.size == 0:
            continue
        order = np.argsort(arr[:, 0])
        ms, me = _merge_intervals_sorted(arr[order, 0], arr[order, 1])
        index[chrom] = {"start": ms, "end": me}
    return index


def compute_coverage_from_index(chrom: str, sv_start: int, sv_end: int, merged_index: Dict) -> float:
    """CNVに対するBED領域のカバレッジ (%) を返す。"""
    sv_len = sv_end - sv_start
    if sv_len <= 0 or chrom not in merged_index:
        return 0.0
    ms = merged_index[chrom]["start"]
    me = merged_index[chrom]["end"]
    if ms.size == 0:
        return 0.0
    idx_start = bisect_right(ms, sv_end) - 1
    idx_end = bisect_left(ms, sv_start)
    cov = 0
    for i in range(max(0, idx_end - 1), min(len(ms), idx_start + 2)):
        ov = max(0, min(sv_end, int(me[i])) - max(sv_start, int(ms[i])))
        cov += ov
    return min(100.0, 100.0 * cov / sv_len)


# ============================================================================
# PHASE 0: Validation (v15互換: *_diffbound.bed.gz ファイル確認)
# ============================================================================
def run_phase0_validate():
    """Phase 0: L2 diffbound BED ファイルと raw h5ad が存在することを確認"""
    log("=" * 70)
    log("PHASE 0: Validation of L2 diffbound BED files and raw h5ad")
    log("=" * 70)

    diff_dir = Path(_L2_DIFFBOUND_DIR)
    if not diff_dir.exists():
        raise FileNotFoundError(f"L2 diffbound directory not found: {_L2_DIFFBOUND_DIR}")
    log(f"  L2 diffbound directory: {_L2_DIFFBOUND_DIR}")

    diff_files = sorted(diff_dir.glob("*_diffbound.bed.gz"))
    if len(diff_files) == 0:
        raise FileNotFoundError(f"No *_diffbound.bed.gz files found in: {diff_dir}")
    for f in diff_files:
        log(f"  Found: {f.name}")

    if not os.path.exists(_RAW_H5AD):
        raise FileNotFoundError(f"Raw h5ad not found: {_RAW_H5AD}")
    log(f"  Raw h5ad: {_RAW_H5AD}")

    log("Phase 0 validation complete.")


# ============================================================================
# PHASE 1: Load CNV + common CNV filtering (platform+DEL/DUP別) + liftOver
# *** v15から完全移植 ***
# ============================================================================
def load_common_cnv_file(filepath, platform_name):
    """Platform別・DEL/DUP別にcommon CNVファイルを解析 (v15互換)"""
    result = {"DEL": [], "DUP": []}
    n_del, n_dup = 0, 0
    if not os.path.exists(filepath):
        log(f"  WARNING: Common CNV file not found: {filepath}")
        return result, 0, 0
    with open(filepath, "r") as f:
        for line_num, line in enumerate(f):
            if line_num < 3:
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            region = parts[0]
            event = parts[3]
            chrom, start, end = parse_region_str(region)
            if chrom is None:
                continue
            chrom_norm = norm_chrom(chrom)
            if "Loss" in event:
                result["DEL"].append((chrom_norm, start, end))
                n_del += 1
            elif "Gain" in event:
                result["DUP"].append((chrom_norm, start, end))
                n_dup += 1
    log(f"  {platform_name} common CNVs: DEL={n_del}, DUP={n_dup}")
    return result, n_del, n_dup


def run_phase1_load_cnv():
    """Phase 1: CNVデータ読み込み + common CNV除外 + liftOver (v15互換)"""
    log("=" * 70)
    log("PHASE 1: Loading CNV data + common CNV filtering + liftOver hg18→hg38")
    log("=" * 70)

    import openpyxl  # availability check

    cnv_raw = pd.read_excel(_CNV_FILE)
    log(f"  Raw CNVs: {len(cnv_raw)}")

    cnv_raw["CNV_type"] = cnv_raw["state"].map(STATE_MAP)
    cnv_raw = cnv_raw[cnv_raw["CNV_type"].notna()].copy()
    log(f"  After state mapping (DEL/DUP): {len(cnv_raw)}")

    # --- Common CNV filtering (platform + DEL/DUP別, v15互換) ---
    log("--- Filtering common CNVs (0.1%% frequency, platform+type specific) ---")

    sample_platform_df = pd.read_excel(_SAMPLE_FILE, usecols=["CNV解析ID", "Data Type"])
    platform_map = dict(zip(sample_platform_df["CNV解析ID"], sample_platform_df["Data Type"]))
    log(f"  Platform map loaded: {len(platform_map)} samples")

    common_cnvs = {"Agilent FE": {"DEL": [], "DUP": []}, "Nimblegen Normalized": {"DEL": [], "DUP": []}}
    agilent_data, _, _ = load_common_cnv_file(_COMMON_CNV_AGILENT, "Agilent")
    nimblegen_data, _, _ = load_common_cnv_file(_COMMON_CNV_NIMBLEGEN, "Nimblegen")
    common_cnvs["Agilent FE"] = agilent_data
    common_cnvs["Nimblegen Normalized"] = nimblegen_data

    cnv_to_remove = set()
    removed_by_platform = {"Agilent FE": 0, "Nimblegen Normalized": 0}
    for idx, row in cnv_raw.iterrows():
        ne_id = row["NE_id"]
        platform = platform_map.get(ne_id, None)
        if platform not in common_cnvs:
            continue
        cnv_chrom = norm_chrom(row["seqnames"])
        cnv_type = row["CNV_type"]
        cnv_start_0based = int(row["start"]) - 1
        cnv_end = int(row["end"])
        if cnv_type in common_cnvs[platform]:
            for common_chrom, common_start, common_end in common_cnvs[platform][cnv_type]:
                if common_chrom == cnv_chrom:
                    ov_frac = one_direction_overlap(cnv_start_0based, cnv_end, common_start, common_end)
                    if ov_frac >= COMMON_CNV_OVERLAP_THR:
                        cnv_to_remove.add(idx)
                        removed_by_platform[platform] += 1
                        break

    n_before_filter = len(cnv_raw)
    cnv_raw = cnv_raw.drop(cnv_to_remove).copy()
    log(f"  Common CNV filtering: Agilent removed={removed_by_platform['Agilent FE']}, "
        f"Nimblegen removed={removed_by_platform['Nimblegen Normalized']}, "
        f"total {n_before_filter} -> {len(cnv_raw)}")

    # --- LiftOver hg18 → hg38 (v15互換) ---
    cnv_raw["SV_chrom_hg18"] = cnv_raw["seqnames"].apply(norm_chrom)
    cnv_raw["sv_start_hg18"] = cnv_raw["start"].astype(int) - 1
    cnv_raw["sv_end_hg18"] = cnv_raw["end"].astype(int)
    cnv_raw["_idx"] = range(len(cnv_raw))

    outdir = Path(_OUTDIR)
    outdir.mkdir(parents=True, exist_ok=True)
    tmp_bed_in = str(outdir / "cnv_hg18.bed")
    tmp_bed_out = str(outdir / "cnv_hg38.bed")
    tmp_bed_unmapped = str(outdir / "cnv_unmapped.bed")

    with open(tmp_bed_in, "w") as f:
        for _, row in cnv_raw.iterrows():
            f.write(f"{row['SV_chrom_hg18']}\t{row['sv_start_hg18']}\t{row['sv_end_hg18']}\t{row['_idx']}\n")

    cmd = [_LIFTOVER_BIN, tmp_bed_in, _CHAIN_FILE, tmp_bed_out, tmp_bed_unmapped, "-minMatch=0.5"]
    log(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log(f"  liftOver stderr: {result.stderr}")
        raise RuntimeError("liftOver failed")

    lifted = pd.read_csv(tmp_bed_out, sep="\t", header=None, names=["chrom_hg38", "start_hg38", "end_hg38", "idx"])
    n_unmapped = len(cnv_raw) - len(lifted)
    log(f"  LiftOver: {len(lifted)}/{len(cnv_raw)} mapped, {n_unmapped} unmapped")

    cnv_raw = cnv_raw.merge(lifted.rename(columns={"idx": "_idx"}), on="_idx", how="inner")
    cnv_raw["SV_chrom"] = cnv_raw["chrom_hg38"].apply(norm_chrom)
    cnv_raw["sv_start"] = cnv_raw["start_hg38"].astype(int)
    cnv_raw["sv_end"] = cnv_raw["end_hg38"].astype(int)
    cnv_raw["abs_len"] = (cnv_raw["sv_end"] - cnv_raw["sv_start"]).abs()
    log(f"  CNVs with hg38 coordinates: {len(cnv_raw)}")
    return cnv_raw, platform_map


# ============================================================================
# PHASE 2: Sample QC + discovery overlap exclusion
# v19: overlap_srWGS列（sample_data_01102023.xlsx内）で直接除外
# ============================================================================
def run_phase2_sample_qc():
    """Phase 2: サンプルQC + discovery重複除外（overlap_srWGS列ベース, ASD + SCZ）"""
    log("=" * 70)
    log("PHASE 2: Sample QC + discovery overlap exclusion (ASD + SCZ)")
    log("  v19: overlap_srWGS列ベースの除外方式")
    log("=" * 70)

    sample_df_raw = pd.read_excel(_SAMPLE_FILE)
    log(f"  Raw samples: {len(sample_df_raw)}")

    # --- overlap_srWGS列の存在確認 (v19) ---
    overlap_col = "overlap_srWGS"
    if overlap_col not in sample_df_raw.columns:
        raise ValueError(f"'{overlap_col}' 列が {_SAMPLE_FILE} に見つかりません。"
                         f"利用可能な列: {list(sample_df_raw.columns)}")

    sample_df_raw = sample_df_raw[sample_df_raw["pass"] == "pass"].copy()
    log(f"  After pass filter: {len(sample_df_raw)}")
    sample_df_raw = sample_df_raw[sample_df_raw["noisy_sample"].isna()].copy()
    log(f"  After noisy exclusion: {len(sample_df_raw)}")
    sample_df_raw = sample_df_raw[sample_df_raw["sex_mismatch"].isna()].copy()
    log(f"  After sex_mismatch exclusion: {len(sample_df_raw)}")
    sample_df_raw = sample_df_raw[sample_df_raw["quality_0.2"].isna()].copy()
    log(f"  After quality_0.2 exclusion: {len(sample_df_raw)}")

    # v16: 3群を対象
    sample_df_raw = sample_df_raw[sample_df_raw["diagnosis"].isin([CASE_LABEL, CTRL_LABEL, SCZ_LABEL])].copy()
    log(f"  ASD + SCZ + CONT: {len(sample_df_raw)}")

    # v16: CNV解析IDの重複検査
    cnv_id_col = "CNV解析ID"
    assert_unique_keys(sample_df_raw, [cnv_id_col], "sample_df_raw after QC")

    for dx in [CASE_LABEL, SCZ_LABEL, CTRL_LABEL]:
        log(f"    {dx}: {len(sample_df_raw[sample_df_raw['diagnosis'] == dx])}")

    # --- Discovery (srWGS) overlap exclusion (v19: overlap_srWGS列ベース) ---
    log("--- Excluding discovery srWGS overlapping samples (overlap_srWGS=Yes) ---")

    # overlap_srWGS列を文字列化して大文字比較
    overlap_flag = sample_df_raw[overlap_col].fillna("").astype(str).str.strip().str.upper()

    # ASD overlap exclusion
    n_asd_before = len(sample_df_raw[sample_df_raw["diagnosis"] == CASE_LABEL])
    overlap_mask_asd = (sample_df_raw["diagnosis"] == CASE_LABEL) & (overlap_flag == "YES")
    n_overlap_asd = int(overlap_mask_asd.sum())
    overlap_asd_ids = set(sample_df_raw.loc[overlap_mask_asd, cnv_id_col].values)
    sample_df_raw = sample_df_raw[~overlap_mask_asd].copy()
    overlap_flag = overlap_flag[~overlap_mask_asd]  # マスクと同期
    n_asd_after = len(sample_df_raw[sample_df_raw["diagnosis"] == CASE_LABEL])
    log(f"  ASD before overlap exclusion: {n_asd_before}")
    log(f"  Overlapping ASD (overlap_srWGS=Yes): {n_overlap_asd}")
    log(f"  ASD after overlap exclusion: {n_asd_after}")

    # SCZ overlap exclusion
    n_scz_before = len(sample_df_raw[sample_df_raw["diagnosis"] == SCZ_LABEL])
    overlap_mask_scz = (sample_df_raw["diagnosis"] == SCZ_LABEL) & (overlap_flag == "YES")
    n_overlap_scz = int(overlap_mask_scz.sum())
    overlap_scz_ids = set(sample_df_raw.loc[overlap_mask_scz, cnv_id_col].values)
    sample_df_raw = sample_df_raw[~overlap_mask_scz].copy()
    n_scz_after = len(sample_df_raw[sample_df_raw["diagnosis"] == SCZ_LABEL])
    log(f"  SCZ before overlap exclusion: {n_scz_before}")
    log(f"  Overlapping SCZ (overlap_srWGS=Yes): {n_overlap_scz}")
    log(f"  SCZ after overlap exclusion: {n_scz_after}")
    log(f"  CONT (unchanged): {len(sample_df_raw[sample_df_raw['diagnosis'] == CTRL_LABEL])}")

    # Save overlap ID lists
    outdir = Path(_OUTDIR)
    outdir.mkdir(parents=True, exist_ok=True)
    for label, ids in [("asd", overlap_asd_ids), ("scz", overlap_scz_ids)]:
        overlap_list_path = str(outdir / f"excluded_overlap_{label}_ids.txt")
        with open(overlap_list_path, "w") as f:
            for sid in sorted(ids):
                f.write(f"{sid}\n")
        log(f"  Saved overlap ID list: {overlap_list_path}")

    # Sex & Platform (v15互換)
    sample_df_raw["Sex_numeric"] = sample_df_raw["性別"].map({"男性": 1.0, "女性": 0.0})
    sample_df_raw = sample_df_raw.dropna(subset=["Sex_numeric"]).copy()
    sample_df_raw["Platform_nimblegen"] = (sample_df_raw["Data Type"] == "Nimblegen Normalized").astype(float)

    sample_info = {}
    for _, row in sample_df_raw.iterrows():
        sid = row[cnv_id_col]
        sample_info[sid] = {
            "diagnosis": row["diagnosis"],
            "sex": row["Sex_numeric"],
            "platform": row["Platform_nimblegen"],
        }

    log(f"  Valid samples: {len(sample_info)}")
    for dx in [CASE_LABEL, SCZ_LABEL, CTRL_LABEL]:
        n_dx = sum(1 for v in sample_info.values() if v["diagnosis"] == dx)
        n_ag = sum(1 for v in sample_info.values() if v["diagnosis"] == dx and v["platform"] == 0.0)
        n_ni = sum(1 for v in sample_info.values() if v["diagnosis"] == dx and v["platform"] == 1.0)
        log(f"  {dx}: total={n_dx}, Agilent={n_ag}, Nimblegen={n_ni}")

    return sample_info, n_overlap_asd, n_overlap_scz


# ============================================================================
# PHASE 3: CNV filtering + NAHR flagging + gene counting
# *** v15から完全移植 ***
# ============================================================================
def run_phase3_filter_cnv(cnv_raw, sample_info):
    """Phase 3: CNVフィルタ + NAHR GD flag + 遺伝子カウント (v15互換)"""
    log("=" * 70)
    log("PHASE 3: CNV filtering + NAHR flagging + gene counting")
    log("=" * 70)

    cnv_raw["sample"] = cnv_raw["NE_id"]
    cnv_raw = cnv_raw[cnv_raw["sample"].isin(sample_info.keys())].copy()
    log(f"  After sample filter: {len(cnv_raw)}")

    autosome_chroms = {f"chr{i}" for i in range(1, 23)}
    cnv_raw = cnv_raw[cnv_raw["SV_chrom"].isin(autosome_chroms)].copy()
    log(f"  After autosome filter: {len(cnv_raw)}")

    cnv_raw = cnv_raw[cnv_raw["abs_len"] >= MIN_SV_LEN].copy()
    log(f"  After size filter (>={MIN_SV_LEN / 1000:.0f}kb): {len(cnv_raw)}")

    cnv_raw = cnv_raw[cnv_raw["CNV_type"].isin(SV_TYPES)].copy()
    log(f"  After DEL/DUP filter: {len(cnv_raw)}")

    # --- v18追加: segdup filter (WGS v18と同一, hg38座標に対して適用) ---
    log(f"--- Segdup filter (< {SEGDUP_MAX_PCT}%) ---")
    log(f"  Reading segdup BED: {_SEGDUP_BED}")
    segdup_raw = load_bed3(Path(_SEGDUP_BED))
    segdup_index = build_merged_interval_index(segdup_raw)
    segdup_pcts = np.array([
        compute_coverage_from_index(c, s, e, segdup_index)
        for c, s, e in zip(cnv_raw["SV_chrom"], cnv_raw["sv_start"], cnv_raw["sv_end"])
    ])
    cnv_raw["segdup_pct"] = segdup_pcts
    n_before_segdup = len(cnv_raw)
    cnv_raw = cnv_raw[cnv_raw["segdup_pct"] < SEGDUP_MAX_PCT].copy()
    log(f"  After segdup filter: {n_before_segdup} -> {len(cnv_raw)} "
        f"({n_before_segdup - len(cnv_raw)} removed)")

    # --- v18追加: exclusion BED filter (WGS v18と同一, hg38座標に対して適用) ---
    log(f"--- Exclusion BED filter (< {EXCLUSION_OVERLAP_THR * 100:.0f}%) ---")
    log(f"  Reading exclusion BED: {_EXCLUSION_BED}")
    excl_raw = load_bed3(Path(_EXCLUSION_BED))
    excl_index = build_merged_interval_index(excl_raw)
    excl_pcts = np.array([
        compute_coverage_from_index(c, s, e, excl_index)
        for c, s, e in zip(cnv_raw["SV_chrom"], cnv_raw["sv_start"], cnv_raw["sv_end"])
    ])
    cnv_raw["exclusion_pct"] = excl_pcts
    n_before_excl = len(cnv_raw)
    cnv_raw = cnv_raw[cnv_raw["exclusion_pct"] < EXCLUSION_OVERLAP_THR * 100.0].copy()
    log(f"  After exclusion BED filter: {n_before_excl} -> {len(cnv_raw)} "
        f"({n_before_excl - len(cnv_raw)} removed)")

    sv = cnv_raw.copy()
    log(f"  Final CNVs for analysis: {len(sv)}")

    # --- NAHR GD flagging (v15互換: nahr列参照 + CNV type照合) ---
    log("--- Loading NAHR GD loci ---")
    gd_df = pd.read_csv(_CURATED_GD_FILE, sep="\t", dtype=str)
    col_map = {c.lower(): c for c in gd_df.columns}
    gd_df["nahr_bool"] = gd_df[col_map["nahr"]].str.strip().str.lower().isin(["true", "1", "yes"])
    gd_nahr = gd_df[gd_df["nahr_bool"]].copy()
    gd_nahr["chrom_norm"] = gd_nahr[col_map["chr"]].apply(norm_chrom)
    gd_nahr["start_int"] = gd_nahr[col_map["start"]].astype(int)
    gd_nahr["end_int"] = gd_nahr[col_map["end"]].astype(int)
    gd_nahr["cnv_type"] = gd_nahr[col_map["cnv"]].str.upper().str.strip()
    gd_nahr["gd_id"] = gd_nahr[col_map["gd_id"]].str.strip()
    log(f"  Loaded {len(gd_nahr)} NAHR GD loci")

    def is_nahr_gd_cnv(chrom, start, end, svtype):
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

    nahr_flags = np.array([is_nahr_gd_cnv(
        sv["SV_chrom"].values[i], sv["sv_start"].values[i], sv["sv_end"].values[i],
        sv["CNV_type"].values[i]) for i in range(len(sv))])
    sv["is_nahr_gd"] = nahr_flags.astype(int)
    log(f"  NAHR GD flagged: {int(sv['is_nahr_gd'].sum())}/{len(sv)}")

    # --- Gene counting (v15互換) ---
    log("--- Building gene interval index from GENCODE ---")
    gene_index = defaultdict(list)
    n_genes = 0
    with gzip.open(_GTF_FILE, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            if 'gene_type "protein_coding"' not in attrs:
                continue
            chrom = norm_chrom(parts[0])
            start = int(parts[3])
            end = int(parts[4])
            m = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_name"):
                    m = attr.split('"')[1]
                    break
            if m:
                gene_index[chrom].append((start, end, m))
                n_genes += 1
    for chrom in gene_index:
        gene_index[chrom].sort()
    log(f"  Loaded {n_genes} protein-coding genes")

    def count_overlapping_genes(chrom, sv_start, sv_end):
        if chrom not in gene_index:
            return 0
        genes = set()
        for gs, ge, gname in gene_index[chrom]:
            if gs >= sv_end:
                break
            if ge > sv_start:
                genes.add(gname)
        return len(genes)

    sv["gene_count"] = [count_overlapping_genes(
        sv["SV_chrom"].values[i], sv["sv_start"].values[i], sv["sv_end"].values[i]
    ) for i in range(len(sv))]
    log(f"  Gene counts computed. Mean={sv['gene_count'].mean():.2f}, Max={sv['gene_count'].max()}")

    return sv, gd_nahr


# ============================================================================
# PHASE 4: Load boundary bins from L2 diffbound BED files
# *** v15から完全移植 (diffbound除外付き static_all) ***
# ============================================================================
def run_phase4_load_bins():
    """Phase 4: L2 diffbound BEDから10 lineage class + static_all のbin setsを構築 (v15互換)"""
    log("=" * 70)
    log("PHASE 4: Loading boundary bins from L2 diffbound BED files")
    log("=" * 70)

    import anndata as ad

    # Read all diffbound BED files (v15互換: *_diffbound.bed.gz glob)
    diff_dir = Path(_L2_DIFFBOUND_DIR)
    diff_files = sorted(diff_dir.glob("*_diffbound.bed.gz"))
    if len(diff_files) == 0:
        raise FileNotFoundError(f"No diffbound files found in: {diff_dir}")

    all_diffbound_bins = set()  # for static exclusion
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

    # Filter by min_bin_threshold (v15互換)
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
                             f"Available: {list(bin_sets.keys())}. Excluded: {excluded_classes}")

    # Static bins from raw h5ad (support=2, EXCLUDE diffbound bins — v15互換)
    log("  Reading static bins from raw h5ad (support=2, excluding diffbound)")
    adata_raw = ad.read_h5ad(_RAW_H5AD)
    raw_support2 = extract_nonzero_sparse_long(adata_raw, "raw", value_filter=2.0)
    raw_support2["chrom"] = raw_support2["chrom"].apply(norm_chrom)
    static_bins = raw_support2.loc[~raw_support2["bin_id"].isin(all_diffbound_bins),
                                    ["chrom", "start0", "end", "bin_id"]].drop_duplicates("bin_id")
    bin_sets["static_all"] = static_bins
    log(f"  static_all: {len(static_bins)} bins (after excluding {len(all_diffbound_bins)} diffbound)")

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
# PHASE 5: Compute unique disrupted bin counts per sample (Pattern A + C)
# *** v15から完全移植 ***
# ============================================================================
def run_phase5_bin_counts(sv, sample_info, bin_index):
    """Phase 5: サンプルごとのunique disrupted bin count (v15互換)"""
    log("=" * 70)
    log("PHASE 5: Computing unique disrupted bin counts (Pattern A + Pattern C)")
    log("=" * 70)

    ALL_BCLASSES = PRIMARY_BOUNDARY_CLASSES + ["static_all"]
    all_samples = sorted(sample_info.keys())
    sample_to_idx = {s: i for i, s in enumerate(all_samples)}
    N_samples = len(all_samples)

    disrupted_bins_a = {}
    for bclass in ALL_BCLASSES:
        disrupted_bins_a[bclass] = {svt: [set() for _ in range(N_samples)] for svt in SV_TYPES}

    disrupted_bins_c = {}
    for bclass in ALL_BCLASSES:
        disrupted_bins_c[bclass] = {svt: [set() for _ in range(N_samples)] for svt in SV_TYPES}

    sample_total_base = {svt: np.zeros(N_samples) for svt in SV_TYPES}
    sample_total_gene = {svt: np.zeros(N_samples) for svt in SV_TYPES}
    sample_total_base_c = {svt: np.zeros(N_samples) for svt in SV_TYPES}
    sample_total_gene_c = {svt: np.zeros(N_samples) for svt in SV_TYPES}

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

        if is_nahr:
            continue
        sample_total_base[svt][si] += abs_len
        sample_total_gene[svt][si] += int(row["gene_count"])

        if abs_len <= MAX_SV_LEN_PATTERN_C:
            sample_total_base_c[svt][si] += abs_len
            sample_total_gene_c[svt][si] += int(row["gene_count"])

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

    # v17: disrupted_bins_a/c の削除を Phase 5B 後に遅延
    # (Phase 5B でファイル出力に使用するため)

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
            sample_total_base_c, sample_total_gene_c,
            disrupted_bins_a, disrupted_bins_c)


# ============================================================================
# PHASE 5B: Per-sample event-bin overlap + sample covariates 出力 (v17新規)
# ============================================================================
def run_phase5b_save_bin_overlap(all_samples, sample_to_idx, sample_info,
                                 sv, bin_index,
                                 disrupted_bins_a, disrupted_bins_c,
                                 sample_total_base, sample_total_gene,
                                 sample_total_base_c, sample_total_gene_c):
    """Phase 5B (v17新規): per-sample event-bin overlap とサンプル共変量をTSV出力。
    factorial v4 スクリプトで arrayCGH factorial 解析に使用する。

    SVを再スキャンし、event_id（各CNV固有のID）付きでbin overlapを記録。
    これによりv4でbin_count, event_count, carrierの3つのexposure定義を計算可能。"""
    log("=" * 70)
    log("PHASE 5B (v17): Saving per-sample event-bin overlap + sample covariates")
    log("=" * 70)

    outdir = Path(_OUTDIR)
    ALL_BCLASSES = PRIMARY_BOUNDARY_CLASSES + ["static_all"]
    N_samples = len(all_samples)

    # --- 1. Per-sample event-bin overlap TSV (SV再スキャン、event_id付き) ---
    overlap_path = outdir / "sample_event_bin_overlap_v22.tsv.gz"
    n_rows = 0
    with gzip.open(overlap_path, "wt") as f:
        f.write("sample_id\tbin_id\tevent_id\tsv_type\tpattern\n")
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

            if is_nahr:
                continue

            sid = all_samples[si]
            event_id = f"{sid}:{chrom}:{sv_s}-{sv_e}:{svt}"

            # Pattern A bins (全SV)
            bins_hit_a = set()
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
                            bins_hit_a.add(bid_arr[j])
            for bid in bins_hit_a:
                f.write(f"{sid}\t{bid}\t{event_id}\t{svt}\tA\n")
                n_rows += 1

            # Pattern C bins (abs_len <= 1MB)
            if abs_len <= MAX_SV_LEN_PATTERN_C:
                bins_hit_c = set()
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
                                bins_hit_c.add(bid_arr[j])
                for bid in bins_hit_c:
                    f.write(f"{sid}\t{bid}\t{event_id}\t{svt}\tC\n")
                    n_rows += 1

    log(f"  Saved: {overlap_path} ({n_rows:,} rows)")

    # --- 2. Sample covariates TSV ---
    covar_path = outdir / "sample_covariates_v22.tsv"
    covar_rows = []
    for si, sid in enumerate(all_samples):
        info = sample_info.get(sid, {})
        covar_rows.append({
            "sample_id": sid,
            "diagnosis": info.get("diagnosis", ""),
            "sex": info.get("sex", ""),
            "platform_nimblegen": info.get("platform", ""),
            "log1p_total_del_bases_A": float(np.log1p(sample_total_base["DEL"][si])),
            "log1p_total_gene_DEL_A": float(np.log1p(sample_total_gene["DEL"][si])),
            "log1p_total_del_bases_C": float(np.log1p(sample_total_base_c["DEL"][si])),
            "log1p_total_gene_DEL_C": float(np.log1p(sample_total_gene_c["DEL"][si])),
            "log1p_total_dup_bases_A": float(np.log1p(sample_total_base["DUP"][si])),
            "log1p_total_gene_DUP_A": float(np.log1p(sample_total_gene["DUP"][si])),
        })
    covar_df = pd.DataFrame(covar_rows)
    covar_df.to_csv(covar_path, sep="\t", index=False)
    log(f"  Saved: {covar_path} ({len(covar_df)} samples)")

    for dx in [CASE_LABEL, SCZ_LABEL, CTRL_LABEL]:
        n_dx = int((covar_df["diagnosis"] == dx).sum())
        log(f"    {dx}: {n_dx}")

    # --- 3. メモリ解放 ---
    del disrupted_bins_a, disrupted_bins_c
    log("  Released disrupted_bins memory")

    return overlap_path, covar_path


# ============================================================================
# PHASE 6 / 6B: Binary logistic regression helper
# ============================================================================
def _run_binary_regression(group_label, case_label, ctrl_label,
                           all_samples, sample_info,
                           bin_counts_a, bin_counts_c,
                           sample_total_base, sample_total_gene,
                           sample_total_base_c, sample_total_gene_c,
                           is_replication=True):
    """Binary logistic regression (shared by Phase 6 and Phase 6B)."""
    log(f"  === {group_label}: {case_label} vs {ctrl_label} ===")

    ALL_BCLASSES = PRIMARY_BOUNDARY_CLASSES + ["static_all"]
    N_samples = len(all_samples)

    rows = []
    for si, sid in enumerate(all_samples):
        info = sample_info[sid]
        rows.append([si, sid, info["diagnosis"], info["sex"], info["platform"]])
    cols = ["sample_idx", "SampleID", "Status", "Sex", "Platform_nimblegen"]
    sample_df = pd.DataFrame(rows, columns=cols)

    case_mask = (sample_df["Status"] == case_label).values
    ctrl_mask = (sample_df["Status"] == ctrl_label).values
    keep = case_mask | ctrl_mask
    y = case_mask[keep].astype(int)
    n_case, n_ctrl = int(y.sum()), int((1 - y).sum())
    log(f"  n_case={n_case} ({case_label}), n_ctrl={n_ctrl} ({ctrl_label})")
    valid_si = sample_df.loc[keep, "sample_idx"].values

    plat_vals_all = sample_df.loc[keep, "Platform_nimblegen"].values
    agilent_mask = (plat_vals_all == 0)
    nimblegen_mask = (plat_vals_all == 1)
    y_ag = y[agilent_mask]
    y_ni = y[nimblegen_mask]
    n_ag = int(agilent_mask.sum())
    n_ni = int(nimblegen_mask.sum())
    log(f"  Platform: Agilent={n_ag}, Nimblegen={n_ni}")

    all_results = []

    for pattern, bin_counts in [("A", bin_counts_a), ("C", bin_counts_c)]:
        if pattern == "C":
            _total_base = sample_total_base_c
            _total_gene = sample_total_gene_c
        else:
            _total_base = sample_total_base
            _total_gene = sample_total_gene

        for svt in SV_TYPES:
            log1p_base = np.log1p(_total_base[svt][valid_si])
            log1p_gene = np.log1p(_total_gene[svt][valid_si])
            sex_vals = sample_df.loc[keep, "Sex"].values
            plat_vals = sample_df.loc[keep, "Platform_nimblegen"].values

            X_null = np.column_stack([np.ones(int(keep.sum())), sex_vals,
                                      log1p_base, log1p_gene, plat_vals])
            null_llf = np.nan
            for m in ["newton", "lbfgs", "bfgs"]:
                try:
                    res = sm.Logit(y, X_null).fit(disp=0, maxiter=200, method=m)
                    if np.isfinite(res.llf):
                        null_llf = float(res.llf)
                        break
                except Exception:
                    continue
            if not np.isfinite(null_llf):
                continue

            X_null_ag = np.column_stack([np.ones(n_ag), sex_vals[agilent_mask],
                                          log1p_base[agilent_mask], log1p_gene[agilent_mask]])
            null_llf_ag = np.nan
            for m in ["newton", "lbfgs", "bfgs"]:
                try:
                    res = sm.Logit(y_ag, X_null_ag).fit(disp=0, maxiter=200, method=m)
                    if np.isfinite(res.llf):
                        null_llf_ag = float(res.llf)
                        break
                except Exception:
                    continue

            X_null_ni = np.column_stack([np.ones(n_ni), sex_vals[nimblegen_mask],
                                          log1p_base[nimblegen_mask], log1p_gene[nimblegen_mask]])
            null_llf_ni = np.nan
            for m in ["newton", "lbfgs", "bfgs"]:
                try:
                    res = sm.Logit(y_ni, X_null_ni).fit(disp=0, maxiter=200, method=m)
                    if np.isfinite(res.llf):
                        null_llf_ni = float(res.llf)
                        break
                except Exception:
                    continue

            for bclass in ALL_BCLASSES:
                is_primary = (pattern == "A" and svt == "DEL" and bclass in PRIMARY_BOUNDARY_CLASSES)
                is_secondary_c = (pattern == "C" and svt == "DEL" and bclass in PRIMARY_BOUNDARY_CLASSES)
                is_negative_ctrl = (svt == "DUP" and bclass in PRIMARY_BOUNDARY_CLASSES) or \
                                   (svt == "DEL" and bclass == "static_all")
                if not is_primary and not is_secondary_c and not is_negative_ctrl:
                    continue

                bsv_col = bin_counts[bclass][svt][valid_si].astype(float)
                n_exposed_case = int(np.sum(bsv_col[y == 1] > 0))
                n_exposed_ctrl = int(np.sum(bsv_col[y == 0] > 0))

                # Pooled regression
                X_full = np.column_stack([X_null[:, 0], bsv_col, X_null[:, 1:]])
                obs_beta = obs_se = obs_or = obs_p_twosided = np.nan
                for m in ["newton", "lbfgs", "bfgs"]:
                    try:
                        res = sm.Logit(y, X_full).fit(disp=0, maxiter=200, method=m)
                        if np.isfinite(res.llf) and np.isfinite(res.params[1]):
                            obs_beta = float(res.params[1])
                            obs_se = float(res.bse[1])
                            obs_or = float(np.exp(res.params[1]))
                            lr = 2.0 * (res.llf - null_llf)
                            obs_p_twosided = float(chi2.sf(max(lr, 0), df=1))
                            break
                    except Exception:
                        continue

                if np.isfinite(obs_p_twosided) and np.isfinite(obs_beta):
                    if is_replication and bclass in DISCOVERY_SIG_CLASSES:
                        obs_p_onesided = obs_p_twosided / 2.0 if obs_beta > 0 else 1.0 - obs_p_twosided / 2.0
                    else:
                        obs_p_onesided = obs_p_twosided
                else:
                    obs_p_onesided = np.nan

                ci_lo = float(np.exp(obs_beta - 1.96 * obs_se)) if np.isfinite(obs_se) else np.nan
                ci_hi = float(np.exp(obs_beta + 1.96 * obs_se)) if np.isfinite(obs_se) else np.nan

                pooled_repl_sig = False
                if is_replication and (is_primary or is_secondary_c):
                    if bclass in DISCOVERY_SIG_CLASSES:
                        pooled_repl_sig = (np.isfinite(obs_p_onesided) and obs_p_onesided < 0.05 and obs_beta > 0)
                    else:
                        pooled_repl_sig = (np.isfinite(obs_p_onesided) and obs_p_onesided < 0.05)

                # Platform-stratified
                bsv_col_ag = bsv_col[agilent_mask]
                agilent_beta = agilent_se = agilent_or = agilent_p_twosided = agilent_p_onesided = np.nan
                if np.isfinite(null_llf_ag) and n_ag > 0:
                    X_full_ag = np.column_stack([X_null_ag[:, 0], bsv_col_ag, X_null_ag[:, 1:]])
                    for m in ["newton", "lbfgs", "bfgs"]:
                        try:
                            res = sm.Logit(y_ag, X_full_ag).fit(disp=0, maxiter=200, method=m)
                            if np.isfinite(res.llf) and np.isfinite(res.params[1]):
                                agilent_beta = float(res.params[1])
                                agilent_se = float(res.bse[1])
                                agilent_or = float(np.exp(res.params[1]))
                                lr = 2.0 * (res.llf - null_llf_ag)
                                agilent_p_twosided = float(chi2.sf(max(lr, 0), df=1))
                                if np.isfinite(agilent_p_twosided):
                                    if is_replication and bclass in DISCOVERY_SIG_CLASSES:
                                        agilent_p_onesided = agilent_p_twosided / 2.0 if agilent_beta > 0 else 1.0 - agilent_p_twosided / 2.0
                                    else:
                                        agilent_p_onesided = agilent_p_twosided
                                break
                        except Exception:
                            continue
                agilent_ci_lo = float(np.exp(agilent_beta - 1.96 * agilent_se)) if np.isfinite(agilent_se) else np.nan
                agilent_ci_hi = float(np.exp(agilent_beta + 1.96 * agilent_se)) if np.isfinite(agilent_se) else np.nan

                bsv_col_ni = bsv_col[nimblegen_mask]
                nimblegen_beta = nimblegen_se = nimblegen_or = nimblegen_p_twosided = nimblegen_p_onesided = np.nan
                if np.isfinite(null_llf_ni) and n_ni > 0:
                    X_full_ni = np.column_stack([X_null_ni[:, 0], bsv_col_ni, X_null_ni[:, 1:]])
                    for m in ["newton", "lbfgs", "bfgs"]:
                        try:
                            res = sm.Logit(y_ni, X_full_ni).fit(disp=0, maxiter=200, method=m)
                            if np.isfinite(res.llf) and np.isfinite(res.params[1]):
                                nimblegen_beta = float(res.params[1])
                                nimblegen_se = float(res.bse[1])
                                nimblegen_or = float(np.exp(res.params[1]))
                                lr = 2.0 * (res.llf - null_llf_ni)
                                nimblegen_p_twosided = float(chi2.sf(max(lr, 0), df=1))
                                if np.isfinite(nimblegen_p_twosided):
                                    if is_replication and bclass in DISCOVERY_SIG_CLASSES:
                                        nimblegen_p_onesided = nimblegen_p_twosided / 2.0 if nimblegen_beta > 0 else 1.0 - nimblegen_p_twosided / 2.0
                                    else:
                                        nimblegen_p_onesided = nimblegen_p_twosided
                                break
                        except Exception:
                            continue
                nimblegen_ci_lo = float(np.exp(nimblegen_beta - 1.96 * nimblegen_se)) if np.isfinite(nimblegen_se) else np.nan
                nimblegen_ci_hi = float(np.exp(nimblegen_beta + 1.96 * nimblegen_se)) if np.isfinite(nimblegen_se) else np.nan

                # IVW Meta-analysis
                platform_meta_beta = platform_meta_se = platform_meta_or = np.nan
                platform_meta_ci_lo = platform_meta_ci_hi = np.nan
                platform_meta_p_main = platform_meta_p_two = platform_meta_p_one_positive = np.nan
                platform_meta_q = platform_meta_q_p = platform_meta_i2 = np.nan
                platform_same_direction = False

                if (np.isfinite(agilent_beta) and np.isfinite(agilent_se) and
                    np.isfinite(nimblegen_beta) and np.isfinite(nimblegen_se) and
                    agilent_se > 0 and nimblegen_se > 0):

                    w_ag = 1.0 / (agilent_se ** 2)
                    w_ni = 1.0 / (nimblegen_se ** 2)
                    sum_w = w_ag + w_ni
                    platform_meta_beta = (w_ag * agilent_beta + w_ni * nimblegen_beta) / sum_w
                    platform_meta_se = float(np.sqrt(1.0 / sum_w))
                    platform_meta_or = float(np.exp(platform_meta_beta))
                    platform_meta_ci_lo = float(np.exp(platform_meta_beta - 1.96 * platform_meta_se))
                    platform_meta_ci_hi = float(np.exp(platform_meta_beta + 1.96 * platform_meta_se))
                    platform_same_direction = (agilent_beta > 0 and nimblegen_beta > 0) or \
                                            (agilent_beta < 0 and nimblegen_beta < 0)
                    z_meta = platform_meta_beta / platform_meta_se if platform_meta_se > 0 else np.nan

                    if is_negative_ctrl or not (is_replication and bclass in DISCOVERY_SIG_CLASSES):
                        if np.isfinite(z_meta):
                            platform_meta_p_main = 2.0 * norm.sf(abs(z_meta))
                            platform_meta_p_two = platform_meta_p_main
                        else:
                            platform_meta_p_main = platform_meta_p_two = np.nan
                    else:
                        if np.isfinite(z_meta):
                            platform_meta_p_main = norm.sf(z_meta)
                            platform_meta_p_one_positive = platform_meta_p_main
                            platform_meta_p_two = 2.0 * norm.sf(abs(z_meta))
                        else:
                            platform_meta_p_main = platform_meta_p_two = np.nan

                    q_ag = w_ag * (agilent_beta - platform_meta_beta) ** 2
                    q_ni = w_ni * (nimblegen_beta - platform_meta_beta) ** 2
                    platform_meta_q = q_ag + q_ni
                    platform_meta_q_p = float(chi2.sf(max(platform_meta_q, 0), df=1))
                    platform_meta_i2 = max(0.0, (platform_meta_q - 1.0) / platform_meta_q) if platform_meta_q > 0 else 0.0

                platform_meta_repl_sig = False
                if is_replication and (is_primary or is_secondary_c):
                    if np.isfinite(platform_meta_p_main):
                        if bclass in DISCOVERY_SIG_CLASSES:
                            platform_meta_repl_sig = (platform_meta_p_main < 0.05 and
                                                     platform_same_direction and
                                                     platform_meta_beta > 0)
                        else:
                            platform_meta_repl_sig = (platform_meta_p_main < 0.05)

                if is_primary:
                    analysis_type = "PRIMARY_replication" if is_replication else "PRIMARY_SZ_burden"
                elif is_secondary_c:
                    analysis_type = "SECONDARY_pattern_c" if is_replication else "SECONDARY_SZ_pattern_c"
                else:
                    analysis_type = "negative_control"
                    obs_p_onesided = np.nan
                    agilent_p_onesided = np.nan
                    nimblegen_p_onesided = np.nan

                formal_test = bool(is_replication and is_primary and bclass in DISCOVERY_SIG_CLASSES)
                exploratory_test = bool(is_replication and is_primary and bclass in DISCOVERY_NONSIG_CLASSES)

                result_row = {
                    "comparison": f"{case_label}_vs_{ctrl_label}",
                    "Pattern": pattern, "Boundary_class": bclass, "SV_type": svt,
                    "N_case": n_case, "N_control": n_ctrl, "N_complete": n_case + n_ctrl,
                    "N_exposed_case": n_exposed_case, "N_exposed_control": n_exposed_ctrl,
                    "Beta": obs_beta, "SE": obs_se, "OR": obs_or,
                    "CI_lo": ci_lo, "CI_hi": ci_hi,
                    "P_twosided": obs_p_twosided, "P_onesided": obs_p_onesided,
                    "direction_match": bool(obs_beta > 0) if np.isfinite(obs_beta) else False,
                    "analysis_type": analysis_type,
                    "is_primary": is_primary, "is_secondary_c": is_secondary_c,
                    "pooled_replicated": pooled_repl_sig,
                    "formal_test": formal_test, "exploratory_test": exploratory_test,
                    "formal_replicated": bool(formal_test and pooled_repl_sig),
                    "exploratory_signal": bool(exploratory_test and pooled_repl_sig),
                    "agilent_beta": agilent_beta, "agilent_se": agilent_se,
                    "agilent_or": agilent_or, "agilent_ci_lo": agilent_ci_lo, "agilent_ci_hi": agilent_ci_hi,
                    "agilent_p_twosided": agilent_p_twosided, "agilent_p_onesided": agilent_p_onesided,
                    "agilent_n_case": int(y_ag.sum()), "agilent_n_control": int((1 - y_ag).sum()),
                    "nimblegen_beta": nimblegen_beta, "nimblegen_se": nimblegen_se,
                    "nimblegen_or": nimblegen_or, "nimblegen_ci_lo": nimblegen_ci_lo, "nimblegen_ci_hi": nimblegen_ci_hi,
                    "nimblegen_p_twosided": nimblegen_p_twosided, "nimblegen_p_onesided": nimblegen_p_onesided,
                    "nimblegen_n_case": int(y_ni.sum()), "nimblegen_n_control": int((1 - y_ni).sum()),
                    "platform_meta_beta": platform_meta_beta, "platform_meta_se": platform_meta_se,
                    "platform_meta_or": platform_meta_or,
                    "platform_meta_ci_lo": platform_meta_ci_lo, "platform_meta_ci_hi": platform_meta_ci_hi,
                    "platform_meta_p_main": platform_meta_p_main,
                    "platform_meta_p_two": platform_meta_p_two,
                    "platform_meta_p_one_positive": platform_meta_p_one_positive,
                    "platform_meta_q": platform_meta_q, "platform_meta_q_p": platform_meta_q_p,
                    "platform_meta_i2": platform_meta_i2,
                    "platform_same_direction": platform_same_direction,
                    "platform_meta_replicated": platform_meta_repl_sig,
                    "formal_meta_replicated": bool(is_replication and is_primary and bclass in DISCOVERY_SIG_CLASSES and platform_meta_repl_sig),
                    "exploratory_meta_signal": bool(is_replication and is_primary and bclass in DISCOVERY_NONSIG_CLASSES and platform_meta_repl_sig),
                }
                all_results.append(result_row)

    return all_results, n_case, n_ctrl


def run_phase6_asd_regression(all_samples, sample_info, bin_counts_a, bin_counts_c,
                              sample_total_base, sample_total_gene,
                              sample_total_base_c, sample_total_gene_c):
    """Phase 6: ASD vs CONT (PRIMARY replication)"""
    log("=" * 70)
    log("PHASE 6: ASD vs CONT Regression (PRIMARY replication)")
    log("=" * 70)
    return _run_binary_regression(
        "Phase 6", CASE_LABEL, CTRL_LABEL, all_samples, sample_info,
        bin_counts_a, bin_counts_c, sample_total_base, sample_total_gene,
        sample_total_base_c, sample_total_gene_c, is_replication=True)


# ============================================================================
# PHASE 6B: SZ vs CONT (v16 新規)
# ============================================================================
def run_phase6b_sz_regression(all_samples, sample_info, bin_counts_a, bin_counts_c,
                              sample_total_base, sample_total_gene,
                              sample_total_base_c, sample_total_gene_c):
    """Phase 6B: SZ vs CONT (external SZ burden assessment)"""
    log("=" * 70)
    log("PHASE 6B: SZ vs CONT Regression (external SZ burden assessment)")
    log("=" * 70)
    return _run_binary_regression(
        "Phase 6B", SCZ_LABEL, CTRL_LABEL, all_samples, sample_info,
        bin_counts_a, bin_counts_c, sample_total_base, sample_total_gene,
        sample_total_base_c, sample_total_gene_c, is_replication=False)


# ============================================================================
# PHASE 6C: MNLogit heterogeneity test (v16 新規)
# + Platform別感度解析 (Agilent-only, Nimblegen-only)
# ============================================================================
def _run_mnlogit_heterogeneity(outcome, X_df, x_colnames, label="pooled"):
    """MNLogit helper: fit model and extract z_het."""
    mnl_model = None
    converged = False
    for fit_method in ['newton', 'bfgs', 'lbfgs']:
        try:
            candidate = sm.MNLogit(outcome, X_df).fit(disp=0, maxiter=200, method=fit_method)
            _ = candidate.cov_params()
            mnl_model = candidate
            converged = True
            break
        except Exception:
            continue

    if mnl_model is None:
        return None

    try:
        beta_asd = float(mnl_model.params.loc["burden", 0])
        beta_sz = float(mnl_model.params.loc["burden", 1])
        cov_full = mnl_model.cov_params()

        try:
            var_b_asd = float(cov_full.loc[('1', 'burden'), ('1', 'burden')])
            var_b_sz = float(cov_full.loc[('2', 'burden'), ('2', 'burden')])
            cov_b_asd_sz = float(cov_full.loc[('1', 'burden'), ('2', 'burden')])
        except KeyError:
            n_params = len(x_colnames)
            var_b_asd = float(cov_full.iloc[1, 1])
            var_b_sz = float(cov_full.iloc[1 + n_params, 1 + n_params])
            cov_b_asd_sz = float(cov_full.iloc[1, 1 + n_params])

        diff_beta = beta_asd - beta_sz
        se_diff = np.sqrt(var_b_asd + var_b_sz - 2 * cov_b_asd_sz)

        if se_diff <= 0 or np.isnan(se_diff):
            return None

        z_het = diff_beta / se_diff
        p_het_twosided = 2 * scipy_stats.norm.sf(abs(z_het))

        return {
            "beta_asd": beta_asd, "beta_sz": beta_sz,
            "z_het": z_het, "se_diff": se_diff,
            "diff_beta": diff_beta, "p_het_twosided": p_het_twosided,
            "cov_b_asd_sz": cov_b_asd_sz, "converged": converged,
        }
    except Exception:
        return None


def run_phase6c_heterogeneity(all_samples, sample_info, bin_counts_a,
                               sample_total_base, sample_total_gene):
    """Phase 6C: ASD-SZ heterogeneity test (Pattern A, DEL only)
    + Platform別感度解析
    """
    log("=" * 70)
    log("PHASE 6C: ASD-SZ Heterogeneity Test (MNLogit)")
    log("  z_het = (beta_ASD - beta_SZ) / sqrt(Var_ASD + Var_SZ - 2*Cov)")
    log("  Global: Correlation-adjusted Stouffer + sign test")
    log("  Sensitivity: Agilent-only + Nimblegen-only heterogeneity")
    log("=" * 70)

    N_samples = len(all_samples)
    sample_to_idx = {s: i for i, s in enumerate(all_samples)}

    rows = []
    for si, sid in enumerate(all_samples):
        info = sample_info[sid]
        rows.append([si, sid, info["diagnosis"], info["sex"], info["platform"]])
    cols = ["sample_idx", "SampleID", "Status", "Sex", "Platform_nimblegen"]
    sample_df = pd.DataFrame(rows, columns=cols)

    asd_mask = (sample_df["Status"] == CASE_LABEL).values
    scz_mask = (sample_df["Status"] == SCZ_LABEL).values
    ctrl_mask = (sample_df["Status"] == CTRL_LABEL).values
    keep_all = asd_mask | scz_mask | ctrl_mask
    valid_si_all = sample_df.loc[keep_all, "sample_idx"].values

    outcome_all = np.zeros(int(keep_all.sum()), dtype=int)
    outcome_all[asd_mask[keep_all]] = 1
    outcome_all[scz_mask[keep_all]] = 2

    n_asd = int(np.sum(outcome_all == 1))
    n_scz = int(np.sum(outcome_all == 2))
    n_ctrl = int(np.sum(outcome_all == 0))
    log(f"  ASD={n_asd}, SCZ={n_scz}, CONT={n_ctrl}")

    sex_all = sample_df.loc[keep_all, "Sex"].values
    plat_all = sample_df.loc[keep_all, "Platform_nimblegen"].values

    # Platform masks within keep_all
    agilent_within = (plat_all == 0)
    nimblegen_within = (plat_all == 1)

    svt = "DEL"
    log1p_base_all = np.log1p(sample_total_base[svt][valid_si_all])
    log1p_gene_all = np.log1p(sample_total_gene[svt][valid_si_all])

    het_results = []
    z_het_list_pooled = []
    z_het_list_agilent = []
    z_het_list_nimblegen = []
    class_names_found = []           # pooled用
    class_names_agilent = []         # v16 fix: agilent用
    class_names_nimblegen = []       # v16 fix: nimblegen用

    for bclass in PRIMARY_BOUNDARY_CLASSES:
        bsv_col_all = bin_counts_a[bclass][svt][valid_si_all].astype(float)

        # === Pooled MNLogit ===
        X_df_pooled = pd.DataFrame({
            "const": np.ones(int(keep_all.sum())),
            "burden": bsv_col_all,
            "sex": sex_all,
            "log1p_base": log1p_base_all,
            "log1p_gene": log1p_gene_all,
            "platform": plat_all,
        })
        x_colnames_pooled = list(X_df_pooled.columns)

        res_pooled = _run_mnlogit_heterogeneity(outcome_all, X_df_pooled, x_colnames_pooled, "pooled")

        if res_pooled is not None:
            z_het_list_pooled.append(res_pooled["z_het"])
            class_names_found.append(bclass)

            log(f"  {bclass:20s}: bASD={res_pooled['beta_asd']:+.4f} bSZ={res_pooled['beta_sz']:+.4f} "
                f"z_het={res_pooled['z_het']:+.3f} P={res_pooled['p_het_twosided']:.2e} "
                f"{'ASD>SZ' if res_pooled['z_het'] > 0 else 'SZ>ASD'}")

            het_results.append({
                "cohort": "ArrayCGH_External",
                "comparison": "ASD_vs_SZ_heterogeneity",
                "pattern": "A", "sv_type": svt,
                "boundary_class": bclass,
                "model": "multinomial_logistic",
                "analysis": "pooled",
                "exposure_type": "heterogeneity_individual",
                "n_asd": n_asd, "n_scz": n_scz, "n_control": n_ctrl,
                "n_complete": n_asd + n_scz + n_ctrl,
                "beta_asd": res_pooled["beta_asd"], "beta_sz": res_pooled["beta_sz"],
                "z_het": res_pooled["z_het"], "se_diff": res_pooled["se_diff"],
                "OR_diff": float(np.exp(res_pooled["diff_beta"])),
                "CI_lower": float(np.exp(res_pooled["diff_beta"] - 1.96 * res_pooled["se_diff"])),
                "CI_upper": float(np.exp(res_pooled["diff_beta"] + 1.96 * res_pooled["se_diff"])),
                "P_het_onesided": float(scipy_stats.norm.sf(res_pooled["z_het"])),
                "P_het_twosided": res_pooled["p_het_twosided"],
                "direction_ASD_gt_SZ": bool(res_pooled["z_het"] > 0),
                "converged": res_pooled["converged"],
            })
        else:
            log(f"  {bclass:20s}: MNLogit FAILED (pooled)")

        # === Agilent-only MNLogit (sensitivity) ===
        ag_sel = agilent_within
        if int(np.sum(ag_sel & (outcome_all == 1))) >= 10 and int(np.sum(ag_sel & (outcome_all == 2))) >= 10:
            X_df_ag = pd.DataFrame({
                "const": np.ones(int(ag_sel.sum())),
                "burden": bsv_col_all[ag_sel],
                "sex": sex_all[ag_sel],
                "log1p_base": log1p_base_all[ag_sel],
                "log1p_gene": log1p_gene_all[ag_sel],
            })
            res_ag = _run_mnlogit_heterogeneity(outcome_all[ag_sel], X_df_ag, list(X_df_ag.columns), "agilent")
            if res_ag is not None:
                z_het_list_agilent.append(res_ag["z_het"])
                class_names_agilent.append(bclass)
                het_results.append({
                    "cohort": "ArrayCGH_External", "comparison": "ASD_vs_SZ_heterogeneity",
                    "pattern": "A", "sv_type": svt, "boundary_class": bclass,
                    "model": "multinomial_logistic", "analysis": "agilent_only",
                    "exposure_type": "heterogeneity_individual",
                    "n_asd": int(np.sum(ag_sel & (outcome_all == 1))),
                    "n_scz": int(np.sum(ag_sel & (outcome_all == 2))),
                    "n_control": int(np.sum(ag_sel & (outcome_all == 0))),
                    "n_complete": int(ag_sel.sum()),
                    "beta_asd": res_ag["beta_asd"], "beta_sz": res_ag["beta_sz"],
                    "z_het": res_ag["z_het"], "se_diff": res_ag["se_diff"],
                    "OR_diff": float(np.exp(res_ag["diff_beta"])),
                    "CI_lower": float(np.exp(res_ag["diff_beta"] - 1.96 * res_ag["se_diff"])),
                    "CI_upper": float(np.exp(res_ag["diff_beta"] + 1.96 * res_ag["se_diff"])),
                    "P_het_onesided": float(scipy_stats.norm.sf(res_ag["z_het"])),
                    "P_het_twosided": res_ag["p_het_twosided"],
                    "direction_ASD_gt_SZ": bool(res_ag["z_het"] > 0),
                    "converged": res_ag["converged"],
                })

        # === Nimblegen-only MNLogit (sensitivity) ===
        ni_sel = nimblegen_within
        if int(np.sum(ni_sel & (outcome_all == 1))) >= 10 and int(np.sum(ni_sel & (outcome_all == 2))) >= 10:
            X_df_ni = pd.DataFrame({
                "const": np.ones(int(ni_sel.sum())),
                "burden": bsv_col_all[ni_sel],
                "sex": sex_all[ni_sel],
                "log1p_base": log1p_base_all[ni_sel],
                "log1p_gene": log1p_gene_all[ni_sel],
            })
            res_ni = _run_mnlogit_heterogeneity(outcome_all[ni_sel], X_df_ni, list(X_df_ni.columns), "nimblegen")
            if res_ni is not None:
                z_het_list_nimblegen.append(res_ni["z_het"])
                class_names_nimblegen.append(bclass)
                het_results.append({
                    "cohort": "ArrayCGH_External", "comparison": "ASD_vs_SZ_heterogeneity",
                    "pattern": "A", "sv_type": svt, "boundary_class": bclass,
                    "model": "multinomial_logistic", "analysis": "nimblegen_only",
                    "exposure_type": "heterogeneity_individual",
                    "n_asd": int(np.sum(ni_sel & (outcome_all == 1))),
                    "n_scz": int(np.sum(ni_sel & (outcome_all == 2))),
                    "n_control": int(np.sum(ni_sel & (outcome_all == 0))),
                    "n_complete": int(ni_sel.sum()),
                    "beta_asd": res_ni["beta_asd"], "beta_sz": res_ni["beta_sz"],
                    "z_het": res_ni["z_het"], "se_diff": res_ni["se_diff"],
                    "OR_diff": float(np.exp(res_ni["diff_beta"])),
                    "CI_lower": float(np.exp(res_ni["diff_beta"] - 1.96 * res_ni["se_diff"])),
                    "CI_upper": float(np.exp(res_ni["diff_beta"] + 1.96 * res_ni["se_diff"])),
                    "P_het_onesided": float(scipy_stats.norm.sf(res_ni["z_het"])),
                    "P_het_twosided": res_ni["p_het_twosided"],
                    "direction_ASD_gt_SZ": bool(res_ni["z_het"] > 0),
                    "converged": res_ni["converged"],
                })

    # === Global tests ===
    # v16 fix: sample_indices引数を追加し、platform別subsetの負荷量相関を使う
    def _global_tests(z_list, analysis_label, class_names, sample_indices):
        """Global Stouffer + sign test.
        sample_indices: burden相関計算に使うサンプルのインデックス配列
                        (pooled=valid_si_all, agilent=valid_si_all[agilent_within], etc.)
        """
        if len(z_list) < 2:
            log(f"  {analysis_label}: Only {len(z_list)} classes — skipping global")
            return []
        z_arr = np.array(z_list)
        M = len(z_arr)
        n_pos = int(np.sum(z_arr > 0))
        sum_z = float(np.sum(z_arr))

        # Correlation-adjusted Stouffer (Brown's variance adjustment)
        var_brown = float(M)
        try:
            burden_matrix = np.column_stack([
                bin_counts_a[cn][svt][sample_indices].astype(float)
                for cn in class_names
            ])
            corr_mat = np.corrcoef(burden_matrix, rowvar=False)
            sum_rho = 0.0
            for ii in range(M):
                for jj in range(ii + 1, M):
                    rr = corr_mat[ii, jj]
                    if not np.isnan(rr):
                        sum_rho += rr
            var_brown = float(M) + 2.0 * sum_rho
        except Exception as e:
            log(f"  WARNING: burden correlation failed: {e}")

        z_brown = sum_z / np.sqrt(var_brown)
        p_brown = float(scipy_stats.norm.sf(z_brown))
        p_sign = float(scipy_stats.binomtest(n_pos, M, 0.5, alternative="greater").pvalue)

        log(f"  {analysis_label} (N_samples={len(sample_indices)}): "
            f"{M} classes, {n_pos}/{M} positive, "
            f"Stouffer Z={z_brown:.3f} P={p_brown:.2e}, sign P={p_sign:.4f}")

        results = []
        results.append({
            "cohort": "ArrayCGH_External", "comparison": "ASD_vs_SZ_heterogeneity",
            "pattern": "A", "sv_type": svt,
            "boundary_class": "GLOBAL_corr_stouffer",
            "model": "multinomial_logistic", "analysis": analysis_label,
            "exposure_type": "heterogeneity_global",
            "n_asd": np.nan, "n_scz": np.nan, "n_control": np.nan,
            "n_complete": M,
            "beta_asd": np.nan, "beta_sz": np.nan,
            "z_het": z_brown, "se_diff": np.sqrt(var_brown),
            "OR_diff": np.nan, "CI_lower": np.nan, "CI_upper": np.nan,
            "P_het_onesided": p_brown,
            "P_het_twosided": np.nan,
            "direction_ASD_gt_SZ": bool(z_brown > 0),
            "converged": True,
        })
        results.append({
            "cohort": "ArrayCGH_External", "comparison": "ASD_vs_SZ_heterogeneity",
            "pattern": "A", "sv_type": svt,
            "boundary_class": "GLOBAL_sign_test",
            "model": "multinomial_logistic", "analysis": analysis_label,
            "exposure_type": "heterogeneity_global",
            "n_asd": np.nan, "n_scz": np.nan, "n_control": np.nan,
            "n_complete": M,
            "beta_asd": np.nan, "beta_sz": np.nan,
            "z_het": float(n_pos), "se_diff": np.nan,
            "OR_diff": np.nan, "CI_lower": np.nan, "CI_upper": np.nan,
            "P_het_onesided": p_sign,
            "P_het_twosided": np.nan,
            "direction_ASD_gt_SZ": bool(n_pos > M / 2),
            "converged": True,
        })
        return results

    # Platform別サンプルインデックス (Brown補正用)
    si_pooled = valid_si_all
    si_agilent = valid_si_all[agilent_within]
    si_nimblegen = valid_si_all[nimblegen_within]

    log("\n  === GLOBAL TESTS ===")
    het_results.extend(_global_tests(z_het_list_pooled, "pooled", class_names_found, si_pooled))
    if z_het_list_agilent:
        het_results.extend(_global_tests(z_het_list_agilent, "agilent_only", class_names_agilent, si_agilent))
    if z_het_list_nimblegen:
        het_results.extend(_global_tests(z_het_list_nimblegen, "nimblegen_only", class_names_nimblegen, si_nimblegen))

    return het_results


# ============================================================================
# PHASE 7: Output
# ============================================================================
def run_phase7_output(asd_results, sz_results, het_results,
                      n_asd, n_ctrl_asd, n_scz, n_ctrl_scz,
                      n_overlap_asd, n_overlap_scz, gd_nahr_count):
    """Phase 7: 全結果出力"""
    log(f"\n{'=' * 70}")
    log("PHASE 7: Output")
    log("=" * 70)

    outdir = Path(_OUTDIR)
    outdir.mkdir(parents=True, exist_ok=True)

    # ASD results
    asd_df = pd.DataFrame(asd_results)
    asd_out_path = str(outdir / "tad_replication_asd_vs_cont_v22.tsv")
    asd_df.to_csv(asd_out_path, sep="\t", index=False, float_format="%.6g")
    log(f"  ASD vs CONT results: {asd_out_path}")

    # SZ results
    sz_df = pd.DataFrame(sz_results)
    sz_out_path = str(outdir / "tad_sz_vs_cont_v22.tsv")
    sz_df.to_csv(sz_out_path, sep="\t", index=False, float_format="%.6g")
    log(f"  SZ vs CONT results: {sz_out_path}")

    # Heterogeneity results
    het_df = pd.DataFrame(het_results)
    het_out_path = str(outdir / "tad_heterogeneity_asd_vs_sz_v22.tsv")
    het_df.to_csv(het_out_path, sep="\t", index=False, float_format="%.6g")
    log(f"  Heterogeneity results: {het_out_path}")

    # Summary
    log("\n=== ASD vs CONT — PRIMARY REPLICATION (Pattern A, DEL) ===")
    primary_asd = asd_df[(asd_df["analysis_type"] == "PRIMARY_replication")].copy()
    for _, r in primary_asd.sort_values("P_onesided").iterrows():
        sig = "**" if r["formal_replicated"] else "  "
        log(f"  {sig} {r['Boundary_class']:20s} OR={r['OR']:.3f} P_1s={r['P_onesided']:.2e}")

    log("\n=== SZ vs CONT — PRIMARY (Pattern A, DEL) ===")
    primary_sz = sz_df[(sz_df["analysis_type"] == "PRIMARY_SZ_burden") &
                       (sz_df["SV_type"] == "DEL") & (sz_df["Pattern"] == "A")].copy()
    for _, r in primary_sz.sort_values("P_twosided").iterrows():
        sig = "**" if np.isfinite(r["P_twosided"]) and r["P_twosided"] < 0.05 else "  "
        log(f"  {sig} {r['Boundary_class']:20s} OR={r['OR']:.3f} P_2s={r['P_twosided']:.2e}")

    log("\n=== HETEROGENEITY — POOLED (Pattern A, DEL) ===")
    pooled_het = het_df[(het_df["exposure_type"] == "heterogeneity_individual") &
                        (het_df["analysis"] == "pooled")]
    for _, r in pooled_het.iterrows():
        dir_str = "ASD>SZ" if r["direction_ASD_gt_SZ"] else "SZ>ASD"
        sig = "**" if r["P_het_twosided"] < 0.05 else "  "
        log(f"  {sig} {r['boundary_class']:20s} z_het={r['z_het']:+.3f} P={r['P_het_twosided']:.2e} {dir_str}")

    global_het = het_df[het_df["exposure_type"] == "heterogeneity_global"]
    for _, r in global_het.iterrows():
        p_col = "P_het_onesided" if "stouffer" in r["boundary_class"].lower() else "P_het_onesided"
        log(f"  [{r['analysis']:15s}] {r['boundary_class']:25s}: P={r[p_col]:.4f}")

    # Run config
    run_config = {
        "script": "tad_replication_arraycgh_v22.py",
        "version": "v22",
        "base_version": "v18 + overlap_srWGS列ベースの除外方式",
        "comparisons": ["ASD_vs_CONT", "SZ_vs_CONT", "ASD_vs_SZ_heterogeneity"],
        "model_pooled": "Logit_LRT_Bprime_covariates",
        "model_platform_meta": "IVW_fixed_effect",
        "model_heterogeneity": "MNLogit (3-level: pooled + platform-stratified sensitivity)",
        "common_cnv_filtering": "platform+DEL/DUP specific, one_direction_overlap >= 0.30",
        "static_all": "support==2 bins EXCLUDING all diffbound bins",
        "n_asd": int(n_asd), "n_ctrl_asd": int(n_ctrl_asd),
        "n_scz": int(n_scz), "n_ctrl_scz": int(n_ctrl_scz),
        "n_overlap_asd_removed": int(n_overlap_asd),
        "n_overlap_scz_removed": int(n_overlap_scz),
        "n_gd_nahr_flagged": int(gd_nahr_count),
        "timestamp": datetime.now().isoformat(),
    }
    config_path = str(outdir / "run_config.json")
    with open(config_path, "w") as f:
        json.dump(run_config, f, indent=2)
    log(f"  Run config: {config_path}")

    return asd_out_path, sz_out_path, het_out_path


# ============================================================================
# Utility: extract nonzero sparse (v15互換)
# ============================================================================
def extract_nonzero_sparse_long(adata, matrix_name: str, value_filter=None) -> pd.DataFrame:
    from scipy import sparse as sp
    X = adata.X.tocoo() if sp.issparse(adata.X) else sp.coo_matrix(adata.X)
    df = pd.DataFrame({"row_idx": X.row, "col_idx": X.col, "value": X.data})
    if value_filter is not None:
        df = df.loc[df["value"] == value_filter].copy()
    var_df = adata.var.reset_index(drop=True).copy()
    required_var_cols = ["chrom", "start", "end"]
    missing = [c for c in required_var_cols if c not in var_df.columns]
    if missing:
        raise ValueError(f"{matrix_name}: var に必要列がありません: {missing}")
    df["chrom"] = var_df.loc[df["col_idx"], "chrom"].to_numpy()
    df["start0"] = var_df.loc[df["col_idx"], "start"].to_numpy()
    df["end"] = var_df.loc[df["col_idx"], "end"].to_numpy()
    df["bin_id"] = (
        df["chrom"].astype(str) + ":" + df["start0"].astype(str) + "-" + df["end"].astype(str)
    )
    return df[["chrom", "start0", "end", "bin_id"]].reset_index(drop=True)


# ============================================================================
# MAIN
# ============================================================================
def main():
    total_start = time.time()
    log("=" * 70)
    log("tad_replication_arraycgh_v22.py")
    log("  Array CGH: ASD replication + SZ burden + ASD-SZ heterogeneity")
    log("  Phase 0-5: v15互換 (common CNV: platform+type, static_all: diffbound除外)")
    log("  Phase 2: overlap_srWGS列ベースの除外方式 (v19変更)")
    log("  Phase 3: segdup/exclusion BED filter追加 (v18新規)")
    log("  Phase 5B: per-sample bin overlap + covariates出力 (v17新規)")
    log("  Phase 6: ASD vs CONT (PRIMARY replication)")
    log("  Phase 6B: SZ vs CONT (external SZ burden)")
    log("  Phase 6C: MNLogit heterogeneity + platform sensitivity")
    log("  Phase 7: output")
    log("=" * 70)

    run_phase0_validate()
    cnv_raw, platform_map = run_phase1_load_cnv()
    sample_info, n_overlap_asd, n_overlap_scz = run_phase2_sample_qc()
    sv, gd_nahr = run_phase3_filter_cnv(cnv_raw, sample_info)
    bin_index = run_phase4_load_bins()
    (all_samples, sample_to_idx, bin_counts_a, bin_counts_c,
     sample_total_base, sample_total_gene,
     sample_total_base_c, sample_total_gene_c,
     disrupted_bins_a, disrupted_bins_c) = run_phase5_bin_counts(sv, sample_info, bin_index)

    # v17新規: Phase 5B (event_id付きbin overlap + sample covariates)
    overlap_path, covar_path = run_phase5b_save_bin_overlap(
        all_samples, sample_to_idx, sample_info,
        sv, bin_index,
        disrupted_bins_a, disrupted_bins_c,
        sample_total_base, sample_total_gene,
        sample_total_base_c, sample_total_gene_c)

    asd_results, n_asd, n_ctrl_asd = run_phase6_asd_regression(
        all_samples, sample_info, bin_counts_a, bin_counts_c,
        sample_total_base, sample_total_gene, sample_total_base_c, sample_total_gene_c)

    sz_results, n_scz, n_ctrl_scz = run_phase6b_sz_regression(
        all_samples, sample_info, bin_counts_a, bin_counts_c,
        sample_total_base, sample_total_gene, sample_total_base_c, sample_total_gene_c)

    het_results = run_phase6c_heterogeneity(
        all_samples, sample_info, bin_counts_a, sample_total_base, sample_total_gene)

    asd_path, sz_path, het_path = run_phase7_output(
        asd_results, sz_results, het_results,
        n_asd, n_ctrl_asd, n_scz, n_ctrl_scz,
        n_overlap_asd, n_overlap_scz, len(gd_nahr))

    total_elapsed = time.time() - total_start
    log(f"\nPipeline complete. Total elapsed: {total_elapsed:.1f}s ({total_elapsed / 60:.1f}min)")
    log(f"Results: ASD={asd_path}, SZ={sz_path}, Het={het_path}")
    log(f"v17 new: overlap={overlap_path}, covariates={covar_path}")


if __name__ == "__main__":
    main()
