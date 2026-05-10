#!/usr/bin/env python3
# 42_compute_wgs_exon_exclusion_v2.py
#
# 処理内容:
#  - WGS discovery cohort (tad04212026 パイプライン) の exon-exclusion sensitivity
#    を計算する。heffel run_exonfree_aggregate_sensitivity_v4.py の tad04212026 移植
#    に Script 11 v4 スタイルの Firth sensitivity + BH/Holm 補正を追加した v2。
#  - Input: tad04212026 v9/v2 の処理済みファイル
#      01_heffel_boundary_master/output_v9/heffel_boundary_master_v9.tsv.gz
#      02_bin_l2_annotation/output_v2/bin_l2_annotation_v2.tsv.gz
#          (10 L2 class membership + bin 座標; hpc_astro は元から除外)
#      04_wgs_sv_boundary_overlap/output_v9/sample_boundary_event_overlap_v9.tsv.gz
#          (sample × event × bin)
#      05_wgs_sample_burden/output_v2/sample_burden_L2_and_specificity_v2.tsv
#          (共変量 + サンプル単位 burden; Diagnosis/Sex_numeric/PC1-10/log1p_* 等)
#      Annotation: /lustre12/home/kushima-pg/annotationInfo/gencode.v46.annotation.gtf.gz
#  - Steps:
#      1. GENCODE v46 protein-coding exon を merge
#      2. bin_l2_annotation の 25-kb bin 座標と exon を交差し、bin ごとに
#         exon_free / exon_overlap を判定。10 L2 class のいずれかに属する bin
#         (any_l2_bins) のみを解析対象とする。
#      3. sample × event × bin overlap を読み、any_l2_bins に絞って:
#           A1/A2 軸 (bin レベル): exonfree / exonoverlap / all
#           A3 軸 (SV レベル)   : SV 全体が exon 非重複かを sv_chr/sv_start0/sv_end
#                                 から判定 → exonfreeCNV
#           A3-strict            : exonfreeCNV かつ hit bin が exon-free
#         サンプル × SV type (DEL/DUP) × 5 stratification の unique bin 数を集計
#      4. burden TSV に 5 exposure × 2 SVT 列を追加し、以下の回帰を実施:
#           Main: statsmodels.Logit (B' 共変量、両側 P、Wald 95% CI)
#           Sensitivity: 自前 Firth (Heinze & Schemper 2002, Newton-Raphson)
#                         QR-with-pivoting で rank-deficient 列を自動 drop
#                         protected_idx = (intercept, exposure[, exposure2])
#      5. 解析グリッド: 2 comparisons × 2 SVT × 4 analyses
#           Comparisons   : ASD_vs_Healthy, SZ_vs_Healthy
#           SVT           : DEL, DUP
#           Analyses      : A1_separate (3 bin_types), A2_joint (2 bin_types),
#                           A3_exonfreeCNV (1), A3_strict (1)
#           → 1 comparison × 1 SVT あたり 7 結果行 × 2 methods (Logit + Firth)
#      6. Multiple testing: ASD_vs_Healthy × DEL discovery-positive 4 行
#         (A1 exonfree, A2 exonfree_in_joint, A3 exonfreeCNV, A3 strict) に対して
#         BH-FDR + Holm を適用 (main script は Logit P と Firth P 両方)
#      7. 出力 TSV:
#           wgs_exon_exclusion_v2.tsv (全結果行)
#           wgs_exon_exclusion_summary_v2.tsv (サマリ: ASD DEL 4 解析の multiple
#                                              testing 込み)
#  - 出力先: /lustre12/home/kushima-pg/tad04212026/12_exon_exclusion_wgs/output_v2/
#  - 処理時間を先頭と末尾で記録
#
# v1 (heffel run_exonfree_aggregate_sensitivity_v4.py) からの変更点:
#  1. PATHS: heffel → tad04212026 (v9 input + v2 output)
#  2. L2 class 判定: heffel v16 boundary master 内 10 列 →
#     tad04212026 v2 bin_l2_annotation との merge で対応
#  3. [追加] 自前 Firth penalized logistic (sensitivity) + QR rank reduction
#            (Script 11 v4 から porting)
#  4. [追加] BH-FDR + Holm multiple testing 補正 (Logit P, Firth P それぞれ)
#  5. [追加] 結果 DataFrame に fit_status / carrier_case / carrier_ctrl / note /
#            fit_status_firth 列を追加
#  6. log 形式を [YYYY-MM-DD HH:MM:SS] Script 11 v4 スタイルに統一
#  7. OUT 出力ファイル名: exonfree_aggregate_sensitivity_patternA_v4_results.tsv
#     → wgs_exon_exclusion_v2.tsv + wgs_exon_exclusion_summary_v2.tsv
#  8. 主要解析ロジック (A1/A2/A3/A3-strict) は heffel v4 と同一 (結果の一致を
#     保証するため、GENCODE intersection アルゴリズムと集計コードは変えない)

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=42_compute_wgs_exon_exclusion_v2_%j.out
#SBATCH --error=42_compute_wgs_exon_exclusion_v2_%j.err

from __future__ import annotations

import bisect
import gzip
import time
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd
import scipy.linalg as sp_linalg
import scipy.stats as sp_stats

import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


# =========================================================
# PATHS (tad04212026 パイプライン)
# =========================================================
BASE_TAD = Path("/lustre12/home/kushima-pg/tad04212026")

GENCODE_GTF = Path(
    "/lustre12/home/kushima-pg/annotationInfo/gencode.v46.annotation.gtf.gz"
)
BOUNDARY_MASTER = (
    BASE_TAD / "01_heffel_boundary_master" / "output_v9"
    / "heffel_boundary_master_v9.tsv.gz"
)
BIN_L2_ANNOT = (
    BASE_TAD / "02_bin_l2_annotation" / "output_v2"
    / "bin_l2_annotation_v2.tsv.gz"
)
EVENT_OVERLAP = (
    BASE_TAD / "04_wgs_sv_boundary_overlap" / "output_v9"
    / "sample_boundary_event_overlap_v9.tsv.gz"
)
BURDEN_TABLE = (
    BASE_TAD / "05_wgs_sample_burden" / "output_v2"
    / "sample_burden_L2_and_specificity_v2.tsv"
)

OUT_DIR = BASE_TAD / "12_exon_exclusion_wgs" / "output_v2"
OUT_TSV = OUT_DIR / "wgs_exon_exclusion_v2.tsv"
OUT_SUMMARY = OUT_DIR / "wgs_exon_exclusion_summary_v2.tsv"


# =========================================================
# CONSTANTS
# =========================================================
L2_MEMBERSHIP_COLS = [
    "membership_HPC_Exc-CA",
    "membership_HPC_Exc-DG",
    "membership_HPC_Exc-ENT",
    "membership_HPC_Inh-CGE",
    "membership_HPC_Inh-MGE",
    "membership_PFC_Astro",
    "membership_PFC_Exc-DL",
    "membership_PFC_Exc-UL",
    "membership_PFC_Inh-CGE",
    "membership_PFC_Inh-MGE",
]

SVT_LIST = ["DEL", "DUP"]
COMPARISONS = [
    ("ASD_vs_Healthy", "ASD", "Healthy"),
    ("SZ_vs_Healthy", "SZ", "Healthy"),
]

# Multiple-testing stratum (discovery-positive, ASD × DEL)
DISCOVERY_POSITIVE_KEYS = [
    ("A1_separate", "exonfree"),
    ("A2_joint", "exonfree_in_joint"),
    ("A3_exonfreeCNV", "exonfreeCNV"),
    ("A3_strict_exonfreeCNV_exonfreeBin", "exonfreeCNV_exonfreeBin"),
]

PC_COLS = [f"PC{i}" for i in range(1, 11)]
MIN_CELL_COUNT = 5

# =========================================================
# Logging
# =========================================================
_T0 = time.time()


def log(msg: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    elapsed = time.time() - _T0
    print(f"[{stamp}] [{elapsed:7.1f}s] {msg}", flush=True)


# =========================================================
# STEP 1: Merge protein-coding exon intervals from GENCODE
# =========================================================
def build_merged_exon_bed(gtf_path: Path) -> Dict[str, List[Tuple[int, int]]]:
    log(f"Reading GENCODE GTF: {gtf_path}")
    exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    n_exon = 0
    with gzip.open(str(gtf_path), "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "exon":
                continue
            if 'gene_type "protein_coding"' not in fields[8]:
                continue
            chrom = fields[0]
            s = int(fields[3]) - 1  # 0-based
            e = int(fields[4])
            exons[chrom].append((s, e))
            n_exon += 1
    log(f"  Raw protein-coding exon intervals: {n_exon}")

    merged: Dict[str, List[Tuple[int, int]]] = {}
    total_merged = 0
    for chrom in sorted(exons.keys()):
        intervals = sorted(exons[chrom])
        m = [intervals[0]]
        for s, e in intervals[1:]:
            if s <= m[-1][1]:
                m[-1] = (m[-1][0], max(m[-1][1], e))
            else:
                m.append((s, e))
        merged[chrom] = m
        total_merged += len(m)
    log(f"  Merged exon intervals: {total_merged}")
    return merged


def overlaps_any_exon(
    chrom: str,
    start: int,
    end: int,
    merged_exons: Dict[str, List[Tuple[int, int]]],
    exon_starts: Dict[str, List[int]],
) -> bool:
    if chrom not in exon_starts:
        return False
    starts = exon_starts[chrom]
    intervals = merged_exons[chrom]
    idx = bisect.bisect_right(starts, start) - 1
    for i in range(max(0, idx), len(intervals)):
        es, ee = intervals[i]
        if es >= end:
            break
        if ee > start:
            return True
    return False


# =========================================================
# STEP 2: Load bin_l2_annotation -> exon-free bins + 10 L2 pool
# =========================================================
def load_bin_l2_annotation(
    path: Path,
    merged_exons: Dict[str, List[Tuple[int, int]]],
    exon_starts: Dict[str, List[int]],
) -> Tuple[Set[str], Set[str], Set[str]]:
    """Return (exon_free_bins, exon_overlap_bins, any_l2_bins)."""
    log(f"Reading bin L2 annotation: {path}")
    exon_free_bins: Set[str] = set()
    exon_overlap_bins: Set[str] = set()
    any_l2_bins: Set[str] = set()

    with gzip.open(str(path), "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        req_cols = ["bin_id", "chrom", "start0", "end"] + L2_MEMBERSHIP_COLS
        missing = [c for c in req_cols if c not in col]
        if missing:
            raise RuntimeError(
                f"bin_l2_annotation missing columns: {missing}"
            )
        idx_bin = col["bin_id"]
        idx_chr = col["chrom"]
        idx_s = col["start0"]
        idx_e = col["end"]
        idx_mem = [col[c] for c in L2_MEMBERSHIP_COLS]

        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < max(idx_mem) + 1:
                continue
            bin_id = fields[idx_bin]
            # L2 membership (union across 10 columns)
            is_l2 = False
            for mi in idx_mem:
                try:
                    v = float(fields[mi]) if fields[mi] else 0.0
                except ValueError:
                    v = 0.0
                if v > 0:
                    is_l2 = True
                    break
            if is_l2:
                any_l2_bins.add(bin_id)

            chrom = fields[idx_chr]
            try:
                bs = int(fields[idx_s])
                be = int(fields[idx_e])
            except ValueError:
                continue
            if overlaps_any_exon(chrom, bs, be, merged_exons, exon_starts):
                exon_overlap_bins.add(bin_id)
            else:
                exon_free_bins.add(bin_id)

    total = len(exon_free_bins) + len(exon_overlap_bins)
    log(
        f"  Total annotated bins: {total}  "
        f"exon-overlap={len(exon_overlap_bins)} "
        f"({100*len(exon_overlap_bins)/max(total,1):.1f}%)  "
        f"exon-free={len(exon_free_bins)} "
        f"({100*len(exon_free_bins)/max(total,1):.1f}%)"
    )
    log(f"  Any-L2 pool (10 classes, hpc_astro excluded): {len(any_l2_bins)}")
    return exon_free_bins, exon_overlap_bins, any_l2_bins


# =========================================================
# STEP 3: Count disrupted bins per sample × sv_type
# =========================================================
def count_disrupted_bins(
    event_path: Path,
    exon_free_bins: Set[str],
    any_l2_bins: Set[str],
    merged_exons: Dict[str, List[Tuple[int, int]]],
    exon_starts: Dict[str, List[int]],
) -> Tuple[dict, dict, dict, dict, dict]:
    """Return 5 defaultdicts:
       sample_bins_exonfree[sid][svt] -> set of bin_id
       sample_bins_exonoverlap[sid][svt]
       sample_bins_all[sid][svt]
       sample_bins_exonfreeCNV[sid][svt]
       sample_bins_exonfreeCNV_exonfreeBin[sid][svt]
    """
    log(f"Reading event overlap: {event_path}")
    s_ef: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    s_eo: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    s_all: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    s_efcnv: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))
    s_efcnv_efbin: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: defaultdict(set))

    sv_exonfree_cache: Dict[str, bool] = {}

    n_rows = 0
    n_kept = 0
    n_ef_bin = 0
    n_a3_rows = 0
    n_a3_strict_rows = 0
    a3_hit_ef = 0
    a3_hit_eo = 0

    with gzip.open(str(event_path), "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        req_cols = [
            "sample_id",
            "event_id",
            "sv_type_norm",
            "sv_chr",
            "sv_start0",
            "sv_end",
            "bin_id",
        ]
        missing = [c for c in req_cols if c not in col]
        if missing:
            raise RuntimeError(
                f"event overlap missing columns: {missing}"
            )
        idx_sid = col["sample_id"]
        idx_eid = col["event_id"]
        idx_svt = col["sv_type_norm"]
        idx_svch = col["sv_chr"]
        idx_svs = col["sv_start0"]
        idx_sve = col["sv_end"]
        idx_bin = col["bin_id"]

        for line in f:
            n_rows += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) < max(idx_bin, idx_sve) + 1:
                continue
            bin_id = fields[idx_bin]
            if bin_id not in any_l2_bins:
                continue
            n_kept += 1

            sid = fields[idx_sid]
            svt = fields[idx_svt]
            if svt not in SVT_LIST:
                continue
            event_id = fields[idx_eid]

            bin_is_ef = bin_id in exon_free_bins

            # A1 / A2: bin-level
            s_all[sid][svt].add(bin_id)
            if bin_is_ef:
                s_ef[sid][svt].add(bin_id)
                n_ef_bin += 1
            else:
                s_eo[sid][svt].add(bin_id)

            # A3: SV-level exon-free judgment (cached)
            if event_id not in sv_exonfree_cache:
                try:
                    ch = fields[idx_svch]
                    ss = int(fields[idx_svs])
                    ee = int(fields[idx_sve])
                    sv_exonfree_cache[event_id] = not overlaps_any_exon(
                        ch, ss, ee, merged_exons, exon_starts
                    )
                except Exception:
                    sv_exonfree_cache[event_id] = False

            if sv_exonfree_cache[event_id]:
                s_efcnv[sid][svt].add(bin_id)
                n_a3_rows += 1
                if bin_is_ef:
                    s_efcnv_efbin[sid][svt].add(bin_id)
                    n_a3_strict_rows += 1
                    a3_hit_ef += 1
                else:
                    a3_hit_eo += 1

    n_sv = len(sv_exonfree_cache)
    n_ef_sv = sum(1 for v in sv_exonfree_cache.values() if v)
    log(f"  Raw event-bin rows: {n_rows}")
    log(f"  Kept (any-L2 bin, DEL/DUP): {n_kept}")
    log(f"  Exon-free bin rows: {n_ef_bin}")
    log(f"  Unique SVs on L2 bins: {n_sv}")
    log(
        f"  Exon-free SVs (axis B): {n_ef_sv} "
        f"({100*n_ef_sv/max(n_sv,1):.1f}%)"
    )
    log(f"  A3 rows (exon-free CNV × any L2 bin): {n_a3_rows}")
    log(f"    hit bin is exon-free: {a3_hit_ef}")
    log(f"    hit bin is exon-overlap: {a3_hit_eo}")
    log(f"  A3-strict rows: {n_a3_strict_rows}")

    return s_ef, s_eo, s_all, s_efcnv, s_efcnv_efbin


# =========================================================
# STEP 4: QR rank reduction + Firth (Script 11 v4 から porting)
# =========================================================
def _reduce_rank_via_qr(
    X: np.ndarray,
    protected_idx: Tuple[int, ...] = (0, 1),
    tol_rel: float = 1e-10,
) -> Tuple[np.ndarray, np.ndarray, bool]:
    """QR with pivoting で rank-deficient 列を drop。
    protected_idx は必ず保持する列 (intercept, exposure 等)。"""
    n, p = X.shape
    col_norms = np.linalg.norm(X, axis=0)
    col_norms_safe = np.where(col_norms > 0, col_norms, 1.0)
    X_scaled = X / col_norms_safe
    try:
        _, R, pivot = sp_linalg.qr(X_scaled, pivoting=True, mode="economic")
    except Exception:
        return X, np.arange(p), False
    diag_R = np.abs(np.diag(R))
    if diag_R.size == 0:
        keep = np.array(sorted(set(protected_idx)))
        return X[:, keep], keep, True
    tol = max(n, p) * np.finfo(float).eps * diag_R[0]
    tol = max(tol, tol_rel * diag_R[0])
    rank = int(np.sum(diag_R > tol))
    keep_set = set(pivot[:rank].tolist())
    for pi in protected_idx:
        if 0 <= pi < p:
            keep_set.add(pi)
    keep_sorted = sorted(keep_set)
    if len(keep_sorted) > rank:
        X_sub = X[:, keep_sorted]
        rank_sub = int(np.linalg.matrix_rank(X_sub, tol=tol))
        non_protected = [c for c in keep_sorted if c not in protected_idx]
        for c in reversed(non_protected):
            if rank_sub == len(keep_sorted):
                break
            keep_sorted.remove(c)
            X_sub = X[:, keep_sorted]
            rank_sub = int(np.linalg.matrix_rank(X_sub, tol=tol))
    keep_idx = np.array(keep_sorted, dtype=int)
    reduced = len(keep_idx) < p
    return X[:, keep_idx], keep_idx, reduced


def firth_logit_fit(
    X: np.ndarray,
    y: np.ndarray,
    maxiter: int = 100,
    tol: float = 1e-6,
    protected_idx: Tuple[int, ...] = (0, 1),
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, bool, str]:
    """Firth penalized logistic (Heinze & Schemper 2002).
    QR-with-pivoting で rank deficient 列を事前 drop (protected_idx は保持)。
    full-size beta/se/p を返却 (drop 列は NaN)。"""
    n, p_full = X.shape
    try:
        X_red, keep_idx, reduced = _reduce_rank_via_qr(
            X, protected_idx=protected_idx
        )
        p = X_red.shape[1]
        for pi in protected_idx:
            if pi not in keep_idx:
                return (
                    np.full(p_full, np.nan),
                    np.full(p_full, np.nan),
                    np.full(p_full, np.nan),
                    False,
                    f"protected_col_dropped:{pi}",
                )

        beta = np.zeros(p)
        converged = False
        status = "not_converged"

        for _ in range(maxiter):
            eta = np.clip(X_red @ beta, -30.0, 30.0)
            mu = 1.0 / (1.0 + np.exp(-eta))
            mu = np.clip(mu, 1e-10, 1.0 - 1e-10)
            W = mu * (1.0 - mu)
            XW = X_red * W[:, None]
            XWX = X_red.T @ XW
            try:
                XWX_inv = np.linalg.inv(XWX)
            except np.linalg.LinAlgError:
                return (
                    np.full(p_full, np.nan),
                    np.full(p_full, np.nan),
                    np.full(p_full, np.nan),
                    False,
                    "singular_hessian_after_qr",
                )
            Wsqrt = np.sqrt(W)
            XWhalf = X_red * Wsqrt[:, None]
            H_diag = np.einsum("ij,jk,ik->i", XWhalf, XWX_inv, XWhalf)
            U = X_red.T @ (y - mu + (H_diag / 2.0) * (1.0 - 2.0 * mu))
            delta = XWX_inv @ U
            beta_new = beta + delta
            if np.max(np.abs(delta)) < tol:
                beta = beta_new
                converged = True
                status = "ok"
                break
            beta = beta_new

        # SE from unpenalized Fisher info at convergence
        eta = np.clip(X_red @ beta, -30.0, 30.0)
        mu = 1.0 / (1.0 + np.exp(-eta))
        mu = np.clip(mu, 1e-10, 1.0 - 1e-10)
        W = mu * (1.0 - mu)
        XWX = X_red.T @ (X_red * W[:, None])
        try:
            cov = np.linalg.inv(XWX)
        except np.linalg.LinAlgError:
            return (
                np.full(p_full, np.nan),
                np.full(p_full, np.nan),
                np.full(p_full, np.nan),
                False,
                "singular_cov_at_end",
            )
        var = np.diag(cov)
        var = np.where(var > 0, var, np.nan)
        se = np.sqrt(var)
        z = np.where(se > 0, beta / se, np.nan)
        p_two = 2.0 * sp_stats.norm.sf(np.abs(z))

        beta_full = np.full(p_full, np.nan)
        se_full = np.full(p_full, np.nan)
        p_full_arr = np.full(p_full, np.nan)
        beta_full[keep_idx] = beta
        se_full[keep_idx] = se
        p_full_arr[keep_idx] = p_two
        if reduced and status == "ok":
            status = f"ok_rank_reduced:{p_full - p}"
        return beta_full, se_full, p_full_arr, converged, status
    except Exception as e:
        return (
            np.full(p_full, np.nan),
            np.full(p_full, np.nan),
            np.full(p_full, np.nan),
            False,
            f"error:{str(e)[:80]}",
        )


# =========================================================
# STEP 5: Standard Logit fit (heffel B' 互換) + carrier 集計
# =========================================================
def _empty_fit_result(status: str, n_case: int, n_ctrl: int,
                       carrier_case: int, carrier_ctrl: int) -> dict:
    return {
        "n_case": int(n_case),
        "n_ctrl": int(n_ctrl),
        "carrier_case": int(carrier_case),
        "carrier_ctrl": int(carrier_ctrl),
        "beta": np.nan,
        "se": np.nan,
        "or": np.nan,
        "or_lo95": np.nan,
        "or_hi95": np.nan,
        "p_two": np.nan,
        "fit_status": status,
    }


def fit_logit_wgs(
    df: pd.DataFrame,
    y_col: str,
    x_cols: List[str],
    covar_cols: List[str],
    min_cell: int = MIN_CELL_COUNT,
) -> List[dict]:
    """statsmodels.Logit 回帰。x_cols は 1 or 2 列 (A2 joint 用)。
    各 exposure について結果 dict を返す。"""
    sub = df[covar_cols + list(x_cols) + [y_col]].copy()
    sub = sub.replace([np.inf, -np.inf], np.nan).dropna()
    n_case = int((sub[y_col] == 1).sum())
    n_ctrl = int((sub[y_col] == 0).sum())
    results = []

    # carrier counts per exposure
    for xc in x_cols:
        cc = int(((sub[y_col] == 1) & (sub[xc] >= 1)).sum())
        cn = int(((sub[y_col] == 0) & (sub[xc] >= 1)).sum())
        results.append(
            {"_xc": xc, "_cc": cc, "_cn": cn}
        )

    # Pre-checks
    if sub[y_col].nunique() < 2:
        return [
            _empty_fit_result("no_outcome_variance", n_case, n_ctrl, r["_cc"], r["_cn"])
            | {"exposure": r["_xc"]}
            for r in results
        ]
    for xc in x_cols:
        if sub[xc].nunique() < 2:
            return [
                _empty_fit_result(
                    f"no_variance_{xc}", n_case, n_ctrl, r["_cc"], r["_cn"]
                )
                | {"exposure": r["_xc"]}
                for r in results
            ]
    total_carriers = sum(r["_cc"] + r["_cn"] for r in results)
    if total_carriers < min_cell:
        return [
            _empty_fit_result(
                "insufficient_carriers", n_case, n_ctrl, r["_cc"], r["_cn"]
            )
            | {"exposure": r["_xc"]}
            for r in results
        ]

    X_df = sub[list(x_cols) + covar_cols].astype(float)
    X_df = sm.add_constant(X_df, has_constant="add")
    y = sub[y_col].astype(int).to_numpy()

    last_err = ""
    fit = None
    for method in ["newton", "lbfgs", "bfgs"]:
        try:
            fit_tmp = sm.Logit(y, X_df).fit(disp=0, maxiter=200, method=method)
            if all(np.isfinite(fit_tmp.params.get(xc, np.nan)) for xc in x_cols):
                fit = fit_tmp
                converged = fit_tmp.mle_retvals.get("converged", False)
                break
        except Exception as e:
            last_err = str(e)[:120]

    out = []
    if fit is None:
        return [
            _empty_fit_result(
                f"glm_error:{last_err}", n_case, n_ctrl, r["_cc"], r["_cn"]
            )
            | {"exposure": r["_xc"]}
            for r in results
        ]

    for r in results:
        xc = r["_xc"]
        try:
            beta = float(fit.params[xc])
            se = float(fit.bse[xc])
            p = float(fit.pvalues[xc])
            out.append(
                {
                    "exposure": xc,
                    "n_case": n_case,
                    "n_ctrl": n_ctrl,
                    "carrier_case": r["_cc"],
                    "carrier_ctrl": r["_cn"],
                    "beta": beta,
                    "se": se,
                    "or": float(np.exp(beta)),
                    "or_lo95": float(np.exp(beta - 1.959964 * se)),
                    "or_hi95": float(np.exp(beta + 1.959964 * se)),
                    "p_two": p,
                    "fit_status": (
                        "ok"
                        if fit.mle_retvals.get("converged", False)
                        else "not_converged"
                    ),
                }
            )
        except Exception as e:
            out.append(
                _empty_fit_result(
                    f"extract_error:{str(e)[:80]}",
                    n_case,
                    n_ctrl,
                    r["_cc"],
                    r["_cn"],
                )
                | {"exposure": xc}
            )
    return out


def fit_firth_wgs(
    df: pd.DataFrame,
    y_col: str,
    x_cols: List[str],
    covar_cols: List[str],
    min_cell: int = MIN_CELL_COUNT,
) -> List[dict]:
    """自前 Firth penalized logistic (QR rank reduction つき)。
    x_cols は 1 or 2 列。protected_idx は (0, 1) または (0, 1, 2)。"""
    sub = df[covar_cols + list(x_cols) + [y_col]].copy()
    sub = sub.replace([np.inf, -np.inf], np.nan).dropna()
    n_case = int((sub[y_col] == 1).sum())
    n_ctrl = int((sub[y_col] == 0).sum())

    carrier_info = []
    for xc in x_cols:
        cc = int(((sub[y_col] == 1) & (sub[xc] >= 1)).sum())
        cn = int(((sub[y_col] == 0) & (sub[xc] >= 1)).sum())
        carrier_info.append({"_xc": xc, "_cc": cc, "_cn": cn})

    if sub[y_col].nunique() < 2:
        return [
            _empty_fit_result("no_outcome_variance", n_case, n_ctrl, r["_cc"], r["_cn"])
            | {"exposure": r["_xc"]}
            for r in carrier_info
        ]
    for xc in x_cols:
        if sub[xc].nunique() < 2:
            return [
                _empty_fit_result(
                    f"no_variance_{xc}", n_case, n_ctrl, r["_cc"], r["_cn"]
                )
                | {"exposure": r["_xc"]}
                for r in carrier_info
            ]
    total_carriers = sum(r["_cc"] + r["_cn"] for r in carrier_info)
    if total_carriers < min_cell:
        return [
            _empty_fit_result(
                "insufficient_carriers", n_case, n_ctrl, r["_cc"], r["_cn"]
            )
            | {"exposure": r["_xc"]}
            for r in carrier_info
        ]

    # Design matrix layout: [intercept, x1(, x2), covariates...]
    X_list = [np.ones(len(sub))]
    X_list.extend(sub[xc].astype(float).to_numpy() for xc in x_cols)
    X_list.extend(sub[c].astype(float).to_numpy() for c in covar_cols)
    X = np.column_stack(X_list)
    y = sub[y_col].astype(int).to_numpy()

    # protected: intercept + all x_cols
    protected = tuple(range(1 + len(x_cols)))  # (0,1) or (0,1,2)
    beta_arr, se_arr, p_arr, converged, status = firth_logit_fit(
        X, y, maxiter=100, tol=1e-6, protected_idx=(0,) + protected[1:]
    )
    out = []
    for i, xc in enumerate(x_cols):
        col_idx = 1 + i
        beta = beta_arr[col_idx]
        se = se_arr[col_idx]
        p = p_arr[col_idx]
        r = carrier_info[i]
        if np.isfinite(beta) and np.isfinite(se):
            out.append(
                {
                    "exposure": xc,
                    "n_case": n_case,
                    "n_ctrl": n_ctrl,
                    "carrier_case": r["_cc"],
                    "carrier_ctrl": r["_cn"],
                    "beta": float(beta),
                    "se": float(se),
                    "or": float(np.exp(beta)),
                    "or_lo95": float(np.exp(beta - 1.959964 * se)),
                    "or_hi95": float(np.exp(beta + 1.959964 * se)),
                    "p_two": float(p),
                    "fit_status": status,
                }
            )
        else:
            out.append(
                _empty_fit_result(
                    status, n_case, n_ctrl, r["_cc"], r["_cn"]
                )
                | {"exposure": xc}
            )
    return out


# =========================================================
# STEP 6: Main orchestration
# =========================================================
def main() -> None:
    log("=" * 72)
    log("Script 12 v2 START: WGS exon-exclusion sensitivity (tad04212026)")
    log("  42_compute_wgs_exon_exclusion_v2.py")
    log("=" * 72)
    log("Analyses: [A1] Separate, [A2] Joint, [A3] exon-free CNV, [A3-strict]")
    log("Comparisons: ASD_vs_Healthy, SZ_vs_Healthy; SVT: DEL, DUP")
    log("Methods: statsmodels.Logit (main) + Firth (sensitivity, QR rank reduction)")
    log("Multiple testing: BH + Holm on ASD×DEL discovery-positive 4 tests")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    log(f"Output dir: {OUT_DIR}")

    # ----- Step 1: GENCODE exon index -----
    merged_exons = build_merged_exon_bed(GENCODE_GTF)
    exon_starts = {c: [iv[0] for iv in merged_exons[c]] for c in merged_exons}

    # ----- Step 2: bin_l2_annotation -> exon_free / exon_overlap / any_l2 -----
    exon_free_bins, exon_overlap_bins, any_l2_bins = load_bin_l2_annotation(
        BIN_L2_ANNOT, merged_exons, exon_starts
    )

    # ----- Step 3: event-bin overlap -> sample counts (5 stratifications) -----
    s_ef, s_eo, s_all, s_efcnv, s_efcnv_efbin = count_disrupted_bins(
        EVENT_OVERLAP, exon_free_bins, any_l2_bins, merged_exons, exon_starts
    )

    # ----- Step 4: burden TSV (covariates + sample_id) -----
    log(f"Reading burden TSV: {BURDEN_TABLE}")
    burden = pd.read_csv(str(BURDEN_TABLE), sep="\t")
    log(f"  Raw burden rows: {len(burden)}")

    essential_cols = (
        ["Sex_numeric", "log1p_total_del_bases", "log1p_total_dup_bases",
         "log1p_total_gene_DEL", "log1p_total_gene_DUP"]
        + PC_COLS
    )
    missing_cov = [c for c in essential_cols + ["Diagnosis"] if c not in burden.columns]
    if missing_cov:
        raise RuntimeError(f"Burden TSV missing columns: {missing_cov}")

    n_before = len(burden)
    burden = burden.dropna(subset=essential_cols).reset_index(drop=True)
    log(f"  After dropna on covariates: {len(burden)} (dropped {n_before - len(burden)})")

    # Attach per-sample 5 exposure columns for each SVT
    def _ncount(d: Dict[str, Dict[str, Set[str]]], sid: str, svt: str) -> int:
        return len(d.get(sid, {}).get(svt, set()))

    for svt in SVT_LIST:
        col_ef = f"n_exonfree_boundary_aggregate_{svt}"
        col_eo = f"n_exonoverlap_boundary_aggregate_{svt}"
        col_all = f"n_all_boundary_aggregate_{svt}"
        col_cnv = f"n_exonfreeCNV_boundary_aggregate_{svt}"
        col_strict = f"n_exonfreeCNV_exonfreeBin_boundary_aggregate_{svt}"

        burden[col_ef] = burden["sample_id"].map(
            lambda sid, s=svt: _ncount(s_ef, sid, s)
        ).astype(int)
        burden[col_eo] = burden["sample_id"].map(
            lambda sid, s=svt: _ncount(s_eo, sid, s)
        ).astype(int)
        burden[col_all] = burden["sample_id"].map(
            lambda sid, s=svt: _ncount(s_all, sid, s)
        ).astype(int)
        burden[col_cnv] = burden["sample_id"].map(
            lambda sid, s=svt: _ncount(s_efcnv, sid, s)
        ).astype(int)
        burden[col_strict] = burden["sample_id"].map(
            lambda sid, s=svt: _ncount(s_efcnv_efbin, sid, s)
        ).astype(int)

        log(
            f"  {svt}: exon-free-bin n_exposed="
            f"{int((burden[col_ef] > 0).sum())}, "
            f"exon-overlap-bin={int((burden[col_eo] > 0).sum())}, "
            f"all-bin={int((burden[col_all] > 0).sum())}, "
            f"exon-free-CNV={int((burden[col_cnv] > 0).sum())}, "
            f"exon-free-CNV+bin={int((burden[col_strict] > 0).sum())}"
        )

    # Diagnosis summary
    for dx in ["ASD", "SZ", "Healthy"]:
        n_dx = int((burden["Diagnosis"] == dx).sum())
        log(f"  Samples Diagnosis={dx}: {n_dx}")

    # ----- Step 5: regression grid -----
    all_rows: List[dict] = []

    for comp_name, case_label, ctrl_label in COMPARISONS:
        sub = burden[burden["Diagnosis"].isin([case_label, ctrl_label])].copy()
        sub["is_case"] = (sub["Diagnosis"] == case_label).astype(int)
        n_case = int(sub["is_case"].sum())
        n_ctrl = int((sub["is_case"] == 0).sum())
        log(f"\n-- {comp_name}: n_case={n_case}, n_control={n_ctrl} --")

        for svt in SVT_LIST:
            svt_lower = svt.lower()
            covar_cols = (
                ["Sex_numeric", f"log1p_total_{svt_lower}_bases"]
                + PC_COLS
                + [f"log1p_total_gene_{svt}"]
            )
            col_ef = f"n_exonfree_boundary_aggregate_{svt}"
            col_eo = f"n_exonoverlap_boundary_aggregate_{svt}"
            col_all = f"n_all_boundary_aggregate_{svt}"
            col_cnv = f"n_exonfreeCNV_boundary_aggregate_{svt}"
            col_strict = f"n_exonfreeCNV_exonfreeBin_boundary_aggregate_{svt}"

            def _run_pair(
                x_cols: List[str],
                analysis: str,
                bin_types: List[str],
            ) -> None:
                logit_results = fit_logit_wgs(sub, "is_case", x_cols, covar_cols)
                firth_results = fit_firth_wgs(sub, "is_case", x_cols, covar_cols)
                # Merge by exposure (x_cols order is preserved)
                for idx, xc in enumerate(x_cols):
                    lr = next(r for r in logit_results if r["exposure"] == xc)
                    fr = next(r for r in firth_results if r["exposure"] == xc)
                    row = {
                        "analysis": analysis,
                        "comparison": comp_name,
                        "sv_type": svt,
                        "bin_type": bin_types[idx],
                        "model": "Logit_main+Firth_sens",
                        "exposure": xc,
                        "n_case": lr["n_case"],
                        "n_ctrl": lr["n_ctrl"],
                        "carrier_case": lr["carrier_case"],
                        "carrier_ctrl": lr["carrier_ctrl"],
                        # Main Logit
                        "beta_logit": lr["beta"],
                        "se_logit": lr["se"],
                        "or_logit": lr["or"],
                        "or_lo95_logit": lr["or_lo95"],
                        "or_hi95_logit": lr["or_hi95"],
                        "p_two_logit": lr["p_two"],
                        "fit_status_logit": lr["fit_status"],
                        # Sensitivity Firth
                        "beta_firth": fr["beta"],
                        "se_firth": fr["se"],
                        "or_firth": fr["or"],
                        "or_lo95_firth": fr["or_lo95"],
                        "or_hi95_firth": fr["or_hi95"],
                        "p_two_firth": fr["p_two"],
                        "fit_status_firth": fr["fit_status"],
                    }
                    all_rows.append(row)
                    if np.isfinite(row["or_logit"]):
                        log(
                            f"    [{analysis} {bin_types[idx]}] Logit "
                            f"OR={row['or_logit']:.3f} "
                            f"[{row['or_lo95_logit']:.3f}-{row['or_hi95_logit']:.3f}] "
                            f"P={row['p_two_logit']:.4g}  "
                            f"Firth OR={row['or_firth']:.3f} "
                            f"P={row['p_two_firth']:.4g}  "
                            f"carriers={row['carrier_case']}/{row['carrier_ctrl']} "
                            f"({row['fit_status_firth']})"
                        )

            log(f"  -- {svt} --")
            # A1_separate: 3 bin_types
            _run_pair([col_ef], "A1_separate", ["exonfree"])
            _run_pair([col_eo], "A1_separate", ["exonoverlap"])
            _run_pair([col_all], "A1_separate", ["all"])
            # A2_joint: 2 bin_types together
            _run_pair(
                [col_ef, col_eo],
                "A2_joint",
                ["exonfree_in_joint", "exonoverlap_in_joint"],
            )
            # A3_exonfreeCNV
            _run_pair([col_cnv], "A3_exonfreeCNV", ["exonfreeCNV"])
            # A3_strict
            _run_pair(
                [col_strict],
                "A3_strict_exonfreeCNV_exonfreeBin",
                ["exonfreeCNV_exonfreeBin"],
            )

    # ----- Step 6: multiple testing (BH + Holm) -----
    df_out = pd.DataFrame(all_rows)

    # Multiple testing stratum: ASD × DEL × 4 discovery-positive keys
    mask_disc = (
        (df_out["comparison"] == "ASD_vs_Healthy")
        & (df_out["sv_type"] == "DEL")
        & df_out.apply(
            lambda r: (r["analysis"], r["bin_type"]) in DISCOVERY_POSITIVE_KEYS,
            axis=1,
        )
    )
    disc_idx = df_out.index[mask_disc].tolist()

    for pcol, bh_col, holm_col in [
        ("p_two_logit", "p_bh_logit", "p_holm_logit"),
        ("p_two_firth", "p_bh_firth", "p_holm_firth"),
    ]:
        df_out[bh_col] = np.nan
        df_out[holm_col] = np.nan
        if disc_idx:
            pvals = df_out.loc[disc_idx, pcol].to_numpy(dtype=float)
            finite_mask = np.isfinite(pvals)
            if finite_mask.any():
                sub_idx = [
                    disc_idx[i] for i, ok in enumerate(finite_mask) if ok
                ]
                sub_p = pvals[finite_mask]
                _, bh_vals, _, _ = multipletests(sub_p, method="fdr_bh")
                _, holm_vals, _, _ = multipletests(sub_p, method="holm")
                df_out.loc[sub_idx, bh_col] = bh_vals
                df_out.loc[sub_idx, holm_col] = holm_vals

    # Column ordering
    col_order = [
        "analysis", "comparison", "sv_type", "bin_type", "model", "exposure",
        "n_case", "n_ctrl", "carrier_case", "carrier_ctrl",
        # Logit main
        "beta_logit", "se_logit", "or_logit", "or_lo95_logit", "or_hi95_logit",
        "p_two_logit", "p_bh_logit", "p_holm_logit", "fit_status_logit",
        # Firth sensitivity
        "beta_firth", "se_firth", "or_firth", "or_lo95_firth", "or_hi95_firth",
        "p_two_firth", "p_bh_firth", "p_holm_firth", "fit_status_firth",
    ]
    col_order = [c for c in col_order if c in df_out.columns]
    df_out = df_out[col_order]

    df_out.to_csv(str(OUT_TSV), sep="\t", index=False)
    log(f"\nSaved main result TSV: {OUT_TSV}  ({len(df_out)} rows)")

    # ----- Step 7: summary TSV (ASD × DEL 4 discovery-positive) -----
    summary = df_out.loc[disc_idx].copy() if disc_idx else pd.DataFrame()
    summary.to_csv(str(OUT_SUMMARY), sep="\t", index=False)
    log(f"Saved summary TSV: {OUT_SUMMARY}  ({len(summary)} rows)")

    # ----- Step 8: preview -----
    log("\n" + "=" * 72)
    log("PREVIEW: ASD_vs_Healthy × DEL × 4 discovery-positive")
    log("=" * 72)
    for ana, btype in DISCOVERY_POSITIVE_KEYS:
        sel = df_out[
            (df_out["analysis"] == ana)
            & (df_out["comparison"] == "ASD_vs_Healthy")
            & (df_out["sv_type"] == "DEL")
            & (df_out["bin_type"] == btype)
        ]
        if len(sel) == 0:
            continue
        r = sel.iloc[0]
        log(
            f"  {ana}/{btype}:  "
            f"Logit OR={r['or_logit']:.3f} P={r['p_two_logit']:.4g} "
            f"p_BH={r['p_bh_logit']:.4g} p_Holm={r['p_holm_logit']:.4g}  "
            f"Firth OR={r['or_firth']:.3f} P={r['p_two_firth']:.4g} "
            f"p_BH={r['p_bh_firth']:.4g} p_Holm={r['p_holm_firth']:.4g}  "
            f"carriers={int(r['carrier_case'])}/{int(r['carrier_ctrl'])}"
        )

    log("\nPREVIEW: SZ_vs_Healthy × DEL (key rows)")
    for ana, btype in DISCOVERY_POSITIVE_KEYS:
        sel = df_out[
            (df_out["analysis"] == ana)
            & (df_out["comparison"] == "SZ_vs_Healthy")
            & (df_out["sv_type"] == "DEL")
            & (df_out["bin_type"] == btype)
        ]
        if len(sel) == 0:
            continue
        r = sel.iloc[0]
        log(
            f"  {ana}/{btype}:  "
            f"Logit OR={r['or_logit']:.3f} P={r['p_two_logit']:.4g}  "
            f"Firth OR={r['or_firth']:.3f} P={r['p_two_firth']:.4g}"
        )

    elapsed_total = time.time() - _T0
    log("\n" + "=" * 72)
    log(
        f"Script 12 v2 DONE: 42_compute_wgs_exon_exclusion_v2.py  "
        f"total elapsed = {elapsed_total:.1f} s"
    )
    log("=" * 72)


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()
