#!/usr/bin/env python3
# 42_compute_mssng_exon_exclusion_v4.py
#
# 処理内容:
#  - MSSNG replication cohort の exon-exclusion sensitivity を計算する
#    (tad04292026 パイプライン)。
#  - [v3 -> v4 変更点]:
#    * [refactor] 共通 helper module common/exon_helpers.py に切り出し:
#       - build_merged_exon_bed() / overlaps_any_exon() / load_bin_l2_annotation()
#       - L2_MEMBERSHIP_COLS / log()
#      これらは Module 12 (WGS exon-exclusion) と完全に同一ロジックだったため、
#      集約してコード重複を削減 (~150 行削減)。helper-internal の数値結果は
#      v3 と bit-identical (関数本体に変更なし)。
#    * [path] 入力 path を common.paths_v1 (envvar override 可、NIG default
#      保持) から import するように切替:
#       - GENCODE_GTF (旧: 直接記載) → common.paths_v1.GENCODE_GTF
#       - BIN_L2_ANNOT (旧: bin_l2_annotation_v2.tsv.gz) → F_02_BIN_L2_ANNOTATION
#         (= bin_l2_annotation_v3.tsv.gz, S2-S5 削除済 v3 出力)
#       - BASE_TAD → PIPELINE_ROOT
#    * 出力 dir / file: output_v3 → output_v4
#       - mssng_exon_exclusion_v3.tsv → mssng_exon_exclusion_v4.tsv
#       - mssng_exon_exclusion_summary_v3.tsv → mssng_exon_exclusion_summary_v4.tsv
#    * 解析ロジック (A1/A2/A3/A3-strict、GEE Binom Indep main + Firth sens、
#      BH/Holm、3 段 pre-check) は v3 と完全一致。
#    * group_primary / S2-S5 は使わない (Module 13 は exon-free 軸のみ)。
#  - [v3 既存の機能 (継承)]: Script 12 v2 (WGS) と同じ MIN_CELL_COUNT=5
#    ガードを GEE/Firth 両方の fit 関数に持つ。
#  - Input:
#      tad04292026/02_bin_l2_annotation/output_v3/bin_l2_annotation_v3.tsv.gz
#          (10 L2 class membership + bin 座標; hpc_astro 元から除外、S2-S5 削除済)
#      noncoding_tad_mssng_03132026/output_mssng_factorial_overlap_v1/
#          mssng_sample_event_bin_overlap_v1.tsv.gz  (sample × event × bin)
#          mssng_sample_covariates_v1.tsv           (Status/Sex/FAMILYID 等)
#      Annotation: GENCODE v46 GTF (paths_v1.GENCODE_GTF)
#  - 前提: MSSNG QC は heffel と tad04292026 Script 09 v18 で一致
#    (n_case=5955, n_control=511, n_families=4978)。factorial overlap は
#    DEL のみ (v1 仕様) のため、本 sensitivity は DEL のみ実施。
#  - Steps:
#      1. GENCODE v46 protein-coding exon を merge (helper)
#      2. bin_l2_annotation の 25-kb bin 座標と exon を交差 (helper)
#      3. sample × event × bin overlap を読み、any_l2_bins に絞って:
#           A1/A2 軸 (bin レベル): exonfree / exonoverlap / all
#           A3 軸 (SV レベル)   : event_id (chrN:start-end) から SV 座標を
#                                 parse し、SV 全体が exon 非重複かで判定
#                                 (parse 失敗数を log に記録)
#           A3-strict            : exonfreeCNV かつ hit bin が exon-free
#         サンプル単位で unique bin 数を集計
#      4. covariate TSV に 5 exposure 列を追加し、以下の回帰を実施:
#           Main: statsmodels GEE Binomial Independence (FAMILYID cluster)
#           Sensitivity: 自前 Firth (Heinze & Schemper 2002, Newton-Raphson,
#                         pooled; cluster は heffel 仕様に従い無視)
#                         QR-with-pivoting で rank-deficient 列を自動 drop
#                         (MSSNG platform one-hot 5 列 + intercept の collinearity)
#                         protected_idx = (intercept, exposure[, exposure2])
#           MIN_CELL_COUNT=5 ガードを fit_gee_mssng / fit_firth_mssng に
#           適用 (outcome variance / per-exposure variance / total_carriers
#           の 3 段 pre-check; A2 joint は total_carriers を合算して判定)
#      5. 解析グリッド: ASD_vs_unaffected_sibling × DEL × 4 analyses
#           A1_separate (3 bin_types), A2_joint (2 bin_types),
#           A3_exonfreeCNV (1), A3_strict (1)  → 7 結果行 × 2 methods
#      6. Multiple testing: 4 discovery-positive
#         (A1 exonfree, A2 exonfree_in_joint, A3 exonfreeCNV, A3 strict)
#         に対して BH-FDR + Holm (GEE P, Firth P それぞれ)
#      7. 出力 TSV:
#           mssng_exon_exclusion_v4.tsv (全 7 行)
#           mssng_exon_exclusion_summary_v4.tsv (4 discovery-positive 行)
#  - 出力先: ${PIPELINE_ROOT}/13_exon_exclusion_mssng/output_v4/
#  - 処理時間を先頭と末尾で記録

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=42_compute_mssng_exon_exclusion_v4_%j.out
#SBATCH --error=42_compute_mssng_exon_exclusion_v4_%j.err

from __future__ import annotations

import gzip
import os
import sys
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
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.genmod.families import Binomial
from statsmodels.genmod.cov_struct import Independence
from statsmodels.stats.multitest import multipletests

# ----- common helpers (v4 で切り出し) -----
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.exon_helpers import (  # noqa: E402
    L2_MEMBERSHIP_COLS,
    build_merged_exon_bed,
    overlaps_any_exon,
    load_bin_l2_annotation,
    log,
)
from common.paths_v1 import (  # noqa: E402
    PIPELINE_ROOT,
    GENCODE_GTF,
    F_02_BIN_L2_ANNOTATION,
)


# =========================================================
# PATHS (v4: PIPELINE_ROOT / F_02_BIN_L2_ANNOTATION / GENCODE_GTF を common.paths_v1 から)
# =========================================================
# MSSNG event overlap + covariates は noncoding_tad_mssng_03132026 配下の read-only
# データ (heffel/tad04292026 共有)。env var BASE_MSSNG_HEFFEL で override 可。
BASE_MSSNG_HEFFEL = Path(os.environ.get(
    "BASE_MSSNG_HEFFEL",
    "/lustre12/home/kushima-pg/noncoding_tad_mssng_03132026",
))  # CONFIGURE

BIN_L2_ANNOT = F_02_BIN_L2_ANNOTATION  # bin_l2_annotation_v3.tsv.gz (S2-S5 削除済)

EVENT_OVERLAP = (
    BASE_MSSNG_HEFFEL / "output_mssng_factorial_overlap_v1"
    / "mssng_sample_event_bin_overlap_v1.tsv.gz"
)
COVAR_TABLE = (
    BASE_MSSNG_HEFFEL / "output_mssng_factorial_overlap_v1"
    / "mssng_sample_covariates_v1.tsv"
)

OUT_DIR = PIPELINE_ROOT / "13_exon_exclusion_mssng" / "output_v4"
OUT_TSV = OUT_DIR / "mssng_exon_exclusion_v4.tsv"
OUT_SUMMARY = OUT_DIR / "mssng_exon_exclusion_summary_v4.tsv"


# =========================================================
# CONSTANTS (Module 13-specific; L2_MEMBERSHIP_COLS は helper から)
# =========================================================
CASE_LABEL = "ASD"
CTRL_LABEL = "unaffected_sibling"
TARGET_SVT = "DEL"  # MSSNG factorial overlap v1 は DEL のみ

# 4 discovery-positive (ASD vs sib × DEL only) for BH/Holm
DISCOVERY_POSITIVE_KEYS = [
    ("A1_separate", "exonfree"),
    ("A2_joint", "exonfree_in_joint"),
    ("A3_exonfreeCNV", "exonfreeCNV"),
    ("A3_strict_exonfreeCNV_exonfreeBin", "exonfreeCNV_exonfreeBin"),
]

# Minimum cell count guard (Script 12 v2 WGS と同値)
MIN_CELL_COUNT = 5


# =========================================================
# MSSNG event_id parser: "chrN:start-end" (Module 13 specific)
# =========================================================
def parse_mssng_event_id(event_id: str) -> Tuple[str, int, int] | None:
    try:
        parts = event_id.split(":")
        if len(parts) != 2:
            return None
        chrom = parts[0]
        s, e = parts[1].split("-")
        return (chrom, int(s), int(e))
    except Exception:
        return None


# =========================================================
# Count disrupted bins per sample (5 stratifications) — Module 13 specific
# =========================================================
def count_disrupted_bins(
    event_overlap_path: Path,
    exon_free_bins: Set[str],
    any_l2_bins: Set[str],
    merged_exons: Dict[str, List[Tuple[int, int]]],
    exon_starts: Dict[str, List[int]],
) -> Tuple[
    Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Set[str]],
    Dict[str, Set[str]], Dict[str, Set[str]]
]:
    """MSSNG sample×event×bin overlap から、sample 単位で 5 種の unique bin 集合を構築。"""
    log(f"Reading event overlap: {event_overlap_path}")
    s_ef: Dict[str, Set[str]] = defaultdict(set)
    s_eo: Dict[str, Set[str]] = defaultdict(set)
    s_all: Dict[str, Set[str]] = defaultdict(set)
    s_efcnv: Dict[str, Set[str]] = defaultdict(set)
    s_efcnv_efbin: Dict[str, Set[str]] = defaultdict(set)

    sv_exonfree_cache: Dict[str, bool] = {}
    sv_parse_failed: Set[str] = set()
    sample_example_failed: List[str] = []
    n_rows = 0
    n_kept = 0
    n_ef_bin_hit = 0
    n_a3 = 0
    n_a3_strict = 0

    with gzip.open(str(event_overlap_path), "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        for c in ["sample_id", "bin_id", "event_id"]:
            if c not in col:
                raise KeyError(f"event overlap: missing column {c}")
        idx_sid = col["sample_id"]
        idx_bin = col["bin_id"]
        idx_eid = col["event_id"]

        for line in f:
            n_rows += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(idx_sid, idx_bin, idx_eid):
                continue
            bin_id = fields[idx_bin]
            if bin_id not in any_l2_bins:
                continue
            n_kept += 1

            sid = fields[idx_sid]
            event_id = fields[idx_eid]
            bin_is_exonfree = bin_id in exon_free_bins

            s_all[sid].add(bin_id)
            if bin_is_exonfree:
                s_ef[sid].add(bin_id)
                n_ef_bin_hit += 1
            else:
                s_eo[sid].add(bin_id)

            if event_id not in sv_exonfree_cache:
                parsed = parse_mssng_event_id(event_id)
                if parsed is None:
                    sv_exonfree_cache[event_id] = False
                    sv_parse_failed.add(event_id)
                    if len(sample_example_failed) < 5:
                        sample_example_failed.append(event_id)
                else:
                    ch, ss, ee = parsed
                    sv_exonfree_cache[event_id] = not overlaps_any_exon(
                        ch, ss, ee, merged_exons, exon_starts)

            if sv_exonfree_cache[event_id]:
                s_efcnv[sid].add(bin_id)
                n_a3 += 1
                if bin_is_exonfree:
                    s_efcnv_efbin[sid].add(bin_id)
                    n_a3_strict += 1

    n_sv = len(sv_exonfree_cache)
    n_efsv = sum(1 for v in sv_exonfree_cache.values() if v)
    n_parse_fail = len(sv_parse_failed)
    log(f"  Raw event-bin rows: {n_rows}")
    log(f"  Kept (any-L2 bin): {n_kept}")
    log(f"  Exon-free bin rows: {n_ef_bin_hit}")
    log(f"  Unique SVs on L2 bins: {n_sv}")
    log(f"  event_id parse failures: {n_parse_fail} "
        f"({100*n_parse_fail/max(n_sv,1):.2f}% of unique SVs)"
        f"  -> treated as non-exon-free (conservative)")
    if sample_example_failed:
        ex_str = ", ".join(repr(e) for e in sample_example_failed)
        log(f"    example failed event_ids: {ex_str}")
    log(f"  Exon-free SVs (axis B): {n_efsv} "
        f"({100*n_efsv/max(n_sv,1):.1f}%)")
    log(f"  A3 rows (exon-free CNV × any L2 bin): {n_a3}")
    log(f"    ... hit bin is exon-free (A3-strict): {n_a3_strict}")

    return s_ef, s_eo, s_all, s_efcnv, s_efcnv_efbin


# =========================================================
# QR-based rank reduction (Module 13 specific — drop_idx 形式)
# =========================================================
def _reduce_rank_via_qr(
    X: np.ndarray,
    protected_idx: Tuple[int, ...] = (0, 1),
    tol_rel: float = 1e-10,
) -> List[int]:
    """X の列の中で rank-deficient な列を検出し、drop すべき列 index リストを返す。
    protected_idx の列は必ず残す (intercept + exposure)。
    """
    n, p = X.shape
    if p == 0:
        return []
    col_norm = np.linalg.norm(X, axis=0)
    col_norm = np.where(col_norm == 0, 1.0, col_norm)
    Xs = X / col_norm

    Q, R, piv = sp_linalg.qr(Xs, mode="economic", pivoting=True)
    diag = np.abs(np.diag(R))
    if len(diag) == 0:
        return []
    thr = diag[0] * tol_rel
    rank = int(np.sum(diag > thr))

    kept_by_pivot = set(piv[:rank].tolist())
    protected = set(protected_idx)
    kept = kept_by_pivot | protected
    if len(kept) > rank:
        non_prot_kept = [c for c in piv[:rank] if c not in protected]
        excess = len(kept) - rank
        for c in reversed(non_prot_kept):
            if excess <= 0:
                break
            if c in kept:
                kept.remove(c)
                excess -= 1

    drop = [c for c in range(p) if c not in kept]
    return drop


# =========================================================
# Firth penalized logistic (Module 13 specific — Dict return)
# =========================================================
def firth_logit_fit(
    X: np.ndarray,
    y: np.ndarray,
    maxiter: int = 100,
    tol: float = 1e-8,
    protected_idx: Tuple[int, ...] = (0, 1),
) -> Dict:
    """Firth penalized logistic regression (Newton-Raphson)。
    QR-with-pivoting で rank-deficient 列を drop し、reduced X で fit。
    Returns full-size beta/se (drop 列は NaN)。
    """
    n, p = X.shape
    drop_idx = _reduce_rank_via_qr(X, protected_idx=protected_idx)
    keep_idx = [c for c in range(p) if c not in drop_idx]
    X_red = X[:, keep_idx]
    pr = X_red.shape[1]
    status_rank_note = f"ok_rank_reduced:{len(drop_idx)}" if drop_idx else "ok"

    beta = np.zeros(pr)
    converged = False
    for it in range(maxiter):
        eta = X_red @ beta
        eta = np.clip(eta, -30.0, 30.0)
        pi = 1.0 / (1.0 + np.exp(-eta))
        w = pi * (1.0 - pi)
        W = np.diag(w)
        XtWX = X_red.T @ W @ X_red
        try:
            inv_XtWX = np.linalg.inv(XtWX)
        except np.linalg.LinAlgError:
            return {
                "beta": np.full(p, np.nan),
                "se": np.full(p, np.nan),
                "converged": False,
                "status": "singular_hessian",
                "iter": it,
            }
        w_sqrt = np.sqrt(w)
        H = (X_red * w_sqrt[:, None]) @ inv_XtWX @ (X_red.T * w_sqrt[None, :])
        h = np.diag(H)
        u_star = X_red.T @ (y - pi + h * (0.5 - pi))
        delta = inv_XtWX @ u_star
        beta_new = beta + delta
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            converged = True
            break
        beta = beta_new

    eta = np.clip(X_red @ beta, -30.0, 30.0)
    pi = 1.0 / (1.0 + np.exp(-eta))
    w = pi * (1.0 - pi)
    W = np.diag(w)
    XtWX = X_red.T @ W @ X_red
    try:
        inv_XtWX = np.linalg.inv(XtWX)
        se_red = np.sqrt(np.diag(inv_XtWX))
    except np.linalg.LinAlgError:
        se_red = np.full(pr, np.nan)

    beta_full = np.full(p, np.nan)
    se_full = np.full(p, np.nan)
    for j, ki in enumerate(keep_idx):
        beta_full[ki] = beta[j]
        se_full[ki] = se_red[j]

    return {
        "beta": beta_full,
        "se": se_full,
        "converged": converged,
        "status": status_rank_note if converged else "not_converged",
        "iter": it + 1,
    }


# =========================================================
# Empty-fit helpers (Module 13 specific)
# =========================================================
def _empty_fit_result_gee(
    status: str, exposure: str,
    n_case: int, n_ctrl: int,
    carrier_case: int, carrier_ctrl: int,
    n_clusters: int,
) -> Dict:
    return {
        "exposure": exposure,
        "beta": np.nan, "se": np.nan,
        "or": np.nan, "ci_lo": np.nan, "ci_hi": np.nan,
        "p_two": np.nan,
        "n_case": int(n_case), "n_ctrl": int(n_ctrl),
        "carrier_case": int(carrier_case),
        "carrier_ctrl": int(carrier_ctrl),
        "n_clusters": int(n_clusters),
        "fit_status": status,
    }


def _empty_fit_result_firth(
    status: str, exposure: str,
    n_case: int, n_ctrl: int,
    carrier_case: int, carrier_ctrl: int,
) -> Dict:
    return {
        "exposure": exposure,
        "beta": np.nan, "se": np.nan,
        "or": np.nan, "ci_lo": np.nan, "ci_hi": np.nan,
        "p_two": np.nan,
        "n_case": int(n_case), "n_ctrl": int(n_ctrl),
        "carrier_case": int(carrier_case),
        "carrier_ctrl": int(carrier_ctrl),
        "fit_status": status,
    }


# =========================================================
# GEE Binomial Independence fit (main method) — Module 13 specific
# =========================================================
def fit_gee_mssng(
    df: pd.DataFrame,
    y_col: str,
    x_cols: List[str],
    covar_cols: List[str],
    cluster_col: str,
    min_cell: int = MIN_CELL_COUNT,
) -> List[Dict]:
    cols_needed = covar_cols + list(x_cols) + [y_col]
    sub = df[cols_needed].copy()
    sub = sub.replace([np.inf, -np.inf], np.nan).dropna()
    groups_all = df.loc[sub.index, cluster_col].values
    n_case = int((sub[y_col] == 1).sum())
    n_ctrl = int((sub[y_col] == 0).sum())
    n_clusters = int(len(np.unique(groups_all)))

    carrier_info: List[Dict] = []
    for xc in x_cols:
        cc = int(((sub[y_col] == 1) & (sub[xc] >= 1)).sum())
        cn = int(((sub[y_col] == 0) & (sub[xc] >= 1)).sum())
        carrier_info.append({"_xc": xc, "_cc": cc, "_cn": cn})

    if sub[y_col].nunique() < 2:
        return [
            _empty_fit_result_gee(
                "no_outcome_variance", r["_xc"],
                n_case, n_ctrl, r["_cc"], r["_cn"], n_clusters,
            )
            for r in carrier_info
        ]

    for xc in x_cols:
        if sub[xc].nunique() < 2:
            return [
                _empty_fit_result_gee(
                    f"no_variance_{xc}", r["_xc"],
                    n_case, n_ctrl, r["_cc"], r["_cn"], n_clusters,
                )
                for r in carrier_info
            ]

    total_carriers = sum(r["_cc"] + r["_cn"] for r in carrier_info)
    if total_carriers < min_cell:
        return [
            _empty_fit_result_gee(
                "insufficient_carriers", r["_xc"],
                n_case, n_ctrl, r["_cc"], r["_cn"], n_clusters,
            )
            for r in carrier_info
        ]

    out: List[Dict] = []
    try:
        cols_order = covar_cols + list(x_cols)
        X_df = sub[cols_order].copy()
        X_df = sm.add_constant(X_df, has_constant="add")
        y = sub[y_col].values.astype(float)
        groups = groups_all

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = GEE(
                endog=y, exog=X_df.values, groups=groups,
                family=Binomial(), cov_struct=Independence(),
            )
            res = model.fit(maxiter=200)

        col_names = list(X_df.columns)
        for r in carrier_info:
            xc = r["_xc"]
            idx = col_names.index(xc)
            beta = float(res.params[idx])
            se = float(res.bse[idx])
            pval = float(res.pvalues[idx])
            or_val = float(np.exp(beta))
            ci_lo = float(np.exp(beta - 1.96 * se))
            ci_hi = float(np.exp(beta + 1.96 * se))
            out.append({
                "exposure": xc,
                "beta": beta, "se": se,
                "or": or_val, "ci_lo": ci_lo, "ci_hi": ci_hi,
                "p_two": pval,
                "n_case": n_case, "n_ctrl": n_ctrl,
                "carrier_case": r["_cc"],
                "carrier_ctrl": r["_cn"],
                "n_clusters": n_clusters,
                "fit_status": "ok" if getattr(res, "converged", True) else "not_converged",
            })
    except Exception as e:
        out = [
            _empty_fit_result_gee(
                f"gee_failed:{str(e)[:60]}", r["_xc"],
                n_case, n_ctrl, r["_cc"], r["_cn"], n_clusters,
            )
            for r in carrier_info
        ]
    return out


# =========================================================
# Firth MSSNG fit (sensitivity; pooled) — Module 13 specific
# =========================================================
def fit_firth_mssng(
    df: pd.DataFrame,
    y_col: str,
    x_cols: List[str],
    covar_cols: List[str],
    min_cell: int = MIN_CELL_COUNT,
) -> List[Dict]:
    cols_needed = covar_cols + list(x_cols) + [y_col]
    sub = df[cols_needed].copy()
    sub = sub.replace([np.inf, -np.inf], np.nan).dropna()
    n_case = int((sub[y_col] == 1).sum())
    n_ctrl = int((sub[y_col] == 0).sum())

    carrier_info: List[Dict] = []
    for xc in x_cols:
        cc = int(((sub[y_col] == 1) & (sub[xc] >= 1)).sum())
        cn = int(((sub[y_col] == 0) & (sub[xc] >= 1)).sum())
        carrier_info.append({"_xc": xc, "_cc": cc, "_cn": cn})

    if sub[y_col].nunique() < 2:
        return [
            _empty_fit_result_firth(
                "no_outcome_variance", r["_xc"],
                n_case, n_ctrl, r["_cc"], r["_cn"],
            )
            for r in carrier_info
        ]

    for xc in x_cols:
        if sub[xc].nunique() < 2:
            return [
                _empty_fit_result_firth(
                    f"no_variance_{xc}", r["_xc"],
                    n_case, n_ctrl, r["_cc"], r["_cn"],
                )
                for r in carrier_info
            ]

    total_carriers = sum(r["_cc"] + r["_cn"] for r in carrier_info)
    if total_carriers < min_cell:
        return [
            _empty_fit_result_firth(
                "insufficient_carriers", r["_xc"],
                n_case, n_ctrl, r["_cc"], r["_cn"],
            )
            for r in carrier_info
        ]

    out: List[Dict] = []
    try:
        cols_order = list(x_cols) + covar_cols
        X_df = sub[cols_order].copy()
        X_df.insert(0, "const", 1.0)
        y = sub[y_col].values.astype(float)
        X = X_df.values.astype(float)
        col_names = list(X_df.columns)

        protected = tuple(range(1, 1 + len(x_cols)))
        protected = (0,) + protected

        res = firth_logit_fit(X, y, maxiter=200, tol=1e-8, protected_idx=protected)
        for r in carrier_info:
            xc = r["_xc"]
            idx = col_names.index(xc)
            beta = float(res["beta"][idx])
            se = float(res["se"][idx])
            if np.isfinite(beta) and np.isfinite(se) and se > 0:
                z = beta / se
                p_two = 2.0 * (1.0 - sp_stats.norm.cdf(abs(z)))
                or_val = float(np.exp(beta))
                ci_lo = float(np.exp(beta - 1.96 * se))
                ci_hi = float(np.exp(beta + 1.96 * se))
            else:
                p_two = np.nan
                or_val = np.nan
                ci_lo = np.nan
                ci_hi = np.nan
            out.append({
                "exposure": xc,
                "beta": beta, "se": se,
                "or": or_val, "ci_lo": ci_lo, "ci_hi": ci_hi,
                "p_two": p_two,
                "n_case": n_case, "n_ctrl": n_ctrl,
                "carrier_case": r["_cc"],
                "carrier_ctrl": r["_cn"],
                "fit_status": res["status"],
            })
    except Exception as e:
        out = [
            _empty_fit_result_firth(
                f"firth_failed:{str(e)[:60]}", r["_xc"],
                n_case, n_ctrl, r["_cc"], r["_cn"],
            )
            for r in carrier_info
        ]
    return out


# =========================================================
# MAIN
# =========================================================
def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    log("=" * 72)
    log("Script 13 v4 START: MSSNG exon-exclusion sensitivity (tad04292026)")
    log("  42_compute_mssng_exon_exclusion_v4.py")
    log("  Helpers: common.exon_helpers (build_merged_exon_bed, load_bin_l2_annotation, etc.)")
    log("  Input: bin_l2_annotation_v3 (S2-S5 削除済) + MSSNG factorial overlap v1")
    log("=" * 72)
    log("Analyses: [A1] Separate, [A2] Joint, [A3] exon-free CNV, [A3-strict]")
    log("Comparison: ASD_vs_unaffected_sibling; SVT: DEL only")
    log("Methods: GEE Binomial Indep. (FAMILYID cluster, main) + "
        "Firth (sensitivity, QR rank reduction)")
    log(f"Guard: MIN_CELL_COUNT = {MIN_CELL_COUNT} "
        "(total_carriers basis; A2 joint では合算で判定)")
    log("Multiple testing: BH + Holm on 4 discovery-positive")
    log(f"Output dir: {OUT_DIR}")

    # -- Exon index --
    merged_exons = build_merged_exon_bed(GENCODE_GTF)
    exon_starts = {ch: [iv[0] for iv in merged_exons[ch]] for ch in merged_exons}

    # -- bin L2 annotation -> exon-free bins + any_l2 pool --
    exon_free, exon_overlap, any_l2 = load_bin_l2_annotation(
        BIN_L2_ANNOT, merged_exons, exon_starts
    )

    # -- Sample-level bin counts (5 stratifications) --
    (s_ef, s_eo, s_all, s_efcnv, s_efcnv_efbin) = count_disrupted_bins(
        EVENT_OVERLAP, exon_free, any_l2, merged_exons, exon_starts
    )

    # -- Covariates --
    log(f"Reading covariates: {COVAR_TABLE}")
    covar = pd.read_csv(str(COVAR_TABLE), sep="\t")
    log(f"  Raw rows: {len(covar)}")

    anc_cols = [c for c in covar.columns if c.startswith("anc_")]
    plat_cols = [c for c in covar.columns if c.startswith("plat_")]
    ess_cols = [
        "Sex_numeric", "log1p_total_del_bases", "log1p_total_gene_DEL",
    ] + anc_cols + plat_cols

    missing = [c for c in ess_cols + ["Status", "FAMILYID", "sample_id"] if c not in covar.columns]
    if missing:
        raise RuntimeError(f"Covariate TSV missing columns: {missing}")

    n_before = len(covar)
    covar = covar.dropna(subset=["Status", "FAMILYID"] + ess_cols).reset_index(drop=True)
    log(f"  After dropna on covariates: {len(covar)} (dropped {n_before - len(covar)})")

    # -- Add exposure columns (sample -> unique bin count) --
    def _n(d: Dict[str, Set[str]], sid: str) -> int:
        return len(d.get(sid, set()))

    col_exonfree = "n_exonfree_boundary_aggregate_DEL"
    col_exonoverlap = "n_exonoverlap_boundary_aggregate_DEL"
    col_all = "n_all_boundary_aggregate_DEL"
    col_efcnv = "n_exonfreeCNV_boundary_aggregate_DEL"
    col_efcnv_efbin = "n_exonfreeCNV_exonfreeBin_boundary_aggregate_DEL"

    covar[col_exonfree] = covar["sample_id"].map(lambda s: _n(s_ef, s)).astype(int)
    covar[col_exonoverlap] = covar["sample_id"].map(lambda s: _n(s_eo, s)).astype(int)
    covar[col_all] = covar["sample_id"].map(lambda s: _n(s_all, s)).astype(int)
    covar[col_efcnv] = covar["sample_id"].map(lambda s: _n(s_efcnv, s)).astype(int)
    covar[col_efcnv_efbin] = covar["sample_id"].map(lambda s: _n(s_efcnv_efbin, s)).astype(int)

    for lab, c in [
        ("exon-free-bin", col_exonfree),
        ("exon-overlap-bin", col_exonoverlap),
        ("all-bin", col_all),
        ("exon-free-CNV", col_efcnv),
        ("exon-free-CNV+bin", col_efcnv_efbin),
    ]:
        n_exp = int((covar[c] > 0).sum())
        log(f"  {lab}: n_exposed (all Status)={n_exp}")

    for st in [CASE_LABEL, CTRL_LABEL]:
        n = int((covar["Status"] == st).sum())
        log(f"  Samples Status={st}: {n}")

    # -- Subset ASD / sibling --
    sub = covar[covar["Status"].isin([CASE_LABEL, CTRL_LABEL])].copy()
    sub["y"] = (sub["Status"] == CASE_LABEL).astype(int)
    log(f"")
    log(f"-- ASD_vs_unaffected_sibling: n_case={int(sub['y'].sum())}, "
        f"n_control={int((sub['y'] == 0).sum())} --")
    log(f"   Families: {sub['FAMILYID'].nunique()}")

    results: List[Dict] = []

    def _run_one(
        analysis: str, bin_type: str, x_cols: List[str], sub_df: pd.DataFrame,
    ) -> None:
        gee_rows = fit_gee_mssng(sub_df, "y", x_cols, ess_cols, "FAMILYID")
        firth_rows = fit_firth_mssng(sub_df, "y", x_cols, ess_cols)
        gee_by = {r["exposure"]: r for r in gee_rows}
        firth_by = {r["exposure"]: r for r in firth_rows}

        for xc in x_cols:
            if len(x_cols) == 1:
                bt = bin_type
            else:
                bt = (
                    "exonfree_in_joint" if "exonfree" in xc and "overlap" not in xc
                    else "exonoverlap_in_joint"
                )

            g = gee_by.get(xc, {})
            f = firth_by.get(xc, {})
            row = {
                "analysis": analysis,
                "comparison": "ASD_vs_unaffected_sibling",
                "sv_type": TARGET_SVT,
                "bin_type": bt,
                "model": "GEE_BinomIndep_main+Firth_sens",
                "exposure": xc,
                "n_case": g.get("n_case", np.nan),
                "n_ctrl": g.get("n_ctrl", np.nan),
                "n_clusters": g.get("n_clusters", np.nan),
                "carrier_case": g.get("carrier_case", np.nan),
                "carrier_ctrl": g.get("carrier_ctrl", np.nan),
                "beta_gee": g.get("beta", np.nan),
                "se_gee": g.get("se", np.nan),
                "or_gee": g.get("or", np.nan),
                "or_lo95_gee": g.get("ci_lo", np.nan),
                "or_hi95_gee": g.get("ci_hi", np.nan),
                "p_two_gee": g.get("p_two", np.nan),
                "fit_status_gee": g.get("fit_status", ""),
                "beta_firth": f.get("beta", np.nan),
                "se_firth": f.get("se", np.nan),
                "or_firth": f.get("or", np.nan),
                "or_lo95_firth": f.get("ci_lo", np.nan),
                "or_hi95_firth": f.get("ci_hi", np.nan),
                "p_two_firth": f.get("p_two", np.nan),
                "fit_status_firth": f.get("fit_status", ""),
            }
            results.append(row)

            or_g = row["or_gee"]
            p_g = row["p_two_gee"]
            or_f = row["or_firth"]
            p_f = row["p_two_firth"]
            ci_lo = row["or_lo95_gee"]
            ci_hi = row["or_hi95_gee"]
            status = row["fit_status_gee"]
            log(f"    [{analysis} {bt}] GEE OR={or_g:.3f} "
                f"[{ci_lo:.3f}-{ci_hi:.3f}] P={p_g:.4g}  "
                f"Firth OR={or_f:.3f} P={p_f:.4g}  "
                f"carriers={row['carrier_case']}/{row['carrier_ctrl']} "
                f"({status})")

    # -- A1 separate --
    log(f"  -- A1 separate --")
    _run_one("A1_separate", "exonfree", [col_exonfree], sub)
    _run_one("A1_separate", "exonoverlap", [col_exonoverlap], sub)
    _run_one("A1_separate", "all", [col_all], sub)

    # -- A2 joint --
    log(f"  -- A2 joint --")
    r_corr = sub[col_exonfree].corr(sub[col_exonoverlap])
    vif = 1.0 / (1.0 - r_corr ** 2) if abs(r_corr) < 1.0 else float("inf")
    log(f"    r(ef, eo)={r_corr:.3f}, VIF={vif:.2f}")
    _run_one("A2_joint", "joint", [col_exonfree, col_exonoverlap], sub)

    # -- A3 exon-free CNV (any bin) --
    log(f"  -- A3 exon-free CNV --")
    _run_one("A3_exonfreeCNV", "exonfreeCNV", [col_efcnv], sub)

    # -- A3-strict --
    log(f"  -- A3 strict (exon-free CNV × exon-free bin) --")
    _run_one("A3_strict_exonfreeCNV_exonfreeBin",
             "exonfreeCNV_exonfreeBin", [col_efcnv_efbin], sub)

    # -- Multiple testing --
    df_all = pd.DataFrame(results)
    for prefix in ["p_two_gee", "p_two_firth"]:
        df_all[f"{prefix}_bh"] = np.nan
        df_all[f"{prefix}_holm"] = np.nan

    keys_set = set(DISCOVERY_POSITIVE_KEYS)
    mask_disc = df_all.apply(
        lambda r: (r["analysis"], r["bin_type"]) in keys_set, axis=1
    )
    disc_rows = df_all[mask_disc].copy()
    if len(disc_rows) != len(DISCOVERY_POSITIVE_KEYS):
        log(f"  [WARN] Expected {len(DISCOVERY_POSITIVE_KEYS)} discovery-positive "
            f"rows but found {len(disc_rows)}")

    if len(disc_rows) > 0:
        for prefix in ["p_two_gee", "p_two_firth"]:
            pvals = disc_rows[prefix].values.astype(float)
            mask = np.isfinite(pvals)
            if mask.sum() > 0:
                _, bh_adj, _, _ = multipletests(pvals[mask], method="fdr_bh")
                _, holm_adj, _, _ = multipletests(pvals[mask], method="holm")
                bh_full = np.full(len(pvals), np.nan)
                holm_full = np.full(len(pvals), np.nan)
                bh_full[mask] = bh_adj
                holm_full[mask] = holm_adj
                disc_rows.loc[:, f"{prefix}_bh"] = bh_full
                disc_rows.loc[:, f"{prefix}_holm"] = holm_full
        df_all.loc[mask_disc, "p_two_gee_bh"] = disc_rows["p_two_gee_bh"].values
        df_all.loc[mask_disc, "p_two_gee_holm"] = disc_rows["p_two_gee_holm"].values
        df_all.loc[mask_disc, "p_two_firth_bh"] = disc_rows["p_two_firth_bh"].values
        df_all.loc[mask_disc, "p_two_firth_holm"] = disc_rows["p_two_firth_holm"].values

    cols_order = [
        "analysis", "comparison", "sv_type", "bin_type", "model", "exposure",
        "n_case", "n_ctrl", "n_clusters", "carrier_case", "carrier_ctrl",
        "beta_gee", "se_gee", "or_gee", "or_lo95_gee", "or_hi95_gee",
        "p_two_gee", "p_two_gee_bh", "p_two_gee_holm", "fit_status_gee",
        "beta_firth", "se_firth", "or_firth", "or_lo95_firth", "or_hi95_firth",
        "p_two_firth", "p_two_firth_bh", "p_two_firth_holm", "fit_status_firth",
    ]
    df_all = df_all[[c for c in cols_order if c in df_all.columns]]
    df_all.to_csv(str(OUT_TSV), sep="\t", index=False)
    log(f"")
    log(f"Saved main result TSV: {OUT_TSV}  ({len(df_all)} rows)")

    df_sum = df_all[mask_disc].copy()
    df_sum.to_csv(str(OUT_SUMMARY), sep="\t", index=False)
    log(f"Saved summary TSV: {OUT_SUMMARY}  ({len(df_sum)} rows)")

    # -- Preview --
    log(f"")
    log("=" * 72)
    log("PREVIEW: ASD_vs_unaffected_sibling × DEL × 4 discovery-positive")
    log("=" * 72)
    for ana, bt in DISCOVERY_POSITIVE_KEYS:
        sel = df_all[(df_all["analysis"] == ana) & (df_all["bin_type"] == bt)]
        if len(sel) == 0:
            log(f"  {ana}/{bt}: NOT FOUND")
            continue
        r = sel.iloc[0]
        log(
            f"  {ana}/{bt}:  "
            f"GEE OR={r['or_gee']:.3f} P={r['p_two_gee']:.4g} "
            f"p_BH={r['p_two_gee_bh']:.4g} p_Holm={r['p_two_gee_holm']:.4g}  "
            f"Firth OR={r['or_firth']:.3f} P={r['p_two_firth']:.4g} "
            f"p_BH={r['p_two_firth_bh']:.4g} p_Holm={r['p_two_firth_holm']:.4g}  "
            f"carriers={int(r['carrier_case'])}/{int(r['carrier_ctrl'])}"
        )

    log(f"")
    log("PREVIEW: other rows (context)")
    for ana, bt in [
        ("A1_separate", "exonoverlap"),
        ("A1_separate", "all"),
        ("A2_joint", "exonoverlap_in_joint"),
    ]:
        sel = df_all[(df_all["analysis"] == ana) & (df_all["bin_type"] == bt)]
        if len(sel) == 0:
            continue
        r = sel.iloc[0]
        log(
            f"  {ana}/{bt}:  "
            f"GEE OR={r['or_gee']:.3f} P={r['p_two_gee']:.4g}  "
            f"Firth OR={r['or_firth']:.3f} P={r['p_two_firth']:.4g}  "
            f"carriers={int(r['carrier_case'])}/{int(r['carrier_ctrl'])}"
        )

    log(f"")
    log("=" * 72)
    log("Script 13 v4 DONE: 42_compute_mssng_exon_exclusion_v4.py")
    log("=" * 72)


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()
