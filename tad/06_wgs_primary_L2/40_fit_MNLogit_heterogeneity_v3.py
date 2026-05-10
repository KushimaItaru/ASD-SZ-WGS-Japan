#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ファイル名: 40_fit_MNLogit_heterogeneity_v3.py
# 処理内容 (v2 -> v3 修正点):
#   - [path fix only] OUTDIR を OUT_06_B_PRIME_L2 (= .../06_wgs_primary_L2/output_v3)
#     の直下 join から OUT_06_B_PRIME_L2.parent (= .../06_wgs_primary_L2) の join に
#     変更。v2 では OUTDIR が `.../06_wgs_primary_L2/output_v3/output_v4/
#     output_MNLogit_v2/` のように output_v3 が二重化した不整合パスに書き出していた
#     (sbatch の -o/-e は `.../output_v4/output_MNLogit_v2/` のため log と結果 TSV が
#     別ディレクトリに分離)。v3 で `.../06_wgs_primary_L2/output_v4/output_MNLogit_v3/`
#     に正しく揃うよう修正。
#   - [naming] 出力ファイル suffix: _v2 -> _v3。OUTDIR 名: output_MNLogit_v2
#     -> output_MNLogit_v3。QC JSON の script/version フィールドも v3 に更新し、
#     v2_changes -> v3_changes に変更。
#   - 統計・モデルロジック (10 L2 × 2 SV = 20 individual fit + Stouffer/Sign × 2)
#     および v1 -> v2 で入れた全てのバグ修正 (sex_col fix / tuple unpack / converged
#     要件 / NA 行保持 / Healthy=0 assertion / Stouffer wording / Cov warn / QC log)
#     は v2 と完全一致。数値結果は v2 とビット単位で同一になる想定。
#
# 処理内容 (v1 -> v2 修正点, 継続記録):
#   - [致命的バグ修正] load_burden_table() の戻り値に sex_col を追加。
#     v1 では sex_col を検出しながら return しておらず、main() が sample_col を
#     sex_col として run_one_class_one_svt() に渡していたため、Sex_numeric では
#     なく sample_id が共変量に入る誤実装だった。v2 で完全修正。
#   - [致命的バグ修正] normalize_l2_burden_columns() の戻り値 interface を修正。
#     common/naming_v1.py の実装は (df, rename_map) のタプルを返す仕様だが、
#     v1 は単一戻り値を仮定していた。v2 で tuple unpack に修正。
#   - [robustness] fit_mnlogit_with_fallback() で res.mle_retvals['converged']
#     == True を要求。未収束 optimizer はスキップし、次の optimizer に fallback。
#     どれも未収束なら (None, 'failed_not_converged') を返す。
#   - [transparency] 失敗した (bclass, svt) ペアも結果 TSV に NA 行として保持。
#     これにより想定 20 individual 行を必ず保証し、"なぜ 22 行?" を避ける。
#   - [assertion] Outcome encoding (Healthy=0, ASD=1, SZ=2) と MNLogit params
#     の column 順 (0=ASD, 1=SZ) を assert + log で明示化。
#   - [wording] docstring / log の "Brown's method" を削除し、
#     "correlation-adjusted Stouffer" に統一。
#   - [warning] Cov(β_ASD, β_SZ) < 0 の場合に WARN log を出す
#     (共有 control の MNLogit では正の共分散が期待される)。
#   - [transparency] 起動時 log に "Step 05 v3 post-QC (top-1% sample-level QC
#     already applied)" を明示。
#   - 入力: Step 05 v3 (Pattern A = main analysis) の post-QC burden table
#     (v1 / v2 / v3 で同一、上流変更なし)
#   - 出力: output_MNLogit_v3/
#       MNLogit_heterogeneity_results_v3.tsv  (20 individual + 4 global = 24 行)
#       MNLogit_heterogeneity_summary_v3.tsv  (sv 別 global 集計)
#       MNLogit_heterogeneity_qc_v3.json      (実行サマリ + 設定 + failed class list)
#   - 実行時間記録あり

import os
import re
import sys
import json
import math
import time
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats as scipy_stats

# v2: tad04212026 pipeline - centralized paths via common.paths_v1
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths_v1 import (
    OUT_06_B_PRIME_L2,
    ensure_output_dirs,
)
from common.naming_v1 import normalize_l2_burden_columns

warnings.filterwarnings("ignore")

# ============================================================
# CONFIG
# ============================================================

# Step 05 v3 (post-QC, Pattern A = main analysis) burden table を
# hardcoded override で参照 (paths_v1.F_05_SAMPLE_BURDEN_L2 は output_v2 を
# 指すため override)。
_BURDEN_TABLE = os.environ.get(
    "BURDEN",
    "/lustre12/home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v3/"
    "sample_burden_L2_and_specificity_v3.tsv"
)

_OUTDIR = os.environ.get(
    "OUTDIR",
    # v3: OUT_06_B_PRIME_L2 は既に ".../06_wgs_primary_L2/output_v3" を指しているため、
    # 直下に "output_v4/..." を join すると ".../output_v3/output_v4/..." の二重化に
    # なる (v2 のバグ)。OUT_06_B_PRIME_L2.parent で ".../06_wgs_primary_L2" に戻し、
    # sbatch の -o/-e と完全一致する ".../06_wgs_primary_L2/output_v4/output_MNLogit_v3/"
    # を生成する。
    str(OUT_06_B_PRIME_L2.parent / "output_v4" / "output_MNLogit_v3")
)

# 10 PRIMARY L2 classes (HPC_Astro 除外, snake_case after normalize_l2_burden_columns)
PRIMARY_BOUNDARY_CLASSES = [
    "hpc_exc_ca", "hpc_exc_dg", "hpc_exc_ent",
    "hpc_inh_cge", "hpc_inh_mge",
    "pfc_astro", "pfc_exc_dl", "pfc_exc_ul",
    "pfc_inh_cge", "pfc_inh_mge",
]

SV_TYPES = ["DEL", "DUP"]

# Sample / outcome encoding (Healthy=0 が MNLogit の reference になることを保証)
DIAGNOSIS_CODES = {"Healthy": 0, "ASD": 1, "SZ": 2}
KEEP_DX = set(DIAGNOSIS_CODES.keys())

# fit min sample guard
MIN_PER_GROUP = 10  # ASD/SZ/Healthy 各群で 10 例以上必要 (sparsity 保護)

# Output filenames (v3 suffix)
_RESULTS_TSV = "MNLogit_heterogeneity_results_v3.tsv"
_SUMMARY_TSV = "MNLogit_heterogeneity_summary_v3.tsv"
_QC_JSON = "MNLogit_heterogeneity_qc_v3.json"

# ============================================================
# UTILITY
# ============================================================

def now_str() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")

def log(msg: str) -> None:
    print(f"[{now_str()}] {msg}", file=sys.stderr, flush=True)

def detect_column(df: pd.DataFrame, candidates: List[str], required: bool = True) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        raise KeyError(f"None of the columns {candidates} found in dataframe (cols: {list(df.columns)[:30]}...)")
    return None


def make_na_individual_row(bclass: str, svt: str, status: str) -> Dict:
    """Failed individual fit 用の NA 行。想定 20 行を保証するために使用。"""
    return {
        "boundary_class": bclass,
        "sv_type": svt,
        "exposure_col": f"n_boundary_{bclass}_{svt}",
        "model": "MNLogit",
        "fit_method": status,         # 'skipped_no_column' / 'skipped_min_n' / 'failed_not_converged' 等
        "n_control": np.nan,
        "n_ASD": np.nan,
        "n_SZ": np.nan,
        "beta_ASD": np.nan,
        "SE_ASD": np.nan,
        "beta_SZ": np.nan,
        "SE_SZ": np.nan,
        "Cov_ASD_SZ": np.nan,
        "diff_beta": np.nan,
        "SE_diff": np.nan,
        "z_het": np.nan,
        "P_het_two_sided": np.nan,
        "OR_diff": np.nan,
        "CI95_lower_OR_diff": np.nan,
        "CI95_upper_OR_diff": np.nan,
        "converged": False,
    }


# ============================================================
# DATA LOADING
# ============================================================

def load_burden_table(
    burden_path: str,
) -> Tuple[pd.DataFrame, str, str, str, List[str]]:
    """Load Step 05 v3 burden table.

    Returns
    -------
    (df_normalized, sample_col, diag_col, sex_col, pc_cols)
        v2 修正: sex_col を戻り値に追加 (v1 の致命的バグ修正)。

    列名は normalize_l2_burden_columns で snake_case 化:
      n_boundary_HPC_Exc-CA_DEL -> n_boundary_hpc_exc_ca_DEL
      n_events_HPC_Exc-CA_DEL   -> n_events_hpc_exc_ca_DEL
      carrier_boundary_HPC_Exc-CA_DEL -> carrier_boundary_hpc_exc_ca_DEL
    """
    log(f"Reading burden table: {burden_path}")
    df = pd.read_csv(burden_path, sep="\t", low_memory=False)
    log(f"Input: {df.shape[0]} rows x {df.shape[1]} cols")

    # v2 修正: normalize_l2_burden_columns() は (df, rename_map) のタプルを返す
    df, rename_map = normalize_l2_burden_columns(df)
    log(f"L2 column normalize: {len(rename_map)} columns renamed (CamelCase-dash -> snake_case)")

    # Detect required columns
    sample_col = detect_column(df, ["sample_id", "Sample_ID", "IID"])
    diag_col = detect_column(df, ["Diagnosis", "diagnosis", "DX"])
    sex_col = detect_column(df, ["Sex_numeric", "sex_numeric"])
    log(f"  sample_col = {sample_col}")
    log(f"  diag_col   = {diag_col}")
    log(f"  sex_col    = {sex_col}")

    pc_cols = sorted(
        [c for c in df.columns if re.match(r"^PC\d+$", c)],
        key=lambda c: int(re.sub(r"^PC", "", c))
    )[:10]
    if len(pc_cols) < 10:
        raise ValueError(f"Need PC1..PC10; found only {len(pc_cols)} PC columns")

    df[sex_col] = pd.to_numeric(df[sex_col], errors="coerce")
    for pc in pc_cols:
        df[pc] = pd.to_numeric(df[pc], errors="coerce")
    df[diag_col] = df[diag_col].astype(str)

    df = df.loc[df[diag_col].isin(KEEP_DX)].copy()
    log(f"After Diagnosis filter ({sorted(KEEP_DX)}):")
    log(f"  Diagnosis counts:\n{df[diag_col].value_counts().to_string()}")

    # B' model required covariates (svt 別)
    required_covars = [
        "log1p_total_del_bases", "log1p_total_dup_bases",
        "log1p_total_gene_DEL", "log1p_total_gene_DUP",
    ]
    missing = [c for c in required_covars if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required covariate columns in burden table: {missing}")

    return df, sample_col, diag_col, sex_col, pc_cols


def get_b_prime_covariates(svt: str, sex_col: str, pc_cols: List[str]) -> List[str]:
    svt_lower = svt.lower()
    return [sex_col, f"log1p_total_{svt_lower}_bases"] + pc_cols + [f"log1p_total_gene_{svt}"]


# ============================================================
# MNLogit FIT
# ============================================================

def fit_mnlogit_with_fallback(y: pd.Series, X: pd.DataFrame) -> Tuple[Optional[object], str]:
    """MNLogit を newton -> bfgs -> lbfgs の順で fit。

    v2 修正: res.mle_retvals['converged'] == True を要求。未収束の場合は次の
    optimizer に fallback。全 optimizer が未収束なら (None, 'failed_not_converged')。
    cov_params() 取得失敗なら (None, 'failed_cov_unavailable')。
    """
    last_status = "failed"
    for method in ("newton", "bfgs", "lbfgs"):
        try:
            res = sm.MNLogit(y, X).fit(disp=0, maxiter=200, method=method)
            converged = bool(res.mle_retvals.get("converged", False))
            if not converged:
                last_status = f"not_converged_{method}"
                continue
            _ = res.cov_params()  # availability check
            return res, method
        except Exception:
            last_status = f"exception_{method}"
            continue
    return None, last_status


def extract_beta_se_cov(res, exp_col: str, X: pd.DataFrame) -> Tuple[float, float, float, float, float]:
    """MNLogit fit 結果から β_ASD, β_SZ, Var(β_ASD), Var(β_SZ),
    Cov(β_ASD, β_SZ) を抽出。

    前提 (load_burden_table + DIAGNOSIS_CODES で保証):
      - Outcome = {Healthy=0 (reference), ASD=1, SZ=2}
      - MNLogit.params の columns は [0, 1] = [vs_0(ASD), vs_0(SZ)]
      - cov_params() は MultiIndex で ('1', col) / ('2', col) でアクセス可能な
        場合と、統計 version によっては positional (int) の場合がある。両者
        フォールバック。
    """
    # baseline / column order の動的 assert
    param_cols = list(res.params.columns)
    if len(param_cols) != 2:
        raise RuntimeError(
            f"MNLogit params has {len(param_cols)} columns (expected 2 for 3-class outcome). "
            f"columns={param_cols}"
        )

    beta_asd = float(res.params.loc[exp_col, param_cols[0]])  # col 0 = ASD (vs Healthy)
    beta_sz = float(res.params.loc[exp_col, param_cols[1]])   # col 1 = SZ  (vs Healthy)

    cov_full = res.cov_params()

    try:
        var_b_asd = float(cov_full.loc[("1", exp_col), ("1", exp_col)])
        var_b_sz = float(cov_full.loc[("2", exp_col), ("2", exp_col)])
        cov_b_asd_sz = float(cov_full.loc[("1", exp_col), ("2", exp_col)])
    except (KeyError, TypeError):
        # Position-based fallback
        n_params = X.shape[1]
        idx_exp = list(X.columns).index(exp_col)
        var_b_asd = float(cov_full.iloc[idx_exp, idx_exp])
        var_b_sz = float(cov_full.iloc[idx_exp + n_params, idx_exp + n_params])
        cov_b_asd_sz = float(cov_full.iloc[idx_exp, idx_exp + n_params])

    return beta_asd, beta_sz, var_b_asd, var_b_sz, cov_b_asd_sz


def run_one_class_one_svt(
    df: pd.DataFrame,
    bclass: str,
    svt: str,
    diag_col: str,
    sex_col: str,
    pc_cols: List[str],
) -> Dict:
    """1 boundary class × 1 SV type の MNLogit fit + heterogeneity test。

    v2 修正: 失敗時も NA 行 (fit_method に理由を格納) を返す (None は返さない)。
    これにより想定 20 行の個別結果を必ず保証する。
    """
    exp_col = f"n_boundary_{bclass}_{svt}"
    if exp_col not in df.columns:
        log(f"    [{bclass} {svt}] SKIP: column {exp_col} not found")
        return make_na_individual_row(bclass, svt, "skipped_no_column")

    covars = get_b_prime_covariates(svt, sex_col, pc_cols)
    needed = [diag_col, exp_col] + covars
    missing = [c for c in needed if c not in df.columns]
    if missing:
        log(f"    [{bclass} {svt}] SKIP: missing columns {missing}")
        return make_na_individual_row(bclass, svt, "skipped_missing_columns")

    sub = df[needed].dropna().copy()

    sub["_Outcome"] = sub[diag_col].map(DIAGNOSIS_CODES)
    sub = sub.loc[sub["_Outcome"].isin([0, 1, 2])].copy()

    n_ctrl = int((sub["_Outcome"] == 0).sum())
    n_asd = int((sub["_Outcome"] == 1).sum())
    n_sz = int((sub["_Outcome"] == 2).sum())

    if n_asd < MIN_PER_GROUP or n_sz < MIN_PER_GROUP or n_ctrl < MIN_PER_GROUP:
        log(f"    [{bclass} {svt}] SKIP: too few samples (ASD={n_asd}, SZ={n_sz}, ctrl={n_ctrl})")
        row = make_na_individual_row(bclass, svt, "skipped_min_n")
        row["n_control"] = n_ctrl
        row["n_ASD"] = n_asd
        row["n_SZ"] = n_sz
        return row

    X = sub[[exp_col] + covars].astype(float)
    X = sm.add_constant(X)
    y = sub["_Outcome"].astype(int)

    # v2 assertion: y が [0,1,2] の 3 値に揃っていること (Healthy=0 が reference)
    y_unique = sorted(y.unique().tolist())
    if y_unique != [0, 1, 2]:
        log(f"    [{bclass} {svt}] FAIL: y_unique={y_unique} (expected [0,1,2])")
        row = make_na_individual_row(bclass, svt, "failed_outcome_encoding")
        row["n_control"] = n_ctrl
        row["n_ASD"] = n_asd
        row["n_SZ"] = n_sz
        return row

    res, method = fit_mnlogit_with_fallback(y, X)
    if res is None:
        log(f"    [{bclass} {svt}] FAIL: MNLogit did not converge ({method})")
        row = make_na_individual_row(bclass, svt, f"failed_{method}")
        row["n_control"] = n_ctrl
        row["n_ASD"] = n_asd
        row["n_SZ"] = n_sz
        return row

    try:
        beta_asd, beta_sz, var_asd, var_sz, cov_asd_sz = extract_beta_se_cov(res, exp_col, X)
    except Exception as e:
        log(f"    [{bclass} {svt}] FAIL: could not extract beta/cov: {e}")
        row = make_na_individual_row(bclass, svt, "failed_extract")
        row["n_control"] = n_ctrl
        row["n_ASD"] = n_asd
        row["n_SZ"] = n_sz
        return row

    # v2 追加: 共有 control の MNLogit は通常 Cov(β_ASD, β_SZ) > 0 が期待される
    if cov_asd_sz < 0:
        log(
            f"    [{bclass} {svt}] WARN: negative Cov(β_ASD, β_SZ) = {cov_asd_sz:.4g} "
            f"(expected > 0 under shared control; possible index mismatch in extract_beta_se_cov)"
        )

    diff = beta_asd - beta_sz
    var_diff = var_asd + var_sz - 2.0 * cov_asd_sz
    if var_diff <= 0 or not np.isfinite(var_diff):
        log(f"    [{bclass} {svt}] FAIL: invalid Var(β_ASD-β_SZ) = {var_diff}")
        row = make_na_individual_row(bclass, svt, "failed_var_diff_nonpositive")
        row["n_control"] = n_ctrl
        row["n_ASD"] = n_asd
        row["n_SZ"] = n_sz
        row["beta_ASD"] = beta_asd
        row["beta_SZ"] = beta_sz
        row["Cov_ASD_SZ"] = cov_asd_sz
        return row

    se_diff = math.sqrt(var_diff)
    z_het = diff / se_diff
    p_het_two = float(2.0 * scipy_stats.norm.sf(abs(z_het)))

    se_asd = math.sqrt(var_asd) if var_asd > 0 else float("nan")
    se_sz = math.sqrt(var_sz) if var_sz > 0 else float("nan")

    return {
        "boundary_class": bclass,
        "sv_type": svt,
        "exposure_col": exp_col,
        "model": "MNLogit",
        "fit_method": method,
        "n_control": n_ctrl,
        "n_ASD": n_asd,
        "n_SZ": n_sz,
        "beta_ASD": beta_asd,
        "SE_ASD": se_asd,
        "beta_SZ": beta_sz,
        "SE_SZ": se_sz,
        "Cov_ASD_SZ": cov_asd_sz,
        "diff_beta": diff,
        "SE_diff": se_diff,
        "z_het": z_het,
        "P_het_two_sided": p_het_two,
        "OR_diff": float(np.exp(diff)),
        "CI95_lower_OR_diff": float(np.exp(diff - 1.96 * se_diff)),
        "CI95_upper_OR_diff": float(np.exp(diff + 1.96 * se_diff)),
        "converged": True,
    }


# ============================================================
# GLOBAL TESTS (per SV type, across 10 classes)
# ============================================================

def run_global_tests(
    individual_rows: List[Dict],
    df: pd.DataFrame,
    svt: str,
) -> Tuple[Optional[Dict], Optional[Dict]]:
    """1 sv 型に対する Global tests:
       (a) Correlation-adjusted Stouffer Z (one-sided, alternative ASD>SZ)
       (b) Sign test (one-sided binomial)
    返り値: (stouffer_dict, sign_dict)

    v2 修正: "Brown's method" の呼称を削除し、
    "correlation-adjusted Stouffer" に統一 (厳密には exposure columns の Pearson
    correlation を用いた Lipták-Stouffer 系の補正であり、Brown 1975 の p-value
    合成とは異なる)。
    """
    # converged 個別 fit だけを対象にする (NA 行は除外)
    rows_svt = [
        r for r in individual_rows
        if r["sv_type"] == svt and bool(r.get("converged", False))
    ]
    if len(rows_svt) < 2:
        log(f"  [global {svt}] SKIP: only {len(rows_svt)} converged class results")
        return None, None

    z_arr = np.array([r["z_het"] for r in rows_svt])
    classes_in = [r["boundary_class"] for r in rows_svt]
    M = len(z_arr)
    sum_z = float(np.sum(z_arr))
    n_positive = int(np.sum(z_arr > 0))

    # ---- Correlation-adjusted Stouffer ----
    # exposure 列の相関を取り、Var = M + 2*sum_{i<j} rho_ij で補正。
    # これは "correlation-adjusted Stouffer" で、厳密な Brown's method とは異なる。
    var_adj = float(M)  # default = independence
    burden_cols = [f"n_boundary_{c}_{svt}" for c in classes_in]
    existing = [c for c in burden_cols if c in df.columns]
    mean_rho = 0.0
    if len(existing) >= 2:
        try:
            corr_mat = df[existing].corr().values
            sum_rho = 0.0
            n_pair = 0
            for i in range(len(existing)):
                for j in range(i + 1, len(existing)):
                    rho = corr_mat[i, j]
                    if np.isfinite(rho):
                        sum_rho += float(rho)
                        n_pair += 1
            if n_pair > 0:
                mean_rho = sum_rho / n_pair
            var_adj = float(M) + 2.0 * sum_rho
            if var_adj <= 0:
                log(f"  [global {svt}] WARNING: Var_adj={var_adj:.3f} <=0, fallback to M={M}")
                var_adj = float(M)
        except Exception as e:
            log(f"  [global {svt}] WARNING: correlation computation failed: {e}")

    z_adj = sum_z / math.sqrt(var_adj)
    p_adj_one = float(scipy_stats.norm.sf(z_adj))  # one-sided, alternative ASD>SZ

    stouffer_row = {
        "boundary_class": "GLOBAL_corr_stouffer",
        "sv_type": svt,
        "exposure_col": f"n_boundary_GLOBAL_{svt}",
        "model": "correlation_adjusted_Stouffer",
        "fit_method": "n/a",
        "n_control": np.nan,
        "n_ASD": np.nan,
        "n_SZ": np.nan,
        "beta_ASD": np.nan,
        "SE_ASD": np.nan,
        "beta_SZ": np.nan,
        "SE_SZ": np.nan,
        "Cov_ASD_SZ": np.nan,
        "diff_beta": np.nan,
        "SE_diff": math.sqrt(var_adj),
        "z_het": z_adj,
        "P_het_two_sided": np.nan,
        "OR_diff": np.nan,
        "CI95_lower_OR_diff": np.nan,
        "CI95_upper_OR_diff": np.nan,
        "converged": True,
        "global_M": M,
        "global_sum_z": sum_z,
        "global_var_adj": var_adj,
        "global_mean_rho": mean_rho,
        "global_n_positive": n_positive,
        "P_global_one_sided": p_adj_one,
    }

    # ---- Sign test (one-sided binomial) ----
    p_sign_one = float(
        scipy_stats.binomtest(n_positive, M, 0.5, alternative="greater").pvalue
    )

    sign_row = {
        "boundary_class": "GLOBAL_sign_test",
        "sv_type": svt,
        "exposure_col": f"n_boundary_GLOBAL_{svt}",
        "model": "binomial_sign_test",
        "fit_method": "n/a",
        "n_control": np.nan,
        "n_ASD": np.nan,
        "n_SZ": np.nan,
        "beta_ASD": np.nan,
        "SE_ASD": np.nan,
        "beta_SZ": np.nan,
        "SE_SZ": np.nan,
        "Cov_ASD_SZ": np.nan,
        "diff_beta": np.nan,
        "SE_diff": np.nan,
        "z_het": np.nan,
        "P_het_two_sided": np.nan,
        "OR_diff": np.nan,
        "CI95_lower_OR_diff": np.nan,
        "CI95_upper_OR_diff": np.nan,
        "converged": True,
        "global_M": M,
        "global_sum_z": np.nan,
        "global_var_adj": np.nan,
        "global_mean_rho": np.nan,
        "global_n_positive": n_positive,
        "P_global_one_sided": p_sign_one,
    }

    log(
        f"  [global {svt}] M={M}, n_positive={n_positive}/{M}, "
        f"sum_z={sum_z:.3f}, Var_adj={var_adj:.2f}, "
        f"Stouffer_Z={z_adj:.3f} P={p_adj_one:.3e}, sign P={p_sign_one:.3e}"
    )
    return stouffer_row, sign_row


# ============================================================
# MAIN
# ============================================================

def main() -> int:
    t0 = time.time()
    log("=" * 60)
    log("40_fit_MNLogit_heterogeneity_v3.py")
    log("  ASD vs SZ heterogeneity test (Multinomial Logistic)")
    log("=" * 60)
    log(f"BURDEN: {_BURDEN_TABLE}")
    log(f"OUTDIR: {_OUTDIR}")
    log("Input burden: Step 05 v3 post-QC (top-1% sample-level CNV-count QC "
        "already applied at Step 04 v10 / Step 05 v3).")

    ensure_output_dirs()

    outdir = Path(_OUTDIR)
    outdir.mkdir(parents=True, exist_ok=True)

    # v2 修正: load_burden_table() の戻り値に sex_col を追加
    df, sample_col, diag_col, sex_col, pc_cols = load_burden_table(_BURDEN_TABLE)
    log(f"PC columns used: {pc_cols}")

    # v2 assertion: DIAGNOSIS_CODES が Healthy=0 を baseline にしているか
    if DIAGNOSIS_CODES.get("Healthy") != 0:
        raise RuntimeError(
            f"DIAGNOSIS_CODES must set Healthy=0 for MNLogit baseline. "
            f"Got: {DIAGNOSIS_CODES}"
        )
    log(f"Outcome encoding: {DIAGNOSIS_CODES} (Healthy=0 is MNLogit reference)")
    log("MNLogit params column convention: col 0 = ASD (vs Healthy), col 1 = SZ (vs Healthy)")

    individual_rows: List[Dict] = []
    log("")
    log("---- Individual class × SV type fits ----")
    for svt in SV_TYPES:
        for bclass in PRIMARY_BOUNDARY_CLASSES:
            # v2 修正: sample_col ではなく sex_col を渡す (致命的バグ修正)
            row = run_one_class_one_svt(
                df, bclass, svt, diag_col, sex_col, pc_cols
            )
            individual_rows.append(row)
            if bool(row.get("converged", False)):
                log(
                    f"  [{bclass} {svt}] z_het={row['z_het']:.3f} "
                    f"P={row['P_het_two_sided']:.3e} (method={row['fit_method']})"
                )
            else:
                log(f"  [{bclass} {svt}] row preserved with fit_method={row['fit_method']}")

    log("")
    log("---- Global tests (correlation-adjusted Stouffer + Sign) per SV type ----")
    global_rows: List[Dict] = []
    summary_rows: List[Dict] = []
    for svt in SV_TYPES:
        st_row, sg_row = run_global_tests(individual_rows, df, svt)
        if st_row is not None:
            global_rows.append(st_row)
            summary_rows.append({
                "sv_type": svt,
                "test": "correlation_adjusted_Stouffer",
                "M": st_row["global_M"],
                "n_positive_z": st_row["global_n_positive"],
                "sum_z": st_row["global_sum_z"],
                "Var_adj": st_row["global_var_adj"],
                "mean_rho": st_row["global_mean_rho"],
                "Z": st_row["z_het"],
                "P_one_sided_ASD_gt_SZ": st_row["P_global_one_sided"],
            })
        if sg_row is not None:
            global_rows.append(sg_row)
            summary_rows.append({
                "sv_type": svt,
                "test": "binomial_sign_test",
                "M": sg_row["global_M"],
                "n_positive_z": sg_row["global_n_positive"],
                "sum_z": np.nan,
                "Var_adj": np.nan,
                "mean_rho": np.nan,
                "Z": np.nan,
                "P_one_sided_ASD_gt_SZ": sg_row["P_global_one_sided"],
            })

    # failed class 情報 (QC json 用)
    failed_classes = [
        {"boundary_class": r["boundary_class"],
         "sv_type": r["sv_type"],
         "fit_method": r["fit_method"]}
        for r in individual_rows if not bool(r.get("converged", False))
    ]
    n_converged = sum(1 for r in individual_rows if bool(r.get("converged", False)))
    log(f"Converged individual fits: {n_converged}/{len(individual_rows)}")

    # Output 1: full results (individual + global)
    out_full = pd.DataFrame(individual_rows + global_rows)
    out_full_path = outdir / _RESULTS_TSV
    out_full.to_csv(out_full_path, sep="\t", index=False)
    log(f"Wrote: {out_full_path}  ({out_full.shape[0]} rows × {out_full.shape[1]} cols)")

    # Output 2: summary (sv 別 global 集計)
    out_summary = pd.DataFrame(summary_rows)
    out_summary_path = outdir / _SUMMARY_TSV
    out_summary.to_csv(out_summary_path, sep="\t", index=False)
    log(f"Wrote: {out_summary_path}  ({out_summary.shape[0]} rows × {out_summary.shape[1]} cols)")

    # Output 3: QC JSON
    elapsed = time.time() - t0
    qc = {
        "script": "40_fit_MNLogit_heterogeneity_v3.py",
        "version": "v3",
        "v3_changes": [
            "PATH FIX ONLY: OUTDIR uses OUT_06_B_PRIME_L2.parent (remove output_v3 duplication)",
            "NAMING: output filename suffix _v2 -> _v3, dir output_MNLogit_v2 -> output_MNLogit_v3",
            "Statistical/model logic identical to v2 (numerical results bit-identical expected)",
        ],
        "v2_changes_inherited": [
            "FIX: load_burden_table returns sex_col (v1 bug: sample_id was used as sex)",
            "FIX: normalize_l2_burden_columns tuple unpack (v1 bug: df became tuple)",
            "ROBUST: fit_mnlogit_with_fallback requires converged==True",
            "TRANSPARENCY: failed individual fits kept as NA rows (n=20 guaranteed)",
            "ASSERTION: Healthy=0 baseline + MNLogit params column order",
            "WORDING: 'Brown's method' -> 'correlation-adjusted Stouffer'",
            "WARN: negative Cov(beta_ASD, beta_SZ)",
            "LOG: top-1% QC already applied upstream",
        ],
        "timestamp": now_str(),
        "elapsed_sec": round(elapsed, 2),
        "elapsed_min": round(elapsed / 60.0, 3),
        "input": {
            "burden_table": _BURDEN_TABLE,
            "note": "Step 05 v3 post-QC (top-1% sample-level CNV-count QC already applied)",
        },
        "config": {
            "primary_boundary_classes": PRIMARY_BOUNDARY_CLASSES,
            "sv_types": SV_TYPES,
            "diagnosis_codes": DIAGNOSIS_CODES,
            "min_per_group": MIN_PER_GROUP,
            "model": "MNLogit (Outcome ~ exposure + Sex + PC1..10 + log1p_total_{sv}_bases + log1p_total_gene_{SV})",
            "fit_methods_attempted": ["newton", "bfgs", "lbfgs"],
            "global_tests": [
                "Correlation-adjusted Stouffer Z (one-sided ASD>SZ)",
                "One-sided binomial sign test",
            ],
        },
        "n_individual_rows": len(individual_rows),
        "n_individual_converged": n_converged,
        "n_individual_failed": len(individual_rows) - n_converged,
        "failed_classes": failed_classes,
        "n_global_rows": len(global_rows),
        "n_total_rows": len(individual_rows) + len(global_rows),
    }
    qc_path = outdir / _QC_JSON
    with open(qc_path, "w") as f:
        json.dump(qc, f, indent=2, ensure_ascii=False)
    log(f"Wrote: {qc_path}")

    log("")
    log(f"Total elapsed: {elapsed:.1f} sec ({elapsed / 60.0:.2f} min)")
    log("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
