#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ファイル名: tad_replication_meta_analysis_v12.py
# 処理内容:
#   - arrayCGH(v22, tad04212026 08_arraycgh_sample_burden/output_v22) と
#     MSSNG(v18, tad04212026 09_mssng_sample_burden/output_v18) の
#     replication 結果 TSV を読み込む
#   - PRIMARY / SECONDARY Pattern C / negative control を cohort 間で対応付け
#   - Main meta: IVW fixed-effect meta-analysis
#   - Sensitivity meta: weighted Stouffer (signed Z)
#   - PRIMARY main / exploratory / negative control ごとに TSV 出力
#   - 実行時間をログ出力
#
# v11 から v12 への変更点 (2026-04-21 tad04212026 パイプライン移植):
#   1. arrayCGH 入力: v21 -> v22 (tad04212026/08_arraycgh_sample_burden/output_v22)
#      解析ロジックは v21 と同一 (filename と path のみ更新)。
#   2. MSSNG 入力: v17 -> v18 (tad04212026/09_mssng_sample_burden/output_v18)
#      解析ロジックは v17 と同一 (filename と path のみ更新)。
#   3. 出力: /lustre12/home/kushima-pg/tad04212026/09_mssng_sample_burden/meta_v12
#      出力ファイル名サフィックスを _v11 -> _v12 に更新。
#   4. cohort ラベルを "arrayCGH_v22" / "MSSNG_v18_fullfamily" に更新。
#   5. メタ解析ロジック (IVW, Stouffer, analysis_set 判定) は v11 から一切変更なし。
#
# v6 からの変更点 (v7):
#   1. arrayCGH 入力: v15 -> v18 (segdup < 50% + exclusion BED < 50% フィルタ追加版)
#   2. v18 出力は ASD/SZ/heterogeneity の3ファイルに分割 → ASD 結果ファイルを読み込み
#   3. detect_column: platform_meta_p_two の候補名追加
#   4. MSSNG 入力は変更なし (v17)
#
# 主な方針 (v5 以降変更なし):
#   1. Discovery 有意 8 クラス (formal sig8):
#        - PRIMARY_sig8_main: IVW one-sided positive (ivw_p_main)
#        - SECONDARY_C_sig8_main: 同上
#   2. Discovery 非有意 2 クラス (exploratory nonsig2):
#        - IVW two-sided (ivw_p_main)
#   3. Negative control:
#        - IVW two-sided (ivw_p_main)
#   4. Sensitivity:
#        - 全セットに対して weighted Stouffer (stouffer_p_sensitivity)
#
# 入力:
#   - tad_replication_asd_vs_cont_v22.tsv (arrayCGH v22 ASD 結果)
#   - tad_replication_mssng_v18.tsv (MSSNG v18 GEE Independence 結果)
#
# 出力:
#   - meta_results_ivw_v12.tsv
#   - meta_results_stouffer_v12.tsv
#   - meta_primary_sig8_main_v12.tsv
#   - meta_run_config_v12.json
#
# 実行例:
#   python tad_replication_meta_analysis_v12.py
#
# 環境変数で上書き可能:
#   ARRAY_TSV=/path/to/tad_replication_asd_vs_cont_v22.tsv
#   MSSNG_TSV=/path/to/tad_replication_mssng_v18.tsv
#   OUTDIR=/path/to/meta_v12

import os
import sys
import json
import time
from math import sqrt
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.stats import norm, chi2

# ============================================================
# CONFIG
# ============================================================

_BASEDIR_ARRAY = "/lustre12/home/kushima-pg/tad04212026/08_arraycgh_sample_burden/output_v22"
_BASEDIR_MSSNG = "/lustre12/home/kushima-pg/tad04212026/09_mssng_sample_burden/output_v18"
_DEFAULT_ARRAY_TSV = f"{_BASEDIR_ARRAY}/tad_replication_asd_vs_cont_v22.tsv"
_DEFAULT_MSSNG_TSV = f"{_BASEDIR_MSSNG}/tad_replication_mssng_v18.tsv"
_DEFAULT_OUTDIR = "/lustre12/home/kushima-pg/tad04212026/09_mssng_sample_burden/meta_v12"

ARRAY_TSV = os.environ.get("ARRAY_TSV", _DEFAULT_ARRAY_TSV)
MSSNG_TSV = os.environ.get("MSSNG_TSV", _DEFAULT_MSSNG_TSV)
OUTDIR = os.environ.get("OUTDIR", _DEFAULT_OUTDIR)

DISCOVERY_SIG_CLASSES = [
    "hpc_exc_dg", "hpc_exc_ent", "hpc_inh_cge", "hpc_inh_mge",
    "pfc_exc_dl", "pfc_exc_ul", "pfc_inh_cge", "pfc_inh_mge",
]

DISCOVERY_NONSIG_CLASSES = [
    "hpc_exc_ca", "pfc_astro",
]

PRIMARY_BOUNDARY_CLASSES = DISCOVERY_SIG_CLASSES + DISCOVERY_NONSIG_CLASSES

# weighted Stouffer の重み: sqrt(N_eff) — 感度解析で使用
# N_eff = 4 / (1/n_case + 1/n_control)
WEIGHT_MODE = "sqrt_neff"

# p 値の数値安定化
P_MIN = 1e-300
P_MAX = 1.0 - 1e-16

# ============================================================
# LOGGING
# ============================================================

_T0 = time.time()


def log(msg: str) -> None:
    elapsed = time.time() - _T0
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


# ============================================================
# UTIL
# ============================================================


def clamp_p(p: float) -> float:
    if pd.isna(p):
        return np.nan
    return float(min(max(p, P_MIN), P_MAX))


def detect_column(df: pd.DataFrame, candidates: List[str], required: bool = True) -> Optional[str]:
    col_map = {str(c).lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in col_map:
            return col_map[cand.lower()]
    if required:
        raise KeyError(f"Missing column. candidates={candidates}, actual={list(df.columns)}")
    return None


def assert_unique_keys(df: pd.DataFrame, key_cols: List[str], label: str, n_show: int = 10) -> None:
    dup_mask = df.duplicated(subset=key_cols, keep=False)
    if dup_mask.any():
        n_dup = int(dup_mask.sum())
        ex = df.loc[dup_mask, key_cols].head(n_show).to_string(index=False)
        raise ValueError(f"[{label}] duplicate rows found: n={n_dup}\n{ex}")


def calc_neff(n_case: float, n_ctrl: float) -> float:
    if pd.isna(n_case) or pd.isna(n_ctrl):
        return np.nan
    n_case = float(n_case)
    n_ctrl = float(n_ctrl)
    if n_case <= 0 or n_ctrl <= 0:
        return np.nan
    return 4.0 / ((1.0 / n_case) + (1.0 / n_ctrl))


def calc_weight(n_case: float, n_ctrl: float, mode: str = "sqrt_neff") -> float:
    neff = calc_neff(n_case, n_ctrl)
    if pd.isna(neff):
        return np.nan
    if mode == "sqrt_neff":
        return sqrt(neff)
    elif mode == "neff":
        return neff
    elif mode == "sqrt_ncomplete":
        return sqrt(float(n_case) + float(n_ctrl))
    else:
        raise ValueError(f"Unknown weight mode: {mode}")


def signed_z_from_one_sided_p(p_one: float) -> float:
    """
    one-sided p が 'expected positive direction' 用に既に計算済みであることを前提とする。
    p が小さいほど正方向支持が強く、p > 0.5 なら負方向寄りになる。
    """
    p_one = clamp_p(p_one)
    if pd.isna(p_one):
        return np.nan
    return float(norm.isf(p_one))


def signed_z_from_two_sided_p(p_two: float, beta: float) -> float:
    """
    two-sided p と beta の符号から signed Z を作る。
    """
    p_two = clamp_p(p_two)
    if pd.isna(p_two) or pd.isna(beta):
        return np.nan
    sign = 1.0 if beta > 0 else (-1.0 if beta < 0 else 0.0)
    if sign == 0.0:
        return 0.0
    z_abs = float(norm.isf(p_two / 2.0))
    return sign * z_abs


def weighted_stouffer(z_list: List[float], w_list: List[float]) -> Dict[str, float]:
    z = np.asarray(z_list, dtype=float)
    w = np.asarray(w_list, dtype=float)
    ok = np.isfinite(z) & np.isfinite(w) & (w > 0)
    if ok.sum() == 0:
        return {"z_meta": np.nan, "p_one_positive": np.nan, "p_two": np.nan}
    z = z[ok]
    w = w[ok]
    z_meta = float(np.sum(w * z) / np.sqrt(np.sum(w ** 2)))
    p_one_positive = float(norm.sf(z_meta))
    p_two = float(2.0 * norm.sf(abs(z_meta)))
    return {"z_meta": z_meta, "p_one_positive": p_one_positive, "p_two": p_two}


def ivw_fixed_effect(beta_list: List[float], se_list: List[float]) -> Dict[str, float]:
    beta = np.asarray(beta_list, dtype=float)
    se = np.asarray(se_list, dtype=float)
    ok = np.isfinite(beta) & np.isfinite(se) & (se > 0)
    if ok.sum() == 0:
        return {
            "beta_meta": np.nan, "se_meta": np.nan, "or_meta": np.nan,
            "ci_lo": np.nan, "ci_hi": np.nan,
            "z_meta": np.nan, "p_two": np.nan, "p_one_positive": np.nan,
            "q": np.nan, "q_p": np.nan, "i2": np.nan
        }

    beta = beta[ok]
    se = se[ok]
    var = se ** 2
    w = 1.0 / var

    beta_meta = float(np.sum(w * beta) / np.sum(w))
    se_meta = float(np.sqrt(1.0 / np.sum(w)))
    z_meta = float(beta_meta / se_meta)
    p_two = float(2.0 * norm.sf(abs(z_meta)))
    p_one_positive = float(norm.sf(z_meta))
    or_meta = float(np.exp(beta_meta))
    ci_lo = float(np.exp(beta_meta - 1.96 * se_meta))
    ci_hi = float(np.exp(beta_meta + 1.96 * se_meta))

    # Heterogeneity
    q = float(np.sum(w * (beta - beta_meta) ** 2))
    df = max(0, len(beta) - 1)
    q_p = float(chi2.sf(q, df=df)) if df > 0 else np.nan
    i2 = float(max(0.0, (q - df) / q) * 100.0) if (df > 0 and q > 0) else 0.0

    return {
        "beta_meta": beta_meta, "se_meta": se_meta, "or_meta": or_meta,
        "ci_lo": ci_lo, "ci_hi": ci_hi,
        "z_meta": z_meta, "p_two": p_two, "p_one_positive": p_one_positive,
        "q": q, "q_p": q_p, "i2": i2
    }


def classify_analysis_set(analysis_type: str, pattern: str, sv_type: str, boundary_class: str) -> str:
    if analysis_type == "PRIMARY_replication" and pattern == "A" and sv_type == "DEL":
        if boundary_class in DISCOVERY_SIG_CLASSES:
            return "PRIMARY_sig8_main"
        elif boundary_class in DISCOVERY_NONSIG_CLASSES:
            return "PRIMARY_nonsig2_exploratory"
    if analysis_type == "SECONDARY_pattern_c" and pattern == "C" and sv_type == "DEL":
        if boundary_class in DISCOVERY_SIG_CLASSES:
            return "SECONDARY_C_sig8_main"
        elif boundary_class in DISCOVERY_NONSIG_CLASSES:
            return "SECONDARY_C_nonsig2_exploratory"
    if analysis_type == "negative_control":
        if sv_type == "DUP" and boundary_class in PRIMARY_BOUNDARY_CLASSES:
            return f"NEGCTRL_DUP_pattern_{pattern}"
        if sv_type == "DEL" and boundary_class == "static_all":
            return f"NEGCTRL_STATIC_pattern_{pattern}"
    return "OTHER"


# ============================================================
# LOAD / STANDARDIZE
# ============================================================


def load_replication_result(path: str, cohort_name: str) -> pd.DataFrame:
    log(f"Reading {cohort_name}: {path}")
    df = pd.read_csv(path, sep="\t", low_memory=False)

    col_pattern = detect_column(df, ["pattern", "Pattern"])
    col_sv_type = detect_column(df, ["sv_type", "SV_type"])
    col_bclass = detect_column(df, ["boundary_class", "Boundary_class"])
    col_analysis = detect_column(df, ["analysis_type"])
    col_n_case = detect_column(df, ["n_case", "N_case"])
    col_n_ctrl = detect_column(df, ["n_control", "N_control"])
    col_beta = detect_column(df, ["beta", "Beta"])
    col_se = detect_column(df, ["SE", "se"])
    col_or = detect_column(df, ["OR", "or"])
    col_ci_lo = detect_column(df, ["CI_lower", "ci_lower", "CI_lo", "ci_lo"], required=False)
    col_ci_hi = detect_column(df, ["CI_upper", "ci_upper", "CI_hi", "ci_hi"], required=False)
    col_p_two = detect_column(df, ["P_twosided", "p_twosided", "P_two", "p_two"])
    col_p_one = detect_column(df, ["P_onesided", "p_onesided", "P_one", "p_one"], required=False)

    # v6+: detect platform_meta columns (arrayCGH v15/v18/v22)
    col_pm_beta = detect_column(df, ["platform_meta_beta"], required=False)
    col_pm_se = detect_column(df, ["platform_meta_SE", "platform_meta_se"], required=False)
    col_pm_p_main = detect_column(df, ["platform_meta_P_main", "platform_meta_p_main"], required=False)
    # v7+: platform_meta_p_two as candidate (v18 以降の column name)
    col_pm_p_two = detect_column(df, ["platform_meta_P_twosided", "platform_meta_p_twosided", "platform_meta_p_two"], required=False)
    col_pm_p_one = detect_column(df, ["platform_meta_P_one_positive", "platform_meta_p_one_positive"], required=False)

    out = pd.DataFrame({
        "cohort": cohort_name,
        "pattern": df[col_pattern].astype(str),
        "sv_type": df[col_sv_type].astype(str),
        "boundary_class": df[col_bclass].astype(str),
        "analysis_type": df[col_analysis].astype(str),
        "n_case": pd.to_numeric(df[col_n_case], errors="coerce"),
        "n_control": pd.to_numeric(df[col_n_ctrl], errors="coerce"),
        "beta": pd.to_numeric(df[col_beta], errors="coerce"),
        "se": pd.to_numeric(df[col_se], errors="coerce"),
        "or": pd.to_numeric(df[col_or], errors="coerce"),
        "ci_lo": pd.to_numeric(df[col_ci_lo], errors="coerce") if col_ci_lo else np.nan,
        "ci_hi": pd.to_numeric(df[col_ci_hi], errors="coerce") if col_ci_hi else np.nan,
        "p_two": pd.to_numeric(df[col_p_two], errors="coerce"),
        "p_one": pd.to_numeric(df[col_p_one], errors="coerce") if col_p_one else np.nan,
    })

    # v6+: add platform_meta columns if available
    if col_pm_beta:
        out["platform_meta_beta"] = pd.to_numeric(df[col_pm_beta], errors="coerce")
        out["platform_meta_se"] = pd.to_numeric(df[col_pm_se], errors="coerce") if col_pm_se else np.nan
        out["platform_meta_p_main"] = pd.to_numeric(df[col_pm_p_main], errors="coerce") if col_pm_p_main else np.nan
        out["platform_meta_p_two"] = pd.to_numeric(df[col_pm_p_two], errors="coerce") if col_pm_p_two else np.nan
        out["platform_meta_p_one"] = pd.to_numeric(df[col_pm_p_one], errors="coerce") if col_pm_p_one else np.nan
        log(f"  {cohort_name}: platform_meta columns detected and loaded")
    else:
        out["platform_meta_beta"] = np.nan
        out["platform_meta_se"] = np.nan
        out["platform_meta_p_main"] = np.nan
        out["platform_meta_p_two"] = np.nan
        out["platform_meta_p_one"] = np.nan

    out["n_eff"] = [calc_neff(a, b) for a, b in zip(out["n_case"], out["n_control"])]
    out["weight"] = [calc_weight(a, b, mode=WEIGHT_MODE) for a, b in zip(out["n_case"], out["n_control"])]
    out["analysis_set"] = [
        classify_analysis_set(a, p, s, b)
        for a, p, s, b in zip(out["analysis_type"], out["pattern"], out["sv_type"], out["boundary_class"])
    ]

    key_cols = ["pattern", "sv_type", "boundary_class", "analysis_type"]
    assert_unique_keys(out, key_cols, f"{cohort_name}_keys")

    log(f"  {cohort_name}: {len(out)} rows")
    return out


# ============================================================
# META
# ============================================================


def build_pair_table(array_df: pd.DataFrame, mssng_df: pd.DataFrame) -> pd.DataFrame:
    key_cols = ["pattern", "sv_type", "boundary_class", "analysis_type", "analysis_set"]

    arr = array_df.copy().rename(columns={
        "n_case": "array_n_case",
        "n_control": "array_n_control",
        "n_eff": "array_n_eff",
        "weight": "array_weight",
        "beta": "array_pooled_beta",
        "se": "array_pooled_se",
        "or": "array_pooled_or",
        "ci_lo": "array_pooled_ci_lo",
        "ci_hi": "array_pooled_ci_hi",
        "p_two": "array_pooled_p_two",
        "p_one": "array_pooled_p_one",
        "platform_meta_beta": "array_beta",
        "platform_meta_se": "array_se",
        "platform_meta_p_main": "array_p_main",
        "platform_meta_p_two": "array_p_two",
        "platform_meta_p_one": "array_p_one",
    }).drop(columns=["cohort"])

    mss = mssng_df.copy().rename(columns={
        "n_case": "mssng_n_case",
        "n_control": "mssng_n_control",
        "n_eff": "mssng_n_eff",
        "weight": "mssng_weight",
        "beta": "mssng_beta",
        "se": "mssng_se",
        "or": "mssng_or",
        "ci_lo": "mssng_ci_lo",
        "ci_hi": "mssng_ci_hi",
        "p_two": "mssng_p_two",
        "p_one": "mssng_p_one",
    }).drop(columns=["cohort"])

    merged = arr.merge(mss, on=key_cols, how="inner", validate="1:1")
    assert_unique_keys(merged, key_cols, "meta_pair_table")

    # 期待行数チェック: 各 cohort は pattern(A,C) × sv_type(DEL,DUP) × 11 classes = 44 行
    # ただし static_all は DEL のみで DUP なし → 42 行
    expected_n = 42
    if len(merged) != expected_n:
        msg = (f"FATAL: Expected {expected_n} paired rows, but got {len(merged)}. "
               f"Check input TSV alignment.")
        log(msg)
        raise ValueError(msg)

    log(f"Paired rows for meta-analysis: {len(merged)}")
    return merged


def do_meta_for_row(row: pd.Series) -> Dict[str, object]:
    analysis_set = row["analysis_set"]

    # =====================================================
    # IVW fixed-effect (主解析 in v5)
    # =====================================================
    ivw = ivw_fixed_effect(
        [row["array_beta"], row["mssng_beta"]],
        [row["array_se"], row["mssng_se"]]
    )

    # IVW 主解析の p 値: sig8 は one-sided positive、exploratory/negctrl は two-sided
    if analysis_set in {"PRIMARY_sig8_main", "SECONDARY_C_sig8_main"}:
        ivw_p_main = ivw["p_one_positive"]
    else:
        ivw_p_main = ivw["p_two"]

    # =====================================================
    # Stouffer (感度解析 in v5)
    # =====================================================
    if analysis_set in {"PRIMARY_sig8_main", "SECONDARY_C_sig8_main"}:
        z_array = signed_z_from_one_sided_p(row["array_p_one"])
        z_mssng = signed_z_from_one_sided_p(row["mssng_p_one"])
        st = weighted_stouffer(
            [z_array, z_mssng],
            [row["array_weight"], row["mssng_weight"]]
        )
        stouffer_p_sensitivity = st["p_one_positive"]
    else:
        z_array = signed_z_from_two_sided_p(row["array_p_two"], row["array_beta"])
        z_mssng = signed_z_from_two_sided_p(row["mssng_p_two"], row["mssng_beta"])
        st = weighted_stouffer(
            [z_array, z_mssng],
            [row["array_weight"], row["mssng_weight"]]
        )
        stouffer_p_sensitivity = st["p_two"]

    # =====================================================
    # Direction / support flags
    # =====================================================
    same_direction = (
        np.isfinite(row["array_beta"]) and np.isfinite(row["mssng_beta"]) and
        ((row["array_beta"] > 0 and row["mssng_beta"] > 0) or (row["array_beta"] < 0 and row["mssng_beta"] < 0))
    )

    formal_meta_supported = (
        analysis_set == "PRIMARY_sig8_main"
        and same_direction
        and ivw["beta_meta"] > 0
        and ivw_p_main < 0.05
    )

    secondary_meta_supported = (
        analysis_set == "SECONDARY_C_sig8_main"
        and same_direction
        and ivw["beta_meta"] > 0
        and ivw_p_main < 0.05
    )

    exploratory_meta_signal = (
        analysis_set in {"PRIMARY_nonsig2_exploratory", "SECONDARY_C_nonsig2_exploratory"}
        and ivw_p_main < 0.05
    )

    negative_control_nominal = (
        analysis_set.startswith("NEGCTRL_")
        and ivw_p_main < 0.05
    )

    out = {
        "analysis_set": analysis_set,
        "pattern": row["pattern"],
        "sv_type": row["sv_type"],
        "boundary_class": row["boundary_class"],
        "analysis_type": row["analysis_type"],

        "array_n_case": row["array_n_case"],
        "array_n_control": row["array_n_control"],
        "array_n_eff": row["array_n_eff"],
        "array_weight": row["array_weight"],
        "array_beta": row["array_beta"],
        "array_se": row["array_se"],
        "array_or": float(np.exp(row["array_beta"])) if np.isfinite(row["array_beta"]) else np.nan,
        "array_p_two": row["array_p_two"],
        "array_p_one": row["array_p_one"],
        "array_p_main": row["array_p_main"],
        "array_pooled_beta": row["array_pooled_beta"],
        "array_pooled_se": row["array_pooled_se"],
        "array_pooled_or": row["array_pooled_or"],

        "mssng_n_case": row["mssng_n_case"],
        "mssng_n_control": row["mssng_n_control"],
        "mssng_n_eff": row["mssng_n_eff"],
        "mssng_weight": row["mssng_weight"],
        "mssng_beta": row["mssng_beta"],
        "mssng_se": row["mssng_se"],
        "mssng_or": row["mssng_or"],
        "mssng_ci_lo": row["mssng_ci_lo"],
        "mssng_ci_hi": row["mssng_ci_hi"],
        "mssng_p_two": row["mssng_p_two"],
        "mssng_p_one": row["mssng_p_one"],

        "same_direction": same_direction,

        # IVW (主解析)
        "ivw_beta_meta": ivw["beta_meta"],
        "ivw_se_meta": ivw["se_meta"],
        "ivw_or_meta": ivw["or_meta"],
        "ivw_ci_lo": ivw["ci_lo"],
        "ivw_ci_hi": ivw["ci_hi"],
        "ivw_z_meta": ivw["z_meta"],
        "ivw_p_two": ivw["p_two"],
        "ivw_p_one_positive": ivw["p_one_positive"],
        "ivw_p_main": ivw_p_main,
        "ivw_q": ivw["q"],
        "ivw_q_p": ivw["q_p"],
        "ivw_i2": ivw["i2"],

        # Stouffer (感度解析)
        "stouffer_z_meta": st["z_meta"],
        "stouffer_p_one_positive": st["p_one_positive"],
        "stouffer_p_two": st["p_two"],
        "stouffer_p_sensitivity": stouffer_p_sensitivity,

        # Support flags
        "formal_meta_supported": formal_meta_supported,
        "secondary_meta_supported": secondary_meta_supported,
        "exploratory_meta_signal": exploratory_meta_signal,
        "negative_control_nominal": negative_control_nominal,
    }
    return out


# ============================================================
# OUTPUT
# ============================================================


def save_outputs(meta_df: pd.DataFrame, outdir: str) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # IVW 主解析 TSV
    ivw_cols = [
        "analysis_set", "pattern", "sv_type", "boundary_class", "analysis_type",
        "same_direction",
        "array_n_case", "array_n_control", "array_n_eff",
        "array_beta", "array_se", "array_or", "array_p_two", "array_p_one", "array_p_main",
        "array_pooled_beta", "array_pooled_se", "array_pooled_or",
        "mssng_n_case", "mssng_n_control", "mssng_n_eff",
        "mssng_beta", "mssng_se", "mssng_or", "mssng_ci_lo", "mssng_ci_hi", "mssng_p_two", "mssng_p_one",
        "ivw_beta_meta", "ivw_se_meta", "ivw_or_meta", "ivw_ci_lo", "ivw_ci_hi",
        "ivw_z_meta", "ivw_p_two", "ivw_p_one_positive", "ivw_p_main",
        "ivw_q", "ivw_q_p", "ivw_i2",
        "formal_meta_supported", "secondary_meta_supported",
        "exploratory_meta_signal", "negative_control_nominal",
    ]

    # Stouffer 感度解析 TSV
    stouffer_cols = [
        "analysis_set", "pattern", "sv_type", "boundary_class", "analysis_type",
        "same_direction",
        "array_n_case", "array_n_control", "array_n_eff", "array_weight",
        "array_beta", "array_se", "array_or", "array_p_two", "array_p_one", "array_p_main",
        "mssng_n_case", "mssng_n_control", "mssng_n_eff", "mssng_weight",
        "mssng_beta", "mssng_se", "mssng_or", "mssng_p_two", "mssng_p_one",
        "stouffer_z_meta", "stouffer_p_one_positive", "stouffer_p_two", "stouffer_p_sensitivity",
    ]

    ivw_path = outdir / "meta_results_ivw_v12.tsv"
    stouffer_path = outdir / "meta_results_stouffer_v12.tsv"
    main_path = outdir / "meta_primary_sig8_main_v12.tsv"
    config_path = outdir / "meta_run_config_v12.json"

    meta_df[ivw_cols].to_csv(ivw_path, sep="\t", index=False, float_format="%.6g")
    meta_df[stouffer_cols].to_csv(stouffer_path, sep="\t", index=False, float_format="%.6g")

    # メイン結果: IVW 基準でソート
    main_df = meta_df.loc[meta_df["analysis_set"] == "PRIMARY_sig8_main"].copy()
    main_df = main_df.sort_values("ivw_p_main", ascending=True)
    main_sig8_cols = [
        "boundary_class",
        "array_or", "array_p_main",
        "mssng_or", "mssng_p_one",
        "ivw_or_meta", "ivw_ci_lo", "ivw_ci_hi", "ivw_p_main",
        "same_direction", "ivw_i2",
        "stouffer_p_sensitivity",
        "formal_meta_supported",
    ]
    main_df[main_sig8_cols].to_csv(main_path, sep="\t", index=False, float_format="%.6g")

    config = {
        "script": "tad_replication_meta_analysis_v12.py",
        "version": "v12",
        "changes_from_v11": [
            "arrayCGH input updated: v21 -> v22 (tad04212026/08_arraycgh_sample_burden/output_v22)",
            "MSSNG input updated: v17 -> v18 (tad04212026/09_mssng_sample_burden/output_v18)",
            "OUTDIR updated: tad04212026/09_mssng_sample_burden/meta_v12",
            "Output TSV suffixes: _v11 -> _v12",
            "Cohort labels: arrayCGH_v22, MSSNG_v18_fullfamily",
            "Meta logic (IVW, Stouffer, analysis_set classification) UNCHANGED from v11",
        ],
        "changes_from_v7": [
            "arrayCGH input updated: v18 -> v19 (overlap_srWGS column-based exclusion replacing SampleID matching)",
            "arrayCGH v19 outputs split into ASD/SZ/het files; reading ASD result only",
            "detect_column: platform_meta_p_two added as candidate for v18 column name",
            "MSSNG v17 unchanged (already has segdup/exclusion BED filters)",
            "All three datasets now share unified QC: segdup < 50% + exclusion BED < 50%",
        ],
        "array_tsv": ARRAY_TSV,
        "mssng_tsv": MSSNG_TSV,
        "outdir": str(outdir),
        "primary_main_method": "inverse_variance_fixed_effect",
        "primary_sensitivity_method": "directional_one_sided_weighted_stouffer",
        "secondary_main_method": "inverse_variance_fixed_effect",
        "negative_control_method": "inverse_variance_fixed_effect_two_sided",
        "weight_mode_stouffer": WEIGHT_MODE,
        "discovery_sig_classes": DISCOVERY_SIG_CLASSES,
        "discovery_nonsig_classes": DISCOVERY_NONSIG_CLASSES,
        "arraycgh_version": "v22",
        "arraycgh_inference": "Platform-stratified Logit+LRT -> IVW fixed-effect meta (v22, tad04212026 pipeline)",
        "mssng_version": "v18_full_family",
        "mssng_inference": "GEE Independence (v18, tad04212026 pipeline)",
        "qc_alignment": "All cohorts apply: segdup < 50%, exclusion BED < 50% (unified with WGS v18)",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)

    log(f"Saved: {ivw_path}")
    log(f"Saved: {stouffer_path}")
    log(f"Saved: {main_path}")
    log(f"Saved: {config_path}")


def print_summary(meta_df: pd.DataFrame) -> None:
    log("=" * 70)
    log("SUMMARY (IVW = primary analysis)")
    log("=" * 70)

    # PRIMARY sig8
    primary_sig8 = meta_df.loc[meta_df["analysis_set"] == "PRIMARY_sig8_main"].copy()
    if len(primary_sig8) > 0:
        n_supported = int(primary_sig8["formal_meta_supported"].sum())
        log(f"PRIMARY_sig8_main (formal replication): {n_supported}/{len(primary_sig8)} meta-supported")
        show = primary_sig8.sort_values("ivw_p_main")[[
            "boundary_class", "array_or", "array_p_main", "mssng_or", "mssng_p_one",
            "ivw_or_meta", "ivw_p_main", "same_direction", "ivw_i2",
            "stouffer_p_sensitivity", "formal_meta_supported"
        ]]
        with pd.option_context("display.max_rows", 50, "display.width", 220, "display.float_format", "{:.4g}".format):
            log("\n" + show.to_string(index=False))

    # PRIMARY nonsig2
    primary_nonsig2 = meta_df.loc[meta_df["analysis_set"] == "PRIMARY_nonsig2_exploratory"].copy()
    if len(primary_nonsig2) > 0:
        n_expl = int(primary_nonsig2["exploratory_meta_signal"].sum())
        log(f"\nPRIMARY_nonsig2_exploratory: {n_expl}/{len(primary_nonsig2)} exploratory signal")
        show = primary_nonsig2.sort_values("ivw_p_main")[[
            "boundary_class", "ivw_or_meta", "ivw_p_main", "stouffer_p_sensitivity",
            "exploratory_meta_signal"
        ]]
        with pd.option_context("display.max_rows", 20, "display.width", 160, "display.float_format", "{:.4g}".format):
            log("\n" + show.to_string(index=False))

    # SECONDARY C sig8
    sec_sig8 = meta_df.loc[meta_df["analysis_set"] == "SECONDARY_C_sig8_main"].copy()
    if len(sec_sig8) > 0:
        n_sec = int(sec_sig8["secondary_meta_supported"].sum())
        log(f"\nSECONDARY_C_sig8_main: {n_sec}/{len(sec_sig8)} meta-supported")
        show = sec_sig8.sort_values("ivw_p_main")[[
            "boundary_class", "array_or", "array_p_main", "mssng_or", "mssng_p_one",
            "ivw_or_meta", "ivw_p_main", "same_direction", "ivw_i2",
            "stouffer_p_sensitivity", "secondary_meta_supported"
        ]]
        with pd.option_context("display.max_rows", 50, "display.width", 220, "display.float_format", "{:.4g}".format):
            log("\n" + show.to_string(index=False))

    # Negative controls
    neg = meta_df.loc[meta_df["analysis_set"].str.startswith("NEGCTRL_")].copy()
    if len(neg) > 0:
        n_neg_sig = int(neg["negative_control_nominal"].sum())
        log(f"\nNEGATIVE_CONTROL: {n_neg_sig}/{len(neg)} rows with IVW meta P < 0.05")
        show = neg.sort_values("ivw_p_main")[[
            "analysis_set", "pattern", "sv_type", "boundary_class",
            "ivw_or_meta", "ivw_p_main", "stouffer_p_sensitivity",
            "negative_control_nominal"
        ]]
        with pd.option_context("display.max_rows", 50, "display.width", 220, "display.float_format", "{:.4g}".format):
            log("\n" + show.to_string(index=False))


# ============================================================
# MAIN
# ============================================================


def main() -> None:
    total_start = time.time()
    log("=" * 70)
    log("tad_replication_meta_analysis_v12.py")
    log("  Main meta: IVW fixed-effect")
    log("  Sensitivity: weighted Stouffer (signed Z)")
    log("  PRIMARY main = IVW one-sided positive for sig8")
    log("  PRIMARY nonsig2 = IVW two-sided for exploratory")
    log("  Pattern C = IVW (same rules as PRIMARY)")
    log("  Negative controls = IVW two-sided")
    log("  v12 change: arrayCGH v21 -> v22, MSSNG v17 -> v18 (tad04212026 pipeline)")
    log("  arrayCGH v22 (platform-meta, segdup/excl filtered) + MSSNG v18 input")
    log("=" * 70)

    if not Path(ARRAY_TSV).exists():
        raise FileNotFoundError(f"ARRAY_TSV not found: {ARRAY_TSV}")
    if not Path(MSSNG_TSV).exists():
        raise FileNotFoundError(f"MSSNG_TSV not found: {MSSNG_TSV}")

    array_df = load_replication_result(ARRAY_TSV, "arrayCGH_v22")
    mssng_df = load_replication_result(MSSNG_TSV, "MSSNG_v18_fullfamily")

    pair_df = build_pair_table(array_df, mssng_df)

    meta_rows = []
    for _, row in pair_df.iterrows():
        meta_rows.append(do_meta_for_row(row))
    meta_df = pd.DataFrame(meta_rows)

    save_outputs(meta_df, OUTDIR)
    print_summary(meta_df)

    elapsed = time.time() - total_start
    log(f"Done. elapsed={elapsed:.1f}s")


if __name__ == "__main__":
    main()
