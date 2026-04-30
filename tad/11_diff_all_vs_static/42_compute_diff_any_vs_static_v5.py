#!/usr/bin/env python3
# 42_compute_diff_any_vs_static_v5.py
#
# 処理内容:
#  - Script 10 v5 (tad04292026 パイプライン、WGS discovery の top-1% sample-level
#    QC 経路を反映) の出力
#       arraycgh_burden_L2_and_specificity_v5.tsv
#       mssng_burden_L2_and_specificity_v5.tsv
#    を読み込み
#  - group_primary の Diff_any_DEL 列を exposure-aware に動的追加:
#      n_boundary, n_events: Diff_specific_n1_DEL + Diff_shared_n2plus_DEL (count sum)
#      carrier_boundary   : (specific>0 OR shared>0) を binary (0/1) に clip (Dynamic_any carrier)
#  - Fit & meta:
#      arrayCGH: B' pooled logistic (statsmodels.Logit) [main]
#      MSSNG   : GEE (Binomial + Independence, groups=FAMILYID) [main]
#      MSSNG   : 自前 Firth penalized logistic [sensitivity, pooled, family 無視]
#      Meta    : IVW fixed-effect + DerSimonian-Laird random-effects
#  - 6 tests (3 exposure × 2 category=Diff_any,Static) に BH + Holm 補正
#  - 出力: meta_diff_any_vs_static_v5.tsv (6 行、列は fixed+RE+Firth+FDR)
#  - 処理時間を先頭と末尾で記録
#
# v4 -> v5 変更点 (path update のみ):
#  1. REPLICATION_OUT_DIR: output_v4 -> output_v5
#  2. 入力 burden: arraycgh_burden_L2_and_specificity_v4.tsv -> v5
#                  mssng_burden_L2_and_specificity_v4.tsv    -> v5
#  3. OUT_DIR: output_v4 -> output_v5
#  4. 出力ファイル名: meta_diff_any_vs_static_v4.tsv -> v5
#  5. 解析ロジック (arrayCGH Logit, MSSNG GEE, Firth rank reduction, IVW/DL meta,
#     BH/Holm 補正, Diff_any 列加算) は v4 と完全一致

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=42_compute_diff_any_vs_static_v5_%j.out
#SBATCH --error=42_compute_diff_any_vs_static_v5_%j.err

from __future__ import annotations

import time
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import scipy.stats as sp_stats
import scipy.linalg as sp_linalg

import statsmodels.api as sm
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.genmod.families import Binomial
from statsmodels.genmod.cov_struct import Independence
from statsmodels.stats.multitest import multipletests


# =========================================================
# PATHS (tad04292026 パイプライン、v5: 10_v5 replication output を参照)
# =========================================================
REPLICATION_OUT_DIR = Path(
    "/lustre12/home/kushima-pg/tad04292026/10_replication_2way_meta/output_v5"
)

ACGH_BURDEN = REPLICATION_OUT_DIR / "arraycgh_burden_L2_and_specificity_v5.tsv"
MSSNG_BURDEN = REPLICATION_OUT_DIR / "mssng_burden_L2_and_specificity_v5.tsv"

OUT_DIR = Path(
    "/lustre12/home/kushima-pg/tad04292026/11_diff_all_vs_static/output_v5"
)
OUT_TSV = OUT_DIR / "meta_diff_any_vs_static_v5.tsv"


# =========================================================
# 定数 (v4/v5 pipeline と一致)
# =========================================================
EXPOSURES = ["n_boundary", "n_events", "carrier_boundary"]
MIN_CELL_COUNT = 5

ACGH_CASE_LABEL = "ASD"
ACGH_CTRL_LABEL = "CONT"
MSSNG_CASE_LABEL = "ASD"
MSSNG_CTRL_LABEL = "unaffected_sibling"

ACGH_COVS = [
    "sex",
    "log1p_total_del_bases_A",
    "log1p_total_gene_DEL_A",
    "platform_nimblegen",
]

MSSNG_COVS_BASE = ["Sex_numeric", "log1p_total_del_bases", "log1p_total_gene_DEL"]
MSSNG_ANCESTRY = ["anc_OTH", "anc_SAS", "anc_EAS", "anc_AMR", "anc_AFR"]
MSSNG_PLATFORM = [
    "plat_NovaSeq",
    "plat_HiSeq",
    "plat_HiSeq2000",
    "plat_HiSeq2500",
    "plat_CG",
]


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


# =========================================================
# Diff_any 列の exposure-aware 追加
# =========================================================
def add_diff_any_columns(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    out = df.copy()
    for exp in EXPOSURES:
        specific_col = f"{exp}_group_primary__Diff_specific_n1_DEL"
        shared_col = f"{exp}_group_primary__Diff_shared_n2plus_DEL"
        new_col = f"{exp}_group_primary__Diff_any_DEL"
        missing = [c for c in [specific_col, shared_col] if c not in out.columns]
        if missing:
            raise KeyError(
                f"[{prefix}] Missing columns in burden file: {missing}"
            )
        spec_vals = pd.to_numeric(out[specific_col], errors="coerce").fillna(0).astype(float)
        shared_vals = pd.to_numeric(out[shared_col], errors="coerce").fillna(0).astype(float)
        if exp == "carrier_boundary":
            out[new_col] = ((spec_vals > 0) | (shared_vals > 0)).astype(float)
            aggregation_type = "binary_OR"
        else:
            out[new_col] = spec_vals + shared_vals
            aggregation_type = "count_sum"
        n_nonzero = (out[new_col] > 0).sum()
        log(
            f"[{prefix}] Added {new_col} [{aggregation_type}]: "
            f"n_nonzero_case_or_ctrl={n_nonzero}, "
            f"max={out[new_col].max()}, mean={out[new_col].mean():.4f}"
        )
    return out


# =========================================================
# Fit helpers (v4/v5 pipeline と同一 logic)
# =========================================================
def _empty_res(status, n_case, n_ctrl, cc, cn):
    return {
        "n_case": int(n_case),
        "n_ctrl": int(n_ctrl),
        "carrier_case": int(cc),
        "carrier_ctrl": int(cn),
        "beta": np.nan,
        "se": np.nan,
        "or": np.nan,
        "or_lo95": np.nan,
        "or_hi95": np.nan,
        "p_value": np.nan,
        "fit_status": status,
    }


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
    needed = [exposure_col] + ACGH_COVS + ["is_case"]
    missing_cov = [c for c in needed if c not in sub.columns]
    if missing_cov:
        return _empty_res(f"missing_cov:{','.join(missing_cov)}", 0, 0, 0, 0)
    sub = (
        sub[needed]
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
        .copy()
    )
    sub[exposure_col] = pd.to_numeric(sub[exposure_col], errors="coerce")
    sub = sub.dropna(subset=[exposure_col])
    n_case = int(sub["is_case"].sum())
    n_ctrl = int(len(sub) - n_case)
    carrier_case = int(((sub["is_case"] == 1) & (sub[exposure_col] >= 1)).sum())
    carrier_ctrl = int(((sub["is_case"] == 0) & (sub[exposure_col] >= 1)).sum())

    if sub[exposure_col].nunique() < 2 or sub["is_case"].nunique() < 2:
        return _empty_res("no_variance", n_case, n_ctrl, carrier_case, carrier_ctrl)
    if (carrier_case + carrier_ctrl) < min_cell:
        return _empty_res(
            "insufficient_carriers", n_case, n_ctrl, carrier_case, carrier_ctrl
        )

    X_cols = [exposure_col] + ACGH_COVS
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
                status = (
                    "ok"
                    if fit.mle_retvals.get("converged", False)
                    else "not_converged"
                )
                return {
                    "n_case": n_case,
                    "n_ctrl": n_ctrl,
                    "carrier_case": carrier_case,
                    "carrier_ctrl": carrier_ctrl,
                    "beta": beta,
                    "se": se,
                    "or": float(np.exp(beta)),
                    "or_lo95": float(np.exp(beta - 1.959964 * se)),
                    "or_hi95": float(np.exp(beta + 1.959964 * se)),
                    "p_value": p,
                    "fit_status": status,
                }
        except Exception as e:
            last_err = str(e)[:120]
    return _empty_res(
        f"glm_error:{last_err}", n_case, n_ctrl, carrier_case, carrier_ctrl
    )


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
    needed = (
        [exposure_col, "is_case", "FAMILYID"]
        + MSSNG_COVS_BASE
        + anc_cols
        + plat_cols
    )
    missing_cov = [c for c in needed if c not in sub.columns]
    if missing_cov:
        d = _empty_res(
            f"missing_cov:{','.join(missing_cov)}", 0, 0, 0, 0
        )
        d["n_families"] = 0
        return d
    sub = (
        sub[needed]
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
        .copy()
    )
    sub[exposure_col] = pd.to_numeric(sub[exposure_col], errors="coerce")
    sub = sub.dropna(subset=[exposure_col])
    n_case = int(sub["is_case"].sum())
    n_ctrl = int(len(sub) - n_case)
    carrier_case = int(((sub["is_case"] == 1) & (sub[exposure_col] >= 1)).sum())
    carrier_ctrl = int(((sub["is_case"] == 0) & (sub[exposure_col] >= 1)).sum())
    n_fam = int(sub["FAMILYID"].nunique())

    if sub[exposure_col].nunique() < 2 or sub["is_case"].nunique() < 2:
        d = _empty_res("no_variance", n_case, n_ctrl, carrier_case, carrier_ctrl)
        d["n_families"] = n_fam
        return d
    if (carrier_case + carrier_ctrl) < min_cell:
        d = _empty_res(
            "insufficient_carriers", n_case, n_ctrl, carrier_case, carrier_ctrl
        )
        d["n_families"] = n_fam
        return d

    X_cols = [exposure_col] + MSSNG_COVS_BASE + anc_cols + plat_cols
    X = sub[X_cols].astype(float).to_numpy()
    X = np.column_stack([np.ones(len(X)), X])
    y = sub["is_case"].astype(int).to_numpy()

    try:
        model = GEE(
            endog=y,
            exog=X,
            groups=sub["FAMILYID"].values,
            family=Binomial(),
            cov_struct=Independence(),
        )
        fit = model.fit(maxiter=200)
        if fit.converged:
            beta = float(fit.params[1])
            se = float(fit.bse[1])
            p = float(fit.pvalues[1])
            return {
                "n_case": n_case,
                "n_ctrl": n_ctrl,
                "carrier_case": carrier_case,
                "carrier_ctrl": carrier_ctrl,
                "n_families": n_fam,
                "beta": beta,
                "se": se,
                "or": float(np.exp(beta)),
                "or_lo95": float(np.exp(beta - 1.959964 * se)),
                "or_hi95": float(np.exp(beta + 1.959964 * se)),
                "p_value": p,
                "fit_status": "ok",
            }
        else:
            d = _empty_res(
                "not_converged", n_case, n_ctrl, carrier_case, carrier_ctrl
            )
            d["n_families"] = n_fam
            return d
    except Exception as e:
        d = _empty_res(
            f"gee_error:{str(e)[:120]}",
            n_case,
            n_ctrl,
            carrier_case,
            carrier_ctrl,
        )
        d["n_families"] = n_fam
        return d


# =========================================================
# Firth penalized logistic regression (自前実装、v4 と同一)
# =========================================================
def _reduce_rank_via_qr(
    X: np.ndarray,
    protected_idx: Tuple[int, ...] = (0, 1),
    tol_rel: float = 1e-10,
) -> Tuple[np.ndarray, np.ndarray, bool]:
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

        for it in range(maxiter):
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

        if not converged:
            pass

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
                converged,
                "singular_hessian_at_se",
            )
        se_red = np.sqrt(np.maximum(np.diag(cov), 0.0))
        z_red = np.divide(
            beta, se_red, out=np.full_like(beta, np.nan), where=(se_red > 0)
        )
        p_red = 2.0 * sp_stats.norm.sf(np.abs(z_red))

        beta_full = np.full(p_full, np.nan)
        se_full = np.full(p_full, np.nan)
        p_full_arr = np.full(p_full, np.nan)
        for i, idx in enumerate(keep_idx):
            beta_full[idx] = beta[i]
            se_full[idx] = se_red[i]
            p_full_arr[idx] = p_red[i]

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


def fit_firth_mssng_pooled(
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
    needed = (
        [exposure_col, "is_case"]
        + MSSNG_COVS_BASE
        + anc_cols
        + plat_cols
    )
    missing_cov = [c for c in needed if c not in sub.columns]
    if missing_cov:
        return _empty_res(
            f"missing_cov:{','.join(missing_cov)}", 0, 0, 0, 0
        )
    sub = (
        sub[needed]
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
        .copy()
    )
    sub[exposure_col] = pd.to_numeric(sub[exposure_col], errors="coerce")
    sub = sub.dropna(subset=[exposure_col])
    n_case = int(sub["is_case"].sum())
    n_ctrl = int(len(sub) - n_case)
    carrier_case = int(((sub["is_case"] == 1) & (sub[exposure_col] >= 1)).sum())
    carrier_ctrl = int(((sub["is_case"] == 0) & (sub[exposure_col] >= 1)).sum())

    if sub[exposure_col].nunique() < 2 or sub["is_case"].nunique() < 2:
        return _empty_res("no_variance", n_case, n_ctrl, carrier_case, carrier_ctrl)
    if (carrier_case + carrier_ctrl) < min_cell:
        return _empty_res(
            "insufficient_carriers", n_case, n_ctrl, carrier_case, carrier_ctrl
        )

    X_cols = [exposure_col] + MSSNG_COVS_BASE + anc_cols + plat_cols
    X = sub[X_cols].astype(float).to_numpy()
    X = np.column_stack([np.ones(len(X)), X])
    y = sub["is_case"].astype(int).to_numpy()

    beta_arr, se_arr, p_arr, converged, status = firth_logit_fit(X, y)
    if not np.isfinite(beta_arr[1]) or not np.isfinite(se_arr[1]):
        return _empty_res(
            f"firth_failed:{status}",
            n_case,
            n_ctrl,
            carrier_case,
            carrier_ctrl,
        )
    beta = float(beta_arr[1])
    se = float(se_arr[1])
    p = float(p_arr[1])
    return {
        "n_case": n_case,
        "n_ctrl": n_ctrl,
        "carrier_case": carrier_case,
        "carrier_ctrl": carrier_ctrl,
        "beta": beta,
        "se": se,
        "or": float(np.exp(beta)),
        "or_lo95": float(np.exp(beta - 1.959964 * se)),
        "or_hi95": float(np.exp(beta + 1.959964 * se)),
        "p_value": p,
        "fit_status": status,
    }


# =========================================================
# IVW fixed-effect 2-way meta
# =========================================================
def ivw_meta(beta_a, se_a, beta_b, se_b) -> dict:
    if any(pd.isna([beta_a, se_a, beta_b, se_b])):
        return {
            "beta_meta": np.nan,
            "se_meta": np.nan,
            "or_meta": np.nan,
            "or_meta_lo95": np.nan,
            "or_meta_hi95": np.nan,
            "p_meta": np.nan,
            "q_het": np.nan,
            "p_het": np.nan,
            "i2_het": np.nan,
            "meta_status": "missing_inputs",
        }
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
        "beta_meta": float(beta_m),
        "se_meta": float(se_m),
        "or_meta": float(np.exp(beta_m)),
        "or_meta_lo95": float(np.exp(beta_m - 1.959964 * se_m)),
        "or_meta_hi95": float(np.exp(beta_m + 1.959964 * se_m)),
        "p_meta": float(p_m),
        "q_het": float(q),
        "p_het": float(p_q),
        "i2_het": float(i2),
        "meta_status": "ok",
    }


# =========================================================
# DerSimonian-Laird random-effects 2-way meta
# =========================================================
def dl_random_effects_meta(beta_a, se_a, beta_b, se_b) -> dict:
    if any(pd.isna([beta_a, se_a, beta_b, se_b])):
        return {
            "beta_meta_re": np.nan,
            "se_meta_re": np.nan,
            "or_meta_re": np.nan,
            "or_meta_re_lo95": np.nan,
            "or_meta_re_hi95": np.nan,
            "p_meta_re": np.nan,
            "tau2_re": np.nan,
            "meta_status_re": "missing_inputs",
        }
    w_a = 1.0 / (se_a ** 2)
    w_b = 1.0 / (se_b ** 2)
    beta_fe = (w_a * beta_a + w_b * beta_b) / (w_a + w_b)
    q = w_a * (beta_a - beta_fe) ** 2 + w_b * (beta_b - beta_fe) ** 2
    C = (w_a + w_b) - (w_a ** 2 + w_b ** 2) / (w_a + w_b)
    tau2 = max(0.0, (q - 1.0) / C) if C > 0 else 0.0
    w_a_re = 1.0 / (se_a ** 2 + tau2)
    w_b_re = 1.0 / (se_b ** 2 + tau2)
    beta_re = (w_a_re * beta_a + w_b_re * beta_b) / (w_a_re + w_b_re)
    se_re = np.sqrt(1.0 / (w_a_re + w_b_re))
    z = beta_re / se_re
    p_re = 2.0 * sp_stats.norm.sf(abs(z))
    return {
        "beta_meta_re": float(beta_re),
        "se_meta_re": float(se_re),
        "or_meta_re": float(np.exp(beta_re)),
        "or_meta_re_lo95": float(np.exp(beta_re - 1.959964 * se_re)),
        "or_meta_re_hi95": float(np.exp(beta_re + 1.959964 * se_re)),
        "p_meta_re": float(p_re),
        "tau2_re": float(tau2),
        "meta_status_re": "ok",
    }


# =========================================================
# Main
# =========================================================
def main() -> None:
    t0 = time.time()
    log("=" * 70)
    log("42_compute_diff_any_vs_static_v5.py START")
    log("=" * 70)
    log(f"  ACGH_BURDEN:  {ACGH_BURDEN}")
    log(f"  MSSNG_BURDEN: {MSSNG_BURDEN}")
    log(f"  OUT_TSV:      {OUT_TSV}")

    # --- Load burdens
    log("Loading arrayCGH burden ...")
    acgh = pd.read_csv(ACGH_BURDEN, sep="\t")
    log(f"  arrayCGH rows={len(acgh)}, cols={acgh.shape[1]}")
    log("Loading MSSNG burden ...")
    mssng = pd.read_csv(MSSNG_BURDEN, sep="\t")
    log(f"  MSSNG rows={len(mssng)}, cols={mssng.shape[1]}")

    # --- Add Diff_any columns (exposure-aware)
    acgh = add_diff_any_columns(acgh, prefix="arrayCGH")
    mssng = add_diff_any_columns(mssng, prefix="MSSNG")

    # --- Fit and meta per exposure, per category (Diff_any, Static)
    CATEGORIES = ["Diff_any", "Static"]
    rows: List[dict] = []

    for exp in EXPOSURES:
        for cat in CATEGORIES:
            col = f"{exp}_group_primary__{cat}_DEL"
            log(f"--- exposure={exp}  category={cat}  column={col} ---")

            r_acgh = fit_bprime_logit_acgh(
                acgh,
                exposure_col=col,
                case_label=ACGH_CASE_LABEL,
                control_label=ACGH_CTRL_LABEL,
                min_cell=MIN_CELL_COUNT,
            )
            log(
                f"  arrayCGH (main): beta={r_acgh['beta']}, OR={r_acgh['or']}, "
                f"p={r_acgh['p_value']}, status={r_acgh['fit_status']}"
            )

            r_mssng = fit_gee_mssng(
                mssng,
                exposure_col=col,
                case_label=MSSNG_CASE_LABEL,
                control_label=MSSNG_CTRL_LABEL,
                min_cell=MIN_CELL_COUNT,
            )
            log(
                f"  MSSNG GEE (main): beta={r_mssng['beta']}, OR={r_mssng['or']}, "
                f"p={r_mssng['p_value']}, status={r_mssng['fit_status']}"
            )

            r_mssng_firth = fit_firth_mssng_pooled(
                mssng,
                exposure_col=col,
                case_label=MSSNG_CASE_LABEL,
                control_label=MSSNG_CTRL_LABEL,
                min_cell=MIN_CELL_COUNT,
            )
            log(
                f"  MSSNG Firth (sens): beta={r_mssng_firth['beta']}, "
                f"OR={r_mssng_firth['or']}, p={r_mssng_firth['p_value']}, "
                f"status={r_mssng_firth['fit_status']}"
            )

            meta_fe = ivw_meta(
                r_acgh["beta"], r_acgh["se"], r_mssng["beta"], r_mssng["se"]
            )
            log(
                f"  META fixed: OR={meta_fe['or_meta']} "
                f"[{meta_fe['or_meta_lo95']}, {meta_fe['or_meta_hi95']}] "
                f"p_meta={meta_fe['p_meta']} I2={meta_fe['i2_het']}"
            )

            meta_re = dl_random_effects_meta(
                r_acgh["beta"], r_acgh["se"], r_mssng["beta"], r_mssng["se"]
            )
            log(
                f"  META random: OR={meta_re['or_meta_re']} "
                f"[{meta_re['or_meta_re_lo95']}, {meta_re['or_meta_re_hi95']}] "
                f"p_meta_re={meta_re['p_meta_re']} tau2={meta_re['tau2_re']}"
            )

            row = {
                "group_scheme": "group_primary",
                "group_label": cat,
                "exposure": exp,
                "sv_type": "DEL",
                "beta_acgh": r_acgh["beta"],
                "se_acgh": r_acgh["se"],
                "n_case_acgh": r_acgh["n_case"],
                "n_ctrl_acgh": r_acgh["n_ctrl"],
                "carrier_case_acgh": r_acgh["carrier_case"],
                "carrier_ctrl_acgh": r_acgh["carrier_ctrl"],
                "or_acgh": r_acgh["or"],
                "or_lo95_acgh": r_acgh["or_lo95"],
                "or_hi95_acgh": r_acgh["or_hi95"],
                "p_acgh": r_acgh["p_value"],
                "fit_status_acgh": r_acgh["fit_status"],
                "beta_mssng": r_mssng["beta"],
                "se_mssng": r_mssng["se"],
                "n_case_mssng": r_mssng["n_case"],
                "n_ctrl_mssng": r_mssng["n_ctrl"],
                "carrier_case_mssng": r_mssng["carrier_case"],
                "carrier_ctrl_mssng": r_mssng["carrier_ctrl"],
                "n_families_mssng": r_mssng.get("n_families", np.nan),
                "or_mssng": r_mssng["or"],
                "or_lo95_mssng": r_mssng["or_lo95"],
                "or_hi95_mssng": r_mssng["or_hi95"],
                "p_mssng": r_mssng["p_value"],
                "fit_status_mssng": r_mssng["fit_status"],
                "beta_mssng_firth": r_mssng_firth["beta"],
                "se_mssng_firth": r_mssng_firth["se"],
                "or_mssng_firth": r_mssng_firth["or"],
                "or_lo95_mssng_firth": r_mssng_firth["or_lo95"],
                "or_hi95_mssng_firth": r_mssng_firth["or_hi95"],
                "p_mssng_firth": r_mssng_firth["p_value"],
                "fit_status_mssng_firth": r_mssng_firth["fit_status"],
                **meta_fe,
                **meta_re,
            }
            rows.append(row)

    out_df = pd.DataFrame(rows)

    for pcol, prefix in [("p_meta", "p_meta"), ("p_meta_re", "p_meta_re")]:
        p_vals = out_df[pcol].values
        valid_mask = ~pd.isna(p_vals)
        if valid_mask.sum() > 0:
            bh = np.full(len(p_vals), np.nan)
            bh[valid_mask] = multipletests(p_vals[valid_mask], method="fdr_bh")[1]
            out_df[f"{prefix}_bh"] = bh
            holm = np.full(len(p_vals), np.nan)
            holm[valid_mask] = multipletests(p_vals[valid_mask], method="holm")[1]
            out_df[f"{prefix}_holm"] = holm
        else:
            out_df[f"{prefix}_bh"] = np.nan
            out_df[f"{prefix}_holm"] = np.nan

    out_df["p_meta_onesided"] = np.where(
        out_df["beta_meta"] > 0,
        out_df["p_meta"] / 2.0,
        1.0 - out_df["p_meta"] / 2.0,
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT_TSV, sep="\t", index=False)
    log(f"Wrote {OUT_TSV}  (rows={len(out_df)}, cols={out_df.shape[1]})")

    log("=" * 70)
    log("=== PREVIEW 1: fixed-effect meta (main) ===")
    prev1 = out_df[
        [
            "exposure",
            "group_label",
            "or_acgh",
            "p_acgh",
            "or_mssng",
            "p_mssng",
            "or_meta",
            "or_meta_lo95",
            "or_meta_hi95",
            "p_meta",
            "i2_het",
            "p_meta_bh",
            "p_meta_holm",
        ]
    ]
    log("\n" + prev1.to_string(index=False))

    log("=" * 70)
    log("=== PREVIEW 2: random-effects meta (DL tau^2) sensitivity ===")
    prev2 = out_df[
        [
            "exposure",
            "group_label",
            "or_meta",
            "p_meta",
            "or_meta_re",
            "or_meta_re_lo95",
            "or_meta_re_hi95",
            "p_meta_re",
            "tau2_re",
            "p_meta_re_bh",
            "p_meta_re_holm",
        ]
    ]
    log("\n" + prev2.to_string(index=False))

    log("=" * 70)
    log("=== PREVIEW 3: MSSNG GEE vs Firth sensitivity ===")
    prev3 = out_df[
        [
            "exposure",
            "group_label",
            "carrier_case_mssng",
            "carrier_ctrl_mssng",
            "or_mssng",
            "p_mssng",
            "or_mssng_firth",
            "or_lo95_mssng_firth",
            "or_hi95_mssng_firth",
            "p_mssng_firth",
            "fit_status_mssng_firth",
        ]
    ]
    log("\n" + prev3.to_string(index=False))

    t1 = time.time()
    log("=" * 70)
    log(f"42_compute_diff_any_vs_static_v5.py DONE in {t1 - t0:.1f} sec")
    log("=" * 70)


if __name__ == "__main__":
    main()
