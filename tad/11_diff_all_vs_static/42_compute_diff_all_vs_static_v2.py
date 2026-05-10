#!/usr/bin/env python3
# 42_compute_diff_all_vs_static_v2.py
#
# 処理内容:
#  - Script 10 v4 (tad04212026 パイプライン) の出力
#       arraycgh_burden_L2_and_specificity_v4.tsv
#       mssng_burden_L2_and_specificity_v4.tsv
#    を読み込み
#  - group_primary の Diff_all_DEL 列を動的に追加
#    ( = Diff_specific_n1_DEL + Diff_shared_n2plus_DEL )
#    対象 exposure: n_boundary, n_events, carrier_boundary
#  - arrayCGH: B' pooled logistic (statsmodels.Logit) で fit (v3/v4 と同一 covariates)
#  - MSSNG   : GEE (Binomial + Independence, groups=FAMILYID) で fit (v3/v4 と同一 covariates)
#  - IVW fixed-effect 2-way meta (arrayCGH ASD_vs_CONT + MSSNG ASD_vs_unaffSib)
#  - 出力: Diff_all / Static の per-bin / per-event / carrier OR, 95% CI, p_meta を TSV
#  - 処理時間を先頭と末尾で記録
#
# v1 -> v2 変更点:
#  - REPLICATION_OUT_DIR: heffel_deep_analysis_03242026/replication_full_pipeline_outputs_v2
#                       -> tad04212026/10_replication_2way_meta/output_v4
#  - 入力 TSV: _v1.tsv -> _v4.tsv
#  - OUT_DIR : tad04212026/11_diff_all_vs_static/output_v2
#  - OUT_TSV : meta_diff_all_vs_static_v1.tsv -> meta_diff_all_vs_static_v2.tsv
#  - 解析ロジック (add_diff_all_columns, fit_bprime_logit_acgh, fit_gee_mssng, ivw_meta) は
#    v1 から一切変更なし (列名が v4 burden TSV と完全一致することを事前確認済み)

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=42_compute_diff_all_vs_static_v2_%j.out
#SBATCH --error=42_compute_diff_all_vs_static_v2_%j.err

from __future__ import annotations

import time
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import scipy.stats as sp_stats

import statsmodels.api as sm
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.genmod.families import Binomial
from statsmodels.genmod.cov_struct import Independence


# =========================================================
# PATHS (tad04212026 パイプライン)
# =========================================================
REPLICATION_OUT_DIR = Path(
    "/lustre12/home/kushima-pg/tad04212026/10_replication_2way_meta/output_v4"
)

ACGH_BURDEN = REPLICATION_OUT_DIR / "arraycgh_burden_L2_and_specificity_v4.tsv"
MSSNG_BURDEN = REPLICATION_OUT_DIR / "mssng_burden_L2_and_specificity_v4.tsv"

OUT_DIR = Path(
    "/lustre12/home/kushima-pg/tad04212026/11_diff_all_vs_static/output_v2"
)
OUT_TSV = OUT_DIR / "meta_diff_all_vs_static_v2.tsv"


# =========================================================
# 定数 (v1/v3/v4 pipeline と一致)
# =========================================================
EXPOSURES = ["n_boundary", "n_events", "carrier_boundary"]
MIN_CELL_COUNT = 5

ACGH_CASE_LABEL = "ASD"
ACGH_CTRL_LABEL = "CONT"
MSSNG_CASE_LABEL = "ASD"
MSSNG_CTRL_LABEL = "unaffected_sibling"

# arrayCGH covariates (v3/v4 と同一)
ACGH_COVS = [
    "sex",
    "log1p_total_del_bases_A",
    "log1p_total_gene_DEL_A",
    "platform_nimblegen",
]

# MSSNG covariates (v3/v4 と同一; anc_*, plat_* は df にあるものだけ使う)
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
# Diff_all 列の追加 (列名を動的に参照)
# =========================================================
def add_diff_all_columns(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """
    df に n_boundary_group_primary__Diff_specific_n1_DEL と
    n_boundary_group_primary__Diff_shared_n2plus_DEL の和である
    Diff_all_DEL 列を追加する。

    exposure ごとに列名が異なるため、列名を動的に生成して存在確認する。
    """
    out = df.copy()
    for exp in EXPOSURES:
        specific_col = f"{exp}_group_primary__Diff_specific_n1_DEL"
        shared_col = f"{exp}_group_primary__Diff_shared_n2plus_DEL"
        new_col = f"{exp}_group_primary__Diff_all_DEL"
        missing = [c for c in [specific_col, shared_col] if c not in out.columns]
        if missing:
            raise KeyError(
                f"[{prefix}] Missing columns in burden file: {missing}"
            )
        out[new_col] = (
            pd.to_numeric(out[specific_col], errors="coerce").fillna(0).astype(float)
            + pd.to_numeric(out[shared_col], errors="coerce").fillna(0).astype(float)
        )
        log(
            f"[{prefix}] Added {new_col}: sum({specific_col} + {shared_col}) "
            f"n_nonzero_case_or_ctrl={(out[new_col] > 0).sum()}"
        )
    return out


# =========================================================
# Fit helpers (v3/v4 pipeline と同一 logic)
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
        d = _empty_res(
            "no_variance", n_case, n_ctrl, carrier_case, carrier_ctrl
        )
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
# IVW fixed-effect 2-way meta (v3/v4 と同一)
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
# Main
# =========================================================
def main() -> None:
    t0 = time.time()
    log("=" * 70)
    log("42_compute_diff_all_vs_static_v2.py START")
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

    # --- Add Diff_all columns (列名で動的に加算)
    acgh = add_diff_all_columns(acgh, prefix="arrayCGH")
    mssng = add_diff_all_columns(mssng, prefix="MSSNG")

    # --- Fit and meta per exposure, per category (Diff_all, Static)
    CATEGORIES = ["Diff_all", "Static"]
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
                f"  arrayCGH: beta={r_acgh['beta']}, se={r_acgh['se']}, "
                f"OR={r_acgh['or']}, p={r_acgh['p_value']}, "
                f"status={r_acgh['fit_status']}"
            )

            r_mssng = fit_gee_mssng(
                mssng,
                exposure_col=col,
                case_label=MSSNG_CASE_LABEL,
                control_label=MSSNG_CTRL_LABEL,
                min_cell=MIN_CELL_COUNT,
            )
            log(
                f"  MSSNG:    beta={r_mssng['beta']}, se={r_mssng['se']}, "
                f"OR={r_mssng['or']}, p={r_mssng['p_value']}, "
                f"status={r_mssng['fit_status']}"
            )

            meta = ivw_meta(
                r_acgh["beta"], r_acgh["se"], r_mssng["beta"], r_mssng["se"]
            )
            log(
                f"  META 2way: OR={meta['or_meta']} "
                f"[{meta['or_meta_lo95']}, {meta['or_meta_hi95']}] "
                f"p_meta={meta['p_meta']} I2={meta['i2_het']}"
            )

            row = {
                "group_scheme": "group_primary",
                "group_label": cat,
                "exposure": exp,
                "sv_type": "DEL",
                # arrayCGH
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
                # MSSNG
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
                # meta
                **meta,
            }
            rows.append(row)

    out_df = pd.DataFrame(rows)
    # p_meta の one-sided も付与 (Diff_all は discovery-positive 側、Static は内部参照)
    out_df["p_meta_onesided"] = np.where(
        out_df["beta_meta"] > 0,
        out_df["p_meta"] / 2.0,
        1.0 - out_df["p_meta"] / 2.0,
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT_TSV, sep="\t", index=False)
    log(f"Wrote {OUT_TSV}  (rows={len(out_df)})")

    # Preview (per-bin focus)
    preview = out_df[out_df["exposure"] == "n_boundary"][
        [
            "group_label",
            "exposure",
            "or_acgh",
            "p_acgh",
            "or_mssng",
            "p_mssng",
            "or_meta",
            "or_meta_lo95",
            "or_meta_hi95",
            "p_meta",
            "p_meta_onesided",
            "i2_het",
        ]
    ]
    log("=== PREVIEW: exposure=n_boundary ===")
    log("\n" + preview.to_string(index=False))

    t1 = time.time()
    log("=" * 70)
    log(f"42_compute_diff_all_vs_static_v2.py DONE in {t1 - t0:.1f} sec")
    log("=" * 70)


if __name__ == "__main__":
    main()
