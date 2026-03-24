#!/usr/bin/env python3
# 17_burden_statistical_test_v20.py
#
# 処理内容:
# - v19 per_sampleファイル（observed_clusters_total列を含む）を読み込み、
#   共変量完備サンプルでburden解析を実施
# - Primary: rare_any（≥1保有の二値）に対するロジスティック回帰 → OR, 95%CI, P
#   ※ v19から変更なし
# - Secondary: rare_outlier_count（カウント）に対するPoisson GLM → RR, 95%CI, P
#   ※ v20修正: offset=log(observed_clusters_total)を追加
#   ※ observed_clusters_total=0のサンプルはPoisson GLMから除外（log(0)未定義のため）
# - 補助: Mann-Whitney U検定（ノンパラメトリック）
# - ASD vs Healthy, SZ vs Healthy の2比較
# - 共変量: Sex_M, Depth, PC1-PC10
# - 実行時間を記録
#
# v19→v20変更点:
#   1. 入力ファイルをv19 per_sample（observed_clusters_total列含む）に変更
#   2. Poisson GLMにoffset=log(observed_clusters_total)を追加
#   3. observed_clusters_total=0サンプルをPoisson GLMから除外、除外IDをログ出力
#   4. 出力TSVにexposure関連サマリー列を追加
#   5. ロジスティック回帰（Primary）は変更なし

import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu
import os
import sys
import time
from pathlib import Path


def main():
    start_time = time.time()
    print("=== Burden Statistical Analysis v20 ===")
    print("=== Poisson GLM now includes offset=log(observed_clusters_total) ===")
    project_root = Path(__file__).resolve().parents[1]

    # v20: v19のper_sampleファイルを入力（observed_clusters_total列を含む）
    input_file = Path(os.environ.get("INPUT_TSV", str(
        project_root / "analysis_results_novel" / "outlier_burden_rare_crossfit_v19.per_sample.tsv")))
    output_file = Path(os.environ.get("OUTPUT_TSV", str(
        project_root / "analysis_results_novel" / "burden_stats_results_v20.tsv")))

    if not input_file.exists():
        sys.exit(f"[ERROR] Input file not found: {input_file}")

    df = pd.read_csv(input_file, sep="\t")
    print(f"[INFO] Loaded {len(df)} samples from {input_file.name}")

    # observed_clusters_total列の存在確認
    if "observed_clusters_total" not in df.columns:
        sys.exit("[ERROR] Column 'observed_clusters_total' not found in input. "
                 "Please use v19 per_sample output as input.")

    # 共変量の定義（列名ベースで動的に選択）
    covariates = ["Sex_M", "Depth"] + [f"PC{i}" for i in range(1, 11)]
    missing_cols = [c for c in covariates if c not in df.columns]
    if missing_cols:
        sys.exit(f"[ERROR] Missing covariate columns: {missing_cols}")

    # 解析用データの準備 (3群のみ)
    df = df[df["Group"].isin(["Healthy", "ASD", "SZ"])].copy()

    # 欠損値除去
    len_before = len(df)
    required_cols = ["rare_outlier_count", "observed_clusters_total"] + covariates
    df = df.dropna(subset=required_cols)
    if len(df) < len_before:
        print(f"[INFO] Dropped {len_before - len(df)} samples due to missing values.")

    # rare_any（二値エンドポイント）を作成
    df["rare_any"] = (df["rare_outlier_count"] >= 1).astype(int)

    # observed_clusters_totalのサマリーを表示
    print(f"\n[INFO] observed_clusters_total summary (all groups):")
    print(f"  mean={df['observed_clusters_total'].mean():.1f}, "
          f"median={df['observed_clusters_total'].median():.1f}, "
          f"min={df['observed_clusters_total'].min()}, "
          f"max={df['observed_clusters_total'].max()}")
    n_zero_all = (df["observed_clusters_total"] == 0).sum()
    print(f"  N with observed_clusters_total=0: {n_zero_all}")

    # exposure=0サンプルの特定（Poisson GLMから除外対象）
    zero_exposure_ids = df.loc[df["observed_clusters_total"] == 0, "SampleID"].tolist()
    if zero_exposure_ids:
        print(f"\n[WARN] {len(zero_exposure_ids)} samples have observed_clusters_total=0 "
              f"and will be EXCLUDED from Poisson GLM (log(0) undefined):")
        for sid in zero_exposure_ids:
            grp = df.loc[df["SampleID"] == sid, "Group"].values[0]
            print(f"  Excluded: {sid} (Group={grp})")

    results = []
    target_groups = ["ASD", "SZ"]

    # Control Group
    df_ctrl = df[df["Group"] == "Healthy"]
    print(f"\nControl (Healthy): N={len(df_ctrl)}, "
          f"rare_any={df_ctrl['rare_any'].mean():.3f}, "
          f"mean_count={df_ctrl['rare_outlier_count'].mean():.4f}, "
          f"mean_exposure={df_ctrl['observed_clusters_total'].mean():.1f}")

    cov_formula = "Sex_M + Depth + " + " + ".join([f"PC{i}" for i in range(1, 11)])

    for case_grp in target_groups:
        df_case = df[df["Group"] == case_grp]
        if len(df_case) == 0:
            continue

        print(f"Case ({case_grp}): N={len(df_case)}, "
              f"rare_any={df_case['rare_any'].mean():.3f}, "
              f"mean_count={df_case['rare_outlier_count'].mean():.4f}, "
              f"mean_exposure={df_case['observed_clusters_total'].mean():.1f}")

        df_sub = pd.concat([df_ctrl, df_case]).copy()
        df_sub["IsCase"] = (df_sub["Group"] == case_grp).astype(int)

        row = {
            "Comparison": f"{case_grp} vs Healthy",
            "N_Case": len(df_case),
            "N_Ctrl": len(df_ctrl),
            "rare_any_Case": f"{df_case['rare_any'].mean():.3f}",
            "rare_any_Ctrl": f"{df_ctrl['rare_any'].mean():.3f}",
            "mean_count_Case": f"{df_case['rare_outlier_count'].mean():.4f}",
            "mean_count_Ctrl": f"{df_ctrl['rare_outlier_count'].mean():.4f}",
        }

        # --- 1. Mann-Whitney U Test (on count) ---
        stat, p_mwu = mannwhitneyu(
            df_case["rare_outlier_count"].values,
            df_ctrl["rare_outlier_count"].values,
            alternative="greater")
        row["MWU_P"] = p_mwu

        # --- 2. Primary: Logistic Regression on rare_any → OR ---
        # ※ v19と同一（offsetなし、全サンプル使用）
        formula_logit = f"rare_any ~ IsCase + {cov_formula}"
        try:
            model_logit = smf.logit(formula=formula_logit, data=df_sub).fit(
                disp=0, method="bfgs", maxiter=200)
            coef = model_logit.params["IsCase"]
            ci = model_logit.conf_int().loc["IsCase"]
            row["Logit_OR"] = np.exp(coef)
            row["Logit_OR_CI_low"] = np.exp(ci[0])
            row["Logit_OR_CI_high"] = np.exp(ci[1])
            row["Logit_P"] = model_logit.pvalues["IsCase"]
            row["Logit_N_used"] = len(df_sub)
        except Exception as e:
            print(f"[WARN] Logistic regression failed for {case_grp}: {e}")
            row["Logit_OR"] = np.nan
            row["Logit_OR_CI_low"] = np.nan
            row["Logit_OR_CI_high"] = np.nan
            row["Logit_P"] = np.nan
            row["Logit_N_used"] = np.nan

        # --- 3. Secondary: Poisson GLM on count with offset → RR ---
        # v20修正: offset=log(observed_clusters_total)を追加
        # observed_clusters_total=0のサンプルを除外
        df_pois = df_sub[df_sub["observed_clusters_total"] > 0].copy()
        n_excluded_case = ((df_sub["Group"] == case_grp) & (df_sub["observed_clusters_total"] == 0)).sum()
        n_excluded_ctrl = ((df_sub["Group"] == "Healthy") & (df_sub["observed_clusters_total"] == 0)).sum()
        n_pois_case = ((df_pois["Group"] == case_grp)).sum()
        n_pois_ctrl = ((df_pois["Group"] == "Healthy")).sum()

        print(f"  Poisson GLM ({case_grp}): N_used={len(df_pois)} "
              f"(excluded {n_excluded_case} case + {n_excluded_ctrl} ctrl with exposure=0)")

        formula_pois = f"rare_outlier_count ~ IsCase + {cov_formula}"
        try:
            # offset = log(observed_clusters_total)
            offset_vals = np.log(df_pois["observed_clusters_total"].astype(float).values)

            model_pois = smf.glm(
                formula=formula_pois, data=df_pois,
                family=sm.families.Poisson(),
                offset=offset_vals).fit()
            coef_p = model_pois.params["IsCase"]
            ci_p = model_pois.conf_int().loc["IsCase"]
            row["Poisson_RR"] = np.exp(coef_p)
            row["Poisson_RR_CI_low"] = np.exp(ci_p[0])
            row["Poisson_RR_CI_high"] = np.exp(ci_p[1])
            row["Poisson_P"] = model_pois.pvalues["IsCase"]
            row["Poisson_N_used"] = len(df_pois)
            row["Poisson_N_excluded_exposure0"] = n_excluded_case + n_excluded_ctrl
        except Exception as e:
            print(f"[WARN] Poisson GLM failed for {case_grp}: {e}")
            row["Poisson_RR"] = np.nan
            row["Poisson_RR_CI_low"] = np.nan
            row["Poisson_RR_CI_high"] = np.nan
            row["Poisson_P"] = np.nan
            row["Poisson_N_used"] = np.nan
            row["Poisson_N_excluded_exposure0"] = n_excluded_case + n_excluded_ctrl

        # Exposure summary columns
        row["Exposure_Column"] = "observed_clusters_total"
        row["mean_exposure_Case"] = f"{df_pois.loc[df_pois['Group'] == case_grp, 'observed_clusters_total'].mean():.1f}"
        row["mean_exposure_Ctrl"] = f"{df_pois.loc[df_pois['Group'] == 'Healthy', 'observed_clusters_total'].mean():.1f}"
        row["N_exposure_positive_Case"] = n_pois_case
        row["N_exposure_positive_Ctrl"] = n_pois_ctrl

        results.append(row)

    # 結果保存
    res_df = pd.DataFrame(results)
    res_df.to_csv(output_file, sep="\t", index=False)

    print("\n" + "=" * 70)
    print("=== Primary: Logistic Regression (rare_any -> OR) [unchanged from v19] ===")
    for _, r in res_df.iterrows():
        print(f"  {r['Comparison']}: OR={r['Logit_OR']:.2f} "
              f"(95%CI {r['Logit_OR_CI_low']:.2f}-{r['Logit_OR_CI_high']:.2f}), "
              f"P={r['Logit_P']:.2e}, N={r['Logit_N_used']}")

    print("\n=== Secondary: Poisson GLM with offset (count -> RR) [v20 UPDATED] ===")
    for _, r in res_df.iterrows():
        print(f"  {r['Comparison']}: RR={r['Poisson_RR']:.2f} "
              f"(95%CI {r['Poisson_RR_CI_low']:.2f}-{r['Poisson_RR_CI_high']:.2f}), "
              f"P={r['Poisson_P']:.2e}, N={r['Poisson_N_used']} "
              f"(excluded {r['Poisson_N_excluded_exposure0']} with exposure=0)")

    print(f"\n=== Mann-Whitney U (supplementary) ===")
    for _, r in res_df.iterrows():
        print(f"  {r['Comparison']}: P={r['MWU_P']:.2e}")

    elapsed = time.time() - start_time
    print(f"\n[INFO] Saved to {output_file}")
    print(f"[DONE] Elapsed: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
