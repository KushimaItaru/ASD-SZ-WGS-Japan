#!/bin/bash
# ファイル名: 10_replication_full_pipeline_v5.sbatch
# 処理内容:
#   - Script 10 v5 (40_run_replication_full_pipeline_v5.py) を ncbn-cpu パーティション
#     で実行する sbatch ラッパー (遺伝研スパコン用)。
#   - 出力: tad04212026/10_replication_2way_meta/output_v5/
#       arraycgh_burden_L2_and_specificity_v5.tsv
#       mssng_burden_L2_and_specificity_v5.tsv
#       replication_burden_summary_v5.tsv
#       arraycgh_L2_classes_results_v5.tsv
#       arraycgh_specificity_groups_results_v5.tsv
#       mssng_L2_classes_results_v5.tsv
#       mssng_specificity_groups_results_v5.tsv
#       meta_L2_classes_results_v5.tsv
#       meta_specificity_groups_results_v5.tsv
#       meta3way_L2_classes_results_v5.tsv        (WGS v4 + arrayCGH v22 + MSSNG)
#       meta3way_specificity_groups_results_v5.tsv
#       replication_fit_summary_v5.tsv
#
#   - v4 -> v5 変更点 (path update のみ; WGS discovery sample-level top-1% CNV count
#     QC の伝播):
#       1. WGS_PRIMARY_OUT_DIR: 06_wgs_primary_L2/output_v3 -> output_v4
#       2. WGS_L2_RESULTS:   B_prime_L2_classes_results_v3.tsv  -> _v4.tsv
#       3. WGS_SPEC_RESULTS: B_prime_specificity_groups_results_v3.tsv -> _v4.tsv
#       4. OUT_DIR: output_v4 -> output_v5
#       5. 全出力 TSV の suffix: _v4 -> _v5
#       6. arrayCGH v22 / MSSNG / bin_l2_annotation_v2 は v4 と完全同一
#       7. 解析ロジック (Step A/B/C/D) は v4 と完全一致
#
#   - 依存:
#       tad04212026/06_wgs_primary_L2/output_v4/B_prime_L2_classes_results_v4.tsv
#       tad04212026/06_wgs_primary_L2/output_v4/B_prime_specificity_groups_results_v4.tsv
#       tad04212026/08_arraycgh_sample_burden/output_v22/ (変更なし)
#       tad04212026/02_bin_l2_annotation/output_v2/bin_l2_annotation_v2.tsv.gz (変更なし)
#       noncoding_tad_specificity_04022026/output_mssng_factorial_overlap_v2/ (変更なし)
#
#   - Step 06 v4 (WGS primary fit) 完了後に実行すること。
#   - fit モデル総数: ACGH (L2 10 + spec 15) × exp 3 × cmp 2
#                   + MSSNG (L2 10 + spec 15) × exp 3 × cmp 1
#                   = 150 + 75 = 225 logistic/GEE fit
#     + 2-way meta (225 ペア) + 3-way meta (75 WGS ペア)
#   - 保守的に 16 CPU / 12G / 2h 割当。

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_rep_v5
#SBATCH --cpus-per-task=16
#SBATCH --mem=12G
#SBATCH --time=2:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/10_replication_2way_meta/tad_rep_v5_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/10_replication_2way_meta/tad_rep_v5_%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04212026/10_replication_2way_meta
python3 40_run_replication_full_pipeline_v5.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
