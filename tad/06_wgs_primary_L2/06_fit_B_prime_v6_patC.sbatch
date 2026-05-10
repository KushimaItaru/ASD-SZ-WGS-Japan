#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH -J B_prime_v6_patC
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o /lustre12/home/kushima-pg/tad04212026/06_wgs_primary_L2/output_v6_patC/39_v6_patC.out
#SBATCH -e /lustre12/home/kushima-pg/tad04212026/06_wgs_primary_L2/output_v6_patC/39_v6_patC.err

# ファイル名: 06_fit_B_prime_v6_patC.sbatch
# 処理内容:
#   - Step 6 v6 (PATTERN=patC) を ncbn-cpu で実行
#   - Pattern C (NAHR-excluded + |SV_length| < 1 MB) の B' logistic regression (感度分析)
#   - v5 -> v6 変更点: path update のみ (統計・モデルは v5 と完全一致)
#       * 入力: 05_v5_{PATTERN} (silent-sample-loss fix 適用済み) の burden を参照
#       * 出力: 06_v6_{PATTERN}/ に書き出し
#   - 入力: /lustre12/home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v5_patC/
#            sample_burden_L2_and_specificity_v5_patC.tsv
#   - 出力: output_v6_patC/
#       B_prime_L2_classes_results_v6_patC.tsv
#       B_prime_specificity_groups_results_v6_patC.tsv
#       B_prime_run_summary_v6_patC.tsv
#   - Size cap (<1MB) は Step 03 v9 で既に適用済
#   - Step 05 v5 patC 完了後に実行
#   - 実行時間を記録

set -euo pipefail
SCRIPT_DIR="/lustre12/home/kushima-pg/tad04212026/06_wgs_primary_L2"
OUT_DIR="${SCRIPT_DIR}/output_v6_patC"
mkdir -p "${OUT_DIR}"

export PATTERN="patC"

echo "Job start: $(date)"
echo "PATTERN: ${PATTERN}"
START=$(date +%s)

cd "${SCRIPT_DIR}"
Rscript "${SCRIPT_DIR}/39_fit_B_prime_L2_and_specificity_v6.R"

END=$(date +%s)
ELAPSED=$((END - START))
echo "Job end: $(date)"
echo "Elapsed: ${ELAPSED}s"
