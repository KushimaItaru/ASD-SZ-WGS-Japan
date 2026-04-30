#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH -J B_prime_v5
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o /lustre12/home/kushima-pg/tad04292026/06_wgs_primary_L2/output_v5/39_v5.out
#SBATCH -e /lustre12/home/kushima-pg/tad04292026/06_wgs_primary_L2/output_v5/39_v5.err

# ファイル名: 06_fit_B_prime_v5.sbatch
# 処理内容:
#   - Step 6 v5 (Pattern A primary; B' logistic regression, S2-S5 sensitivity
#     grouping 削除済) を ncbn-cpu で実行。
#   - v4 との差分は GROUP_COLUMNS = c("group_primary") のみ、specificity_groups
#     fit が group_primary だけに縮小。10 L2 classes B' fit (paper Fig 3a) は
#     v4 と逐字維持。
#   - 入力: 05_v4 burden (sample_burden_L2_and_specificity_v4.tsv)
#   - 出力: output_v5/B_prime_L2_classes_results_v5.tsv + group_primary 結果 + summary
#   - 前提: Step 05 v4 完了済
#   - 実行時間を記録

set -euo pipefail
SCRIPT_DIR="/lustre12/home/kushima-pg/tad04292026/06_wgs_primary_L2"
OUT_DIR="${SCRIPT_DIR}/output_v5"
mkdir -p "${OUT_DIR}"
START=$(date +%s)
echo "Start: $(date)"
echo "Script: 06_fit_B_prime_v5.sbatch -> 39_fit_B_prime_L2_and_specificity_v5.R"
echo "Output: ${OUT_DIR}"
Rscript "${SCRIPT_DIR}/39_fit_B_prime_L2_and_specificity_v5.R"
END=$(date +%s)
echo "End  : $(date)"
echo "Elapsed: $((END-START))s"
