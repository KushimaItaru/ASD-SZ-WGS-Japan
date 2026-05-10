#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH -J B_prime_v4
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o /lustre12/home/kushima-pg/tad04212026/06_wgs_primary_L2/output_v4/39_v4.out
#SBATCH -e /lustre12/home/kushima-pg/tad04212026/06_wgs_primary_L2/output_v4/39_v4.err

# ファイル名: 06_fit_B_prime_v4.sbatch
# 処理内容:
#   - Step 6 v4 (B' logistic regression, input=05_v3) を ncbn-cpu で実行
#   - v3 との差分は input burden が 05_v2 → 05_v3 (sample-level top-1% QC 適用済み) になった点のみ
#   - 出力: output_v4/
#   - Steps 05_v3 完了後に実行
#   - 実行時間を記録

set -euo pipefail
SCRIPT_DIR="/lustre12/home/kushima-pg/tad04212026/06_wgs_primary_L2"
OUT_DIR="${SCRIPT_DIR}/output_v4"
mkdir -p "${OUT_DIR}"
START=$(date +%s)
echo "Start: $(date)"
Rscript "${SCRIPT_DIR}/39_fit_B_prime_L2_and_specificity_v4.R"
END=$(date +%s)
echo "End  : $(date)"
echo "Elapsed: $((END-START))s"
