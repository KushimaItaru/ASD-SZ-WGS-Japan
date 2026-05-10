#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_05_burden_v3
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v3/05_burden_v3_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v3/05_burden_v3_%j.err

# ファイル名: 05_sample_burden_v3.sbatch
# 処理内容:
#   - Step 5 v3 (sample burden L2 + specificity; input=04_v10) を ncbn-cpu で実行
#   - v2 との差分は入力が 04_v9 → 04_v10 になった点のみ (sample-level top-1% QC 適用済み)
#   - 出力: output_v3/
#   - Steps 2 AND 4_v10 完了後に実行
#   - 実行時間を記録

set -euo pipefail

SCRIPT_DIR="/lustre12/home/kushima-pg/tad04212026/05_wgs_sample_burden"
OUT_DIR="${SCRIPT_DIR}/output_v3"
mkdir -p "${OUT_DIR}"

echo "Job start: $(date -Iseconds)"
START_TS=$(date +%s)

cd "${SCRIPT_DIR}"
python3 38_compute_sample_burden_L2_and_specificity_v3.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
