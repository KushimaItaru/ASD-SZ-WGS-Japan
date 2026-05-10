#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_05_burden
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v2/05_burden_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v2/05_burden_%j.err

# ファイル名: 05_sample_burden_v2.sbatch
# 処理内容:
#   - Step 5 (sample burden L2 + specificity) を ncbn-cpu で実行
#   - Steps 2 AND 4 完了後に実行
#   - 実行時間を記録

set -euo pipefail

SCRIPT_DIR="/lustre12/home/kushima-pg/tad04212026/05_wgs_sample_burden"
mkdir -p "${SCRIPT_DIR}/output_v2"

echo "Job start: $(date -Iseconds)"
START_TS=$(date +%s)

cd "${SCRIPT_DIR}"
python3 38_compute_sample_burden_L2_and_specificity_v2.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
