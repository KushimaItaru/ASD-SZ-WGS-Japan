#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_02_binl2
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/02_bin_l2_annotation/output_v2/02_binl2_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/02_bin_l2_annotation/output_v2/02_binl2_%j.err

# ファイル名: 02_bin_l2_annotation_v2.sbatch
# 処理内容:
#   - Step 2 (bin-level L2 annotation) を ncbn-cpu で実行 (Step 1 完了後)
#   - 実行時間を記録

set -euo pipefail

SCRIPT_DIR="/lustre12/home/kushima-pg/tad04212026/02_bin_l2_annotation"
mkdir -p "${SCRIPT_DIR}/output_v2"

echo "Job start: $(date -Iseconds)"
START_TS=$(date +%s)

cd "${SCRIPT_DIR}"
python3 37_compute_bin_l2_annotation_v2.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
