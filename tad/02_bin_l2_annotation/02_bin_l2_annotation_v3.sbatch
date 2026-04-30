#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_02_binl2_v3
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04292026/02_bin_l2_annotation/output_v3/02_binl2_v3_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04292026/02_bin_l2_annotation/output_v3/02_binl2_v3_%j.err

# ファイル名: 02_bin_l2_annotation_v3.sbatch
# 処理内容:
#   - Step 2 v3 (bin-level L2 annotation, S2-S5 sensitivity grouping 削除済) を ncbn-cpu で実行
#   - v2 との差分は GROUP_COLUMNS から S2-S5 を削除した点のみ (group_primary は逐字維持)
#   - 出力: output_v3/bin_l2_annotation_v3.tsv.gz + bin_l2_annotation_summary_v3.tsv
#   - 実行時間を記録

set -euo pipefail

SCRIPT_DIR="/lustre12/home/kushima-pg/tad04292026/02_bin_l2_annotation"
OUT_DIR="${SCRIPT_DIR}/output_v3"
mkdir -p "${OUT_DIR}"

echo "Job start: $(date -Iseconds)"
echo "Script: 02_bin_l2_annotation_v3.sbatch -> 37_compute_bin_l2_annotation_v3.py"
echo "Output: ${OUT_DIR}"
START_TS=$(date +%s)

cd "${SCRIPT_DIR}"
python3 37_compute_bin_l2_annotation_v3.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
