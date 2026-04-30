#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_05_burden_v4
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04292026/05_wgs_sample_burden/output_v4/05_burden_v4_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04292026/05_wgs_sample_burden/output_v4/05_burden_v4_%j.err

# ファイル名: 05_sample_burden_v4.sbatch
# 処理内容:
#   - Step 5 v4 (Pattern A primary; sample burden L2 + group_primary, S2-S5 削除済)
#     を ncbn-cpu で実行。
#   - v3 との差分は GROUP_COLUMNS から S2-S5 を削除し、bin annotation 入力を v2 → v3 に
#     更新した点のみ (group_primary は逐字維持; 数値結果は v3 と bit-identical)。
#   - 入力: 04_v10 event overlap + 02_v3 bin annotation
#   - 出力: output_v4/sample_burden_L2_and_specificity_v4.tsv + summary
#   - 前提: Step 02 v3 と Step 04 v10 が完了済
#   - 実行時間を記録

set -euo pipefail

SCRIPT_DIR="/lustre12/home/kushima-pg/tad04292026/05_wgs_sample_burden"
OUT_DIR="${SCRIPT_DIR}/output_v4"
mkdir -p "${OUT_DIR}"

echo "Job start: $(date -Iseconds)"
echo "Script: 05_sample_burden_v4.sbatch -> 38_compute_sample_burden_L2_and_specificity_v4.py"
echo "Output: ${OUT_DIR}"
START_TS=$(date +%s)

cd "${SCRIPT_DIR}"
python3 38_compute_sample_burden_L2_and_specificity_v4.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
