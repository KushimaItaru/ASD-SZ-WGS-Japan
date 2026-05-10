#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_03_sv
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/03_wgs_sv_events/output_v8/03_sv_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/03_wgs_sv_events/output_v8/03_sv_%j.err

# ファイル名: 03_extract_sv_v8.sbatch
# 処理内容:
#   - Step 3 (rare SV events extraction from AnnotSV) を ncbn-cpu で実行
#   - Step 1/2 と独立に並列実行可
#   - 実行時間を記録

set -euo pipefail

SCRIPT_DIR="/lustre12/home/kushima-pg/tad04212026/03_wgs_sv_events"
mkdir -p "${SCRIPT_DIR}/output_v8"

echo "Job start: $(date -Iseconds)"
START_TS=$(date +%s)

cd "${SCRIPT_DIR}"
python3 34a_extract_rare_sv_events_v8.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
