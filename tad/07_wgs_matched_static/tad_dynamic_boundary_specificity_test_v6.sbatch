#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=spec_v6
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/07_wgs_matched_static/spec_v6_L%e_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/07_wgs_matched_static/spec_v6_L%e_%j.err

# ファイル名: tad_dynamic_boundary_specificity_test_v6.sbatch
# 処理内容:
#   - Script 7 (matched-static specificity test) を ncbn-cpu で実行
#   - 引数なし起動、MATCH_LEVEL/N_RESAMPLES は環境変数で上書き可能
#   - 実行時間を記録

set -euo pipefail

export MATCH_LEVEL="${MATCH_LEVEL:-2}"
export N_RESAMPLES="${N_RESAMPLES:-1000}"
export RANDOM_SEED="${RANDOM_SEED:-20260402}"

echo "Job start: $(date -Iseconds)"
echo "MATCH_LEVEL=${MATCH_LEVEL}"
echo "N_RESAMPLES=${N_RESAMPLES}"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04212026/07_wgs_matched_static
python3 tad_dynamic_boundary_specificity_test_v6.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
