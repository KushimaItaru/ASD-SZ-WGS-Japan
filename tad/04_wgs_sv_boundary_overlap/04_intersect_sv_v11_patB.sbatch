#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_04_isec_v11_patB
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04292026/04_wgs_sv_boundary_overlap/output_v11_patB/04_isec_v11_patB_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04292026/04_wgs_sv_boundary_overlap/output_v11_patB/04_isec_v11_patB_%j.err

# ファイル名: 04_intersect_sv_v11_patB.sbatch
# 処理内容:
#   - Step 4 v11 (PATTERN=patB) を ncbn-cpu で実行
#   - Pattern B (NAHR-included) の感度分析 overlap + burden 計算
#   - 入力: /lustre12/home/kushima-pg/tad04292026/03_wgs_sv_events/output_v9/
#            wgs_rare_sv_events_v9_patB.tsv.gz
#   - 出力: output_v11_patB/
#       sample_boundary_event_overlap_v11_patB.tsv.gz
#       sample_mechanistic_burden_v11_patB.tsv
#       sample_mechanistic_burden_summary_v11_patB.tsv
#   - 解析ロジックは v10 (Pattern A) と完全一致
#   - v10 と同一リソース見積
#   - Step 03 v9 完了後に実行
#   - 実行時間を記録

set -euo pipefail

SCRIPT_DIR="/lustre12/home/kushima-pg/tad04292026/04_wgs_sv_boundary_overlap"
OUT_DIR="${SCRIPT_DIR}/output_v11_patB"
mkdir -p "${OUT_DIR}"

export PATTERN="patB"

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "PATTERN: ${PATTERN}"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: ${SCRIPT_DIR}"
START_TS=$(date +%s)

cd "${SCRIPT_DIR}"
python3 34_intersect_sv_with_heffel_master_v11.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
