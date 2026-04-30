#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_01_master
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04292026/01_heffel_boundary_master/output_v9/01_master_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04292026/01_heffel_boundary_master/output_v9/01_master_%j.err

# ファイル名: 01_build_master_v9.sbatch
# 処理内容:
#   - Step 1 (Heffel 2024 boundary master 構築) を ncbn-cpu で実行
#   - 実行時間を記録

set -euo pipefail

SCRIPT_DIR="/lustre12/home/kushima-pg/tad04292026/01_heffel_boundary_master"
mkdir -p "${SCRIPT_DIR}/output_v9"

echo "Job start: $(date -Iseconds)"
START_TS=$(date +%s)

cd "${SCRIPT_DIR}"
python3 33_build_heffel_boundary_master_v9.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
