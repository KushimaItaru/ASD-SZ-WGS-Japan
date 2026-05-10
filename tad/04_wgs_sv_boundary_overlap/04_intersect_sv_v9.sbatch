#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_04_isec
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/04_wgs_sv_boundary_overlap/output_v9/04_isec_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/04_wgs_sv_boundary_overlap/output_v9/04_isec_%j.err

# ファイル名: 04_intersect_sv_v9.sbatch
# 処理内容:
#   - Step 4 (SV × boundary intersection + GTF-based gene covariates) を ncbn-cpu で実行
#   - Steps 1 AND 3 完了後に実行
#   - 実行時間を記録

set -euo pipefail

SCRIPT_DIR="/lustre12/home/kushima-pg/tad04212026/04_wgs_sv_boundary_overlap"
mkdir -p "${SCRIPT_DIR}/output_v9"

echo "Job start: $(date -Iseconds)"
START_TS=$(date +%s)

cd "${SCRIPT_DIR}"
python3 34_intersect_sv_with_heffel_master_v9.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
