#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_04_isec_v10
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04292026/04_wgs_sv_boundary_overlap/output_v10/04_isec_v10_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04292026/04_wgs_sv_boundary_overlap/output_v10/04_isec_v10_%j.err

# ファイル名: 04_intersect_sv_v10.sbatch
# 処理内容:
#   - Step 4 v10 (SV × boundary intersection + GTF-based gene covariates + sample-level top-1% CNV QC) を ncbn-cpu で実行
#   - v9 からの変更点: sample-level top-1% CNV count exclusion を追加
#     (CNV_SAMPLE_COUNTS TSV を参照し、cnv_count が 99 percentile を超えるサンプルを universe から除外)
#   - 出力ディレクトリ: output_v10/
#   - Steps 1 AND 3 完了後に実行
#   - 実行時間を記録

set -euo pipefail

SCRIPT_DIR="/lustre12/home/kushima-pg/tad04292026/04_wgs_sv_boundary_overlap"
OUT_DIR="${SCRIPT_DIR}/output_v10"
mkdir -p "${OUT_DIR}"

echo "Job start: $(date -Iseconds)"
START_TS=$(date +%s)

cd "${SCRIPT_DIR}"
python3 34_intersect_sv_with_heffel_master_v10.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
