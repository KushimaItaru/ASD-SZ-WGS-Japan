#!/bin/bash
# ファイル名: 08_arraycgh_replication_v22.sbatch
# 処理内容:
#   - Script 8 (arrayCGH replication v22) を ncbn-cpu パーティションで実行する
#     sbatch ラッパー。
#   - 遺伝研スパコン用。実行時間を記録。
#   - 出力は tad04212026/08_arraycgh_sample_burden/ 直下に log / err / output_v22/ を生成。

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_acgh_v22
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/08_arraycgh_sample_burden/tad_acgh_v22_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/08_arraycgh_sample_burden/tad_acgh_v22_%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04212026/08_arraycgh_sample_burden
python3 tad_replication_arraycgh_v22.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
