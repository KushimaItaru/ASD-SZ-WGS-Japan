#!/bin/bash
# ファイル名: 09a_mssng_replication_v18.sbatch
# 処理内容:
#   - Script 9A (MSSNG replication v18) を ncbn-cpu パーティションで実行する
#     sbatch ラッパー。
#   - 遺伝研スパコン用。実行時間を記録。
#   - 出力は tad04212026/09_mssng_sample_burden/ 直下に log / err / output_v18/ を生成。
#   - v17 の実行実績 (GEE Independence, n=6466) は数分〜数十分レンジのため
#     8 CPU / 32G / 4h を保守的上限として割当。

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_mssng_v18
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/09_mssng_sample_burden/tad_mssng_v18_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/09_mssng_sample_burden/tad_mssng_v18_%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04212026/09_mssng_sample_burden
python3 tad_replication_mssng_v18.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
