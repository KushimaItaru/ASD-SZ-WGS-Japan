#!/bin/bash
# ファイル名: 09b_meta_analysis_v12.sbatch
# 処理内容:
#   - Script 9B (arrayCGH v22 + MSSNG v18 meta-analysis v12) を ncbn-cpu パーティションで実行する
#     sbatch ラッパー。
#   - 遺伝研スパコン用。実行時間を記録。
#   - 出力は tad04212026/09_mssng_sample_burden/meta_v12/ に生成。
#   - meta 解析は numpy/pandas のみの軽い処理のため最小リソース割当。
#   - 依存: 09_mssng_sample_burden/output_v18/tad_replication_mssng_v18.tsv と
#     08_arraycgh_sample_burden/output_v22/tad_replication_asd_vs_cont_v22.tsv が必要。

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_meta_v12
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=0:30:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/09_mssng_sample_burden/tad_meta_v12_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/09_mssng_sample_burden/tad_meta_v12_%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04212026/09_mssng_sample_burden
python3 tad_replication_meta_analysis_v12.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
