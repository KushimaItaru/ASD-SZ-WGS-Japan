#!/bin/bash
# ファイル名: 11_diff_all_vs_static_v2.sbatch
# 処理内容:
#   - Script 11 (42_compute_diff_all_vs_static_v2.py) を ncbn-cpu パーティションで実行
#     する sbatch ラッパー。
#   - 遺伝研スパコン用。実行時間を記録。
#   - 出力は tad04212026/11_diff_all_vs_static/output_v2/ に生成:
#       meta_diff_all_vs_static_v2.tsv (6 行: 3 exposures × {Diff_all, Static})
#   - Diff_all = Diff_specific_n1 + Diff_shared_n2plus を動的に列加算し、
#     Static と並列に B' Logit (arrayCGH) + GEE Indep (MSSNG) + IVW 2-way meta を fit。
#   - 依存: Script 10 v4 の出力
#       tad04212026/10_replication_2way_meta/output_v4/arraycgh_burden_L2_and_specificity_v4.tsv
#       tad04212026/10_replication_2way_meta/output_v4/mssng_burden_L2_and_specificity_v4.tsv
#   - fit モデル総数: (arrayCGH + MSSNG) × 3 exposures × 2 categories = 12 logistic/GEE fit
#     + 2-way meta 6 ペア
#   - 軽量処理のため 2 CPU / 8G / 15min 割当 (想定所要時間 1-2 分)

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_diffall_v2
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:15:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/11_diff_all_vs_static/tad_diffall_v2_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/11_diff_all_vs_static/tad_diffall_v2_%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04212026/11_diff_all_vs_static
python3 42_compute_diff_all_vs_static_v2.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
