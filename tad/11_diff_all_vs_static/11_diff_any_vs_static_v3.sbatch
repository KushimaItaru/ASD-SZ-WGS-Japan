#!/bin/bash
# ファイル名: 11_diff_any_vs_static_v3.sbatch
# 処理内容:
#   - Script 11 v3 (42_compute_diff_any_vs_static_v3.py) を ncbn-cpu パーティション
#     で実行する sbatch ラッパー (遺伝研スパコン用)。
#   - 出力: tad04212026/11_diff_all_vs_static/output_v3/
#       meta_diff_any_vs_static_v3.tsv (6 行: 3 exposures × {Diff_any, Static})
#   - v2 -> v3 変更点 (ChatGPT レビュー反映):
#       1. Diff_all -> Diff_any に命名統一 (3 exposure)
#       2. carrier_boundary は binary OR (0/1) に clip
#       3. DerSimonian-Laird random-effects meta を追加
#       4. BH + Holm multiple testing 補正を追加 (6 tests 一律)
#       5. MSSNG Firth penalized logistic sensitivity を追加 (pooled, 自前実装)
#       6. OUT_DIR: output_v2 -> output_v3
#       7. 出力ファイル名: meta_diff_all_vs_static_v2.tsv -> meta_diff_any_vs_static_v3.tsv
#   - 依存: Script 10 v4 の出力
#       tad04212026/10_replication_2way_meta/output_v4/arraycgh_burden_L2_and_specificity_v4.tsv
#       tad04212026/10_replication_2way_meta/output_v4/mssng_burden_L2_and_specificity_v4.tsv
#   - fit モデル総数: (arrayCGH Logit + MSSNG GEE + MSSNG Firth) × 3 exposures × 2 categories
#     = 18 fit + 2-way meta 6 ペア × 2 (fixed + random-effects)
#   - 軽量処理のため 2 CPU / 8G / 15min 割当 (想定所要時間 1-2 分、Firth 分で若干増加)

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_diffany_v3
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:15:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/11_diff_all_vs_static/tad_diffany_v3_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/11_diff_all_vs_static/tad_diffany_v3_%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04212026/11_diff_all_vs_static
python3 42_compute_diff_any_vs_static_v3.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
