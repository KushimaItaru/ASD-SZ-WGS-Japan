#!/bin/bash
# ファイル名: 11_diff_any_vs_static_v4.sbatch
# 処理内容:
#   - Script 11 v4 (42_compute_diff_any_vs_static_v4.py) を ncbn-cpu パーティション
#     で実行する sbatch ラッパー (遺伝研スパコン用)。
#   - 出力: tad04212026/11_diff_all_vs_static/output_v4/
#       meta_diff_any_vs_static_v4.tsv (6 行: 3 exposures × {Diff_any, Static})
#   - v3 -> v4 変更点 (Firth fix):
#       1. Firth に QR-with-pivoting ベースの rank reduction を追加
#          (v3 では platform one-hot (5 列) + intercept の完全共線性で全 6 行
#           singular_hessian fail していたのを修正)
#       2. 新規: _reduce_rank_via_qr ヘルパー
#       3. firth_logit_fit に protected_idx=(0, 1) を追加 (intercept + exposure 保持)
#       4. OUT_DIR: output_v3 -> output_v4
#       5. 出力ファイル名: meta_diff_any_vs_static_v3.tsv -> _v4.tsv
#       6. Main 解析ロジック (arrayCGH Logit, MSSNG GEE, fixed/RE meta, BH/Holm,
#          Diff_any 列加算) は v3 と完全一致 (main 数値は不変)
#   - 依存: Script 10 v4 の出力
#       tad04212026/10_replication_2way_meta/output_v4/arraycgh_burden_L2_and_specificity_v4.tsv
#       tad04212026/10_replication_2way_meta/output_v4/mssng_burden_L2_and_specificity_v4.tsv
#   - fit モデル総数: (arrayCGH Logit + MSSNG GEE + MSSNG Firth) × 3 exposures × 2 categories
#     = 18 fit + 2-way meta 6 ペア × 2 (fixed + random-effects)
#   - 軽量処理のため 2 CPU / 8G / 15min 割当 (想定所要時間 1-2 分)

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_diffany_v4
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:15:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/11_diff_all_vs_static/tad_diffany_v4_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/11_diff_all_vs_static/tad_diffany_v4_%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04212026/11_diff_all_vs_static
python3 42_compute_diff_any_vs_static_v4.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
