#!/bin/bash
# ファイル名: 11_diff_any_vs_static_v5.sbatch
# 処理内容:
#   - Script 11 v5 (42_compute_diff_any_vs_static_v5.py) を ncbn-cpu パーティション
#     で実行する sbatch ラッパー (遺伝研スパコン用)。
#   - 出力: tad04292026/11_diff_all_vs_static/output_v5/
#       meta_diff_any_vs_static_v5.tsv (6 行: 3 exposures × {Diff_any, Static})
#   - v4 -> v5 変更点 (path update のみ):
#       1. REPLICATION_OUT_DIR: output_v4 -> output_v5
#          (WGS discovery に sample-level top-1% QC を適用した 10_v5 の出力を参照)
#       2. 入力 burden: arraycgh_burden_..._v4.tsv -> _v5.tsv
#                       mssng_burden_..._v4.tsv    -> _v5.tsv
#       3. OUT_DIR: output_v4 -> output_v5
#       4. 出力ファイル名: meta_diff_any_vs_static_v4.tsv -> v5
#       5. 解析ロジック (arrayCGH Logit, MSSNG GEE, Firth rank reduction, IVW/DL,
#          BH/Holm, Diff_any 列加算) は v4 と完全一致
#   - 依存: Script 10 v5 の出力
#       tad04292026/10_replication_2way_meta/output_v5/arraycgh_burden_L2_and_specificity_v5.tsv
#       tad04292026/10_replication_2way_meta/output_v5/mssng_burden_L2_and_specificity_v5.tsv

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_diffany_v5
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:15:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04292026/11_diff_all_vs_static/tad_diffany_v5_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04292026/11_diff_all_vs_static/tad_diffany_v5_%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04292026/11_diff_all_vs_static
mkdir -p output_v5
python3 42_compute_diff_any_vs_static_v5.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
