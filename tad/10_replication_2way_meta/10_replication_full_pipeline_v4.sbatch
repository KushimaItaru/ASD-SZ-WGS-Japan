#!/bin/bash
# ファイル名: 10_replication_full_pipeline_v4.sbatch
# 処理内容:
#   - Script 10 (40_run_replication_full_pipeline_v4.py) を ncbn-cpu パーティションで実行
#     する sbatch ラッパー。
#   - 遺伝研スパコン用。実行時間を記録。
#   - 出力は tad04212026/10_replication_2way_meta/output_v4/ に生成。
#   - Step A (burden: arrayCGH + MSSNG)
#     Step B (fit: arrayCGH Logit + MSSNG GEE Independence)
#     Step C (2-way meta: arrayCGH + MSSNG)
#     Step D (3-way meta: WGS v3 + arrayCGH v22 + MSSNG)
#   - 依存: arrayCGH v22 output (tad04212026/08_arraycgh_sample_burden/output_v22/)
#           bin_l2_annotation v2 (tad04212026/02_bin_l2_annotation/output_v2/)
#           WGS primary v3 (tad04212026/06_wgs_primary_L2/output_v3/)
#           MSSNG overlap v1 (noncoding_tad_specificity_04022026/output_mssng_factorial_overlap_v2/)
#   - fit モデル総数: ACGH (L2 10 + spec 15) × exp 3 × cmp 2 + MSSNG (L2 10 + spec 15) × exp 3 × cmp 1
#                   = 150 + 75 = 225 logistic/GEE fit
#     + 2-way meta (225 ペア) + 3-way meta (75 WGS ペア)
#   - 保守的に 16 CPU / 64G / 2h 割当。

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_rep_v4
#SBATCH --cpus-per-task=16
#SBATCH --mem=12G
#SBATCH --time=2:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/10_replication_2way_meta/tad_rep_v4_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/10_replication_2way_meta/tad_rep_v4_%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04212026/10_replication_2way_meta
python3 40_run_replication_full_pipeline_v4.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
