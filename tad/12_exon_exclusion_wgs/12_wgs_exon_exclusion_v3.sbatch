#!/bin/bash
# ファイル名: 12_wgs_exon_exclusion_v3.sbatch
# 処理内容:
#   - Script 12 v3 (42_compute_wgs_exon_exclusion_v3.py) を ncbn-cpu パーティション
#     で実行する sbatch ラッパー (遺伝研スパコン用)。
#   - 出力: tad04212026/12_exon_exclusion_wgs/output_v3/
#       wgs_exon_exclusion_v3.tsv         (全結果行; 2 comparisons × 2 SVT × 7 行)
#       wgs_exon_exclusion_summary_v3.tsv (ASD×DEL 4 discovery-positive 行サマリ)
#   - v2 -> v3 変更点 (path update のみ):
#       1. EVENT_OVERLAP: 04_wgs_sv_boundary_overlap/output_v9/
#                         sample_boundary_event_overlap_v9.tsv.gz
#          -> 04_wgs_sv_boundary_overlap/output_v10/
#             sample_boundary_event_overlap_v10.tsv.gz
#          (WGS discovery に sample-level top-1% CNV count QC 適用済み)
#       2. BURDEN_TABLE: 05_wgs_sample_burden/output_v2/
#                        sample_burden_L2_and_specificity_v2.tsv
#          -> 05_wgs_sample_burden/output_v3/
#             sample_burden_L2_and_specificity_v3.tsv
#       3. OUT_DIR: output_v2 -> output_v3
#       4. 出力ファイル名: wgs_exon_exclusion_v2.tsv -> v3
#                          wgs_exon_exclusion_summary_v2.tsv -> v3
#       5. 解析ロジック (GENCODE exon merge, A1/A2/A3/A3-strict, B' Logit main +
#          Firth sensitivity, BH/Holm on ASD×DEL discovery-positive) は v2 と完全一致
#   - 依存: Script 04 v10 + Script 05 v3 の出力
#       tad04212026/04_wgs_sv_boundary_overlap/output_v10/sample_boundary_event_overlap_v10.tsv.gz
#       tad04212026/05_wgs_sample_burden/output_v3/sample_burden_L2_and_specificity_v3.tsv
#       tad04212026/02_bin_l2_annotation/output_v2/bin_l2_annotation_v2.tsv.gz (変更なし)
#       /lustre12/home/kushima-pg/annotationInfo/gencode.v46.annotation.gtf.gz (変更なし)
#
#   - fit モデル総数: 2 comparisons × 2 SVT × (A1 3 bin + A2 2 bin + A3 1 + A3s 1)
#                    × 2 methods (Logit + Firth) = 56 fit (v2 と同一)
#
#   - GENCODE intersection + bin L2 交差 + sample aggregation が主負荷のため
#     4 CPU / 32G / 2h を割当 (想定所要時間 20-60 分)

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_wgs_exonexcl_v3
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/12_exon_exclusion_wgs/tad_wgs_exonexcl_v3_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/12_exon_exclusion_wgs/tad_wgs_exonexcl_v3_%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04212026/12_exon_exclusion_wgs
mkdir -p output_v3
python3 42_compute_wgs_exon_exclusion_v3.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
