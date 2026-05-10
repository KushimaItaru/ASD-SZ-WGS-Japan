#!/bin/bash
# ファイル名: 12_wgs_exon_exclusion_v2.sbatch
# 処理内容:
#   - Script 12 v2 (42_compute_wgs_exon_exclusion_v2.py) を ncbn-cpu パーティション
#     で実行する sbatch ラッパー (遺伝研スパコン用)。
#   - 出力: tad04212026/12_exon_exclusion_wgs/output_v2/
#       wgs_exon_exclusion_v2.tsv         (全結果行; 2 comparisons × 2 SVT × 7 行)
#       wgs_exon_exclusion_summary_v2.tsv (ASD×DEL 4 discovery-positive 行サマリ)
#   - 主処理:
#       1. GENCODE v46 protein-coding exon を merge
#       2. bin_l2_annotation v2 と exon を交差し exon_free/exon_overlap 判定
#          (any_l2_bins のみ解析対象)
#       3. sample_boundary_event_overlap v9 を読み、SV 全体 (sv_chr/start0/end)
#          と exon の交差から exonfreeCNV を判定 (A3 軸)
#       4. 5 stratification × 2 SVT の bin 数を sample × SV type で unique 集計
#       5. burden v2 に exposure 列を追加 → B' Logit (main) + Firth (sensitivity)
#          回帰を 2 comparisons × 2 SVT × 4 analyses のグリッドで実施
#       6. Multiple testing: ASD × DEL × 4 discovery-positive 行 (A1 exonfree,
#          A2 exonfree_in_joint, A3 exonfreeCNV, A3-strict) に BH-FDR + Holm 適用
#
#   - v1 (heffel run_exonfree_aggregate_sensitivity_v4.py) -> v2 変更点:
#       1. PATHS: heffel heffel_deep_analysis_03242026/ -> tad04212026/
#          (v9 boundary/event + v2 bin_l2_annotation + v2 sample burden)
#       2. L2 判定: heffel v16 boundary master 内 10 列 ->
#          tad04212026 v2 bin_l2_annotation (独立ファイル) との merge
#       3. [追加] 自前 Firth penalized logistic (Script 11 v4 ポーティング)
#          + QR-with-pivoting rank reduction (intercept + exposure 保護)
#       4. [追加] BH-FDR + Holm 多重検定補正 (Logit P, Firth P それぞれ)
#       5. [追加] carrier_case / carrier_ctrl / fit_status / fit_status_firth 列
#       6. log 形式を Script 11 v4 スタイル [YYYY-MM-DD HH:MM:SS] に統一
#       7. 出力ファイル名: exonfree_aggregate_sensitivity_patternA_v4_results.tsv
#          -> wgs_exon_exclusion_v2.tsv + wgs_exon_exclusion_summary_v2.tsv
#       8. 主要解析ロジック (A1/A2/A3/A3-strict, B' Logit main) は heffel v4 と同一
#          (数値一致を担保)
#
#   - 依存入力:
#       tad04212026/02_bin_l2_annotation/output_v2/bin_l2_annotation_v2.tsv.gz
#       tad04212026/04_wgs_sv_boundary_overlap/output_v9/sample_boundary_event_overlap_v9.tsv.gz
#       tad04212026/05_wgs_sample_burden/output_v2/sample_burden_L2_and_specificity_v2.tsv
#       /lustre12/home/kushima-pg/annotationInfo/gencode.v46.annotation.gtf.gz
#
#   - fit モデル総数: 2 comparisons × 2 SVT × (A1 3 bin + A2 2 bin + A3 1 + A3s 1)
#                    × 2 methods (Logit + Firth) = 56 fit
#
#   - GENCODE intersection + bin L2 交差 + sample aggregation が主負荷のため
#     4 CPU / 32G / 2h を割当 (想定所要時間 20-60 分)

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_wgs_exonexcl_v2
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/12_exon_exclusion_wgs/tad_wgs_exonexcl_v2_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/12_exon_exclusion_wgs/tad_wgs_exonexcl_v2_%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04212026/12_exon_exclusion_wgs
python3 42_compute_wgs_exon_exclusion_v2.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
