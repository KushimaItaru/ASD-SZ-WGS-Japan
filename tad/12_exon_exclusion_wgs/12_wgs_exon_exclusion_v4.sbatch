#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_wgs_exonexcl_v4
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04292026/12_exon_exclusion_wgs/tad_wgs_exonexcl_v4_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04292026/12_exon_exclusion_wgs/tad_wgs_exonexcl_v4_%j.err

# ファイル名: 12_wgs_exon_exclusion_v4.sbatch
# 処理内容:
#   - Module 12 v4 (42_compute_wgs_exon_exclusion_v4.py) を ncbn-cpu で実行する
#     sbatch ラッパー (tad04292026 配下、Q2 helper-based refactor 後)。
#   - 出力: tad04292026/12_exon_exclusion_wgs/output_v4/
#       wgs_exon_exclusion_v4.tsv         (全結果行; 2 comparisons × 2 SVT × 7 行)
#       wgs_exon_exclusion_summary_v4.tsv (ASD×DEL 4 discovery-positive 行サマリ)
#   - v3 -> v4 変更点:
#       1. [refactor] common/exon_helpers.py に重複コード切り出し
#          (build_merged_exon_bed / overlaps_any_exon / load_bin_l2_annotation /
#           L2_MEMBERSHIP_COLS / log; 関数本体は v3 と同一で数値結果は bit-identical)
#       2. [path] common.paths_v1 から PIPELINE_ROOT / GENCODE_GTF /
#          F_02_BIN_L2_ANNOTATION を import (envvar override 可、NIG default 保持)
#       3. [path] BIN_L2_ANNOT: bin_l2_annotation_v2.tsv.gz -> v3 (S2-S5 削除済;
#          ただし Module 12 は L2_MEMBERSHIP_COLS のみ使うので数値影響なし)
#       4. [path] BURDEN_TABLE: output_v3/sample_burden_L2_and_specificity_v3.tsv
#          -> output_v4 (Q1 Pattern A primary; S2-S5 列削除済だが Module 12 は
#          使わないため数値影響なし)
#       5. 出力 dir / file: output_v3 -> output_v4
#       6. 解析ロジック (A1/A2/A3/A3-strict、Logit main + Firth sensitivity、
#          BH/Holm) は v3 と完全一致
#   - 期待結果: tad04212026/output_v3 と数値 bit-identical
#   - 依存: Module 02 v3 + Module 04 v10 + Module 05 v4 の出力 (NIG 上で完了済)
#       tad04292026/02_bin_l2_annotation/output_v3/bin_l2_annotation_v3.tsv.gz
#       tad04212026/04_wgs_sv_boundary_overlap/output_v10/sample_boundary_event_overlap_v10.tsv.gz
#         (※ Module 04 は tad04212026 から copy back 済)
#       tad04292026/05_wgs_sample_burden/output_v4/sample_burden_L2_and_specificity_v4.tsv
#       /lustre12/home/kushima-pg/annotationInfo/gencode.v46.annotation.gtf.gz
#   - fit モデル総数: 56 fit (v3 と同一)
#   - 想定所要時間: 20-60 分

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04292026/12_exon_exclusion_wgs
mkdir -p output_v4
python3 42_compute_wgs_exon_exclusion_v4.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
