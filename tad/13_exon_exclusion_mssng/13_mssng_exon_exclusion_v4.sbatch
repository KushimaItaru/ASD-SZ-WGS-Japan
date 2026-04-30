#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_mssng_exonexcl_v4
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04292026/13_exon_exclusion_mssng/tad_mssng_exonexcl_v4_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04292026/13_exon_exclusion_mssng/tad_mssng_exonexcl_v4_%j.err

# ファイル名: 13_mssng_exon_exclusion_v4.sbatch
# 処理内容:
#   - Module 13 v4 (42_compute_mssng_exon_exclusion_v4.py) を ncbn-cpu で実行する
#     sbatch ラッパー (tad04292026 配下、Q2 helper-based refactor 後)。
#   - 出力: tad04292026/13_exon_exclusion_mssng/output_v4/
#       mssng_exon_exclusion_v4.tsv         (全 7 行: A1 3 + A2 2 + A3 1 + A3s 1)
#       mssng_exon_exclusion_summary_v4.tsv (ASD vs sibling × DEL 4 discovery-positive)
#   - v3 -> v4 変更点:
#       1. [refactor] common/exon_helpers.py に重複コード切り出し
#          (build_merged_exon_bed / overlaps_any_exon / load_bin_l2_annotation /
#           L2_MEMBERSHIP_COLS / log; 関数本体は v3 と同一で数値結果は bit-identical)
#       2. [path] common.paths_v1 から PIPELINE_ROOT / GENCODE_GTF /
#          F_02_BIN_L2_ANNOTATION を import (envvar override 可、NIG default 保持)
#       3. [path] BIN_L2_ANNOT: bin_l2_annotation_v2.tsv.gz -> v3 (S2-S5 削除済;
#          Module 13 は L2_MEMBERSHIP_COLS のみ使うので数値影響なし)
#       4. [path] BASE_MSSNG_HEFFEL は os.environ で override 可
#          (default = /lustre12/home/kushima-pg/noncoding_tad_mssng_03132026)
#       5. 出力 dir / file: output_v3 -> output_v4
#       6. 解析ロジック (A1/A2/A3/A3-strict、GEE Binom Indep main + Firth sens、
#          BH/Holm、3 段 pre-check / MIN_CELL_COUNT=5) は v3 と完全一致
#   - 期待結果: tad04292026/output_v3 と数値 bit-identical
#   - 依存入力:
#       tad04292026/02_bin_l2_annotation/output_v3/bin_l2_annotation_v3.tsv.gz
#       /lustre12/home/kushima-pg/noncoding_tad_mssng_03132026/output_mssng_factorial_overlap_v1/
#           mssng_sample_event_bin_overlap_v1.tsv.gz
#           mssng_sample_covariates_v1.tsv
#       /lustre12/home/kushima-pg/annotationInfo/gencode.v46.annotation.gtf.gz
#   - fit モデル総数: 14 fit (7 bin_type × 2 methods、v3 と同一)
#   - 想定所要時間: 5-20 分 (v2 では 13.4 秒)

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04292026/13_exon_exclusion_mssng
mkdir -p output_v4
python3 42_compute_mssng_exon_exclusion_v4.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
