#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_03_sv_v9
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04292026/03_wgs_sv_events/output_v9/03_sv_v9_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04292026/03_wgs_sv_events/output_v9/03_sv_v9_%j.err

# ファイル名: 03_extract_sv_v9.sbatch
# 処理内容:
#   - Step 3 v9 (34a_extract_rare_sv_events_v9.py) を ncbn-cpu で実行
#   - 1 run で Pattern A/B/C の rare SV event table を生成する感度分析
#     (v20 unified_heffel_tad_pipeline_v20 の Pattern B/C 対応を新パイプに移植)
#
#   - v8 -> v9 変更点:
#       1. 1 run で 3 Pattern の event table を生成
#          patA: NAHR-excluded (main analysis, v8 と挙動同一)
#          patB: NAHR-included (sensitivity)
#          patC: NAHR-excluded + |SV_length| < 1,000,000 bp (sensitivity)
#       2. 出力ディレクトリ: output_v8 -> output_v9
#       3. 出力ファイル suffix: _v8 -> _v9_patA / _v9_patB / _v9_patC
#       4. NAHR annotation は patB でも保持 (is_nahr_gd_cnv 列)
#       5. per-sample burden は Pattern ごとに再計算 (event set が異なるため)
#       6. common QC step は共通 TSV (overall_qc.tsv) にログ
#
#   - 出力:
#       output_v9/
#         wgs_rare_sv_events_v9_patA.tsv.gz / _patA_summary.tsv
#         wgs_rare_sv_events_v9_patB.tsv.gz / _patB_summary.tsv
#         wgs_rare_sv_events_v9_patC.tsv.gz / _patC_summary.tsv
#         wgs_rare_sv_events_v9_overall_qc.tsv
#         wgs_rare_sv_events_v9_summary.json
#
#   - リソース想定:
#       - AnnotSV 読み込み + exon index 構築 + 3 Pattern の per-sample burden 再計算
#       - v8 (1 Pattern) で ~10-15 分 → v9 (3 Pattern) で ~20-30 分目安
#       - メモリ: AnnotSV raw + exon index + sv_df コピー → 64GB で十分
#       - CPU: 逐次処理中心 (8 CPU 割当)
#
#   - 依存:
#       paths_v1.ANNOTSV_TABLE / SAMPLE_INFO / COMMON_CNV_TABLE / SEGDUP_BED /
#       EXCLUSION_BED / PCA_EIGENVEC / CNV_SAMPLE_COUNTS / GENCODE_GTF /
#       CURATED_GD_FILE
#
#   - 実行時間を記録

set -euo pipefail

SCRIPT_DIR="/lustre12/home/kushima-pg/tad04292026/03_wgs_sv_events"
OUT_DIR="${SCRIPT_DIR}/output_v9"
mkdir -p "${OUT_DIR}"

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: ${SCRIPT_DIR}"
START_TS=$(date +%s)

cd "${SCRIPT_DIR}"
python3 "${SCRIPT_DIR}/34a_extract_rare_sv_events_v9.py"

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
