#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH -J sb_L2_v5_patC
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 4:00:00
#SBATCH -o /lustre12/home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v5_patC/38_v5_patC.out
#SBATCH -e /lustre12/home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v5_patC/38_v5_patC.err

# ファイル名: 05_sample_burden_v5_patC.sbatch
# 処理内容:
#   - Step 5 v5 (PATTERN=patC) を ncbn-cpu で実行
#   - Pattern C (NAHR-excluded + |SV_length| < 1 MB) の sample x (L2 class +
#     specificity group) burden 計算 (感度分析)
#   - v5 変更点 (v4 からの Option II 改善):
#       * event_overlap に居るが cov に居ない sample を chunk 段で明示的に
#         WARN ログ + 除外 (silent sample loss 防止)
#       * cov に居るが event に 1 件も出現しない sample 件数を INFO ログ + summary 記録
#       * v4 の dead code (missing_from_cov = in_bur - in_cov) を削除し
#         assert (in_bur == in_cov) に置換
#       * sample_id で sort して deterministic 出力
#       * summary に QC counter (n_samples_cov_total, n_samples_cov_zero_events,
#         n_samples_event_only_dropped, n_events_kept, n_raw_rows) を追加
#       * 統計・カラム生成ロジックは v4 と完全一致
#   - 入力: Step 04 v11 patC 出力
#       /lustre12/home/kushima-pg/tad04212026/04_wgs_sv_boundary_overlap/output_v11_patC/
#         sample_boundary_event_overlap_v11_patC.tsv.gz
#         sample_mechanistic_burden_v11_patC.tsv
#   - 出力: output_v5_patC/
#       sample_burden_L2_and_specificity_v5_patC.tsv
#       sample_burden_L2_and_specificity_summary_v5_patC.tsv
#   - Size cap (<1MB) は Step 03 v9 で既に適用済
#   - Step 04 v11 patC 完了後に実行
#   - 実行時間を記録

set -euo pipefail
SCRIPT_DIR="/lustre12/home/kushima-pg/tad04212026/05_wgs_sample_burden"
OUT_DIR="${SCRIPT_DIR}/output_v5_patC"
mkdir -p "${OUT_DIR}"

export PATTERN="patC"

echo "Job start: $(date)"
echo "PATTERN: ${PATTERN}"
START=$(date +%s)

cd "${SCRIPT_DIR}"
python3 "${SCRIPT_DIR}/38_compute_sample_burden_L2_and_specificity_v5.py"

END=$(date +%s)
ELAPSED=$((END - START))
echo "Job end: $(date)"
echo "Elapsed: ${ELAPSED}s"
