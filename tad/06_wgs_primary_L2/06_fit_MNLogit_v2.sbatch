#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH -J MNLogit_v2
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /lustre12/home/kushima-pg/tad04212026/06_wgs_primary_L2/output_v4/output_MNLogit_v2/MNLogit_v2.out
#SBATCH -e /lustre12/home/kushima-pg/tad04212026/06_wgs_primary_L2/output_v4/output_MNLogit_v2/MNLogit_v2.err

# ファイル名: 06_fit_MNLogit_v2.sbatch
# 処理内容:
#   - Step 06 supplement v2: 40_fit_MNLogit_heterogeneity_v2.py を ncbn-cpu で実行
#   - v1 からの重要修正 (Option II critical fix):
#       * [CRITICAL] load_burden_table の return に sex_col を追加 (v1 では
#         sample_col が sex_col として誤って main() → run_one_class_one_svt に
#         渡されていたバグを修正)
#       * [CRITICAL] normalize_l2_burden_columns の戻り値 (df, rename_map)
#         を正しく tuple unpack (v1 では tuple を df に代入しており後段で落ちる)
#       * [ROBUST] fit_mnlogit_with_fallback で converged=True を要件化
#         (収束しない fit の結果を棄却、失敗理由を status で記録)
#       * [TRANSPARENCY] fit 失敗クラスを NA row で保存 (20 行は常に出力)
#       * [ASSERTION] Healthy=0 baseline と params 列順 (ASD, SZ) を実行時確認
#       * [WORDING] "Brown's method" → "correlation-adjusted Stouffer" に訂正
#       * [LOG] Cov(β_ASD, β_SZ) < 0 に WARN ログ
#   - ASD vs SZ heterogeneity test (MNLogit, 10 class × 2 SV = 20 fit + 4 global)
#   - 入力: Step 05 v3 (sample-level top-1% QC 適用済み, Pattern A=NAHR-excluded)
#       /lustre12/home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v3/
#         sample_burden_L2_and_specificity_v3.tsv
#   - 出力: tad04212026/06_wgs_primary_L2/output_v4/output_MNLogit_v2/
#       MNLogit_heterogeneity_results_v2.tsv
#       MNLogit_heterogeneity_summary_v2.tsv
#       MNLogit_heterogeneity_qc_v2.json
#   - 想定リソース: 24 fit + correlation matrix のみ (~5-15 分)
#       8 CPU / 32 GB / 2h で十分
#   - Step 05 v3 完了後に実行
#   - 実行時間を記録

set -euo pipefail

SCRIPT_DIR="/lustre12/home/kushima-pg/tad04212026/06_wgs_primary_L2"
OUT_DIR="${SCRIPT_DIR}/output_v4/output_MNLogit_v2"
mkdir -p "${OUT_DIR}"

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: ${SCRIPT_DIR}"
START=$(date +%s)

cd "${SCRIPT_DIR}"
python3 "${SCRIPT_DIR}/40_fit_MNLogit_heterogeneity_v2.py"

END=$(date +%s)
ELAPSED=$((END - START))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
