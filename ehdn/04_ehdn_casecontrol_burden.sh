#!/usr/bin/env bash
# 04_ehdn_casecontrol_burden.sh
# ファイル名: 04_ehdn_casecontrol_burden.sh
# 処理内容:
#   - EHdn burden 解析の入口ラッパー
#   - Step 1: 5-fold cross-fitted rare TRE burden 解析
#   - Step 2: burden statistical test
#             （logistic OR, Poisson RR with offset=log(observed_clusters_total),
#               Mann–Whitney U test）
#   - 各ステップは sbatch dependency で連鎖
#   - 投入した Job ID を manifest ファイルに記録
#
# 前提:
#   - 03_ehdn_profile_and_merge.sh が完了していること
#     （novel_loci_genic_norm.tsv が存在すること）
#
# 使い方:
#   bash 04_ehdn_casecontrol_burden.sh
#
# 注意:
#   - 既存 helper script の処理内容は一切変更しない
#   - 15_annotate_genes / 16_gene_significance_test は主解析ラインに含めない
#     （exploratory / legacy utility）

set -euo pipefail

SCRIPT_START=$(date +%s)
TS(){ date "+%Y-%m-%d %H:%M:%S"; }

# ---- Paths ----
WRAPPER_ROOT="/home/kushima-pg/str_03242026"
HELPER_DIR="/home/kushima-pg/str_12282025/ehdn"

LOG_DIR="${WRAPPER_ROOT}/ehdn/logs"
mkdir -p "${LOG_DIR}"

MANIFEST="${LOG_DIR}/04_burden_manifest_$(date +%Y%m%d_%H%M%S).txt"
echo "[$(TS)] === 04_ehdn_casecontrol_burden.sh ===" | tee "${MANIFEST}"
echo "[$(TS)] WRAPPER_ROOT=${WRAPPER_ROOT}" | tee -a "${MANIFEST}"
echo "[$(TS)] HELPER_DIR=${HELPER_DIR}" | tee -a "${MANIFEST}"

# ---- Step 1: burden analysis (single job) ----
JOB1=$(sbatch --parsable \
    --partition=ncbn-cpu --account=ncbn-cpu \
    --cpus-per-task=8 --mem=128G --time=06:00:00 \
    --job-name=ehdn_burden_v19 \
    --output="${LOG_DIR}/burden_v19_%j.out" \
    --error="${LOG_DIR}/burden_v19_%j.err" \
    --wrap="bash -lc 'source ~/.bashrc && conda activate ngs && cd ${HELPER_DIR} && python3 14_outlier_burden_rare_casecontrol_crossfit_v19.py'")
echo "[$(TS)] [Step1] burden v19 submitted: JobID=${JOB1}" | tee -a "${MANIFEST}"

# ---- Step 2: burden statistical test (single job, after burden) ----
JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} \
    --partition=ncbn-cpu --account=ncbn-cpu \
    --cpus-per-task=4 --mem=64G --time=02:00:00 \
    --job-name=ehdn_burden_stats_v20 \
    --output="${LOG_DIR}/burden_stats_v20_%j.out" \
    --error="${LOG_DIR}/burden_stats_v20_%j.err" \
    --wrap="bash -lc 'source ~/.bashrc && conda activate ngs && cd ${HELPER_DIR} && python3 17_burden_statistical_test_v20.py'")
echo "[$(TS)] [Step2] burden stats v20 submitted: JobID=${JOB2} (afterok:${JOB1})" | tee -a "${MANIFEST}"

# ---- Summary ----
SCRIPT_END=$(date +%s)
ELAPSED=$(( SCRIPT_END - SCRIPT_START ))
echo "" | tee -a "${MANIFEST}"
echo "[$(TS)] All steps submitted." | tee -a "${MANIFEST}"
echo "[$(TS)] Final job: ${JOB2}" | tee -a "${MANIFEST}"
echo "[$(TS)] Submission elapsed: ${ELAPSED}s" | tee -a "${MANIFEST}"
echo "[$(TS)] Manifest: ${MANIFEST}" | tee -a "${MANIFEST}"
