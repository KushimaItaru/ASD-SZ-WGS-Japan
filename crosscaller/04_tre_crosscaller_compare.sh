#!/usr/bin/env bash
# 04_tre_crosscaller_compare.sh
# ファイル名: 04_tre_crosscaller_compare.sh
# 処理内容:
#   - EHdn / STRling の case-case 比較（ASD vs SZ）の入口ラッパー
#   - Step 1: case-case comparison（EHdn v19 + STRling v9 の per_sample を読み込み、
#             Logistic + Poisson + heterogeneity test を実施）
#   - 投入した Job ID を manifest ファイルに記録
#
# 前提:
#   - EHdn burden (v19) と STRling burden (v9) が完了していること
#     （per_sample.tsv が存在すること）
#
# 使い方:
#   bash 04_tre_crosscaller_compare.sh

set -euo pipefail

SCRIPT_START=$(date +%s)
TS(){ date "+%Y-%m-%d %H:%M:%S"; }

# ---- Paths ----
WRAPPER_ROOT="/home/kushima-pg/str_03242026"
CROSSCALLER_DIR="/home/kushima-pg/str_12282025/ehdn-strling_01022026"

LOG_DIR="${WRAPPER_ROOT}/crosscaller/logs"
mkdir -p "${LOG_DIR}"

MANIFEST="${LOG_DIR}/04_crosscaller_manifest_$(date +%Y%m%d_%H%M%S).txt"
echo "[$(TS)] === 04_tre_crosscaller_compare.sh ===" | tee "${MANIFEST}"
echo "[$(TS)] WRAPPER_ROOT=${WRAPPER_ROOT}" | tee -a "${MANIFEST}"
echo "[$(TS)] CROSSCALLER_DIR=${CROSSCALLER_DIR}" | tee -a "${MANIFEST}"

# ---- Step 1: case-case comparison (single job) ----
JOB1=$(sbatch --parsable \
    --partition=ncbn-cpu --account=ncbn-cpu \
    --cpus-per-task=4 --mem=32G --time=00:30:00 \
    --job-name=tre_case_case_v3 \
    --output="${LOG_DIR}/case_case_v3_%j.out" \
    --error="${LOG_DIR}/case_case_v3_%j.err" \
    --wrap="bash -lc 'source ~/.bashrc && conda activate ngs && cd ${CROSSCALLER_DIR} && python3 20_tre_case_case_comparison_v3.py'" )
echo "[$(TS)] [Step1] case-case v3 submitted: JobID=${JOB1}" | tee -a "${MANIFEST}"

# ---- Summary ----
SCRIPT_END=$(date +%s)
ELAPSED=$(( SCRIPT_END - SCRIPT_START ))
echo "" | tee -a "${MANIFEST}"
echo "[$(TS)] All steps submitted." | tee -a "${MANIFEST}"
echo "[$(TS)] Case-case job: ${JOB1}" | tee -a "${MANIFEST}"
echo "[$(TS)] Submission elapsed: ${ELAPSED}s" | tee -a "${MANIFEST}"
echo "[$(TS)] Manifest: ${MANIFEST}" | tee -a "${MANIFEST}"
