#!/usr/bin/env bash
# 04_tre_crosscaller_compare.sh
# Filename: 04_tre_crosscaller_compare.sh
# Description:
#   - Entry-point wrapper for EHdn / STRling case-case comparison (ASD vs SZ)
#   - Step 1: case-case comparison (load EHdn v19 + STRling v9 per_sample,
#             run Logistic + Poisson + heterogeneity test)
#   - Record submitted Job IDs in manifest file
#
# Prerequisites:
#   - EHdn burden (v19) and STRling burden (v9) must have completed
#     (per_sample.tsv must exist)
#
# Usage:
#   bash 04_tre_crosscaller_compare.sh

set -euo pipefail

SCRIPT_START=$(date +%s)
TS(){ date "+%Y-%m-%d %H:%M:%S"; }

# ---- Paths (configurable via environment variables) ----
WRAPPER_ROOT="${WRAPPER_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
CROSSCALLER_DIR="${HELPER_DIR_CROSSCALLER:-${WRAPPER_ROOT}/helpers/crosscaller}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-tre-burden}"

LOG_DIR="${WRAPPER_ROOT}/crosscaller/logs"
mkdir -p "${LOG_DIR}"

MANIFEST="${LOG_DIR}/04_crosscaller_manifest_$(date +%Y%m%d_%H%M%S).txt"
echo "[$(TS)] === 04_tre_crosscaller_compare.sh ===" | tee "${MANIFEST}"
echo "[$(TS)] WRAPPER_ROOT=${WRAPPER_ROOT}" | tee -a "${MANIFEST}"
echo "[$(TS)] CROSSCALLER_DIR=${CROSSCALLER_DIR}" | tee -a "${MANIFEST}"

# ---- Step 1: case-case comparison (single job) ----
JOB1=$(sbatch --parsable \
    --partition="${SLURM_PARTITION:-ncbn-cpu}" --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
    --cpus-per-task=4 --mem=32G --time=00:30:00 \
    --job-name=tre_case_case_v3 \
    --output="${LOG_DIR}/case_case_v3_%j.out" \
    --error="${LOG_DIR}/case_case_v3_%j.err" \
    --wrap="bash -lc 'source ~/.bashrc && conda activate ${CONDA_ENV_NAME} && python3 ${CROSSCALLER_DIR}/20_tre_case_case_comparison_v3.py'")
echo "[$(TS)] [Step1] case-case v3 submitted: JobID=${JOB1}" | tee -a "${MANIFEST}"

# ---- Summary ----
SCRIPT_END=$(date +%s)
ELAPSED=$(( SCRIPT_END - SCRIPT_START ))
echo "" | tee -a "${MANIFEST}"
echo "[$(TS)] All steps submitted." | tee -a "${MANIFEST}"
echo "[$(TS)] Case-case job: ${JOB1}" | tee -a "${MANIFEST}"
echo "[$(TS)] Submission elapsed: ${ELAPSED}s" | tee -a "${MANIFEST}"
echo "[$(TS)] Manifest: ${MANIFEST}" | tee -a "${MANIFEST}"
