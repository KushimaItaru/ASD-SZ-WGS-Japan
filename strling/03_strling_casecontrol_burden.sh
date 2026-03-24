#!/usr/bin/env bash
# 03_strling_casecontrol_burden.sh
# Filename: 03_strling_casecontrol_burden.sh
# Description:
#   - Entry-point wrapper for STRling burden + QC/sensitivity analysis
#   - Step 1: burden analysis (5-fold cross-fit, Logistic + Poisson GLM)
#   - Step 2: QC/sensitivity (per-group summary, SMD, winsorize/trim sensitivity)
#   - Steps are chained via sbatch dependency
#   - Record submitted Job IDs in manifest file
#
# Prerequisites:
#   - 02_strling_casecontrol_call_and_outliers.sh must have completed (STRs.tsv must exist)
#
# Usage:
#   bash 03_strling_casecontrol_burden.sh

set -euo pipefail

SCRIPT_START=$(date +%s)
TS(){ date "+%Y-%m-%d %H:%M:%S"; }

# ---- Paths (configurable via environment variables) ----
WRAPPER_ROOT="${WRAPPER_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
HELPER_DIR="${HELPER_DIR_STRLING:-${WRAPPER_ROOT}/helpers/strling}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-tre-burden}"
CONFIG="${HELPER_DIR}/00_strling_genomewide_config_v1.sh"

source "${CONFIG}"

LOG_DIR="${WRAPPER_ROOT}/strling/logs"
mkdir -p "${LOG_DIR}"

MANIFEST="${LOG_DIR}/03_burden_manifest_$(date +%Y%m%d_%H%M%S).txt"
echo "[$(TS)] === 03_strling_casecontrol_burden.sh ===" | tee "${MANIFEST}"
echo "[$(TS)] WRAPPER_ROOT=${WRAPPER_ROOT}" | tee -a "${MANIFEST}"
echo "[$(TS)] HELPER_DIR=${HELPER_DIR}" | tee -a "${MANIFEST}"

# ---- Step 1: burden analysis (single job) ----
JOB1=$(sbatch --parsable \
    --partition="${SLURM_PARTITION:-ncbn-cpu}" --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
    --cpus-per-task=8 --mem=64G --time=01:00:00 \
    --job-name=strling_burden_v9 \
    --output="${LOG_DIR}/burden_v9_%j.out" \
    --error="${LOG_DIR}/burden_v9_%j.err" \
    --wrap="bash -lc 'source ~/.bashrc && conda activate ${CONDA_ENV_NAME} && python3 ${HELPER_DIR}/08_strling_outlier_burden_rare_casecontrol_crossfit_v9.py --merge_dist 1000'")
echo "[$(TS)] [Step1] burden v9 submitted: JobID=${JOB1}" | tee -a "${MANIFEST}"

# ---- Step 2: QC / sensitivity (single job, after burden) ----
JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} \
    --partition="${SLURM_PARTITION:-ncbn-cpu}" --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
    --cpus-per-task=4 --mem=32G --time=00:30:00 \
    --job-name=strling_qc_v4 \
    --output="${LOG_DIR}/qc_v4_%j.out" \
    --error="${LOG_DIR}/qc_v4_%j.err" \
    --wrap="bash -lc 'source ~/.bashrc && conda activate ${CONDA_ENV_NAME} && python3 ${HELPER_DIR}/09_strling_qc_sensitivity_rare_inbounds_v4.py'")
echo "[$(TS)] [Step2] QC/sensitivity v4 submitted: JobID=${JOB2} (afterok:${JOB1})" | tee -a "${MANIFEST}"

# ---- Summary ----
SCRIPT_END=$(date +%s)
ELAPSED=$(( SCRIPT_END - SCRIPT_START ))
echo "" | tee -a "${MANIFEST}"
echo "[$(TS)] All steps submitted." | tee -a "${MANIFEST}"
echo "[$(TS)] Final QC job: ${JOB2}" | tee -a "${MANIFEST}"
echo "[$(TS)] Submission elapsed: ${ELAPSED}s" | tee -a "${MANIFEST}"
echo "[$(TS)] Manifest: ${MANIFEST}" | tee -a "${MANIFEST}"
