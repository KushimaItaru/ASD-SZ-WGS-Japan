#!/usr/bin/env bash
# 02_strling_casecontrol_call_and_outliers.sh
# Filename: 02_strling_casecontrol_call_and_outliers.sh
# Description:
#   - Entry-point wrapper for STRling genotyping -> in-bounds filter -> outlier detection
#   - Step 1: call (STRling call on genic panel, array job)
#   - Step 2: in-bounds filter (filter call results by genic bounds)
#   - Step 3: outlier detection (compute z-scores via strling-outliers)
#   - Steps are chained via sbatch dependency
#   - Record submitted Job IDs in manifest file
#
# Prerequisites:
#   - 01_strling_build_panel.sh must have completed (genic bounds must exist)
#
# Usage:
#   bash 02_strling_casecontrol_call_and_outliers.sh

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

MANIFEST="${LOG_DIR}/02_call_outliers_manifest_$(date +%Y%m%d_%H%M%S).txt"
echo "[$(TS)] === 02_strling_casecontrol_call_and_outliers.sh ===" | tee "${MANIFEST}"
echo "[$(TS)] WRAPPER_ROOT=${WRAPPER_ROOT}" | tee -a "${MANIFEST}"
echo "[$(TS)] HELPER_DIR=${HELPER_DIR}" | tee -a "${MANIFEST}"

# ---- Step 1: call (array job) ----
# 06_call reads casecontrol_samples.tsv, so sample count is determined from it
CASECONTROL_TSV="${SAMPLE_LIST_CASECONTROL:-${WRAPPER_ROOT}/sample_lists/casecontrol_samples.tsv}"
N_SAMPLES=$(awk 'END{print NR-1}' "${CASECONTROL_TSV}")
echo "[$(TS)] [Step1] call: N_SAMPLES=${N_SAMPLES} (from casecontrol_samples.tsv)" | tee -a "${MANIFEST}"

# 06_call has ARRAY_SIZE placeholder → create temp with actual size
TMP_CALL=""
trap '[[ -n "${TMP_CALL:-}" ]] && rm -f "${TMP_CALL}"' EXIT
TMP_CALL=$(mktemp "${LOG_DIR}/tmp_call_XXXXXX.sh")
sed "s/ARRAY_SIZE/${N_SAMPLES}/g" "${HELPER_DIR}/06_strling_call_array_genic_v1.sh" > "${TMP_CALL}"

JOB1=$(sbatch --parsable \
  --partition="${SLURM_PARTITION:-ncbn-cpu}" \
  --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
  --output="${LOG_DIR}/call_%A_%a.out" \
  --error="${LOG_DIR}/call_%A_%a.err" \
  "${TMP_CALL}")
echo "[$(TS)] [Step1] call submitted: JobID=${JOB1}" | tee -a "${MANIFEST}"

# ---- Step 2: in-bounds filter (single job, after call) ----
JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} \
    --partition="${SLURM_PARTITION:-ncbn-cpu}" --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
    --cpus-per-task=1 --mem=32G --time=00:30:00 \
    --job-name=strling_inbounds \
    --output="${LOG_DIR}/inbounds_%j.out" \
    --error="${LOG_DIR}/inbounds_%j.err" \
    --wrap="bash -lc 'source ~/.bashrc && conda activate ${CONDA_ENV_NAME} && python3 ${HELPER_DIR}/10_make_calls_genic_inbounds_v1.py'")
echo "[$(TS)] [Step2] in-bounds filter submitted: JobID=${JOB2} (afterok:${JOB1})" | tee -a "${MANIFEST}"

# ---- Step 3: outlier detection (single job, after in-bounds) ----
JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} \
    --partition="${SLURM_PARTITION:-ncbn-cpu}" --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
    --output="${LOG_DIR}/outliers_%j.out" \
    --error="${LOG_DIR}/outliers_%j.err" \
    "${HELPER_DIR}/11_strling_outliers_casecontrol_inbounds_v1.sh")
echo "[$(TS)] [Step3] outlier detection submitted: JobID=${JOB3} (afterok:${JOB2})" | tee -a "${MANIFEST}"

# ---- Summary ----
SCRIPT_END=$(date +%s)
ELAPSED=$(( SCRIPT_END - SCRIPT_START ))
echo "" | tee -a "${MANIFEST}"
echo "[$(TS)] All steps submitted." | tee -a "${MANIFEST}"
echo "[$(TS)] Final outlier job: ${JOB3}" | tee -a "${MANIFEST}"
echo "[$(TS)] Submission elapsed: ${ELAPSED}s" | tee -a "${MANIFEST}"
echo "[$(TS)] Manifest: ${MANIFEST}" | tee -a "${MANIFEST}"
