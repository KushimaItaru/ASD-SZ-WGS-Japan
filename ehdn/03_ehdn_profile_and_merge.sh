#!/usr/bin/env bash
# 03_ehdn_profile_and_merge.sh
# Filename: 03_ehdn_profile_and_merge.sh
# Description:
#   - Entry-point wrapper for EHdn profiling + merge
#   - Step 1: EHdn profile (run EHdn per CRAM, array job)
#   - Step 2: merge (merge and normalize locus.tsv from all samples)
#   - Steps are chained via sbatch dependency
#   - Record submitted Job IDs in manifest file
#
# Prerequisites:
#   - 02_ehdn_depth.sh must have completed (depths_all.tsv must exist)
#
# Usage:
#   bash 03_ehdn_profile_and_merge.sh
#
# Notes:
#   - Analytical logic in helper scripts is not modified
#   - Submit temp file with ARRAY_SIZE placeholder replaced via sed

set -euo pipefail

SCRIPT_START=$(date +%s)
TS(){ date "+%Y-%m-%d %H:%M:%S"; }

# ---- Paths (configurable via environment variables) ----
WRAPPER_ROOT="${WRAPPER_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
HELPER_DIR="${HELPER_DIR_EHDN:-${WRAPPER_ROOT}/helpers/ehdn}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-tre-burden}"
SAMPLE_LIST="${SAMPLE_LIST_EHDN:-${WRAPPER_ROOT}/sample_lists/ehdn_all_samples.tsv}"

LOG_DIR="${WRAPPER_ROOT}/ehdn/logs"
mkdir -p "${LOG_DIR}"

MANIFEST="${LOG_DIR}/03_profile_merge_manifest_$(date +%Y%m%d_%H%M%S).txt"
echo "[$(TS)] === 03_ehdn_profile_and_merge.sh ===" | tee "${MANIFEST}"
echo "[$(TS)] WRAPPER_ROOT=${WRAPPER_ROOT}" | tee -a "${MANIFEST}"
echo "[$(TS)] HELPER_DIR=${HELPER_DIR}" | tee -a "${MANIFEST}"

# ---- Step 1: EHdn profile (array job) ----
N_SAMPLES=$(tail -n +2 "${SAMPLE_LIST}" | wc -l | tr -d ' ')
echo "[$(TS)] [Step1] EHdn profile: N_SAMPLES=${N_SAMPLES}" | tee -a "${MANIFEST}"

# 02_run_ehdn has ARRAY_SIZE placeholder → create temp with actual size
TMP_EHDN=""
trap '[[ -n "${TMP_EHDN:-}" ]] && rm -f "${TMP_EHDN}"' EXIT
TMP_EHDN=$(mktemp "${LOG_DIR}/tmp_ehdn_XXXXXX.sh")
sed "s/ARRAY_SIZE/${N_SAMPLES}/g" "${HELPER_DIR}/02_run_ehdn_array_v2.sh" > "${TMP_EHDN}"

JOB1=$(sbatch --parsable \
  --partition="${SLURM_PARTITION:-ncbn-cpu}" \
  --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
  --output="${LOG_DIR}/ehdn_%A_%a.out" \
  --error="${LOG_DIR}/ehdn_%A_%a.err" \
  "${TMP_EHDN}")
echo "[$(TS)] [Step1] EHdn profile submitted: JobID=${JOB1}" | tee -a "${MANIFEST}"

# ---- Step 2: merge (single job, after profile) ----
JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} \
    --partition="${SLURM_PARTITION:-ncbn-cpu}" --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
    --cpus-per-task=4 --mem=64G --time=02:00:00 \
    --job-name=ehdn_merge \
    --output="${LOG_DIR}/merge_%j.out" \
    --error="${LOG_DIR}/merge_%j.err" \
    --wrap="bash -lc 'source ~/.bashrc && conda activate ${CONDA_ENV_NAME} && python3 ${HELPER_DIR}/04_merge_ehdn_novel_norm_v2.py'")
echo "[$(TS)] [Step2] merge submitted: JobID=${JOB2} (afterok:${JOB1})" | tee -a "${MANIFEST}"

# ---- Summary ----
SCRIPT_END=$(date +%s)
ELAPSED=$(( SCRIPT_END - SCRIPT_START ))
echo "" | tee -a "${MANIFEST}"
echo "[$(TS)] All steps submitted." | tee -a "${MANIFEST}"
echo "[$(TS)] Final job: ${JOB2}" | tee -a "${MANIFEST}"
echo "[$(TS)] Submission elapsed: ${ELAPSED}s" | tee -a "${MANIFEST}"
echo "[$(TS)] Manifest: ${MANIFEST}" | tee -a "${MANIFEST}"
