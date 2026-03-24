#!/usr/bin/env bash
# 02_ehdn_depth.sh
# Filename: 02_ehdn_depth.sh
# Description:
#   - Entry-point wrapper for EHdn / STRling shared depth computation pipeline
#   - Step 1: depth calculation (samtools depth, array job)
#   - Step 2: depth aggregation (consolidate all sample depths into one file)
#   - Steps are chained via sbatch dependency
#   - Record submitted Job IDs in manifest file
#
# Prerequisites:
#   - 01_ehdn_setup_and_sample_lists.sh must have completed
#     (ehdn_all_samples.tsv must exist)
#
# Usage:
#   bash 02_ehdn_depth.sh
#
# Notes:
#   - Analytical logic in helper scripts is not modified
#   - Submit temp file with ARRAY_SIZE placeholder replaced via sed

set -euo pipefail

SCRIPT_START=$(date +%s)
TS(){ date "+%Y-%m-%d %H:%M:%S"; }

# ---- Paths (configurable via environment variables) ----
WRAPPER_ROOT="${WRAPPER_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
HELPER_DIR="${HELPER_DIR_STRLING:-${WRAPPER_ROOT}/helpers/strling}"
SAMPLE_LIST="${SAMPLE_LIST_EHDN:-${WRAPPER_ROOT}/sample_lists/ehdn_all_samples.tsv}"

LOG_DIR="${WRAPPER_ROOT}/ehdn/logs"
mkdir -p "${LOG_DIR}"

MANIFEST="${LOG_DIR}/02_depth_manifest_$(date +%Y%m%d_%H%M%S).txt"
echo "[$(TS)] === 02_ehdn_depth.sh ===" | tee "${MANIFEST}"
echo "[$(TS)] WRAPPER_ROOT=${WRAPPER_ROOT}" | tee -a "${MANIFEST}"
echo "[$(TS)] HELPER_DIR=${HELPER_DIR}" | tee -a "${MANIFEST}"

# ---- Step 1: depth array job ----
N_SAMPLES=$(tail -n +2 "${SAMPLE_LIST}" | wc -l | tr -d ' ')
echo "[$(TS)] [Step1] depth: N_SAMPLES=${N_SAMPLES}" | tee -a "${MANIFEST}"

# 01_calc_depth has ARRAY_SIZE placeholder → create temp with actual size
TMP_DEPTH=""
trap '[[ -n "${TMP_DEPTH:-}" ]] && rm -f "${TMP_DEPTH}"' EXIT
TMP_DEPTH=$(mktemp "${LOG_DIR}/tmp_depth_XXXXXX.sh")
sed "s/ARRAY_SIZE/${N_SAMPLES}/g" "${HELPER_DIR}/01_calc_depth_array_fast_v2.sh" > "${TMP_DEPTH}"

JOB1=$(sbatch --parsable \
  --partition="${SLURM_PARTITION:-ncbn-cpu}" \
  --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
  --output="${LOG_DIR}/depth_%A_%a.out" \
  --error="${LOG_DIR}/depth_%A_%a.err" \
  "${TMP_DEPTH}")
echo "[$(TS)] [Step1] depth submitted: JobID=${JOB1}" | tee -a "${MANIFEST}"

# ---- Step 2: collect depths (single job, after depth array) ----
JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} \
    --partition="${SLURM_PARTITION:-ncbn-cpu}" --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
    --cpus-per-task=1 --mem=8G --time=00:30:00 \
    --job-name=ehdn_collect_depths \
    --output="${LOG_DIR}/collect_depths_%j.out" \
    --error="${LOG_DIR}/collect_depths_%j.err" \
    --wrap="bash -lc 'source ~/.bashrc && conda activate ngs && python3 ${HELPER_DIR}/03_collect_depths_v1.py'")
echo "[$(TS)] [Step2] collect_depths submitted: JobID=${JOB2} (afterok:${JOB1})" | tee -a "${MANIFEST}"

# ---- Summary ----
SCRIPT_END=$(date +%s)
ELAPSED=$(( SCRIPT_END - SCRIPT_START ))
echo "" | tee -a "${MANIFEST}"
echo "[$(TS)] All steps submitted." | tee -a "${MANIFEST}"
echo "[$(TS)] Final job: ${JOB2}" | tee -a "${MANIFEST}"
echo "[$(TS)] Submission elapsed: ${ELAPSED}s" | tee -a "${MANIFEST}"
echo "[$(TS)] Manifest: ${MANIFEST}" | tee -a "${MANIFEST}"
