#!/bin/bash
# ============================================================================
# NOTE: This script contains hardcoded paths specific to the original analysis
# environment (NIG supercomputer). To run in a different environment, update
# the paths marked with "# CONFIGURE" below, or set the corresponding
# environment variables defined in config_v1.sh.
# ============================================================================
# strling/01_strling_extract_array_genomewide_v3.sh
# - (1) Run STRling extract on ehdn_all_samples.tsv (ASD/SZ/Healthy/family_member)
# - (2) Source config using absolute path (works even under Slurm spool)
# - (3) Output: strling_output_genomewide/bins/<SampleID>.bin
# - (4) Log: strling_output_genomewide/logs (Slurm out/err + per-sample extract log)
# - (5) Record execution time

#SBATCH -J strling_extract_gw
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --array=1-ARRAY_SIZE%50
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH -t 12:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/extract_%A_%a.out
#SBATCH --error=logs/extract_%A_%a.err

set -euo pipefail

START=$(date +%s)
TS(){ date '+%Y-%m-%d %H:%M:%S'; }

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
CFG="${PROJECT_ROOT}/helpers/strling/00_strling_genomewide_config_v1.sh"

if [ ! -f "${CFG}" ]; then
  echo "[$(TS)] [ERROR] Config not found: ${CFG}" >&2
  exit 2
fi
source "${CFG}"

# Force output directory to match this project
OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/strling_output_genomewide}"
BINS_DIR="${OUT_ROOT}/bins"
LOG_DIR="${OUT_ROOT}/logs"

# Samples same as EHdn (including family_member)
SAMPLES_TSV="${PROJECT_ROOT}/sample_lists/ehdn_all_samples.tsv"

# Fallback if not defined in config
STRLING_BIN="${STRLING_BIN:-strling}"
REFERENCE_FASTA="${REFERENCE_FASTA:-/lustre12/home/grifinpd-pg/resource_2020Aug/Homo_sapiens_assembly38.fasta}"  # CONFIGURE
PROPORTION_REPEAT="${PROPORTION_REPEAT:-0.8}"
MIN_MAPQ="${MIN_MAPQ:-40}"
STRLING_GENOME_INDEX="${STRLING_GENOME_INDEX:-${OUT_ROOT}/hg38.strling.bed}"

mkdir -p "${BINS_DIR}" "${LOG_DIR}"

echo "[$(TS)] [INFO] START extract"
echo "[$(TS)] [INFO] JOB=${SLURM_JOB_ID:-NA} TASK=${SLURM_ARRAY_TASK_ID:-NA} HOST=$(hostname)"
echo "[$(TS)] [INFO] SAMPLES_TSV=${SAMPLES_TSV}"
echo "[$(TS)] [INFO] OUT_ROOT=${OUT_ROOT}"
echo "[$(TS)] [INFO] BINS_DIR=${BINS_DIR}"

if [ ! -f "${SAMPLES_TSV}" ]; then
  echo "[$(TS)] [ERROR] Sample list not found: ${SAMPLES_TSV}" >&2
  exit 3
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:?}"
LINE=$((TASK_ID + 1))
REC=$(sed -n "${LINE}p" "${SAMPLES_TSV}" || true)
if [ -z "${REC}" ]; then
  echo "[$(TS)] [INFO] No record for task=${TASK_ID}. Exit."
  exit 0
fi

# ehdn_all_samples.tsv format: SampleID  Group  CRAM_Path
SAMPLE_ID=$(echo "${REC}" | cut -f1)
GROUP=$(echo "${REC}" | cut -f2)
CRAM_PATH=$(echo "${REC}" | cut -f3)

echo "[$(TS)] [INFO] SampleID=${SAMPLE_ID} Group=${GROUP}"
echo "[$(TS)] [INFO] CRAM=${CRAM_PATH}"

if [ ! -f "${CRAM_PATH}" ]; then
  echo "[$(TS)] [ERROR] CRAM not found: ${CRAM_PATH}" >&2
  exit 4
fi

# Verify CRAI
if [ ! -f "${CRAM_PATH}.crai" ]; then
  echo "[$(TS)] [INFO] Creating CRAM index..."
  samtools index "${CRAM_PATH}"
fi

OUT_BIN="${BINS_DIR}/${SAMPLE_ID}.bin"
SAMPLE_LOG="${LOG_DIR}/${SAMPLE_ID}_extract.log"

# Skip if bin already exists (>1000B as sanity check)
if [ -f "${OUT_BIN}" ]; then
  FILE_SIZE=$(stat -c%s "${OUT_BIN}" 2>/dev/null || wc -c < "${OUT_BIN}" 2>/dev/null || echo 0)
  if [ "${FILE_SIZE}" -gt 1000 ]; then
    echo "[$(TS)] [INFO] Already exists: ${OUT_BIN} (${FILE_SIZE} bytes) -> skip"
    END=$(date +%s)
    echo "[$(TS)] [DONE] Elapsed=$((END-START))s"
    exit 0
  fi
fi

echo "[$(TS)] [INFO] Running strling extract..."
EX_START=$(date +%s)

if [ -f "${STRLING_GENOME_INDEX}" ]; then
  "${STRLING_BIN}" extract \
    -f "${REFERENCE_FASTA}" \
    -g "${STRLING_GENOME_INDEX}" \
    -p "${PROPORTION_REPEAT}" \
    -q "${MIN_MAPQ}" \
    "${CRAM_PATH}" \
    "${OUT_BIN}" 2>&1 | tee "${SAMPLE_LOG}"
else
  "${STRLING_BIN}" extract \
    -f "${REFERENCE_FASTA}" \
    -p "${PROPORTION_REPEAT}" \
    -q "${MIN_MAPQ}" \
    "${CRAM_PATH}" \
    "${OUT_BIN}" 2>&1 | tee "${SAMPLE_LOG}"
fi

EX_END=$(date +%s)
echo "[$(TS)] [INFO] Extract time=${EX_END}-${EX_START} sec"

if [ ! -s "${OUT_BIN}" ]; then
  echo "[$(TS)] [ERROR] Output bin not created or empty: ${OUT_BIN}" >&2
  exit 5
fi

END=$(date +%s)
echo "[$(TS)] [DONE] Sample ${SAMPLE_ID}. Elapsed=$((END-START))s"
