#!/bin/bash
# ============================================================================
# NOTE: This script contains hardcoded paths specific to the original analysis
# environment (NIG supercomputer). To run in a different environment, update
# the paths marked with "# CONFIGURE" below, or set the corresponding
# environment variables defined in config_v1.sh.
# ============================================================================
# 03_strling_merge_by_chrom_genomewide_v3.sh
# - Use bins/*.bin from STRling extract to run strling merge per chromosome (SLURM array)
# - Source config using absolute path to avoid failure under Slurm spool
# - Handle CHROMOSOMES safely even if undefined/non-array (v3 fix)
# - Skip if output already exists and is non-empty (resumable)
# - Record execution time (/usr/bin/time -v logged)

#SBATCH -J strling_merge_gw
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --array=1-24
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --mem=128G
#SBATCH --output=logs/merge_chr%A_%a.out
#SBATCH --error=logs/merge_chr%A_%a.err

set -euo pipefail
shopt -s nullglob

START=$(date +%s)
TS(){ date '+%Y-%m-%d %H:%M:%S'; }

# ★Important: Source config using absolute path (reliable under Slurm spool)
CONFIG="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}/helpers/strling/00_strling_genomewide_config_v1.sh"
if [ ! -f "${CONFIG}" ]; then
  echo "[$(TS)] [ERROR] Config not found: ${CONFIG}" >&2
  exit 1
fi
# shellcheck source=/dev/null
source "${CONFIG}"

# Safety fallback in case config leaves it undefined
OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/strling_output_genomewide}"
BINS_DIR="${BINS_DIR:-${OUT_ROOT}/bins}"
STR_RES_DIR="${STR_RES_DIR:-${OUT_ROOT}/str-results}"
LOG_DIR="${LOG_DIR:-${OUT_ROOT}/logs}"

mkdir -p "${STR_RES_DIR}" "${LOG_DIR}"

# ---- Safely determine chromosome list (v3 fix) ----
# 1) Create default array if CHROMOSOMES is undefined
# 2) If defined as string (not array), split by whitespace into array
if declare -p CHROMOSOMES >/dev/null 2>&1; then
  # defined
  if ! declare -p CHROMOSOMES 2>/dev/null | grep -q 'declare \-a'; then
    # not an array -> convert string to array (split on whitespace)
    _tmp="${CHROMOSOMES}"
    unset CHROMOSOMES
    # shellcheck disable=SC2206
    CHROMOSOMES=(${_tmp})
  fi
else
  CHROMOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
               chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
               chr21 chr22 chrX chrY)
fi

# Use default if array is empty
if [ "${#CHROMOSOMES[@]}" -eq 0 ]; then
  CHROMOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
               chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
               chr21 chr22 chrX chrY)
fi
# --------------------------------------------

TASK_ID="${SLURM_ARRAY_TASK_ID:?}"
IDX=$((TASK_ID - 1))

if [ "${IDX}" -lt 0 ] || [ "${IDX}" -ge "${#CHROMOSOMES[@]}" ]; then
  echo "[$(TS)] [INFO] Task ${TASK_ID} out of range for CHROMOSOMES (n=${#CHROMOSOMES[@]}). Exit 0."
  exit 0
fi

CHR="${CHROMOSOMES[$IDX]}"

echo "[$(TS)] [INFO] Start STRling merge genomewide v3"
echo "[$(TS)] [INFO] JobID=${SLURM_JOB_ID:-NA} Task=${SLURM_ARRAY_TASK_ID:-NA} Host=$(hostname)"
echo "[$(TS)] [INFO] CHR=${CHR}"
echo "[$(TS)] [INFO] OUT_ROOT=${OUT_ROOT}"
echo "[$(TS)] [INFO] BINS_DIR=${BINS_DIR}"
echo "[$(TS)] [INFO] STR_RES_DIR=${STR_RES_DIR}"
echo "[$(TS)] [INFO] LOG_DIR=${LOG_DIR}"
echo "[$(TS)] [INFO] Reference=${REFERENCE_FASTA:-NA}"
echo "[$(TS)] [INFO] STRLING_BIN=${STRLING_BIN:-strling}"
echo "[$(TS)] [INFO] Params: WINDOW=${WINDOW:--1} MIN_SUPPORT=${MIN_SUPPORT:-5} MIN_CLIP_TOTAL=${MIN_CLIP_TOTAL:-0} MIN_MAPQ=${MIN_MAPQ:-40}"

# Verify input bins
if [ ! -d "${BINS_DIR}" ]; then
  echo "[$(TS)] [ERROR] BINS_DIR not found: ${BINS_DIR}" >&2
  exit 2
fi

BIN_FILES=("${BINS_DIR}"/*.bin)
N_BIN=${#BIN_FILES[@]}
if (( N_BIN == 0 )); then
  echo "[$(TS)] [ERROR] No .bin files found in ${BINS_DIR}" >&2
  exit 3
fi
echo "[$(TS)] [INFO] Found ${N_BIN} bin files."

# Output (per chromosome)
OUT_PREFIX="${STR_RES_DIR}/${CHR}"
OUT_BOUNDS="${OUT_PREFIX}-bounds.txt"
TIME_LOG="${LOG_DIR}/${CHR}.merge.time.txt"

# Skip if output exists (resumable)
if [ -s "${OUT_BOUNDS}" ]; then
  N_LOCI=$(tail -n +2 "${OUT_BOUNDS}" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
  if [ "${N_LOCI}" -gt 0 ]; then
    echo "[$(TS)] [INFO] Already done: ${OUT_BOUNDS} (loci=${N_LOCI}) -> skip"
    END=$(date +%s)
    echo "[$(TS)] [DONE] Elapsed=$((END-START))s"
    exit 0
  fi
fi

# Execute merge
echo "[$(TS)] [INFO] Running: strling merge --chromosome ${CHR}"
set +e
/usr/bin/time -v -o "${TIME_LOG}" \
"${STRLING_BIN:-strling}" merge \
  -f "${REFERENCE_FASTA}" \
  --output-prefix "${OUT_PREFIX}" \
  --chromosome "${CHR}" \
  -w "${WINDOW:--1}" \
  -m "${MIN_SUPPORT:-5}" \
  -t "${MIN_CLIP_TOTAL:-0}" \
  -q "${MIN_MAPQ:-40}" \
  "${BIN_FILES[@]}"
RET=$?
set -e

if (( RET != 0 )); then
  echo "[$(TS)] [ERROR] strling merge failed (CHR=${CHR}) exit=${RET}" >&2
  exit "${RET}"
fi

# Quick output check
if [ -s "${OUT_BOUNDS}" ]; then
  N_LOCI=$(tail -n +2 "${OUT_BOUNDS}" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
  echo "[$(TS)] [DONE] ${OUT_BOUNDS} (loci=${N_LOCI})"
else
  echo "[$(TS)] [WARN] Output bounds missing or empty: ${OUT_BOUNDS}"
fi

END=$(date +%s)
echo "[$(TS)] [DONE] Finished CHR=${CHR}. Elapsed=$((END-START))s"
