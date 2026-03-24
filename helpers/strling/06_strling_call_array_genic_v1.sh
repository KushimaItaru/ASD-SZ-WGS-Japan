#!/bin/bash
# ============================================================================
# NOTE: This script contains hardcoded paths specific to the original analysis
# environment (NIG supercomputer). To run in a different environment, update
# the paths marked with "# CONFIGURE" below, or set the corresponding
# environment variables defined in config_v1.sh.
# ============================================================================
# 06_strling_call_array_genic_v1.sh
# - Description:
#   - Run STRling call on casecontrol_samples.tsv (Healthy/ASD/SZ) using genic+repeatunit(3-8bp) joint-bounds
#   - Input: CRAM + bins/<SampleID>.bin + joint-bounds.genic_1kbpad.len3_8.txt
#   - Output: calls_genic/<SampleID>-genotype.txt, calls_genic/<SampleID>-unplaced.txt, calls_genic/<SampleID>-bounds.txt
#   - Skip if genotype/unplaced already exist (resumable)
#   - Save /usr/bin/time -v to ${LOG_DIR}/${SampleID}.call_genic.time.txt (for later MaxRSS analysis)
#   - Log execution time (seconds)

#SBATCH -J strling_call_genic
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --array=1-ARRAY_SIZE%50
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH -t 0:20:00
#SBATCH --mem=32G
#SBATCH --output=logs/call_genic_%A_%a.out
#SBATCH --error=logs/call_genic_%A_%a.err

set -euo pipefail

START=$(date +%s)
TS(){ date '+%Y-%m-%d %H:%M:%S'; }

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
CFG="${PROJECT_ROOT}/helpers/strling/00_strling_genomewide_config_v1.sh"

if [ ! -f "${CFG}" ]; then
  echo "[$(TS)] [ERROR] Config not found: ${CFG}" >&2
  exit 2
fi
# shellcheck source=/dev/null
source "${CFG}"

OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/strling_output_genomewide}"
BINS_DIR="${BINS_DIR:-${OUT_ROOT}/bins}"
STR_RES_DIR="${STR_RES_DIR:-${OUT_ROOT}/str-results}"
CALL_DIR="${OUT_ROOT}/calls_genic"
LOG_DIR="${LOG_DIR:-${OUT_ROOT}/logs}"

# Restrict to case-control only (EHdn verification)
SAMPLES_TSV="${PROJECT_ROOT}/sample_lists/casecontrol_samples.tsv"

# Fallback if not defined in config
STRLING_BIN="${STRLING_BIN:-strling}"
REFERENCE_FASTA="${REFERENCE_FASTA:-/lustre12/home/grifinpd-pg/resource_2020Aug/Homo_sapiens_assembly38.fasta}"  # CONFIGURE

mkdir -p "${CALL_DIR}" "${LOG_DIR}"

# ★Use genic joint-bounds restricted to repeat unit 3-8 bp (matching EHdn)
JOINT_BOUNDS="${STR_RES_DIR}/joint-bounds.genic_1kbpad.len3_8.txt"
if [ ! -s "${JOINT_BOUNDS}" ]; then
  echo "[$(TS)] [ERROR] genic len3_8 joint-bounds not found/empty: ${JOINT_BOUNDS}" >&2
  echo "[$(TS)] [HINT] Create it from joint-bounds.txt with bedtools intersect + awk length(3..8)." >&2
  exit 3
fi

if [ ! -f "${SAMPLES_TSV}" ]; then
  echo "[$(TS)] [ERROR] Sample list not found: ${SAMPLES_TSV}" >&2
  exit 4
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:?}"
LINE=$((TASK_ID + 1))

# Auto-detect SampleID / Group / CRAM columns by header name and extract target row
REC=$(awk -F'\t' -v line="${LINE}" '
  NR==1{
    for(i=1;i<=NF;i++){
      c=tolower($i)
      if($i=="SampleID" || c=="sampleid" || c=="sample") sid=i
      if($i=="Group" || c=="group" || c=="diagnosis") grp=i
      if(c ~ /cram/) cram=i
    }
    next
  }
  NR==line{
    if(!sid || !grp || !cram){exit 10}
    print $sid"\t"$grp"\t"$cram
  }
' "${SAMPLES_TSV}" || true)

if [ -z "${REC}" ]; then
  echo "[$(TS)] [INFO] No record for task=${TASK_ID}. Exit."
  exit 0
fi

SAMPLE_ID=$(echo "${REC}" | cut -f1)
GROUP=$(echo "${REC}" | cut -f2)
CRAM_PATH=$(echo "${REC}" | cut -f3)

BIN="${BINS_DIR}/${SAMPLE_ID}.bin"
OUT_PREFIX="${CALL_DIR}/${SAMPLE_ID}"
OUT_GENOTYPE="${OUT_PREFIX}-genotype.txt"
OUT_UNPLACED="${OUT_PREFIX}-unplaced.txt"
OUT_BOUNDS="${OUT_PREFIX}-bounds.txt"
TIME_LOG="${LOG_DIR}/${SAMPLE_ID}.call_genic.time.txt"
SAMPLE_LOG="${LOG_DIR}/${SAMPLE_ID}.call_genic.log"

echo "[$(TS)] [INFO] START call (genic len3_8)"
echo "[$(TS)] [INFO] JOB=${SLURM_JOB_ID:-NA} TASK=${SLURM_ARRAY_TASK_ID:-NA} HOST=$(hostname)"
echo "[$(TS)] [INFO] SampleID=${SAMPLE_ID} Group=${GROUP}"
echo "[$(TS)] [INFO] CRAM=${CRAM_PATH}"
echo "[$(TS)] [INFO] BIN=${BIN}"
echo "[$(TS)] [INFO] JOINT_BOUNDS=${JOINT_BOUNDS}"
echo "[$(TS)] [INFO] OUT_PREFIX=${OUT_PREFIX}"

# Skip if genotype/unplaced already exist (resumable)
if [ -s "${OUT_GENOTYPE}" ] && [ -s "${OUT_UNPLACED}" ]; then
  echo "[$(TS)] [INFO] Already exists -> skip: ${OUT_GENOTYPE}"
  END=$(date +%s)
  echo "[$(TS)] [DONE] Elapsed=$((END-START))s"
  exit 0
fi

if [ ! -f "${CRAM_PATH}" ]; then
  echo "[$(TS)] [ERROR] CRAM not found: ${CRAM_PATH}" >&2
  exit 5
fi

# Verify CRAI
if [ ! -f "${CRAM_PATH}.crai" ]; then
  echo "[$(TS)] [INFO] Creating CRAM index..."
  samtools index "${CRAM_PATH}"
fi

if [ ! -s "${BIN}" ]; then
  echo "[$(TS)] [ERROR] BIN not found or empty: ${BIN}" >&2
  exit 6
fi

# Remove leftover 0-byte files from previous runs (avoid confusion on restart)
rm -f "${OUT_GENOTYPE}" "${OUT_UNPLACED}" "${OUT_BOUNDS}" "${TIME_LOG}"

echo "[$(TS)] [INFO] Running: strling call (genic len3_8 bounds)"
set +e
/usr/bin/time -v -o "${TIME_LOG}" \
"${STRLING_BIN}" call \
  --output-prefix "${OUT_PREFIX}" \
  -b "${JOINT_BOUNDS}" \
  -f "${REFERENCE_FASTA}" \
  "${CRAM_PATH}" \
  "${BIN}" 2>&1 | tee "${SAMPLE_LOG}"
RET=$?
set -e

if (( RET != 0 )); then
  echo "[$(TS)] [ERROR] strling call failed: exit=${RET} sample=${SAMPLE_ID}" >&2
  exit "${RET}"
fi

# Output check
if [ ! -s "${OUT_GENOTYPE}" ]; then
  echo "[$(TS)] [ERROR] genotype not generated or empty: ${OUT_GENOTYPE}" >&2
  exit 7
fi
# unplaced may be empty, but warn if file is missing entirely
if [ ! -f "${OUT_UNPLACED}" ]; then
  echo "[$(TS)] [WARN] unplaced missing: ${OUT_UNPLACED}" >&2
fi
# bounds normally produced; warn if missing
if [ ! -f "${OUT_BOUNDS}" ]; then
  echo "[$(TS)] [WARN] bounds missing: ${OUT_BOUNDS}" >&2
fi
# time log is required (used for memory optimization)
if [ ! -s "${TIME_LOG}" ]; then
  echo "[$(TS)] [WARN] time log missing/empty: ${TIME_LOG}" >&2
fi

END=$(date +%s)
echo "[$(TS)] [DONE] Sample ${SAMPLE_ID}. Elapsed=$((END-START))s"

