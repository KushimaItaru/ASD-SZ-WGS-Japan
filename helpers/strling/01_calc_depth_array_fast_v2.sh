#!/bin/bash
# ============================================================================
# NOTE: This script contains hardcoded paths specific to the original analysis
# environment (NIG supercomputer). To run in a different environment, update
# the paths marked with "# CONFIGURE" below, or set the corresponding
# environment variables defined in config_v1.sh.
# ============================================================================
#SBATCH -J calc_depth_str1228
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --array=1-ARRAY_SIZE%100
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH -t 1:00:00
#SBATCH --mem=1G
#SBATCH --output=logs/depth/depth_%A_%a.out
#SBATCH --error=logs/depth/depth_%A_%a.err

# strling/01_calc_depth_array_fast_v2.sh (FIXED)
# - Compute depth on sample_lists/ehdn_all_samples.tsv (including family_member)
# - Set PROJECT_ROOT as absolute path to ensure config_v1.sh is readable in Slurm environment
# - Output: depth/depth_parts/<SampleID>.tsv
# - Record execution time

set -euo pipefail

START=$(date +%s)
TS(){ date '+%Y-%m-%d %H:%M:%S'; }

# IMPORTANT: Fix project root as absolute path for reliable access under Slurm
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
source "${PROJECT_ROOT}/config_v1.sh"

mkdir -p "${LOG_DIR}/depth" "${DEPTH_PARTS_DIR}"

SAMPLE_LIST="${SAMPLE_LIST_DIR}/ehdn_all_samples.tsv"

LINE=$((SLURM_ARRAY_TASK_ID + 1))
REC=$(sed -n "${LINE}p" "${SAMPLE_LIST}" || true)
if [ -z "${REC}" ]; then
  exit 0
fi

SAMPLE_ID=$(echo "${REC}" | cut -f1)
CRAM_PATH=$(echo "${REC}" | cut -f3)
OUT="${DEPTH_PARTS_DIR}/${SAMPLE_ID}.tsv"

# Skip if already exists
if [ -s "${OUT}" ]; then
  exit 0
fi

if [ ! -f "${CRAM_PATH}" ]; then
  echo -e "${SAMPLE_ID}\t0.0" > "${OUT}"
  exit 0
fi

samtools idxstats "${CRAM_PATH}" | \
awk -v id="${SAMPLE_ID}" -v rl="${READ_LENGTH_BP}" -v gs="${GENOME_SIZE_BP}" '
  {sum += $3}
  END{
    depth = (sum * rl) / gs;
    printf "%s\t%.4f\n", id, depth
  }
' > "${OUT}"

END=$(date +%s)
echo "[$(TS)] [DONE] ${SAMPLE_ID} depth computed. Elapsed=$((END-START))s"
