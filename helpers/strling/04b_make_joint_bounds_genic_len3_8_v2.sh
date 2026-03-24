#!/usr/bin/env bash
# ============================================================================
# NOTE: This script contains hardcoded paths specific to the original analysis
# environment (NIG supercomputer). To run in a different environment, update
# the paths marked with "# CONFIGURE" below, or set the corresponding
# environment variables defined in config_v1.sh.
# ============================================================================
# 04b_make_joint_bounds_genic_len3_8_v2.sh
# - Description:
#   - Extract loci overlapping gene_regions_1kb_pad.bed from genome-wide joint-bounds.txt (genic 1kb padded)
#   - Restrict repeat (4th column) length to 3-8 bp (matching EHdn min/max unit length)
#   - Stabilize output via sort/uniq (for reproducibility)
#   - Output: joint-bounds.genic_1kbpad.len3_8.txt
#   - Log loci counts (GW / GENIC_LEN3_8)
#   - Record execution time
#
# Usage:
#   Run via the top-level wrapper: strling/01_strling_build_panel.sh
#
# Overridable via environment variables:
#   JOINT=/path/to/joint-bounds.txt
#   GENE=/path/to/gene_regions_1kb_pad.bed
#   OUT=/path/to/joint-bounds.genic_1kbpad.len3_8.txt
#   LOG=/path/to/logfile.log

set -euo pipefail

START=$(date +%s)
TS(){ date '+%Y-%m-%d %H:%M:%S'; }

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/strling_output_genomewide}"
STR_RES="${OUT_ROOT}/str-results"
RES_DIR="${RES_DIR:-${PROJECT_ROOT}/resources}"
LOG_DIR="${OUT_ROOT}/logs"

mkdir -p "${STR_RES}" "${LOG_DIR}"

JOINT="${JOINT:-${STR_RES}/joint-bounds.txt}"
GENE="${GENE:-${RES_DIR}/gene_regions_1kb_pad.bed}"
OUT="${OUT:-${STR_RES}/joint-bounds.genic_1kbpad.len3_8.txt}"
LOG="${LOG:-${LOG_DIR}/make_joint_bounds_genic_len3_8_$(date +%Y%m%d_%H%M%S).log}"

{
  echo "[$(TS)] [INFO] Start 04b_make_joint_bounds_genic_len3_8_v2.sh"
  echo "[$(TS)] [INFO] PROJECT_ROOT=${PROJECT_ROOT}"
  echo "[$(TS)] [INFO] OUT_ROOT=${OUT_ROOT}"
  echo "[$(TS)] [INFO] STR_RES=${STR_RES}"
  echo "[$(TS)] [INFO] JOINT=${JOINT}"
  echo "[$(TS)] [INFO] GENE=${GENE}"
  echo "[$(TS)] [INFO] OUT=${OUT}"
  echo "[$(TS)] [INFO] LOG=${LOG}"
} | tee "${LOG}"

# checks
command -v bedtools >/dev/null 2>&1 || { echo "[$(TS)] [ERROR] bedtools not found in PATH" | tee -a "${LOG}"; exit 2; }
command -v awk >/dev/null 2>&1 || { echo "[$(TS)] [ERROR] awk not found" | tee -a "${LOG}"; exit 2; }
command -v sort >/dev/null 2>&1 || { echo "[$(TS)] [ERROR] sort not found" | tee -a "${LOG}"; exit 2; }

if [ ! -s "${JOINT}" ]; then
  echo "[$(TS)] [ERROR] joint-bounds not found/empty: ${JOINT}" | tee -a "${LOG}"
  exit 3
fi
if [ ! -s "${GENE}" ]; then
  echo "[$(TS)] [ERROR] gene bed not found/empty: ${GENE}" | tee -a "${LOG}"
  exit 3
fi

# Header sanity check (expected: #chrom left right repeat ...)
HDR=$(head -n 1 "${JOINT}" || true)
echo "[$(TS)] [INFO] JOINT header: ${HDR}" | tee -a "${LOG}"
if ! echo "${HDR}" | awk -F'\t' '($1 ~ /chrom/ || $1 ~ /#chrom/) && $2=="left" && $3=="right" {exit 0} {exit 1}'; then
  echo "[$(TS)] [WARN] Header is not the expected format (#chrom,left,right,...). Script assumes columns 1-4 are chrom/left/right/repeat." | tee -a "${LOG}"
fi

TMP="${OUT}.tmp.$$"
trap 'rm -f "${TMP}"' EXIT

echo "[$(TS)] [INFO] Building genic len3_8 joint-bounds (with sort/uniq)..." | tee -a "${LOG}"

{
  head -n 1 "${JOINT}"
  tail -n +2 "${JOINT}" \
    | bedtools intersect -a - -b "${GENE}" -u \
    | awk -F'\t' 'length($4)>=3 && length($4)<=8' \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -k4,4 \
    | uniq
} > "${TMP}"

GW_LOCI=$(tail -n +2 "${JOINT}" | wc -l | tr -d ' ')
OUT_LOCI=$(tail -n +2 "${TMP}" | wc -l | tr -d ' ')

if [ "${OUT_LOCI}" -le 0 ]; then
  echo "[$(TS)] [ERROR] OUT has 0 loci. Check overlap between JOINT and GENE." | tee -a "${LOG}"
  echo "[$(TS)] [HINT] Example:" | tee -a "${LOG}"
  echo "  tail -n +2 ${JOINT} | head -n 5 | bedtools intersect -a - -b ${GENE} -u" | tee -a "${LOG}"
  exit 4
fi

mv -f "${TMP}" "${OUT}"

{
  echo "[$(TS)] [DONE] Wrote: ${OUT}"
  echo "[$(TS)] [INFO] GW_LOCI=${GW_LOCI}"
  echo "[$(TS)] [INFO] GENIC_LEN3_8_LOCI=${OUT_LOCI}"
  echo "[$(TS)] [INFO] Head:"
  head -n 2 "${OUT}" || true
  END=$(date +%s)
  echo "[$(TS)] [DONE] Elapsed=$((END-START))s"
} | tee -a "${LOG}"
