#!/usr/bin/env bash
# ============================================================================
# NOTE: This script contains hardcoded paths specific to the original analysis
# environment (NIG supercomputer). To run in a different environment, update
# the paths marked with "# CONFIGURE" below, or set the corresponding
# environment variables defined in config_v1.sh.
# ============================================================================
# 04_strling_make_joint_bounds_genomewide_v1.sh
# - Description:
#   - Concatenate per-chromosome *-bounds.txt into joint-bounds.txt
#   - Retain header and remove duplicate rows via sort|uniq
#   - Log line count (number of loci) of generated joint-bounds.txt
# - Record execution time

set -euo pipefail
START=$(date +%s)
TS(){ date '+%Y-%m-%d %H:%M:%S'; }

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
CFG="${PROJECT_ROOT}/helpers/strling/00_strling_genomewide_config_v1.sh"
if [ ! -f "${CFG}" ]; then
  echo "[$(TS)] [ERROR] Config not found: ${CFG}" >&2
  exit 1
fi
# shellcheck source=/dev/null
source "${CFG}"

OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/strling_output_genomewide}"
STR_RES_DIR="${STR_RES_DIR:-${OUT_ROOT}/str-results}"
LOG_DIR="${LOG_DIR:-${OUT_ROOT}/logs}"
mkdir -p "${STR_RES_DIR}" "${LOG_DIR}"

echo "[$(TS)] [INFO] OUT_ROOT=${OUT_ROOT}"
echo "[$(TS)] [INFO] STR_RES_DIR=${STR_RES_DIR}"

# Expected per-chromosome bounds filename pattern: ${STR_RES_DIR}/chr1-bounds.txt ... chrY-bounds.txt
shopt -s nullglob
BOUNDS_FILES=("${STR_RES_DIR}"/chr*-bounds.txt)

if [ "${#BOUNDS_FILES[@]}" -eq 0 ]; then
  echo "[$(TS)] [ERROR] No chromosome bounds found: ${STR_RES_DIR}/chr*-bounds.txt" >&2
  exit 2
fi

echo "[$(TS)] [INFO] Found bounds files: ${#BOUNDS_FILES[@]}"
echo "[$(TS)] [INFO] Example: ${BOUNDS_FILES[0]}"

JOINT_BOUNDS="${STR_RES_DIR}/joint-bounds.txt"
TMP_BODY="${STR_RES_DIR}/.tmp_joint_bounds_body.$$"
TMP_UNIQ="${STR_RES_DIR}/.tmp_joint_bounds_uniq.$$"

# Header is the first line of the first file
HEADER=$(head -n 1 "${BOUNDS_FILES[0]}" || true)
if [ -z "${HEADER}" ]; then
  echo "[$(TS)] [ERROR] Header empty: ${BOUNDS_FILES[0]}" >&2
  exit 3
fi

: > "${TMP_BODY}"
for f in "${BOUNDS_FILES[@]}"; do
  if [ ! -s "${f}" ]; then
    echo "[$(TS)] [WARN] Empty bounds file (skip): ${f}" >&2
    continue
  fi
  # Append body (line 2 onward)
  tail -n +2 "${f}" >> "${TMP_BODY}" || true
done

if [ ! -s "${TMP_BODY}" ]; then
  echo "[$(TS)] [ERROR] No body lines collected from chr bounds." >&2
  exit 4
fi

# Sort and deduplicate (bounds format: chrom left right repeatunit ...)
# Treat first 4 columns as primary key for uniq (exact match)
sort -k1,1 -k2,2n -k3,3n -k4,4 "${TMP_BODY}" | uniq > "${TMP_UNIQ}"

{
  echo "${HEADER}"
  cat "${TMP_UNIQ}"
} > "${JOINT_BOUNDS}"

rm -f "${TMP_BODY}" "${TMP_UNIQ}"

N_LOCI=$(tail -n +2 "${JOINT_BOUNDS}" | wc -l | tr -d ' ' || echo 0)
echo "[$(TS)] [DONE] Wrote joint bounds: ${JOINT_BOUNDS} (loci=${N_LOCI})"

END=$(date +%s)
echo "[$(TS)] [DONE] Elapsed=$((END-START))s"
