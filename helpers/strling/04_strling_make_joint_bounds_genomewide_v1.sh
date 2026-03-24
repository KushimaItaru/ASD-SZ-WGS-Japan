#!/usr/bin/env bash
# ============================================================================
# NOTE: This script contains hardcoded paths specific to the original analysis
# environment (NIG supercomputer). To run in a different environment, update
# the paths marked with "# CONFIGURE" below, or set the corresponding
# environment variables defined in config_v1.sh.
# ============================================================================
# 04_strling_make_joint_bounds_genomewide_v1.sh
# - 処理内容:
#   - 染色体ごとに作成済みの *-bounds.txt を結合して joint-bounds.txt を作成
#   - ヘッダを保持し、重複行は sort|uniq で除去
#   - 生成された joint-bounds.txt の行数（loci数）をログ出力
# - 実行時間を記録

set -euo pipefail
START=$(date +%s)
TS(){ date '+%Y-%m-%d %H:%M:%S'; }

PROJECT_ROOT="/lustre12/home/kushima-pg/str_12282025"  # CONFIGURE
CFG="${PROJECT_ROOT}/strling/00_strling_genomewide_config_v1.sh"
if [ ! -f "${CFG}" ]; then
  echo "[$(TS)] [ERROR] Config not found: ${CFG}" >&2
  exit 1
fi
# shellcheck source=/dev/null
source "${CFG}"

OUT_ROOT="${OUT_ROOT:-/lustre12/home/kushima-pg/str_12282025/strling_output_genomewide}"  # CONFIGURE
STR_RES_DIR="${STR_RES_DIR:-${OUT_ROOT}/str-results}"
LOG_DIR="${LOG_DIR:-${OUT_ROOT}/logs}"
mkdir -p "${STR_RES_DIR}" "${LOG_DIR}"

echo "[$(TS)] [INFO] OUT_ROOT=${OUT_ROOT}"
echo "[$(TS)] [INFO] STR_RES_DIR=${STR_RES_DIR}"

# chr別 bounds の想定ファイル名: ${STR_RES_DIR}/chr1-bounds.txt ... chrY-bounds.txt
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

# ヘッダは最初のファイルの1行目
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
  # 2行目以降（ボディ）を追記
  tail -n +2 "${f}" >> "${TMP_BODY}" || true
done

if [ ! -s "${TMP_BODY}" ]; then
  echo "[$(TS)] [ERROR] No body lines collected from chr bounds." >&2
  exit 4
fi

# ソート・重複排除（boundsは概ね: chrom left right repeatunit ...）
# 4列までで主キーとみなして uniq（完全一致でuniq）
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
