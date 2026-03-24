#!/bin/bash
# ============================================================================
# NOTE: This script contains hardcoded paths specific to the original analysis
# environment (NIG supercomputer). To run in a different environment, update
# the paths marked with "# CONFIGURE" below, or set the corresponding
# environment variables defined in config_v1.sh.
# ============================================================================
# 03_strling_merge_by_chrom_genomewide_v3.sh
# - STRling extractで作成した bins/*.bin を使い、染色体ごとに strling merge を実行（SLURM array）
# - Slurm spool 実行でも落ちないよう、config を絶対パスで source
# - CHROMOSOMES が未定義/非配列でも落ちないように安全に扱う（v3修正点）
# - 既に出力が存在し中身がある場合はスキップ（再開可能）
# - 実行時間を記録（/usr/bin/time -v をログへ）

#SBATCH -J strling_merge_gw
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --array=1-24
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --mem=128G
#SBATCH --output=/lustre12/home/kushima-pg/str_12282025/strling_output_genomewide/logs/merge_chr%A_%a.out
#SBATCH --error=/lustre12/home/kushima-pg/str_12282025/strling_output_genomewide/logs/merge_chr%A_%a.err

set -euo pipefail
shopt -s nullglob

START=$(date +%s)
TS(){ date '+%Y-%m-%d %H:%M:%S'; }

# ★重要: config を絶対パスで source（Slurm spool でも確実に見つかる）
CONFIG="/lustre12/home/kushima-pg/str_12282025/strling/00_strling_genomewide_config_v1.sh"  # CONFIGURE
if [ ! -f "${CONFIG}" ]; then
  echo "[$(TS)] [ERROR] Config not found: ${CONFIG}" >&2
  exit 1
fi
# shellcheck source=/dev/null
source "${CONFIG}"

# config 側で未定義でも動くように保険
OUT_ROOT="${OUT_ROOT:-/lustre12/home/kushima-pg/str_12282025/strling_output_genomewide}"  # CONFIGURE
BINS_DIR="${BINS_DIR:-${OUT_ROOT}/bins}"
STR_RES_DIR="${STR_RES_DIR:-${OUT_ROOT}/str-results}"
LOG_DIR="${LOG_DIR:-${OUT_ROOT}/logs}"

mkdir -p "${STR_RES_DIR}" "${LOG_DIR}"

# ---- 染色体リストを安全に確定（v3修正点） ----
# 1) CHROMOSOMES が未定義ならデフォルト配列を作る
# 2) 定義されているが配列でない（文字列）なら空白区切りで配列化する
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

# 空配列ならデフォルト
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

# 入力 bin 確認
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

# 出力（chrごと）
OUT_PREFIX="${STR_RES_DIR}/${CHR}"
OUT_BOUNDS="${OUT_PREFIX}-bounds.txt"
TIME_LOG="${LOG_DIR}/${CHR}.merge.time.txt"

# 既に出力があればスキップ（再開可能）
if [ -s "${OUT_BOUNDS}" ]; then
  N_LOCI=$(tail -n +2 "${OUT_BOUNDS}" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
  if [ "${N_LOCI}" -gt 0 ]; then
    echo "[$(TS)] [INFO] Already done: ${OUT_BOUNDS} (loci=${N_LOCI}) -> skip"
    END=$(date +%s)
    echo "[$(TS)] [DONE] Elapsed=$((END-START))s"
    exit 0
  fi
fi

# merge 実行
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

# 出力簡易チェック
if [ -s "${OUT_BOUNDS}" ]; then
  N_LOCI=$(tail -n +2 "${OUT_BOUNDS}" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
  echo "[$(TS)] [DONE] ${OUT_BOUNDS} (loci=${N_LOCI})"
else
  echo "[$(TS)] [WARN] Output bounds missing or empty: ${OUT_BOUNDS}"
fi

END=$(date +%s)
echo "[$(TS)] [DONE] Finished CHR=${CHR}. Elapsed=$((END-START))s"
