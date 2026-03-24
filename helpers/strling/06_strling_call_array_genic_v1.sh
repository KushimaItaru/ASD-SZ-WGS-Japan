#!/bin/bash
# ============================================================================
# NOTE: This script contains hardcoded paths specific to the original analysis
# environment (NIG supercomputer). To run in a different environment, update
# the paths marked with "# CONFIGURE" below, or set the corresponding
# environment variables defined in config_v1.sh.
# ============================================================================
# 06_strling_call_array_genic_v1.sh
# - 処理内容:
#   - casecontrol_samples.tsv（Healthy/ASD/SZ）を対象に、genic+repeatunit(3–8bp)の joint-bounds を用いて strling call を実行
#   - 入力: CRAM + bins/<SampleID>.bin + joint-bounds.genic_1kbpad.len3_8.txt
#   - 出力: calls_genic/<SampleID>-genotype.txt, calls_genic/<SampleID>-unplaced.txt, calls_genic/<SampleID>-bounds.txt
#   - 既に genotype/unplaced があればスキップ（再開可能）
#   - /usr/bin/time -v を ${LOG_DIR}/${SampleID}.call_genic.time.txt に保存（MaxRSS 等を後で集計可能）
#   - 実行時間（秒）をログ出力

#SBATCH -J strling_call_genic
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --array=1-ARRAY_SIZE%50
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH -t 0:20:00
#SBATCH --mem=32G
#SBATCH --output=/lustre12/home/kushima-pg/str_12282025/strling_output_genomewide/logs/call_genic_%A_%a.out
#SBATCH --error=/lustre12/home/kushima-pg/str_12282025/strling_output_genomewide/logs/call_genic_%A_%a.err

set -euo pipefail

START=$(date +%s)
TS(){ date '+%Y-%m-%d %H:%M:%S'; }

PROJECT_ROOT="/lustre12/home/kushima-pg/str_12282025"  # CONFIGURE
CFG="${PROJECT_ROOT}/strling/00_strling_genomewide_config_v1.sh"

if [ ! -f "${CFG}" ]; then
  echo "[$(TS)] [ERROR] Config not found: ${CFG}" >&2
  exit 2
fi
# shellcheck source=/dev/null
source "${CFG}"

OUT_ROOT="${OUT_ROOT:-/lustre12/home/kushima-pg/str_12282025/strling_output_genomewide}"  # CONFIGURE
BINS_DIR="${BINS_DIR:-${OUT_ROOT}/bins}"
STR_RES_DIR="${STR_RES_DIR:-${OUT_ROOT}/str-results}"
CALL_DIR="${OUT_ROOT}/calls_genic"
LOG_DIR="${LOG_DIR:-${OUT_ROOT}/logs}"

# case-control のみに限定（EHdn再確認目的）
SAMPLES_TSV="${PROJECT_ROOT}/sample_lists/casecontrol_samples.tsv"

# configに無い場合のフォールバック
STRLING_BIN="${STRLING_BIN:-strling}"
REFERENCE_FASTA="${REFERENCE_FASTA:-/lustre12/home/grifinpd-pg/resource_2020Aug/Homo_sapiens_assembly38.fasta}"  # CONFIGURE

mkdir -p "${CALL_DIR}" "${LOG_DIR}"

# ★EHdnに合わせて repeat unit 3–8bp のみに絞った genic joint-bounds を使用
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

# 列名ベースで SampleID / Group / CRAM列を自動検出して該当行を抽出
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

# 既に genotype/unplaced があればスキップ（再開可能）
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

# CRAI確認
if [ ! -f "${CRAM_PATH}.crai" ]; then
  echo "[$(TS)] [INFO] Creating CRAM index..."
  samtools index "${CRAM_PATH}"
fi

if [ ! -s "${BIN}" ]; then
  echo "[$(TS)] [ERROR] BIN not found or empty: ${BIN}" >&2
  exit 6
fi

# 前回途中で 0 byte が残っている場合に備え、開始前に消す（再開時の混乱回避）
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

# 出力チェック
if [ ! -s "${OUT_GENOTYPE}" ]; then
  echo "[$(TS)] [ERROR] genotype not generated or empty: ${OUT_GENOTYPE}" >&2
  exit 7
fi
# unplaced は空のことがあり得るが、ファイル自体が無ければ警告
if [ ! -f "${OUT_UNPLACED}" ]; then
  echo "[$(TS)] [WARN] unplaced missing: ${OUT_UNPLACED}" >&2
fi
# bounds は通常出る。無ければ警告
if [ ! -f "${OUT_BOUNDS}" ]; then
  echo "[$(TS)] [WARN] bounds missing: ${OUT_BOUNDS}" >&2
fi
# time log は必須（mem最適化に使う）
if [ ! -s "${TIME_LOG}" ]; then
  echo "[$(TS)] [WARN] time log missing/empty: ${TIME_LOG}" >&2
fi

END=$(date +%s)
echo "[$(TS)] [DONE] Sample ${SAMPLE_ID}. Elapsed=$((END-START))s"

