#!/bin/bash
# ============================================================================
# NOTE: This script contains hardcoded paths specific to the original analysis
# environment (NIG supercomputer). To run in a different environment, update
# the paths marked with "# CONFIGURE" below, or set the corresponding
# environment variables defined in config_v1.sh.
# ============================================================================
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=ehdn_str1228
#SBATCH --output=logs/ehdn/ehdn_%A_%a.out
#SBATCH --error=logs/ehdn/ehdn_%A_%a.err
#SBATCH --array=1-ARRAY_SIZE%40
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00

# ehdn/02_run_ehdn_array_v2.sh (FIXED)
# - sample_lists/ehdn_all_samples.tsv（ASD/SZ/Healthy/family_member）を対象に EHdn profile を実行
# - Slurm実行環境でも config_v1.sh を確実に読むため、PROJECT_ROOT を絶対パスで指定
# - 出力: ehdn_output/<SampleID>/<SampleID>.locus.tsv 等
# - 実行時間を記録

set -euo pipefail

START=$(date +%s)
TS(){ date '+%Y-%m-%d %H:%M:%S'; }

# ★重要：Slurm上で確実に見つかるようプロジェクトルートを絶対パスで固定
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
source "${PROJECT_ROOT}/config_v1.sh"

# Slurmがログを開けずに落ちるのを避ける（念のため）
mkdir -p "${LOG_DIR}/ehdn" "${EHDN_OUT_DIR}"

echo "[$(TS)] [INFO] EHdn array start: job=${SLURM_JOB_ID:-NA} task=${SLURM_ARRAY_TASK_ID:-NA}"
echo "[$(TS)] [INFO] PROJECT_ROOT=${PROJECT_ROOT}"

SAMPLE_LIST="${SAMPLE_LIST_DIR}/ehdn_all_samples.tsv"

LINE=$((SLURM_ARRAY_TASK_ID + 1))
REC=$(sed -n "${LINE}p" "${SAMPLE_LIST}" || true)
if [ -z "${REC}" ]; then
  echo "[$(TS)] [INFO] No record for task=${SLURM_ARRAY_TASK_ID}. Exit."
  exit 0
fi

SAMPLE_ID=$(echo "${REC}" | cut -f1)
GROUP=$(echo "${REC}" | cut -f2)
CRAM_PATH=$(echo "${REC}" | cut -f3)

echo "[$(TS)] [INFO] SampleID=${SAMPLE_ID} Group=${GROUP}"
echo "[$(TS)] [INFO] CRAM=${CRAM_PATH}"

if [ ! -f "${CRAM_PATH}" ]; then
  echo "[$(TS)] [ERROR] CRAM not found: ${CRAM_PATH}" >&2
  exit 1
fi

if [ ! -f "${CRAM_PATH}.crai" ]; then
  echo "[$(TS)] [INFO] Indexing CRAM..."
  samtools index "${CRAM_PATH}"
fi

OUTDIR="${EHDN_OUT_DIR}/${SAMPLE_ID}"
mkdir -p "${OUTDIR}"
OUTPREFIX="${OUTDIR}/${SAMPLE_ID}"

# 既に locus.tsv があればスキップ（再実行時の高速化）
if [ -s "${OUTPREFIX}.locus.tsv" ]; then
  echo "[$(TS)] [INFO] Already done: ${OUTPREFIX}.locus.tsv -> skip"
  END=$(date +%s)
  echo "[$(TS)] [DONE] Elapsed=$((END-START))s"
  exit 0
fi

echo "[$(TS)] [INFO] Running EHdn profile..."
"${EHDN_BIN}" profile \
  --reads "${CRAM_PATH}" \
  --reference "${REFERENCE_FASTA}" \
  --output-prefix "${OUTPREFIX}" \
  --min-anchor-mapq "${EHDN_MIN_ANCHOR_MAPQ}" \
  --max-irr-mapq "${EHDN_MAX_IRR_MAPQ}" \
  --min-unit-len "${EHDN_MIN_UNIT_LEN}" \
  --max-unit-len "${EHDN_MAX_UNIT_LEN}"

if [ ! -f "${OUTPREFIX}.locus.tsv" ]; then
  echo "[$(TS)] [ERROR] locus.tsv not generated: ${OUTPREFIX}.locus.tsv" >&2
  exit 2
fi

END=$(date +%s)
echo "[$(TS)] [DONE] EHdn completed for ${SAMPLE_ID}. Elapsed=$((END-START))s"
