#!/usr/bin/env bash
# 01_strling_build_panel.sh
# ファイル名: 01_strling_build_panel.sh
# 処理内容:
#   - STRling パネル構築パイプラインの入口ラッパー
#   - Step 1: extract（CRAMから.bin抽出、アレイジョブ）
#   - Step 2: merge（染色体ごとにmerge、アレイジョブ 1-24）
#   - Step 3: joint-bounds（全染色体のboundsを結合）
#   - Step 4: genic filter（genic 3-8bp lociに絞り込み）
#   - 各ステップは sbatch dependency で連鎖し、前段完了後に次段を開始
#   - 投入した Job ID を manifest ファイルに記録
#
# 使い方:
#   bash 01_strling_build_panel.sh
#
# 注意:
#   - 既存 helper script の処理内容は一切変更しない
#   - extract/call のアレイサイズはサンプルリストから動的に決定

set -euo pipefail

SCRIPT_START=$(date +%s)
TS(){ date "+%Y-%m-%d %H:%M:%S"; }

# ---- Paths (configurable via environment variables) ----
WRAPPER_ROOT="${WRAPPER_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
HELPER_DIR="${HELPER_DIR_STRLING:-${WRAPPER_ROOT}/helpers/strling}"
CONFIG="${HELPER_DIR}/00_strling_genomewide_config_v1.sh"

# source config to get STRLING_SAMPLES_TSV, STRLING_OUT_ROOT, etc.
source "${CONFIG}"

LOG_DIR="${WRAPPER_ROOT}/strling/logs"
mkdir -p "${LOG_DIR}"

MANIFEST="${LOG_DIR}/01_build_panel_manifest_$(date +%Y%m%d_%H%M%S).txt"
echo "[$(TS)] === 01_strling_build_panel.sh ===" | tee "${MANIFEST}"
echo "[$(TS)] WRAPPER_ROOT=${WRAPPER_ROOT}" | tee -a "${MANIFEST}"
echo "[$(TS)] HELPER_DIR=${HELPER_DIR}" | tee -a "${MANIFEST}"

# ---- Step 1: extract (array job) ----
N_SAMPLES=$(tail -n +2 "${STRLING_SAMPLES_TSV}" | wc -l)
echo "[$(TS)] [Step1] extract: N_SAMPLES=${N_SAMPLES}" | tee -a "${MANIFEST}"

# 01_extract has ARRAY_SIZE placeholder → create temp with actual size
TMP_EXTRACT=""
trap '[[ -n "${TMP_EXTRACT:-}" ]] && rm -f "${TMP_EXTRACT}"' EXIT
TMP_EXTRACT=$(mktemp "${LOG_DIR}/tmp_extract_XXXXXX.sh")
sed "s/ARRAY_SIZE/${N_SAMPLES}/g" "${HELPER_DIR}/01_strling_extract_array_genomewide_v3.sh" > "${TMP_EXTRACT}"

JOB1=$(sbatch --parsable \
  --partition="${SLURM_PARTITION:-ncbn-cpu}" \
  --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
  --output="${LOG_DIR}/extract_%A_%a.out" \
  --error="${LOG_DIR}/extract_%A_%a.err" \
  "${TMP_EXTRACT}")
echo "[$(TS)] [Step1] extract submitted: JobID=${JOB1}" | tee -a "${MANIFEST}"

# ---- Step 2: merge (array 1-24, after extract) ----
JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} \
    --partition="${SLURM_PARTITION:-ncbn-cpu}" --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
    --output="${LOG_DIR}/merge_%A_%a.out" \
    --error="${LOG_DIR}/merge_%A_%a.err" \
    "${HELPER_DIR}/03_strling_merge_by_chrom_genomewide_v3.sh")
echo "[$(TS)] [Step2] merge submitted: JobID=${JOB2} (afterok:${JOB1})" | tee -a "${MANIFEST}"

# ---- Step 3: joint-bounds (single job, after merge) ----
JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} \
    --partition="${SLURM_PARTITION:-ncbn-cpu}" --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
    --cpus-per-task=1 --mem=64G --time=01:00:00 \
    --job-name=strling_joint_bounds \
    --output="${LOG_DIR}/joint_bounds_%j.out" \
    --error="${LOG_DIR}/joint_bounds_%j.err" \
    --wrap="bash ${HELPER_DIR}/04_strling_make_joint_bounds_genomewide_v1.sh")
echo "[$(TS)] [Step3] joint-bounds submitted: JobID=${JOB3} (afterok:${JOB2})" | tee -a "${MANIFEST}"

# ---- Step 4: genic filter (single job, after joint-bounds) ----
JOB4=$(sbatch --parsable --dependency=afterok:${JOB3} \
    --partition="${SLURM_PARTITION:-ncbn-cpu}" --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
    --cpus-per-task=1 --mem=64G --time=00:30:00 \
    --job-name=strling_genic_filter \
    --output="${LOG_DIR}/genic_filter_%j.out" \
    --error="${LOG_DIR}/genic_filter_%j.err" \
    --wrap="bash ${HELPER_DIR}/04b_make_joint_bounds_genic_len3_8_v2.sh")
echo "[$(TS)] [Step4] genic filter submitted: JobID=${JOB4} (afterok:${JOB3})" | tee -a "${MANIFEST}"

# ---- Summary ----
SCRIPT_END=$(date +%s)
ELAPSED=$(( SCRIPT_END - SCRIPT_START ))
echo "" | tee -a "${MANIFEST}"
echo "[$(TS)] All steps submitted." | tee -a "${MANIFEST}"
echo "[$(TS)] Final panel job: ${JOB4}" | tee -a "${MANIFEST}"
echo "[$(TS)] Submission elapsed: ${ELAPSED}s" | tee -a "${MANIFEST}"
echo "[$(TS)] Manifest: ${MANIFEST}" | tee -a "${MANIFEST}"
