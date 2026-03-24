#!/usr/bin/env bash
# 03_ehdn_profile_and_merge.sh
# ファイル名: 03_ehdn_profile_and_merge.sh
# 処理内容:
#   - EHdn プロファイリング + マージの入口ラッパー
#   - Step 1: EHdn profile（CRAMごとに EHdn 実行、アレイジョブ）
#   - Step 2: merge（全サンプルの locus.tsv をマージ＆正規化）
#   - 各ステップは sbatch dependency で連鎖
#   - 投入した Job ID を manifest ファイルに記録
#
# 前提:
#   - 02_ehdn_depth.sh が完了していること（depths_all.tsv が存在すること）
#
# 使い方:
#   bash 03_ehdn_profile_and_merge.sh
#
# 注意:
#   - 既存 helper script の処理内容は一切変更しない
#   - ARRAY_SIZE placeholder を sed で置換した temp file を submit する

set -euo pipefail

SCRIPT_START=$(date +%s)
TS(){ date "+%Y-%m-%d %H:%M:%S"; }

# ---- Paths ----
WRAPPER_ROOT="/home/kushima-pg/str_03242026"
HELPER_DIR="/home/kushima-pg/str_12282025/ehdn"
SAMPLE_LIST="/home/kushima-pg/str_12282025/sample_lists/ehdn_all_samples.tsv"

LOG_DIR="${WRAPPER_ROOT}/ehdn/logs"
mkdir -p "${LOG_DIR}"

MANIFEST="${LOG_DIR}/03_profile_merge_manifest_$(date +%Y%m%d_%H%M%S).txt"
echo "[$(TS)] === 03_ehdn_profile_and_merge.sh ===" | tee "${MANIFEST}"
echo "[$(TS)] WRAPPER_ROOT=${WRAPPER_ROOT}" | tee -a "${MANIFEST}"
echo "[$(TS)] HELPER_DIR=${HELPER_DIR}" | tee -a "${MANIFEST}"

# ---- Step 1: EHdn profile (array job) ----
N_SAMPLES=$(tail -n +2 "${SAMPLE_LIST}" | wc -l | tr -d ' ')
echo "[$(TS)] [Step1] EHdn profile: N_SAMPLES=${N_SAMPLES}" | tee -a "${MANIFEST}"

# 02_run_ehdn has ARRAY_SIZE placeholder → create temp with actual size
TMP_EHDN=$(mktemp "${LOG_DIR}/tmp_ehdn_XXXXXX.sh")
trap 'rm -f "${TMP_EHDN}"' EXIT
sed "s/ARRAY_SIZE/${N_SAMPLES}/g" "${HELPER_DIR}/02_run_ehdn_array_v2.sh" > "${TMP_EHDN}"

JOB1=$(sbatch --parsable "${TMP_EHDN}")
echo "[$(TS)] [Step1] EHdn profile submitted: JobID=${JOB1}" | tee -a "${MANIFEST}"

# ---- Step 2: merge (single job, after profile) ----
JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} \
    --partition=ncbn-cpu --account=ncbn-cpu \
    --cpus-per-task=4 --mem=64G --time=02:00:00 \
    --job-name=ehdn_merge \
    --output="${LOG_DIR}/merge_%j.out" \
    --error="${LOG_DIR}/merge_%j.err" \
    --wrap="bash -lc 'source ~/.bashrc && conda activate ngs && cd ${HELPER_DIR} && python3 04_merge_ehdn_novel_norm_v2.py'")
echo "[$(TS)] [Step2] merge submitted: JobID=${JOB2} (afterok:${JOB1})" | tee -a "${MANIFEST}"

# ---- Summary ----
SCRIPT_END=$(date +%s)
ELAPSED=$(( SCRIPT_END - SCRIPT_START ))
echo "" | tee -a "${MANIFEST}"
echo "[$(TS)] All steps submitted." | tee -a "${MANIFEST}"
echo "[$(TS)] Final job: ${JOB2}" | tee -a "${MANIFEST}"
echo "[$(TS)] Submission elapsed: ${ELAPSED}s" | tee -a "${MANIFEST}"
echo "[$(TS)] Manifest: ${MANIFEST}" | tee -a "${MANIFEST}"
