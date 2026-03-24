#!/usr/bin/env bash
# 01_ehdn_setup_and_sample_lists.sh
# ファイル名: 01_ehdn_setup_and_sample_lists.sh
# 処理内容:
#   - EHdn パイプラインの初期セットアップ入口ラッパー
#   - Step 1: setup（ディレクトリ作成、リソース確認）
#   - Step 2: sample list 作成（SampleInfo → 用途別サンプルリスト分離）
#   - 各ステップは sbatch dependency で連鎖し、前段完了後に次段を開始
#   - 投入した Job ID を manifest ファイルに記録
#
# 使い方:
#   bash 01_ehdn_setup_and_sample_lists.sh
#
# 注意:
#   - 既存 helper script の処理内容は一切変更しない

set -euo pipefail

SCRIPT_START=$(date +%s)
TS(){ date "+%Y-%m-%d %H:%M:%S"; }

# ---- Paths (configurable via environment variables) ----
WRAPPER_ROOT="${WRAPPER_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
HELPER_ROOT="${HELPER_ROOT:-${WRAPPER_ROOT}/../str_12282025}"

LOG_DIR="${WRAPPER_ROOT}/ehdn/logs"
mkdir -p "${LOG_DIR}"

MANIFEST="${LOG_DIR}/01_setup_samples_manifest_$(date +%Y%m%d_%H%M%S).txt"
echo "[$(TS)] === 01_ehdn_setup_and_sample_lists.sh ===" | tee "${MANIFEST}"
echo "[$(TS)] WRAPPER_ROOT=${WRAPPER_ROOT}" | tee -a "${MANIFEST}"
echo "[$(TS)] HELPER_ROOT=${HELPER_ROOT}" | tee -a "${MANIFEST}"

# ---- Step 1: setup (single job) ----
JOB1=$(sbatch --parsable \
    --partition="${SLURM_PARTITION:-ncbn-cpu}" --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
    --cpus-per-task=1 --mem=4G --time=00:10:00 \
    --job-name=ehdn_setup \
    --output="${LOG_DIR}/setup_%j.out" \
    --error="${LOG_DIR}/setup_%j.err" \
    --wrap="bash ${HELPER_ROOT}/00_setup_project_v4.sh")
echo "[$(TS)] [Step1] setup submitted: JobID=${JOB1}" | tee -a "${MANIFEST}"

# ---- Step 2: sample list creation (single job, after setup) ----
JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} \
    --partition="${SLURM_PARTITION:-ncbn-cpu}" --account="${SLURM_ACCOUNT:-ncbn-cpu}" \
    --cpus-per-task=1 --mem=8G --time=00:10:00 \
    --job-name=ehdn_sample_lists \
    --output="${LOG_DIR}/sample_lists_%j.out" \
    --error="${LOG_DIR}/sample_lists_%j.err" \
    --wrap="bash -lc 'source ~/.bashrc && conda activate ngs && cd ${HELPER_ROOT}/ehdn && python3 01_prepare_sample_lists_v2.py'")
echo "[$(TS)] [Step2] sample lists submitted: JobID=${JOB2} (afterok:${JOB1})" | tee -a "${MANIFEST}"

# ---- Summary ----
SCRIPT_END=$(date +%s)
ELAPSED=$(( SCRIPT_END - SCRIPT_START ))
echo "" | tee -a "${MANIFEST}"
echo "[$(TS)] All steps submitted." | tee -a "${MANIFEST}"
echo "[$(TS)] Final job: ${JOB2}" | tee -a "${MANIFEST}"
echo "[$(TS)] Submission elapsed: ${ELAPSED}s" | tee -a "${MANIFEST}"
echo "[$(TS)] Manifest: ${MANIFEST}" | tee -a "${MANIFEST}"
