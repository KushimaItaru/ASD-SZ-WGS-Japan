#!/bin/bash
# ファイル名: 09a_mssng_replication_v19.sbatch
# 処理内容:
#   - Script 9A (MSSNG replication v19; v18 + patch v3 dual dump) を
#     ncbn-cpu パーティションで実行する sbatch ラッパー。
#   - v19 の追加処理: PHASE 4 末尾で constraint enrichment 用 dual dump
#     (mssng_event_bins_dumped_v1.tsv.gz, mssng_sample_covariates_dumped_v1.tsv.gz)
#     を /lustre12/home/kushima-pg/tad04292026/14_constraint_enrichment_v1/output_v1/ に出力。
#   - 既存の v18 ロジック (Pattern A + Pattern C, GEE Independence) は変更なし、
#     v18 の結果を完全再現する。
#   - 遺伝研スパコン用。実行時間を記録。
#   - 出力は tad04292026/09_mssng_sample_burden/ 直下に log / err / output_v19/ を生成。
#   - v18 の実行実績 (GEE Independence, n=6466) は数分〜数十分レンジのため
#     8 CPU / 32G / 4h を保守的上限として割当。

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_mssng_v19
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04292026/09_mssng_sample_burden/tad_mssng_v19_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04292026/09_mssng_sample_burden/tad_mssng_v19_%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04292026/09_mssng_sample_burden
python3 tad_replication_mssng_v19.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
