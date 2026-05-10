#!/bin/bash
# ファイル名: 13_mssng_exon_exclusion_v2.sbatch
# 処理内容:
#   - Script 13 v2 (42_compute_mssng_exon_exclusion_v2.py) を ncbn-cpu で実行
#     する sbatch ラッパー (遺伝研スパコン用)。
#   - 出力: tad04212026/13_exon_exclusion_mssng/output_v2/
#       mssng_exon_exclusion_v2.tsv         (全 7 行: A1 3 + A2 2 + A3 1 + A3s 1)
#       mssng_exon_exclusion_summary_v2.tsv (ASD vs sibling × DEL 4 discovery-positive)
#
#   - 主処理:
#       1. GENCODE v46 protein-coding exon を merge
#       2. tad04212026 v2 bin_l2_annotation と exon を交差 (any_l2_bins のみ解析対象)
#       3. MSSNG heffel factorial overlap v1 の sample × event × bin overlap を読み、
#          event_id (chrN:start-end) から SV 座標を parse し exonfreeCNV を判定
#       4. 5 stratification の sample 単位 bin 数を集計
#       5. covariates v1 に 5 exposure 列を追加 → GEE Binomial Indep. (FAMILYID cluster)
#          + Firth sensitivity (QR rank reduction) を A1/A2/A3/A3-strict で実施
#       6. Multiple testing: 4 discovery-positive 行に BH + Holm 適用
#
#   - v1 (heffel mssng_exonfree_sensitivity_v1.py) -> v2 変更点:
#       1. PATHS: bin annotation は heffel v16 boundary master →
#          tad04212026 v2 bin_l2_annotation (event overlap と covariates は heffel 共有)
#       2. L2 列名: "is_hpc_exc_*" → "membership_HPC_*" 系
#       3. [追加] 自前 Firth penalized logistic sensitivity (QR rank reduction 付き)
#          MSSNG の platform one-hot 5 列 + intercept は常に collinear のため
#          QR で自動的に 1 列 drop される (Script 11 v4 と同一症状)
#       4. [追加] BH-FDR + Holm 多重検定補正 (GEE P, Firth P それぞれ)
#       5. [追加] carrier_case / carrier_ctrl / fit_status / fit_status_firth 列
#       6. log 形式を Script 11 v4 スタイル [YYYY-MM-DD HH:MM:SS] に統一
#       7. 出力ファイル名: mssng_exonfree_sensitivity_v1_results.tsv
#          -> mssng_exon_exclusion_v2.tsv + mssng_exon_exclusion_summary_v2.tsv
#       8. 主要解析ロジック (A1/A2/A3/A3-strict, GEE Binomial Indep + FAMILYID)
#          は heffel v1 と同一 (数値一致を担保)
#
#   - 依存入力:
#       tad04212026/02_bin_l2_annotation/output_v2/bin_l2_annotation_v2.tsv.gz
#       noncoding_tad_mssng_03132026/output_mssng_factorial_overlap_v1/
#           mssng_sample_event_bin_overlap_v1.tsv.gz
#           mssng_sample_covariates_v1.tsv
#       /lustre12/home/kushima-pg/annotationInfo/gencode.v46.annotation.gtf.gz
#
#   - fit モデル総数: 7 bin_type × 2 methods (GEE + Firth) = 14 fit
#     (A1: 3, A2 joint: 2 in one multi-x fit, A3: 1, A3-strict: 1)
#
#   - GENCODE exon parse + bin 交差 + sample aggregation + 14 fit のため
#     4 CPU / 32G / 2h を割当 (想定所要時間 5-20 分)

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_mssng_exonexcl_v2
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/13_exon_exclusion_mssng/tad_mssng_exonexcl_v2_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/13_exon_exclusion_mssng/tad_mssng_exonexcl_v2_%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04212026/13_exon_exclusion_mssng
python3 42_compute_mssng_exon_exclusion_v2.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
