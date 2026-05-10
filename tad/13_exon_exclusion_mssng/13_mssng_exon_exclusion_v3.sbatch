#!/bin/bash
# ファイル名: 13_mssng_exon_exclusion_v3.sbatch
# 処理内容:
#   - Script 13 v3 (42_compute_mssng_exon_exclusion_v3.py) を ncbn-cpu で実行
#     する sbatch ラッパー (遺伝研スパコン用)。
#   - 出力: tad04212026/13_exon_exclusion_mssng/output_v3/
#       mssng_exon_exclusion_v3.tsv         (全 7 行: A1 3 + A2 2 + A3 1 + A3s 1)
#       mssng_exon_exclusion_summary_v3.tsv (ASD vs sibling × DEL 4 discovery-positive)
#
#   - 主処理 (v2 と共通):
#       1. GENCODE v46 protein-coding exon を merge
#       2. tad04212026 v2 bin_l2_annotation と exon を交差 (any_l2_bins のみ解析対象)
#       3. MSSNG heffel factorial overlap v1 の sample × event × bin overlap を読み、
#          event_id (chrN:start-end) から SV 座標を parse し exonfreeCNV を判定
#          (v3: parse 失敗 event_id を集計し log 表示)
#       4. 5 stratification の sample 単位 bin 数を集計
#       5. covariates v1 に 5 exposure 列を追加 → GEE Binomial Indep. (FAMILYID cluster)
#          + Firth sensitivity (QR rank reduction) を A1/A2/A3/A3-strict で実施
#          (v3: MIN_CELL_COUNT=5 pre-check を GEE/Firth 両方に追加)
#       6. Multiple testing: 4 discovery-positive 行に BH + Holm 適用
#
#   - v2 (42_compute_mssng_exon_exclusion_v2.py) -> v3 変更点:
#       1. [追加] MIN_CELL_COUNT = 5 モジュール定数 (Script 12 v2 と同値)
#       2. [追加] _empty_fit_result_gee() / _empty_fit_result_firth() ヘルパー
#       3. [追加] fit_gee_mssng() に 3 段 pre-check (outcome variance /
#          per-exposure variance / total_carriers < MIN_CELL_COUNT)
#          - A2 joint は x_cols 合算の total_carriers で判定 (Script 12 v2 仕様)
#       4. [追加] fit_firth_mssng() にも同じ 3 段 pre-check
#       5. [追加] parse_mssng_event_id() の失敗件数を count_disrupted_bins() で
#          集計し log に件数 + 代表例 5 件を表示 (サイレント failure 防止)
#       6. 出力ディレクトリ: output_v2 -> output_v3
#          出力ファイル名: mssng_exon_exclusion_v2.tsv -> _v3.tsv, summary_v3.tsv
#       7. 主要解析ロジック (A1/A2/A3/A3-strict, GEE Indep + FAMILYID cluster,
#          Firth + QR rank reduction) は v2 と同一 (MIN_CELL_COUNT に抵触しない
#          7 行は v2 と bit-identical な数値が出る)
#
#   - 依存入力 (v2 と共通):
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
#     4 CPU / 32G / 2h を割当 (想定所要時間 5-20 分; v2 では 13.4 秒)

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_mssng_exonexcl_v3
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/13_exon_exclusion_mssng/tad_mssng_exonexcl_v3_%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/13_exon_exclusion_mssng/tad_mssng_exonexcl_v3_%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

cd /lustre12/home/kushima-pg/tad04212026/13_exon_exclusion_mssng
python3 42_compute_mssng_exon_exclusion_v3.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
