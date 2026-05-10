#!/bin/bash
# ファイル名: tad_dynamic_boundary_specificity_test_v7.sbatch
# 処理内容:
#   - Step 07 v7 (tad_dynamic_boundary_specificity_test_v7.py) を ncbn-cpu パーティション
#     で実行する sbatch ラッパー (遺伝研スパコン用)。
#   - 出力: tad04212026/07_wgs_matched_static/output_tad_dynamic_boundary_specificity_v7_L${MATCH_LEVEL}/
#       specificity_test_results_main_v7_L${MATCH_LEVEL}.tsv
#       specificity_test_per_sample_v7_L${MATCH_LEVEL}.tsv
#       specificity_test_figure_summary_v7_L${MATCH_LEVEL}.tsv
#       specificity_test_qc_v7_L${MATCH_LEVEL}.json
#       specificity_test_event_level_v7_L${MATCH_LEVEL}.tsv
#
#   - v6 -> v7 変更点 (path override only; WGS discovery sample-level top-1%
#     CNV count QC の伝播):
#       1. _OVERLAP_TABLE: paths_v1.F_04_EVENT_OVERLAP (output_v9) ->
#          /lustre12/.../04_wgs_sv_boundary_overlap/output_v10/
#          sample_boundary_event_overlap_v10.tsv.gz
#       2. _BURDEN_TABLE: paths_v1.F_05_SAMPLE_BURDEN_L2 (output_v2) ->
#          /lustre12/.../05_wgs_sample_burden/output_v3/
#          sample_burden_L2_and_specificity_v3.tsv
#       3. 出力ファイル名サフィックス: _v6_L{MATCH_LEVEL} -> _v7_L{MATCH_LEVEL}
#       4. 出力先: output_tad_dynamic_boundary_specificity_v7_L{MATCH_LEVEL}
#       5. paths_v1.py 自体は変更しない (中央管理ファイルへの侵襲を避ける)。
#       6. 解析ロジック (matched-static null / Logit / event-level OR /
#          flanking control / class-wise + 10-class union, MATCH_LEVEL ladder)
#          は v6 と完全一致。
#       7. arrayCGH / MSSNG / GENCODE / SEGDUP / Heffel boundary などの
#          非 WGS 入力は v6 と完全同一。
#
#   - 依存:
#       Step 04 v10:
#         tad04212026/04_wgs_sv_boundary_overlap/output_v10/
#           sample_boundary_event_overlap_v10.tsv.gz
#       Step 05 v3:
#         tad04212026/05_wgs_sample_burden/output_v3/
#           sample_burden_L2_and_specificity_v3.tsv
#       paths_v1 経由 (変更なし):
#         Heffel L2 diffbound BED / BrainDev_raw.boundary.h5ad
#         GENCODE v46 GTF / segdup BED
#         OUT_07_MATCHED_STATIC
#
#   - Step 04 v10 + Step 05 v3 完了後に実行すること。
#
#   - 環境変数 (v6 と同一):
#       MATCH_LEVEL=2     (1=chrom only, 2=chrom+gene+segdup, 3=all 5 vars)
#       N_RESAMPLES=1000  (Level 2 推奨)
#       RANDOM_SEED=20260402
#       TARGET_CLASS=     (空=全クラス実行, 指定=1クラスのみ)
#       OUTDIR=           (デフォルト: paths_v1.OUT_07_MATCHED_STATIC / output_..._v7_L${MATCH_LEVEL})
#       RUN_EVENT_LEVEL=1
#       RUN_FLANKING_CONTROL=0
#
#   - 想定所要時間: 1000 resamples × 10 classes × MATCH_LEVEL=2 で 12-36h
#   - リソース: v6 sbatch と同一 (16 CPU / 128GB / 72h)

#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --job-name=tad_spec_v7
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/07_wgs_matched_static/tad_spec_v7_L%j.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/07_wgs_matched_static/tad_spec_v7_L%j.err

set -euo pipefail

echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "Working dir: $(pwd)"
START_TS=$(date +%s)

# v6 と同一のデフォルト env
export MATCH_LEVEL="${MATCH_LEVEL:-2}"
export N_RESAMPLES="${N_RESAMPLES:-1000}"
export RANDOM_SEED="${RANDOM_SEED:-20260402}"
export RUN_EVENT_LEVEL="${RUN_EVENT_LEVEL:-1}"
export RUN_FLANKING_CONTROL="${RUN_FLANKING_CONTROL:-0}"

echo "MATCH_LEVEL=${MATCH_LEVEL}"
echo "N_RESAMPLES=${N_RESAMPLES}"
echo "RANDOM_SEED=${RANDOM_SEED}"
echo "RUN_EVENT_LEVEL=${RUN_EVENT_LEVEL}"
echo "RUN_FLANKING_CONTROL=${RUN_FLANKING_CONTROL}"

cd /lustre12/home/kushima-pg/tad04212026/07_wgs_matched_static
mkdir -p "output_tad_dynamic_boundary_specificity_v7_L${MATCH_LEVEL}"
python3 tad_dynamic_boundary_specificity_test_v7.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "Job end: $(date -Iseconds)"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
