#!/bin/bash
# ファイル名: tad_dynamic_boundary_specificity_test_v8.sbatch
# 処理内容:
#   - Step 07 v7 (tad_dynamic_boundary_specificity_test_v7.py) を ncbn-cpu で実行する
#     新しい sbatch ラッパー (v7 sbatch を upgrade)。
#   - v7 sbatch との変更点 (Python スクリプト本体は v7.py のまま、変更不要):
#       1. N_RESAMPLES のデフォルトを 1,000 → 10,000 に変更
#          (ChatGPT review 指摘: 1,000 resamples では Monte Carlo error が大きく、
#           empirical P=0.032 のような値が ±0.005 のばらつきを持ち得る。
#           10,000 で min P_emp = 1/10001 ≈ 9.99e-5、解像度が 10倍向上。)
#       2. SLURM array job を導入し、11 クラス (10 primary + dynamic_union) を
#          並列実行することで、N_RESAMPLES 10倍化の追加 wall-time を相殺。
#          各 array task で TARGET_CLASS を 1 クラスのみ指定。
#          ncbn-cpu (576 CPU 共有) では 11 array × 16 CPU = 176 CPU で
#          ユーザー preference の 20-40 ジョブ保守的範囲内。
#       3. 出力ファイル名サフィックスは v7 同様: _v7_L${MATCH_LEVEL}_${TARGET_CLASS}
#          (v7.py 自体が _v7 サフィックスを書き出す。Python スクリプトは v7 のまま。)
#       4. 出力先ディレクトリ: paths_v1.OUT_07_MATCHED_STATIC /
#          output_tad_dynamic_boundary_specificity_v7_L${MATCH_LEVEL}/
#          (v7 と完全同一。10,000 resamples の結果と 1,000 resamples の旧結果は
#           ファイル名が同一になるため、旧結果を別場所に退避してから実行のこと。)
#       5. wall-time を 72h → 60h に短縮 (1 クラスあたり 12-36h 想定、十分余裕)。
#          全クラス並列のため total wall-time も 60h 以内で完了見込み。
#
#   - 依存 (v7 と同一):
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
#   - 環境変数 (v7 と同一構造、N_RESAMPLES デフォルトのみ変更):
#       MATCH_LEVEL=2     (1=chrom only, 2=chrom+gene+segdup, 3=all 5 vars)
#       N_RESAMPLES=10000 (v8 で default を 1000 → 10000 に変更)
#       RANDOM_SEED=20260402
#       TARGET_CLASS=     (array task ごとに自動設定)
#       OUTDIR=           (デフォルト: paths_v1.OUT_07_MATCHED_STATIC / output_..._v7_L${MATCH_LEVEL})
#       RUN_EVENT_LEVEL=1
#       RUN_FLANKING_CONTROL=0
#
#   - 想定所要時間:
#       v7 (1,000 resamples × 10 classes 直列): 12-36h
#       v8 (10,000 resamples × 11 classes 並列, array): 12-36h (各 task ≈ 同じ)
#         理由: 1 class あたりの計算時間 (resample × 検定) が線形に増えるが、
#               array で 11 task 並列実行のため total wall-time は維持。
#
#   - 投入手順:
#       (1) 既存 output_tad_dynamic_boundary_specificity_v7_L2/ を待避:
#             cd /lustre12/home/kushima-pg/tad04212026/07_wgs_matched_static
#             mv output_tad_dynamic_boundary_specificity_v7_L2 output_tad_dynamic_boundary_specificity_v7_L2_1000resamples_BAK
#       (2) v8 sbatch を tad04212026/07_wgs_matched_static/ に転送し sbatch:
#             sbatch tad_dynamic_boundary_specificity_test_v8.sbatch
#       (3) 11 array task すべて完了後、出力ファイルが揃っているか確認:
#             ls output_tad_dynamic_boundary_specificity_v7_L2/specificity_test_main_v7_L2_*.tsv | wc -l   # → 11
#       (4) 統合スクリプト (paper-side aggregation) で 11 クラスを集約。

#SBATCH -J tad_spec_v8
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH --array=1-11
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=60:00:00
#SBATCH --output=/lustre12/home/kushima-pg/tad04212026/07_wgs_matched_static/tad_spec_v8_array%A_%a.log
#SBATCH --error=/lustre12/home/kushima-pg/tad04212026/07_wgs_matched_static/tad_spec_v8_array%A_%a.err

set -euo pipefail

echo "==============================================="
echo "Job start: $(date -Iseconds)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-local}"
echo "SLURM_ARRAY_JOB_ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "Working dir: $(pwd)"
echo "==============================================="
START_TS=$(date +%s)

# --- v8 default env ---
export MATCH_LEVEL="${MATCH_LEVEL:-2}"
export N_RESAMPLES="${N_RESAMPLES:-10000}"   # v8: 1000 → 10000
export RANDOM_SEED="${RANDOM_SEED:-20260402}"
export RUN_EVENT_LEVEL="${RUN_EVENT_LEVEL:-1}"
export RUN_FLANKING_CONTROL="${RUN_FLANKING_CONTROL:-0}"

# --- v8 array task → TARGET_CLASS mapping ---
# v7.py の PRIMARY_BOUNDARY_CLASSES + dynamic_union と完全一致 (順序固定で再現性保証)
CLASSES=(
  "hpc_exc_ca"     # 1
  "hpc_exc_dg"     # 2
  "hpc_exc_ent"    # 3
  "hpc_inh_cge"    # 4
  "hpc_inh_mge"    # 5
  "pfc_astro"      # 6
  "pfc_exc_dl"     # 7
  "pfc_exc_ul"     # 8
  "pfc_inh_cge"    # 9
  "pfc_inh_mge"    # 10
  "dynamic_union"  # 11
)

TASK_IDX=$((SLURM_ARRAY_TASK_ID - 1))
if [ "${TASK_IDX}" -lt 0 ] || [ "${TASK_IDX}" -ge "${#CLASSES[@]}" ]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} out of range 1-${#CLASSES[@]}"
  exit 2
fi
export TARGET_CLASS="${CLASSES[${TASK_IDX}]}"

echo "MATCH_LEVEL=${MATCH_LEVEL}"
echo "N_RESAMPLES=${N_RESAMPLES}"
echo "RANDOM_SEED=${RANDOM_SEED}"
echo "TARGET_CLASS=${TARGET_CLASS}"
echo "RUN_EVENT_LEVEL=${RUN_EVENT_LEVEL}"
echo "RUN_FLANKING_CONTROL=${RUN_FLANKING_CONTROL}"
echo "==============================================="

cd /lustre12/home/kushima-pg/tad04212026/07_wgs_matched_static
mkdir -p "output_tad_dynamic_boundary_specificity_v7_L${MATCH_LEVEL}"

python3 tad_dynamic_boundary_specificity_test_v7.py

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))
echo "==============================================="
echo "Job end: $(date -Iseconds)"
echo "TARGET_CLASS=${TARGET_CLASS}"
echo "N_RESAMPLES=${N_RESAMPLES}"
echo "Elapsed seconds: ${ELAPSED}"
printf "Elapsed: %02d:%02d:%02d\n" $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60))
echo "==============================================="
