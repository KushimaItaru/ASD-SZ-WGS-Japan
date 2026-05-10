#!/bin/bash
# ファイル名: run_phase_E_v1.sh
# 処理内容:
#   - tad04212026 パイプライン Scripts 1-7 を依存関係付きで一括 submit
#   - Step 1 と Step 3 を並列 submit
#   - Step 2 は Step 1 完了後 (afterok)
#   - Step 4 は Step 1 AND Step 3 完了後
#   - Step 5 は Step 2 AND Step 4 完了後
#   - Step 6 は Step 5 完了後
#   - Step 7 は Step 5 完了後 (Step 6 と並列可)
#   - 各 JOBID をログに記録

set -euo pipefail

ROOT="/lustre12/home/kushima-pg/tad04212026"
LOG_FILE="${ROOT}/phase_E_submit_$(date +%Y%m%d_%H%M%S).log"

{
echo "======================================================================"
echo "tad04212026 Phase E submit run"
echo "date: $(date -Iseconds)"
echo "ROOT: ${ROOT}"
echo "======================================================================"

# -----------------------------------------------------------
# Step 1 (no deps) and Step 3 (no deps) — 並列
# -----------------------------------------------------------
JID_1=$(sbatch --parsable "${ROOT}/01_heffel_boundary_master/01_build_master_v9.sbatch")
echo "Step 1 submitted: JobID=${JID_1}"

JID_3=$(sbatch --parsable "${ROOT}/03_wgs_sv_events/03_extract_sv_v8.sbatch")
echo "Step 3 submitted: JobID=${JID_3}"

# -----------------------------------------------------------
# Step 2 depends on Step 1
# -----------------------------------------------------------
JID_2=$(sbatch --parsable --dependency=afterok:${JID_1} \
  "${ROOT}/02_bin_l2_annotation/02_bin_l2_annotation_v2.sbatch")
echo "Step 2 submitted: JobID=${JID_2} (after Step 1=${JID_1})"

# -----------------------------------------------------------
# Step 4 depends on Step 1 AND Step 3
# -----------------------------------------------------------
JID_4=$(sbatch --parsable --dependency=afterok:${JID_1}:${JID_3} \
  "${ROOT}/04_wgs_sv_boundary_overlap/04_intersect_sv_v9.sbatch")
echo "Step 4 submitted: JobID=${JID_4} (after Step 1=${JID_1} AND Step 3=${JID_3})"

# -----------------------------------------------------------
# Step 5 depends on Step 2 AND Step 4
# -----------------------------------------------------------
JID_5=$(sbatch --parsable --dependency=afterok:${JID_2}:${JID_4} \
  "${ROOT}/05_wgs_sample_burden/05_sample_burden_v2.sbatch")
echo "Step 5 submitted: JobID=${JID_5} (after Step 2=${JID_2} AND Step 4=${JID_4})"

# -----------------------------------------------------------
# Step 6 depends on Step 5
# -----------------------------------------------------------
JID_6=$(sbatch --parsable --dependency=afterok:${JID_5} \
  "${ROOT}/06_wgs_primary_L2/39_fit_B_prime_L2_and_specificity_v3.sbatch")
echo "Step 6 submitted: JobID=${JID_6} (after Step 5=${JID_5})"

# -----------------------------------------------------------
# Step 7 depends on Step 5 (parallel with Step 6)
# -----------------------------------------------------------
JID_7=$(sbatch --parsable --dependency=afterok:${JID_5} \
  "${ROOT}/07_wgs_matched_static/tad_dynamic_boundary_specificity_test_v6.sbatch")
echo "Step 7 submitted: JobID=${JID_7} (after Step 5=${JID_5}, parallel with Step 6)"

echo "======================================================================"
echo "All 7 jobs submitted. JobIDs:"
echo "  Step 1: ${JID_1}"
echo "  Step 2: ${JID_2}  (afterok ${JID_1})"
echo "  Step 3: ${JID_3}"
echo "  Step 4: ${JID_4}  (afterok ${JID_1}:${JID_3})"
echo "  Step 5: ${JID_5}  (afterok ${JID_2}:${JID_4})"
echo "  Step 6: ${JID_6}  (afterok ${JID_5})"
echo "  Step 7: ${JID_7}  (afterok ${JID_5})"
echo ""
echo "Monitor with:"
echo "  squeue -u kushima-pg"
echo "  tail -F ${ROOT}/01_heffel_boundary_master/output_v9/01_master_${JID_1}.log"
echo "======================================================================"
} | tee "${LOG_FILE}"

echo ""
echo "Submit log saved: ${LOG_FILE}"
