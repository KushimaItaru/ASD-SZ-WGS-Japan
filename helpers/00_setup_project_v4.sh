#!/bin/bash
# ============================================================================
# NOTE: This script contains hardcoded paths specific to the original analysis
# environment (NIG supercomputer). To run in a different environment, update
# the paths marked with "# CONFIGURE" below, or set the corresponding
# environment variables defined in config_v1.sh.
# ============================================================================
# 00_setup_project_v4.sh
# - Create required directories (sample_lists, logs, outputs, resources, etc.)
# - Copy/symlink gene_regions_1kb_pad.bed from existing location (update path if needed)
# - Record execution time
# v2->v3: Updated Next steps for array job guidance (ARRAY_SIZE replacement note)
# v3→v4 Changes: Updated sed examples to match actual helper script placeholders

set -euo pipefail

START=$(date +%s)
TS(){ date '+%Y-%m-%d %H:%M:%S'; }

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
source "${PROJECT_ROOT}/config_v1.sh"

echo "[$(TS)] [INFO] Start setup: ${PROJECT_ROOT}"

mkdir -p "${RES_DIR}" "${WORK_DIR}" "${LOG_DIR}"
mkdir -p "${SAMPLE_LIST_DIR}" "${EHDN_OUT_DIR}" "${MERGED_NOVEL_DIR}" "${ANALYSIS_NOVEL_DIR}"
mkdir -p "${DEPTH_DIR}" "${DEPTH_PARTS_DIR}"
mkdir -p "${LOG_DIR}/ehdn" "${LOG_DIR}/depth" "${LOG_DIR}/merge" "${LOG_DIR}/burden"

# Prepare gene regions BED
if [ -f "${GENE_REGIONS_BED}" ]; then
  echo "[$(TS)] [INFO] Gene BED exists: ${GENE_REGIONS_BED}"
else
  # Copy from existing project (edit here if needed)
  SRC1="/lustre12/home/kushima-pg/ehdn_knownSTR_10302025/gene_regions_1kb_pad.bed"  # CONFIGURE
  SRC2="/lustre12/home/kushima-pg/strling_knownSTR_11222025/gene_regions_1kb_pad.bed"  # CONFIGURE
  if [ -f "${SRC1}" ]; then
    ln -s "${SRC1}" "${GENE_REGIONS_BED}"
    echo "[$(TS)] [INFO] Linked gene BED: ${GENE_REGIONS_BED} -> ${SRC1}"
  elif [ -f "${SRC2}" ]; then
    ln -s "${SRC2}" "${GENE_REGIONS_BED}"
    echo "[$(TS)] [INFO] Linked gene BED: ${GENE_REGIONS_BED} -> ${SRC2}"
  else
    echo "[$(TS)] [ERROR] gene_regions_1kb_pad.bed not found in default locations." >&2
    echo "[$(TS)] [ERROR] Please place or symlink it to: ${GENE_REGIONS_BED}" >&2
    exit 1
  fi
fi

# Tool checks (existence verification only)
command -v samtools >/dev/null 2>&1 || { echo "[$(TS)] [ERROR] samtools not found in PATH"; exit 2; }
command -v bedtools >/dev/null 2>&1 || { echo "[$(TS)] [ERROR] bedtools not found in PATH"; exit 2; }

# EHdn binary check
if [ ! -x "${EHDN_BIN}" ]; then
  echo "[$(TS)] [ERROR] EHDN_BIN not executable: ${EHDN_BIN}" >&2
  exit 2
fi
if [ ! -f "${REFERENCE_FASTA}" ]; then
  echo "[$(TS)] [ERROR] Reference FASTA not found: ${REFERENCE_FASTA}" >&2
  exit 2
fi
if [ ! -f "${SAMPLE_INFO}" ]; then
  echo "[$(TS)] [ERROR] SampleInfo not found: ${SAMPLE_INFO}" >&2
  exit 2
fi

END=$(date +%s)
ELAPSED=$((END-START))
echo "[$(TS)] [DONE] Setup completed. Elapsed=${ELAPSED}s"
echo ""
echo "Next steps (execute in order):"
echo "  1) python3 helpers/ehdn/01_prepare_sample_lists_v2.py"
echo "  2) helpers/strling/01_calc_depth_array_fast_v2.sh  (array job)"
echo "     NOTE: ARRAY_SIZE placeholder must be replaced before sbatch."
echo "     Use:  N=\$(tail -n +2 sample_lists/ehdn_all_samples.tsv | wc -l)"
echo "           sed \"s/ARRAY_SIZE/\${N}/\" helpers/strling/01_calc_depth_array_fast_v2.sh | sbatch"
echo "  3) python3 helpers/strling/03_collect_depths_v1.py   (after depth jobs finish)"
echo "  4) helpers/ehdn/02_run_ehdn_array_v2.sh  (array job)"
echo "     NOTE: Same ARRAY_SIZE replacement needed."
echo "           sed \"s/ARRAY_SIZE/\${N}/\" helpers/ehdn/02_run_ehdn_array_v2.sh | sbatch"
echo "  5) python3 helpers/ehdn/04_merge_ehdn_novel_norm_v2.py  (after EHdn jobs finish)"
echo "  6) python3 helpers/ehdn/14_outlier_burden_rare_casecontrol_crossfit_v19.py"
echo "  7) python3 helpers/ehdn/17_burden_statistical_test_v20.py"
echo ""
echo "Or use the wrapper scripts in the top-level directories (recommended)."
