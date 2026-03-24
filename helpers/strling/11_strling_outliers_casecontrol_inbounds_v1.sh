#!/bin/bash
# ============================================================================
# NOTE: This script contains hardcoded paths specific to the original analysis
# environment (NIG supercomputer). To run in a different environment, update
# the paths marked with "# CONFIGURE" below, or set the corresponding
# environment variables defined in config_v1.sh.
# ============================================================================
# 11_strling_outliers_casecontrol_inbounds_v1.sh
# - Description:
#   - Aggregate genotypes from calls_genic_inbounds to outliers working directory for casecontrol_samples.tsv
#   - Retrieve unplaced from calls_genic, filter by in-bounds key (3471 loci) if possible (empty file if missing)
#   - Run strling-outliers.py to generate STRs.tsv and control-file.tsv
#   - Record execution log, processed sample count, and missing count; Record execution time
#
# Usage:
#   Run via the top-level wrapper: strling/02_strling_casecontrol_call_and_outliers.sh
#
# Output:
#   <repo_root>/strling_output_genomewide/outliers_casecontrol_inbounds_v1/links/STRs.tsv
#   <repo_root>/strling_output_genomewide/outliers_casecontrol_inbounds_v1/links/control-file.tsv
#   <repo_root>/strling_output_genomewide/outliers_casecontrol_inbounds_v1/links/outliers_run.log
#   <repo_root>/strling_output_genomewide/outliers_casecontrol_inbounds_v1/summary.txt
#
# Notes:
# - Genotypes are pre-filtered to in-bounds only (calls_genic_inbounds); OUT_OF_BOUNDS excluded from main analysis
# - unplaced is also filtered by in-bounds if possible (copy as-is or use empty file if format differs)

#SBATCH -J strling_outliers_inb
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --mem=256G
#SBATCH --output=logs/outliers_inbounds_%A.out
#SBATCH --error=logs/outliers_inbounds_%A.err

set -euo pipefail
START=$(date +%s)
TS(){ date '+%Y-%m-%d %H:%M:%S'; }

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/strling_output_genomewide}"

LOG_DIR="${OUT_ROOT}/logs"
WORK_DIR="${OUT_ROOT}/outliers_casecontrol_inbounds_v1"
LINK_DIR="${WORK_DIR}/links"

# Input
CASECONTROL_TSV="${PROJECT_ROOT}/sample_lists/casecontrol_samples.tsv"
GENO_IN_DIR="${OUT_ROOT}/calls_genic_inbounds"     # in-bounds genotype
UNPLACED_DIR="${OUT_ROOT}/calls_genic"             # Original unplaced (filtered as needed)
BOUNDS="${OUT_ROOT}/str-results/joint-bounds.genic_1kbpad.len3_8.txt"

mkdir -p "${LOG_DIR}" "${WORK_DIR}" "${LINK_DIR}"

echo "[$(TS)] [INFO] START 11_strling_outliers_casecontrol_inbounds_v1.sh"
echo "[$(TS)] [INFO] Host=$(hostname) JobID=${SLURM_JOB_ID:-NA}"
echo "[$(TS)] [INFO] WORK_DIR=${WORK_DIR}"
echo "[$(TS)] [INFO] LINK_DIR=${LINK_DIR}"
echo "[$(TS)] [INFO] CASECONTROL_TSV=${CASECONTROL_TSV}"
echo "[$(TS)] [INFO] GENO_IN_DIR=${GENO_IN_DIR}"
echo "[$(TS)] [INFO] UNPLACED_DIR=${UNPLACED_DIR}"
echo "[$(TS)] [INFO] BOUNDS=${BOUNDS}"

# tool check
command -v strling-outliers.py >/dev/null 2>&1 || { echo "[$(TS)] [ERROR] strling-outliers.py not found in PATH"; exit 2; }
command -v awk >/dev/null 2>&1 || { echo "[$(TS)] [ERROR] awk not found"; exit 2; }

if [ ! -f "${CASECONTROL_TSV}" ]; then
  echo "[$(TS)] [ERROR] casecontrol sample list not found: ${CASECONTROL_TSV}" >&2
  exit 3
fi
if [ ! -d "${GENO_IN_DIR}" ]; then
  echo "[$(TS)] [ERROR] genotype dir not found: ${GENO_IN_DIR}" >&2
  exit 3
fi
if [ ! -f "${BOUNDS}" ]; then
  echo "[$(TS)] [ERROR] bounds file not found: ${BOUNDS}" >&2
  exit 3
fi

# Create in-bounds key (chr,left,right,repeatunit)
KEYS="${WORK_DIR}/inbounds_3471.keys.tsv"
tail -n +2 "${BOUNDS}" | awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' > "${KEYS}"
NKEYS=$(wc -l < "${KEYS}" | tr -d ' ')
echo "[$(TS)] [INFO] IN_BOUNDS keys: ${NKEYS}"

# Clean working directory (avoid confusion on re-run)
# Large files should be Missing, but delete only the contents since there are many links
rm -f "${LINK_DIR}"/*-genotype.txt "${LINK_DIR}"/*-unplaced.txt "${LINK_DIR}"/STRs.tsv "${LINK_DIR}"/control-file.tsv "${LINK_DIR}"/outliers_run.log 2>/dev/null || true

# Extract SampleID from casecontrol list
SIDS=$(awk -F'\t' 'NR>1{print $1}' "${CASECONTROL_TSV}" | sort -u)

N_TOTAL=0
N_GENO_OK=0
N_GENO_MISSING=0
N_UNPLACED_OK=0
N_UNPLACED_EMPTY=0
N_UNPLACED_MISSING=0

echo "[$(TS)] [INFO] Linking/filtering files into ${LINK_DIR} ..."

for sid in ${SIDS}; do
  N_TOTAL=$((N_TOTAL+1))

  g_in="${GENO_IN_DIR}/${sid}-genotype.txt"
  g_link="${LINK_DIR}/${sid}-genotype.txt"

  if [ -s "${g_in}" ]; then
    ln -sf "${g_in}" "${g_link}"
    N_GENO_OK=$((N_GENO_OK+1))
  else
    # If genotype is Missing/empty, exclude this sample from outliers input (preserve analysis population)
    N_GENO_MISSING=$((N_GENO_MISSING+1))
    continue
  fi

  # Handle unplaced:
  # - If source file exists, attempt in-bounds key filter (copy if format differs)
  # - Create empty file if missing (prevent outliers from failing)
  u_in="${UNPLACED_DIR}/${sid}-unplaced.txt"
  u_out="${LINK_DIR}/${sid}-unplaced.txt"

  if [ -s "${u_in}" ]; then
    # If header starts with #chrom, filter by 4-column key for consistency
    first=$(head -n 1 "${u_in}" || true)
    if echo "${first}" | grep -q '^#chrom'; then
      awk -F'\t' '
        BEGIN{OFS="\t"}
        NR==FNR{key[$1 FS $2 FS $3 FS $4]=1; next}
        NR==1{print; next}
        (($1 FS $2 FS $3 FS $4) in key){print}
      ' "${KEYS}" "${u_in}" > "${u_out}"
      N_UNPLACED_OK=$((N_UNPLACED_OK+1))
    else
      # Copy as-is if format is unknown (safe fallback)
      cp -f "${u_in}" "${u_out}"
      N_UNPLACED_OK=$((N_UNPLACED_OK+1))
    fi
  elif [ -f "${u_in}" ]; then
    # File exists but is empty
    : > "${u_out}"
    N_UNPLACED_EMPTY=$((N_UNPLACED_EMPTY+1))
  else
    # Missing
    : > "${u_out}"
    N_UNPLACED_MISSING=$((N_UNPLACED_MISSING+1))
  fi
done

echo "[$(TS)] [INFO] N_TOTAL(casecontrol unique)=${N_TOTAL}"
echo "[$(TS)] [INFO] genotype: OK=${N_GENO_OK} MISSING/EMPTY=${N_GENO_MISSING}"
echo "[$(TS)] [INFO] unplaced: OK=${N_UNPLACED_OK} EMPTY=${N_UNPLACED_EMPTY} MISSING=${N_UNPLACED_MISSING}"

# Run outliers
cd "${LINK_DIR}"

# Abort if input count is too low
N_GENO_FILES=$(ls -1 *-genotype.txt 2>/dev/null | wc -l | tr -d ' ')
if [ "${N_GENO_FILES}" -lt 100 ]; then
  echo "[$(TS)] [ERROR] Too few genotype files in ${LINK_DIR}: ${N_GENO_FILES}" >&2
  exit 4
fi

echo "[$(TS)] [INFO] Running strling-outliers.py (n_genotype=${N_GENO_FILES}) ..."
set +e
strling-outliers.py --genotypes "*-genotype.txt" --unplaced "*-unplaced.txt" --emit control-file.tsv \
  > outliers_run.log 2>&1
RET=$?
set -e

if (( RET != 0 )); then
  echo "[$(TS)] [ERROR] strling-outliers.py failed (exit=${RET}). See: ${LINK_DIR}/outliers_run.log" >&2
  exit "${RET}"
fi

if [ ! -s "STRs.tsv" ]; then
  echo "[$(TS)] [ERROR] STRs.tsv not generated or empty. See: ${LINK_DIR}/outliers_run.log" >&2
  exit 5
fi
if [ ! -s "control-file.tsv" ]; then
  echo "[$(TS)] [WARN] control-file.tsv not generated or empty (may be OK depending on strling version)." >&2
fi

# summary
SUMMARY="${WORK_DIR}/summary.txt"
{
  echo "=== STRling outliers (IN_BOUNDS only) summary ==="
  echo "Created: $(TS)"
  echo "Project: ${PROJECT_ROOT}"
  echo "Bounds : ${BOUNDS}"
  echo "Keys   : ${KEYS} (n=${NKEYS})"
  echo "Casecontrol unique: ${N_TOTAL}"
  echo "Genotype files used: ${N_GENO_FILES}"
  echo "Missing/empty genotype (skipped): ${N_GENO_MISSING}"
  echo "Unplaced OK/EMPTY/MISSING: ${N_UNPLACED_OK}/${N_UNPLACED_EMPTY}/${N_UNPLACED_MISSING}"
  echo "Outputs:"
  echo "  ${LINK_DIR}/STRs.tsv"
  echo "  ${LINK_DIR}/control-file.tsv"
  echo "  ${LINK_DIR}/outliers_run.log"
  echo
  echo "[HEAD STRs.tsv]"
  head -n 3 "${LINK_DIR}/STRs.tsv" || true
} > "${SUMMARY}"

END=$(date +%s)
echo "[$(TS)] [DONE] Outliers completed. Elapsed=$((END-START))s"
echo "[$(TS)] [DONE] Summary: ${SUMMARY}"
