#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH -J B_prime_v3
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o /lustre12/home/kushima-pg/tad04212026/06_wgs_primary_L2/output_v3/39_v3.out
#SBATCH -e /lustre12/home/kushima-pg/tad04212026/06_wgs_primary_L2/output_v3/39_v3.err

set -euo pipefail
SCRIPT_DIR="/lustre12/home/kushima-pg/tad04212026/06_wgs_primary_L2"
mkdir -p "${SCRIPT_DIR}/output_v3"
START=$(date +%s)
echo "Start: $(date)"
Rscript "${SCRIPT_DIR}/39_fit_B_prime_L2_and_specificity_v3.R"
END=$(date +%s)
echo "End  : $(date)"
echo "Elapsed: $((END-START))s"
