#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH -J constraint_enrich_v5
#SBATCH -o 14_v5_%j.log
#SBATCH -e 14_v5_%j.err
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 06:00:00

# 14_run_pipeline_v5.sbatch
# v5: ChatGPT v4 review safety patch
# - Priority 1-1: primary model failure を hard fail (require_fit_ok)
# - Priority 1-2: MSSNG required cov NaN を hard fail
# - Priority 1-3: GEE matrix API + zero-variance covariate drop + dummy 名 sanitize
# - Priority 1-4: arrayCGH stratified IVW を primary に切替
# - Priority 1-5: contrast external IVW meta supporting row 追加
# Step 1 (TAD domain): v3 を再利用 (anchor logic 不変)
# Step 2 (gene annotation): v2 を再利用 (LOEUF logic 不変)
# Step 3-4: v5 で実行

set -e  # any command failure → abort
set -o pipefail

START_TIME=$(date +%s)
echo "=== Pipeline start (v5): $(date) ==="

cd /lustre12/home/kushima-pg/tad04212026/14_constraint_enrichment_v1/

# Required Python modules
python -c "import pandas, numpy, scipy, statsmodels" || {
    echo "ERROR: pandas / numpy / scipy / statsmodels が必要"
    echo "Run: pip install --break-system-packages pandas numpy scipy statsmodels"
    exit 1
}

# Verify MSSNG v19 patch v3 outputs exist (hard fail if missing)
MSSNG_EVBIN="output_v1/mssng_event_bins_dumped_v1.tsv.gz"
MSSNG_COV="output_v1/mssng_sample_covariates_dumped_v1.tsv.gz"
if [ ! -f "$MSSNG_EVBIN" ] || [ ! -f "$MSSNG_COV" ]; then
    echo "ERROR: MSSNG patch v3 outputs not found:"
    echo "  Expected: $MSSNG_EVBIN"
    echo "  Expected: $MSSNG_COV"
    echo ""
    echo "  Run the v19 patched MSSNG pipeline first (use PATCH_v18_to_v19_v3.md):"
    echo "    cd /lustre12/home/kushima-pg/tad04212026/09_mssng_sample_burden/"
    echo "    sbatch 09a_mssng_replication_v19.sbatch"
    exit 1
fi
echo "MSSNG patch outputs verified: OK"

# Verify v3 TAD domain outputs exist (Step 1 reuses v3)
DOM_BED="output_v1/tad_domains_v3.bed"
BIN_DOM="output_v1/bin_to_flanking_domains_v3.tsv"
if [ ! -f "$DOM_BED" ] || [ ! -f "$BIN_DOM" ]; then
    echo "[Step 1/4 needed] TAD domain reconstruction (v3, anchor unchanged)"
    python 14_compute_tad_domains_v3.py
else
    echo "[Step 1/4 skipped] TAD domain v3 outputs exist"
fi

# Verify v2 gene annotation exists (Step 2 reuses v2)
GENES_PER_DOM="output_v1/genes_per_domain_v2.tsv.gz"
if [ ! -f "$GENES_PER_DOM" ]; then
    echo "[Step 2/4 needed] Gene annotation (v2, LOEUF logic unchanged)"
    python 14_annotate_genes_per_domain_v2.py
else
    echo "[Step 2/4 skipped] Gene annotation v2 output exists"
fi

# Step 3: Per-sample exposure (v5)
echo ""
echo "=== [Step 3/4] Per-sample constraint burden (v5: hard fail on missing cov) ==="
python 14_compute_per_sample_burden_v5.py

# Step 4: Statistical models + IVW meta (v5)
echo ""
echo "=== [Step 4/4] Statistical models + IVW meta (v5: hard fail on primary failure) ==="
python 14_fit_constraint_models_v5.py

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo ""
echo "=== Pipeline complete (v5): $(date) ==="
echo "=== Total elapsed: ${ELAPSED}s ($(($ELAPSED / 60))m $(($ELAPSED % 60))s) ==="

echo ""
echo "=== Output files (output_v1/) ==="
ls -la output_v1/ | grep -E '(_v5|_v3\.bed|_v2\.tsv)'
