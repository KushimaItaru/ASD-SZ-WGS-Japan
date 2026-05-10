#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH -J constraint_enrich_v8
#SBATCH -o 14_v8_%j.log
#SBATCH -e 14_v8_%j.err
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 06:00:00

# 14_run_pipeline_v8.sbatch
# v7 → v8 変更点 (ChatGPT v7 review Block 1, 2 fix):
# - Step 4 を v8 fit script に切替
# - Block 1 fix: require_fit_ok() を IVW meta (beta_meta/se_meta) に対応
# - Block 2 fix: primary arrayCGH IVW / MSSNG GEE / Test 4 meta 不成立を hard fail
# - Step 1-3 は v7 を継承 (compute_diffany_bounded_domains_v7 は DIFF_ANY_LABELS
#   に 'Diff_any' 追加で robustness 強化、burden_v7 は docstring + proxy comment 修正)
# - 出力: constraint_enrichment_results_v8.tsv (v8 で primary が確実に hard-fail
#   検証された結果)

set -e
set -o pipefail

START_TIME=$(date +%s)
echo "=== Pipeline start (v8): $(date) ==="

cd /lustre12/home/kushima-pg/tad04212026/14_constraint_enrichment_v1/

python -c "import pandas, numpy, scipy, statsmodels" || {
    echo "ERROR: pandas / numpy / scipy / statsmodels が必要"
    exit 1
}

# Verify MSSNG v19 patch v3 (Sibling-corrected) outputs
MSSNG_EVBIN="output_v1/mssng_event_bins_dumped_v1.tsv.gz"
MSSNG_COV="output_v1/mssng_sample_covariates_dumped_v1.tsv.gz"
if [ ! -f "$MSSNG_EVBIN" ] || [ ! -f "$MSSNG_COV" ]; then
    echo "ERROR: MSSNG v19 patch outputs not found"
    exit 1
fi

N_SIB=$(zcat "$MSSNG_COV" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Diagnosis") c=i; next} $c=="Sibling"' | wc -l)
echo "MSSNG dump Sibling count: ${N_SIB}"
if [ "$N_SIB" -lt "100" ]; then
    echo "ERROR: MSSNG Sibling count too low (${N_SIB}). Re-run v19 with patched _diag_label()."
    exit 1
fi

# Step 1: Diff_any-bounded TAD-scale domain
DOM_BED="output_v1/diffany_bounded_domains_v7.bed"
BIN_DOM="output_v1/bin_to_diffany_bounded_domains_v7.tsv"
if [ ! -f "$DOM_BED" ] || [ ! -f "$BIN_DOM" ]; then
    echo ""
    echo "=== [Step 1/4] Diff_any-bounded TAD-scale domain reconstruction (v7) ==="
    python 14_compute_diffany_bounded_domains_v7.py
else
    echo "[Step 1/4 skipped] v7 domain outputs exist"
fi

# Step 2: Gene annotation v4 (input bed v7)
GENES_PER_DOM="output_v1/genes_per_domain_v4.tsv.gz"
if [ ! -f "$GENES_PER_DOM" ]; then
    echo ""
    echo "=== [Step 2/4] Gene annotation (v4: input diffany_bounded_domains_v7.bed) ==="
    python 14_annotate_genes_per_domain_v4.py
else
    echo "[Step 2/4 skipped] genes_per_domain_v4 exists"
fi

# Step 3: Per-sample burden v7 (5 variants in single file per cohort)
echo ""
echo "=== [Step 3/4] Per-sample constraint burden v7 (5 variants × 3 cohorts) ==="
python 14_compute_per_sample_burden_v7.py

# Step 4: Statistical models v8 (Primary + Sensitivity all variants)
# v8: Block 1 (require_fit_ok meta対応) + Block 2 (primary hard fail) fix
echo ""
echo "=== [Step 4/4] Statistical models + IVW meta (v8) ==="
python 14_fit_constraint_models_v8.py

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo ""
echo "=== Pipeline complete (v8): $(date) ==="
echo "=== Total elapsed: ${ELAPSED}s ($(($ELAPSED / 60))m $(($ELAPSED % 60))s) ==="

echo ""
echo "=== Output files (output_v1/) ==="
ls -la output_v1/ | grep -E '(_v7|_v8|diffany_bounded|_v4\.tsv|_v4\.tsv\.gz)'
