#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH -J constraint_enrich_v6
#SBATCH -o 14_v6_%j.log
#SBATCH -e 14_v6_%j.err
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 06:00:00

# 14_run_pipeline_v6.sbatch
# v6: v5 round 実行で判明した 2 つの根本問題を fix
# - bin_l2_annotation の group_primary は Diff_any ラベルを持たず
#   Diff_specific_n1 + Diff_shared_n2plus + Static の 3 値
#   → 14_compute_tad_domains_v4.py で正しく anchor 化し
#     output で Diff_any/Static に二値化
# - v19 patch v3 の _diag_label() で 'unaffected_sibling' を Sibling に追加
#   (v19.py 上書き済み、v19 sbatch 再投入で MSSNG dump regenerate 必要)
#
# Step 1: TAD domain v4 (anchor schema fix)
# Step 2: Gene annotation v2 (再利用、LOEUF logic 不変)
# Step 3: Per-sample burden v6 (BIN_DOM input path を _v4 に変更)
# Step 4: Statistical models v6 (input path のみ _v6 に変更、統計設計は v5 不変)

set -e
set -o pipefail

START_TIME=$(date +%s)
echo "=== Pipeline start (v6): $(date) ==="

cd /lustre12/home/kushima-pg/tad04212026/14_constraint_enrichment_v1/

# Required Python modules
python -c "import pandas, numpy, scipy, statsmodels" || {
    echo "ERROR: pandas / numpy / scipy / statsmodels が必要"
    echo "Run: pip install --break-system-packages pandas numpy scipy statsmodels"
    exit 1
}

# Verify MSSNG v19 patch v3 (post bug-fix) outputs exist
MSSNG_EVBIN="output_v1/mssng_event_bins_dumped_v1.tsv.gz"
MSSNG_COV="output_v1/mssng_sample_covariates_dumped_v1.tsv.gz"
if [ ! -f "$MSSNG_EVBIN" ] || [ ! -f "$MSSNG_COV" ]; then
    echo "ERROR: MSSNG patch v3 outputs not found:"
    echo "  Expected: $MSSNG_EVBIN"
    echo "  Expected: $MSSNG_COV"
    echo ""
    echo "  Re-run v19 with the patched _diag_label() (unaffected_sibling -> Sibling):"
    echo "    cd /lustre12/home/kushima-pg/tad04212026/09_mssng_sample_burden/"
    echo "    sbatch 09a_mssng_replication_v19.sbatch"
    exit 1
fi
echo "MSSNG patch outputs verified: OK"

# Sanity check: MSSNG dump must contain both ASD and Sibling rows
N_ASD=$(zcat "$MSSNG_COV" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Diagnosis") c=i; next} $c=="ASD"' | wc -l)
N_SIB=$(zcat "$MSSNG_COV" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Diagnosis") c=i; next} $c=="Sibling"' | wc -l)
echo "MSSNG dump: ASD=${N_ASD}, Sibling=${N_SIB}"
if [ "$N_SIB" -lt "100" ]; then
    echo "ERROR: MSSNG Sibling count too low (${N_SIB}). Re-run v19 with the patched"
    echo "       _diag_label() that maps 'unaffected_sibling' -> 'Sibling'."
    exit 1
fi

# Step 1: TAD domain v4 (NEW; anchor schema fix)
DOM_BED="output_v1/tad_domains_v4.bed"
BIN_DOM="output_v1/bin_to_flanking_domains_v4.tsv"
if [ ! -f "$DOM_BED" ] || [ ! -f "$BIN_DOM" ]; then
    echo ""
    echo "=== [Step 1/4] TAD domain reconstruction (v4: anchor schema fix) ==="
    python 14_compute_tad_domains_v4.py
else
    echo "[Step 1/4 skipped] TAD domain v4 outputs exist"
fi

# Step 2: Gene annotation v3 (input TAD bed を v4 に切替, LOEUF logic は v2 不変)
GENES_PER_DOM="output_v1/genes_per_domain_v3.tsv.gz"
if [ ! -f "$GENES_PER_DOM" ]; then
    echo ""
    echo "=== [Step 2/4] Gene annotation (v3: input tad_domains_v4.bed) ==="
    python 14_annotate_genes_per_domain_v3.py
else
    echo "[Step 2/4 skipped] Gene annotation v3 output exists"
fi

# Step 3: Per-sample burden v6 (BIN_DOM input path を _v4 に変更)
echo ""
echo "=== [Step 3/4] Per-sample constraint burden (v6: input bin_to_flanking_domains_v4) ==="
python 14_compute_per_sample_burden_v6.py

# Step 4: Statistical models v6
echo ""
echo "=== [Step 4/4] Statistical models + IVW meta (v6) ==="
python 14_fit_constraint_models_v6.py

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo ""
echo "=== Pipeline complete (v6): $(date) ==="
echo "=== Total elapsed: ${ELAPSED}s ($(($ELAPSED / 60))m $(($ELAPSED % 60))s) ==="

echo ""
echo "=== Output files (output_v1/) ==="
ls -la output_v1/ | grep -E '(_v6|_v4\.bed|_v4\.tsv|_v2\.tsv)'
