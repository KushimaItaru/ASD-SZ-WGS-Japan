#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH -J constraint_enrich_v9
#SBATCH -o 14_v9_%j.log
#SBATCH -e 14_v9_%j.err
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 06:00:00

# 14_run_pipeline_v9.sbatch
# v8 → v9 変更点:
# - direct-overlap exclusion を proxy → exact に変更
# - WGS event_bin: sv_chr/sv_start0/sv_end 直接使用
# - arrayCGH event_bin: event_id parse して sv_chr/sv_start/sv_end 抽出
# - MSSNG event_bin: v19 patch v4 で追加した sv_chr/sv_start/sv_end 列を使用
# - Step 3 を v9 burden script に切替、Step 4 を v9 fit script に切替
# - 統計設計・統計family・Holm 補正は v8 から不変
#
# 前提: v19 patch v4 が NIG で適用済み (mssng_event_bins_dumped_v1.tsv.gz が
#       sv_chr/sv_start/sv_end 列を含むこと)

set -e
# v9: set -o pipefail は zcat | head -1 で SIGPIPE 起こして set -e と組合せで
# 早期 exit するため使わない (v8 では使っていたが、v9 では schema check で
# zcat | head | awk のチェーンが多いため除外)

START_TIME=$(date +%s)
echo "=== Pipeline start (v9): $(date) ==="

cd /lustre12/home/kushima-pg/tad04292026/14_constraint_enrichment_v1/

python -c "import pandas, numpy, scipy, statsmodels" || {
    echo "ERROR: pandas / numpy / scipy / statsmodels が必要"
    exit 1
}

# Verify MSSNG v19 patch v4 (SV-interval-aware) outputs
MSSNG_EVBIN="output_v1/mssng_event_bins_dumped_v1.tsv.gz"
MSSNG_COV="output_v1/mssng_sample_covariates_dumped_v1.tsv.gz"
if [ ! -f "$MSSNG_EVBIN" ] || [ ! -f "$MSSNG_COV" ]; then
    echo "ERROR: MSSNG v19 dumps not found"
    exit 1
fi

# v9: require sv_chr/sv_start/sv_end columns (PATCH v4)
# Dump header to temp file first to avoid SIGPIPE issues with `zcat | head`
HEADER_TMP=$(mktemp)
zcat "$MSSNG_EVBIN" 2>/dev/null | head -1 > "$HEADER_TMP" || true
N_COL=$(awk -F'\t' '{print NF}' "$HEADER_TMP")
HAS_SV_INTERVAL=$(awk -F'\t' '{
    has=0; for(i=1;i<=NF;i++) if($i=="sv_chr" || $i=="sv_start" || $i=="sv_end") has++;
    print has
}' "$HEADER_TMP")
rm -f "$HEADER_TMP"
echo "MSSNG event-bin: ${N_COL} columns, SV interval columns found: ${HAS_SV_INTERVAL}"
if [ "${HAS_SV_INTERVAL:-0}" -lt "3" ]; then
    echo "ERROR: MSSNG event-bin dump missing SV interval columns (sv_chr/sv_start/sv_end)."
    echo "       Re-run v19 with PATCH_v18_to_v19_v4.md applied."
    exit 1
fi

# Sibling count check (v3 patch ガード継承; SIGPIPE-safe)
N_SIB=$(zcat "$MSSNG_COV" 2>/dev/null | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Diagnosis") c=i; next} $c=="Sibling"' | wc -l)
echo "MSSNG dump Sibling count: ${N_SIB}"
if [ "${N_SIB:-0}" -lt "100" ]; then
    echo "ERROR: MSSNG Sibling count too low (${N_SIB})."
    exit 1
fi

# Step 1: Diff_any-bounded TAD-scale domain (v7 を再利用、不変)
DOM_BED="output_v1/diffany_bounded_domains_v7.bed"
BIN_DOM="output_v1/bin_to_diffany_bounded_domains_v7.tsv"
if [ ! -f "$DOM_BED" ] || [ ! -f "$BIN_DOM" ]; then
    echo ""
    echo "=== [Step 1/4] Diff_any-bounded TAD-scale domain reconstruction (v7) ==="
    python 14_compute_diffany_bounded_domains_v7.py
else
    echo "[Step 1/4 skipped] v7 domain outputs exist"
fi

# Step 2: Gene annotation v4 (input bed v7、不変)
GENES_PER_DOM="output_v1/genes_per_domain_v4.tsv.gz"
if [ ! -f "$GENES_PER_DOM" ]; then
    echo ""
    echo "=== [Step 2/4] Gene annotation (v4) ==="
    python 14_annotate_genes_per_domain_v4.py
else
    echo "[Step 2/4 skipped] genes_per_domain_v4 exists"
fi

# Step 3: Per-sample burden v9 (exact direct-overlap exclusion)
echo ""
echo "=== [Step 3/4] Per-sample constraint burden v9 (exact direct-overlap) ==="
python 14_compute_per_sample_burden_v9.py

# Step 4: Statistical models v9
echo ""
echo "=== [Step 4/4] Statistical models + IVW meta (v9) ==="
python 14_fit_constraint_models_v9.py

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo ""
echo "=== Pipeline complete (v9): $(date) ==="
echo "=== Total elapsed: ${ELAPSED}s ($(($ELAPSED / 60))m $(($ELAPSED % 60))s) ==="

echo ""
echo "=== Output files (output_v1/) ==="
ls -la output_v1/ | grep -E '(_v9|diffany_bounded|_v4\.tsv|_v4\.tsv\.gz)'
