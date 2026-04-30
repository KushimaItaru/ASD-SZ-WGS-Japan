#!/bin/bash
#SBATCH -p ncbn-cpu
#SBATCH --account=ncbn-cpu
#SBATCH -J static_v1_pipeline
#SBATCH -o /lustre12/home/kushima-pg/tad04292026/14_constraint_enrichment_v1/log/static_v1_%j.out
#SBATCH -e /lustre12/home/kushima-pg/tad04292026/14_constraint_enrichment_v1/log/static_v1_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

# 14_run_pipeline_static_v1.sbatch
#
# 処理内容 (v231 round, M3 対応 — Static-bounded null control pipeline):
# - Step 1: Static-bounded TAD-scale domain construction
# - Step 2: Per-domain gene annotation (LOEUF / pLI)
# - Step 3: Per-sample 5-variant burden (3 cohort: WGS / arrayCGH / MSSNG)
# - Step 4: 3-cohort meta-analysis with B' framework + Holm primary 2-test family
#
# 期待結果:
# - Static-bounded で OR ≈ 1, P > 0.1 (null) → reviewer M3 反論への直接対応
# - 同じ B' framework / 同じ cohort / 同じ pre-registered design 下で
#   Diff_any anchor 由来の signal が essential であることを示せる

set -e

cd /lustre12/home/kushima-pg/tad04292026/14_constraint_enrichment_v1

mkdir -p log
mkdir -p output_v1

T0=$(date +%s)
echo "=== Static-bounded null control pipeline (v231 M3) start: $(date) ==="

# Step 1: Static-bounded TAD domain construction
echo ""
echo "[Step 1/4] Computing Static-bounded TAD-scale domains ..."
python 14_compute_static_bounded_domains_v1.py
echo "  Step 1 done at $(date)"

# Step 2: Gene annotation per Static-bounded domain
echo ""
echo "[Step 2/4] Annotating genes per Static-bounded domain ..."
python 14_annotate_genes_per_static_domain_v1.py
echo "  Step 2 done at $(date)"

# Step 3: Per-sample burden (3 cohort, 5 variants)
echo ""
echo "[Step 3/4] Computing per-sample burden (3 cohort) ..."
python 14_compute_per_sample_burden_static_v1.py
echo "  Step 3 done at $(date)"

# Step 4: 3-cohort fit + IVW meta + Holm
echo ""
echo "[Step 4/4] Fitting 3-cohort constraint models ..."
python 14_fit_constraint_models_static_v1.py
echo "  Step 4 done at $(date)"

T1=$(date +%s)
ELAPSED=$((T1 - T0))
echo ""
echo "=== Pipeline finished: $(date) (elapsed ${ELAPSED}s) ==="
echo ""
echo "Key output files:"
echo "  output_v1/static_bounded_domains_v1.bed"
echo "  output_v1/bin_to_static_bounded_domains_v1.tsv"
echo "  output_v1/genes_per_static_domain_v1.tsv.gz"
echo "  output_v1/static_domain_constraint_summary_v1.tsv"
echo "  output_v1/sample_constraint_burden_discovery_static_v1.tsv"
echo "  output_v1/sample_constraint_burden_arraycgh_static_v1.tsv"
echo "  output_v1/sample_constraint_burden_mssng_static_v1.tsv"
echo "  output_v1/constraint_enrichment_results_static_v1.tsv"
