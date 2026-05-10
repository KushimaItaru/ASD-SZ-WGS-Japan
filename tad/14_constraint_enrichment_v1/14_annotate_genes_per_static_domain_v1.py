#!/usr/bin/env python3
"""
14_annotate_genes_per_static_domain_v1.py

処理内容 (v231 round, M3 対応 — Static-bounded null control):
- Static-bounded TAD-scale regulatory domain (v1) に gene を midpoint-based assign。
- v4 (Diff_any-bounded gene annotation) と同一ロジックで input bed のみ差し替え。

入力:
- /lustre12/home/kushima-pg/annotationInfo/gencode.v46.annotation.gtf.gz
- /lustre12/home/kushima-pg/annotationInfo/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
- output_v1/static_bounded_domains_v1.bed (Static-only anchor based)

出力:
- output_v1/genes_per_static_domain_v1.tsv.gz
- output_v1/static_domain_constraint_summary_v1.tsv

実行時間記録あり
"""

import gzip
import time
import pandas as pd
import numpy as np
from pathlib import Path

t0 = time.time()

ANNOT_DIR = Path('/lustre12/home/kushima-pg/annotationInfo')
GENCODE_GTF = ANNOT_DIR / 'gencode.v46.annotation.gtf.gz'
GNOMAD_LOEUF = ANNOT_DIR / 'gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz'
BASE_DIR = Path('/lustre12/home/kushima-pg/tad04212026')
DOM_BED = BASE_DIR / '14_constraint_enrichment_v1/output_v1/static_bounded_domains_v1.bed'
OUT_DIR = BASE_DIR / '14_constraint_enrichment_v1/output_v1'

LOEUF_STRICT = 0.35
LOEUF_RELAXED = 0.6
PLI_CUTOFF = 0.9

# ------------------ gnomAD ------------------
print(f"[{time.time()-t0:.1f}s] Loading gnomAD v2.1.1 ...")
gnomad_df = pd.read_csv(GNOMAD_LOEUF, sep='\t', compression='gzip',
                         usecols=['gene', 'transcript', 'oe_lof_upper', 'pLI'])
gnomad_df = gnomad_df.dropna(subset=['gene'])
gnomad_df = gnomad_df.sort_values('oe_lof_upper').drop_duplicates(subset='gene', keep='first')
gnomad_df = gnomad_df.rename(columns={'oe_lof_upper': 'LOEUF'})
n_with_loeuf = gnomad_df['LOEUF'].notna().sum()
print(f"  gnomAD genes total: {len(gnomad_df)}, with LOEUF: {n_with_loeuf}")

# ------------------ GENCODE ------------------
print(f"[{time.time()-t0:.1f}s] Parsing GENCODE v46 ...")
genes = []
with gzip.open(GENCODE_GTF, 'rt') as fh:
    for ln in fh:
        if ln.startswith('#'):
            continue
        parts = ln.rstrip('\n').split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue
        chrom = parts[0]
        start = int(parts[3]) - 1
        end = int(parts[4])
        strand = parts[6]
        attrs = {}
        for kv in parts[8].split(';'):
            kv = kv.strip()
            if not kv:
                continue
            sp = kv.split(' ', 1)
            if len(sp) == 2:
                attrs[sp[0]] = sp[1].strip('"')
        if attrs.get('gene_type', '') != 'protein_coding':
            continue
        genes.append({
            'gene_chr': chrom, 'gene_start': start, 'gene_end': end,
            'gene_strand': strand, 'gene_symbol': attrs.get('gene_name', ''),
            'gene_id': attrs.get('gene_id', ''),
        })
gene_df = pd.DataFrame(genes)
print(f"  GENCODE v46 protein-coding genes: {len(gene_df)}")

gene_df = gene_df.merge(gnomad_df[['gene', 'LOEUF', 'pLI']],
                          left_on='gene_symbol', right_on='gene', how='left').drop(columns=['gene'])

# Constraint flags
gene_df['is_constraint_known'] = gene_df['LOEUF'].notna().astype(int)
gene_df['is_constrained_strict'] = ((gene_df['LOEUF'] < LOEUF_STRICT) & gene_df['is_constraint_known'].astype(bool)).fillna(False).astype(int)
gene_df['is_unconstrained_strict_known'] = ((gene_df['LOEUF'] >= LOEUF_STRICT) & gene_df['is_constraint_known'].astype(bool)).fillna(False).astype(int)
gene_df['is_constrained_relaxed'] = ((gene_df['LOEUF'] < LOEUF_RELAXED) & gene_df['is_constraint_known'].astype(bool)).fillna(False).astype(int)
gene_df['is_unconstrained_relaxed_known'] = ((gene_df['LOEUF'] >= LOEUF_RELAXED) & gene_df['is_constraint_known'].astype(bool)).fillna(False).astype(int)
pli_known = gene_df['pLI'].notna()
gene_df['is_pLI_known'] = pli_known.astype(int)
gene_df['is_pLI_constrained'] = ((gene_df['pLI'] > PLI_CUTOFF) & pli_known).fillna(False).astype(int)
gene_df['is_pLI_unconstrained_known'] = ((gene_df['pLI'] <= PLI_CUTOFF) & pli_known).fillna(False).astype(int)

print(f"  constrained_strict (LOEUF<0.35): {gene_df['is_constrained_strict'].sum()}")
print(f"  unconstrained_strict_known: {gene_df['is_unconstrained_strict_known'].sum()}")
print(f"  constraint unknown: {(gene_df['is_constraint_known'] == 0).sum()}")

# ------------------ Static-bounded domain assign ------------------
print(f"[{time.time()-t0:.1f}s] Loading Static-bounded domains ...")
dom_df = pd.read_csv(DOM_BED, sep='\t')
print(f"  Domains: {len(dom_df)}")

gene_df['gene_mid'] = (gene_df['gene_start'] + gene_df['gene_end']) // 2

print(f"[{time.time()-t0:.1f}s] Assigning genes to Static-bounded domains (midpoint) ...")
domain_records = []
for chrom, dom_sub in dom_df.groupby('chrom'):
    gene_sub = gene_df[gene_df['gene_chr'] == chrom].copy()
    if gene_sub.empty:
        continue
    dom_sub = dom_sub.sort_values('start').reset_index(drop=True)
    gene_sub = gene_sub.sort_values('gene_mid').reset_index(drop=True)
    merged = pd.merge_asof(
        gene_sub[['gene_mid', 'gene_chr', 'gene_start', 'gene_end',
                   'gene_symbol', 'gene_id', 'gene_strand', 'LOEUF', 'pLI',
                   'is_constraint_known', 'is_constrained_strict',
                   'is_unconstrained_strict_known', 'is_constrained_relaxed',
                   'is_unconstrained_relaxed_known', 'is_pLI_known',
                   'is_pLI_constrained', 'is_pLI_unconstrained_known']],
        dom_sub[['start', 'end', 'domain_id']],
        left_on='gene_mid', right_on='start', direction='backward',
    )
    merged = merged[merged['gene_mid'] < merged['end']]
    merged = merged.drop(columns=['gene_mid', 'start', 'end'])
    domain_records.append(merged)

genes_per_dom = pd.concat(domain_records, ignore_index=True)
print(f"  Total gene-domain assignments: {len(genes_per_dom)}")

out_genes = OUT_DIR / 'genes_per_static_domain_v1.tsv.gz'
genes_per_dom.to_csv(out_genes, sep='\t', index=False, compression='gzip')
print(f"  Saved {out_genes}")

# ------------------ Per-domain summary ------------------
dom_summary = genes_per_dom.groupby('domain_id').agg(
    n_genes_total=('gene_symbol', 'count'),
    n_constraint_known=('is_constraint_known', 'sum'),
    n_constrained_strict=('is_constrained_strict', 'sum'),
    n_unconstrained_strict_known=('is_unconstrained_strict_known', 'sum'),
    n_constrained_relaxed=('is_constrained_relaxed', 'sum'),
    n_unconstrained_relaxed_known=('is_unconstrained_relaxed_known', 'sum'),
    n_pLI_known=('is_pLI_known', 'sum'),
    n_pLI_constrained=('is_pLI_constrained', 'sum'),
    n_pLI_unconstrained_known=('is_pLI_unconstrained_known', 'sum'),
).reset_index()

all_dom_ids = dom_df['domain_id'].tolist()
zero_doms = pd.DataFrame({
    'domain_id': [d for d in all_dom_ids if d not in dom_summary['domain_id'].values],
})
for c in dom_summary.columns:
    if c != 'domain_id':
        zero_doms[c] = 0
dom_summary = pd.concat([dom_summary, zero_doms], ignore_index=True)

out_summary = OUT_DIR / 'static_domain_constraint_summary_v1.tsv'
dom_summary.to_csv(out_summary, sep='\t', index=False)
print(f"  Saved {out_summary}")

print(f"\n[Done in {time.time()-t0:.1f}s]")
