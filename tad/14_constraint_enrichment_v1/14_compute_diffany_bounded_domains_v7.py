#!/usr/bin/env python3
"""
14_compute_diffany_bounded_domains_v7.py

処理内容 (v6 → v7 設計変更):
- TAD domain 定義を「Diff_any anchor 間 interval」に変更 (Static は無視)
- v6 の all-boundary reconstruction (Diff_any + Static, 67,130 anchor → mean 86 kb の
  小 domain) では Diff_any flanking domain × constrained gene を持つ ASD sample が
  WGS=0/arrayCGH=1/MSSNG=1 と原理的 power 不足。
- v7 では Diff_any anchor (~4,980) のみで TAD-scale regulatory domain を再構築。
  期待 domain 数 ~2,400, mean ~600 kb (Hi-C TAD scale 200 kb-1 Mb と整合)。
- 命名: "Diff_any-bounded TAD-scale regulatory domains"
  Static boundary は domain subdivision に使わない (但し biologically irrelevant
  とは仮定しない; sensitivity QC として保持可)

入力:
- /lustre12/home/kushima-pg/tad04212026/02_bin_l2_annotation/output_v2/bin_l2_annotation_v2.tsv.gz

出力:
- output_v1/diffany_bounded_domains_v7.bed
- output_v1/bin_to_diffany_bounded_domains_v7.tsv

実行時間記録あり
"""

import time
import pandas as pd
from collections import defaultdict
from pathlib import Path

t0 = time.time()

BASE_DIR = Path('/lustre12/home/kushima-pg/tad04212026')
BIN_ANNOT = BASE_DIR / '02_bin_l2_annotation/output_v2/bin_l2_annotation_v2.tsv.gz'
OUT_DIR = BASE_DIR / '14_constraint_enrichment_v1/output_v1'
OUT_DIR.mkdir(parents=True, exist_ok=True)

# v7: Diff_any anchor のみ
# v7 (post ChatGPT review): 'Diff_any' を追加し、bin_l2_annotation_v2 の schema が
# 将来的に collapse 済み 'Diff_any' label に差し替わっても anchor 数 0 で fail
# しないよう保険を入れる。実 schema (Diff_specific_n1, Diff_shared_n2plus) は
# collapse 前なので両方を anchor として認識する。
DIFF_ANY_LABELS = ['Diff_any', 'Diff_specific_n1', 'Diff_shared_n2plus']
ANCHOR_LABELS = DIFF_ANY_LABELS  # Static は anchor に含めない


def collapse_to_diff_any(label):
    if label in DIFF_ANY_LABELS:
        return 'Diff_any'
    return label


# ============================================================
# Step 1: Load bin annotation
# ============================================================
print(f"[{time.time()-t0:.1f}s] Loading bin_l2_annotation_v2 ...")
bin_df = pd.read_csv(BIN_ANNOT, sep='\t', compression='gzip')
required = {'bin_id', 'chrom', 'start0', 'end', 'group_primary'}
missing = required - set(bin_df.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

print(f"  Total bins: {len(bin_df)}")
print(f"  group_primary distribution: {bin_df['group_primary'].value_counts().to_dict()}")

bin_df = bin_df.sort_values(['chrom', 'start0']).reset_index(drop=True)

# ============================================================
# Step 2: Define Diff_any-only anchors
# ============================================================
bin_df['is_boundary_anchor'] = bin_df['group_primary'].isin(ANCHOR_LABELS).astype(int)
bin_df['group_primary_collapsed'] = bin_df['group_primary'].apply(collapse_to_diff_any)

n_anchor = bin_df['is_boundary_anchor'].sum()
print(f"  Diff_any-only anchors: {n_anchor}")
print(f"    Diff_specific_n1: {(bin_df['group_primary'] == 'Diff_specific_n1').sum()}")
print(f"    Diff_shared_n2plus: {(bin_df['group_primary'] == 'Diff_shared_n2plus').sum()}")

# ============================================================
# Step 3: Reconstruct Diff_any-bounded TAD-scale domains
# ============================================================
print(f"[{time.time()-t0:.1f}s] Reconstructing Diff_any-bounded TAD-scale domains ...")
domains = []
for chrom, sub in bin_df.groupby('chrom', sort=False):
    sub = sub.sort_values('start0').reset_index(drop=True)
    anchor_idx = sub.index[sub['is_boundary_anchor'] == 1].tolist()
    if len(anchor_idx) < 2:
        continue
    for i in range(len(anchor_idx) - 1):
        left_anchor = sub.iloc[anchor_idx[i]]
        right_anchor = sub.iloc[anchor_idx[i + 1]]
        domain_start = int(left_anchor['end'])
        domain_end = int(right_anchor['start0'])
        if domain_end <= domain_start:
            continue
        domains.append({
            'chrom': chrom,
            'start': domain_start,
            'end': domain_end,
            'left_bin_id': left_anchor['bin_id'],
            'right_bin_id': right_anchor['bin_id'],
            'left_group_primary': 'Diff_any',
            'right_group_primary': 'Diff_any',
            'domain_size_bp': domain_end - domain_start,
        })

dom_df = pd.DataFrame(domains)
dom_df['domain_id'] = ['DOM_' + str(i+1).zfill(7) for i in range(len(dom_df))]
print(f"  Reconstructed Diff_any-bounded domains: {len(dom_df)}")
print(f"  Mean domain size: {dom_df['domain_size_bp'].mean()/1000:.1f} kb")
print(f"  Median domain size: {dom_df['domain_size_bp'].median()/1000:.1f} kb")
print(f"  Domain size distribution (kb):")
size_kb = dom_df['domain_size_bp'] / 1000
print(f"    min: {size_kb.min():.1f}, P25: {size_kb.quantile(0.25):.1f}, "
      f"P75: {size_kb.quantile(0.75):.1f}, max: {size_kb.max():.1f}")

# All flanking pairs are Diff_any-Diff_any by construction
n_diff_diff = len(dom_df)
print(f"  QC: All domain flanks are Diff_any-Diff_any: {n_diff_diff}")

out_bed = OUT_DIR / 'diffany_bounded_domains_v7.bed'
dom_df_bed = dom_df[['chrom', 'start', 'end', 'domain_id', 'left_bin_id',
                      'right_bin_id', 'left_group_primary', 'right_group_primary',
                      'domain_size_bp']]
dom_df_bed.to_csv(out_bed, sep='\t', index=False)
print(f"  Saved {out_bed}")

# ============================================================
# Step 4: Map each Diff_any boundary bin to its flanking domains
# ============================================================
print(f"[{time.time()-t0:.1f}s] Mapping Diff_any bins to flanking Diff_any-bounded domains ...")
chrom_dom_map = defaultdict(list)
for _, r in dom_df.iterrows():
    chrom_dom_map[r['chrom']].append((r['start'], r['end'], r['domain_id']))

bin_to_dom = []
for _, b in bin_df[bin_df['is_boundary_anchor'] == 1].iterrows():
    chrom = b['chrom']
    bin_start = int(b['start0'])
    bin_end = int(b['end'])
    domains_on_chr = chrom_dom_map.get(chrom, [])
    left_dom_id = None
    right_dom_id = None
    for ds, de, did in domains_on_chr:
        if de == bin_start:
            left_dom_id = did
        if ds == bin_end:
            right_dom_id = did
        if left_dom_id and right_dom_id:
            break
    bin_to_dom.append({
        'bin_id': b['bin_id'],
        'chrom': chrom,
        'start0': bin_start,
        'end': bin_end,
        'group_primary': 'Diff_any',  # All bins are Diff_any by construction
        'left_domain_id': left_dom_id if left_dom_id else '',
        'right_domain_id': right_dom_id if right_dom_id else '',
        'has_one_flank_only': int(not (left_dom_id and right_dom_id)),
    })

bin_dom_df = pd.DataFrame(bin_to_dom)
n_one_flank = bin_dom_df['has_one_flank_only'].sum()
print(f"  Mapped Diff_any bins: {len(bin_dom_df)}")
print(f"  One-flank-only bins (chromosome ends): {n_one_flank}")

out_bin_dom = OUT_DIR / 'bin_to_diffany_bounded_domains_v7.tsv'
bin_dom_df.to_csv(out_bin_dom, sep='\t', index=False)
print(f"  Saved {out_bin_dom}")

print(f"\n[Done in {time.time()-t0:.1f}s]")
