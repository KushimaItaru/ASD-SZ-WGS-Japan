#!/usr/bin/env python3
"""
14_compute_static_bounded_domains_v1.py

処理内容 (v231 round, M3 対応 — circular reasoning 反論回避):
- Diff_any-bounded TAD-scale regulatory domain (v7) の対称的 null control として、
  Static anchor のみで TAD domain を再構築。
- Reviewer M3 反論「Diff_any anchor → Diff_any-bounded domain → enrichment は circular」
  への直接対応。Diff_any-bounded で OR=1.36 が出た同じ B' framework + cohort 構成下で、
  Static-bounded では effect が attenuate (OR ≈ 1) することを示すのが目標。
- v7 の DIFF_ANY_LABELS = ['Diff_any', 'Diff_specific_n1', 'Diff_shared_n2plus'] を
  STATIC_LABELS = ['Static'] に置き換え。それ以外のロジック (anchor 間 interval 構築、
  flanking domain mapping) は v7 と同一。

入力:
- /lustre12/home/kushima-pg/tad04212026/02_bin_l2_annotation/output_v2/bin_l2_annotation_v2.tsv.gz

出力:
- output_v1/static_bounded_domains_v1.bed
- output_v1/bin_to_static_bounded_domains_v1.tsv

実行時間記録あり

予想される構造:
- Static anchor 数: 62,150 (v6 round で確認済み)
- Static-bounded domain: ~62,000 (Static anchor は密に並ぶため domain mean size は小さい)
- これは "fragmented" に見えるが、reviewer M3 への対応として「Diff_any のみで signal が出る」
  ことを示すのが目的なので、Static-bounded で signal が attenuate するか null になれば成功。
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

# v1 (Static-bounded null control): Static anchor のみ
STATIC_LABELS = ['Static']
ANCHOR_LABELS = STATIC_LABELS

DIFF_ANY_RAW_LABELS = ['Diff_any', 'Diff_specific_n1', 'Diff_shared_n2plus']


def collapse_label(label):
    """For output: Static は Static のまま、Diff_any 系は 'Diff_any' に collapse (但し
    Static-only domain 構築では Static anchor のみ使うため、output bed の anchor は
    すべて 'Static')."""
    if label in STATIC_LABELS:
        return 'Static'
    if label in DIFF_ANY_RAW_LABELS:
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
# Step 2: Define Static-only anchors (M3 null control)
# ============================================================
bin_df['is_boundary_anchor'] = bin_df['group_primary'].isin(ANCHOR_LABELS).astype(int)
bin_df['group_primary_collapsed'] = bin_df['group_primary'].apply(collapse_label)

n_anchor = bin_df['is_boundary_anchor'].sum()
print(f"  Static-only anchors: {n_anchor}")

# ============================================================
# Step 3: Reconstruct Static-bounded TAD-scale domains
# ============================================================
print(f"[{time.time()-t0:.1f}s] Reconstructing Static-bounded TAD-scale domains (M3 null control) ...")
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
            'left_group_primary': 'Static',
            'right_group_primary': 'Static',
            'domain_size_bp': domain_end - domain_start,
        })

dom_df = pd.DataFrame(domains)
dom_df['domain_id'] = ['SDOM_' + str(i+1).zfill(7) for i in range(len(dom_df))]
print(f"  Reconstructed Static-bounded domains: {len(dom_df)}")
print(f"  Mean domain size: {dom_df['domain_size_bp'].mean()/1000:.1f} kb")
print(f"  Median domain size: {dom_df['domain_size_bp'].median()/1000:.1f} kb")
size_kb = dom_df['domain_size_bp'] / 1000
print(f"  Domain size distribution (kb):")
print(f"    min: {size_kb.min():.1f}, P25: {size_kb.quantile(0.25):.1f}, "
      f"P75: {size_kb.quantile(0.75):.1f}, max: {size_kb.max():.1f}")

n_static_static = len(dom_df)
print(f"  QC: All domain flanks are Static-Static: {n_static_static}")

out_bed = OUT_DIR / 'static_bounded_domains_v1.bed'
dom_df_bed = dom_df[['chrom', 'start', 'end', 'domain_id', 'left_bin_id',
                      'right_bin_id', 'left_group_primary', 'right_group_primary',
                      'domain_size_bp']]
dom_df_bed.to_csv(out_bed, sep='\t', index=False)
print(f"  Saved {out_bed}")

# ============================================================
# Step 4: Map each Static boundary bin to its flanking Static-bounded domains
# ============================================================
print(f"[{time.time()-t0:.1f}s] Mapping Static bins to flanking Static-bounded domains ...")
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
        'group_primary': 'Static',
        'left_domain_id': left_dom_id if left_dom_id else '',
        'right_domain_id': right_dom_id if right_dom_id else '',
        'has_one_flank_only': int(not (left_dom_id and right_dom_id)),
    })

bin_dom_df = pd.DataFrame(bin_to_dom)
n_one_flank = bin_dom_df['has_one_flank_only'].sum()
print(f"  Mapped Static bins: {len(bin_dom_df)}")
print(f"  One-flank-only bins (chromosome ends): {n_one_flank}")

out_bin_dom = OUT_DIR / 'bin_to_static_bounded_domains_v1.tsv'
bin_dom_df.to_csv(out_bin_dom, sep='\t', index=False)
print(f"  Saved {out_bin_dom}")

print(f"\n[Done in {time.time()-t0:.1f}s]")
