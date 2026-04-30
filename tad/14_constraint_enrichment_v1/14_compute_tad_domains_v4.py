#!/usr/bin/env python3
"""
14_compute_tad_domains_v4.py

処理内容 (v3 → v4 変更点; v5 round 実行で判明した schema mismatch fix):
- v3 の前提だった "group_primary in [Diff_any, Static]" は誤り。
  実際の bin_l2_annotation_v2 schema は 'Static', 'Diff_specific_n1',
  'Diff_shared_n2plus' の 3 値。Diff_any はラベルではなく
  Diff_specific_n1 + Diff_shared_n2plus の **集合概念**。
- 修正: anchor を ['Diff_specific_n1', 'Diff_shared_n2plus', 'Static'] に
  変更し、出力 bed/tsv の group_primary 列は downstream (burden_v6) との
  互換性のため 'Diff_any' / 'Static' に二値化して保存する。
- 既存修正の継承 (v3 から):
  - HEFFEL_MASTER / L2_CLASSES 削除
  - 隣接 anchor 間 interval を TAD domain 化
  - boundary bin に left/right flanking domain を mapping

入力:
- /lustre12/home/kushima-pg/tad04292026/02_bin_l2_annotation/output_v2/bin_l2_annotation_v2.tsv.gz

出力:
- output_v1/tad_domains_v4.bed
- output_v1/bin_to_flanking_domains_v4.tsv

実行時間記録あり
"""

import time
import pandas as pd
from collections import defaultdict
from pathlib import Path

t0 = time.time()

BASE_DIR = Path('/lustre12/home/kushima-pg/tad04292026')
BIN_ANNOT = BASE_DIR / '02_bin_l2_annotation/output_v2/bin_l2_annotation_v2.tsv.gz'
OUT_DIR = BASE_DIR / '14_constraint_enrichment_v1/output_v1'
OUT_DIR.mkdir(parents=True, exist_ok=True)

# v4: actual labels in bin_l2_annotation_v2.group_primary
DIFF_ANY_LABELS = ['Diff_specific_n1', 'Diff_shared_n2plus']
ANCHOR_LABELS = DIFF_ANY_LABELS + ['Static']


def collapse_to_diff_any(label):
    """Map Diff_specific_n1 / Diff_shared_n2plus -> Diff_any; Static -> Static."""
    if label in DIFF_ANY_LABELS:
        return 'Diff_any'
    return label  # 'Static' のまま、それ以外も pass-through


# ============================================================
# Step 1: Load bin annotation and validate required columns
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
# Step 2: Define boundary anchors (FIXED in v4)
# ============================================================
# v4: use actual labels (Diff_specific_n1, Diff_shared_n2plus, Static).
bin_df['is_boundary_anchor'] = bin_df['group_primary'].isin(ANCHOR_LABELS).astype(int)
bin_df['group_primary_collapsed'] = bin_df['group_primary'].apply(collapse_to_diff_any)

n_anchor = bin_df['is_boundary_anchor'].sum()
n_diff_any = (bin_df['group_primary_collapsed'] == 'Diff_any').sum()
n_static = (bin_df['group_primary_collapsed'] == 'Static').sum()
print(f"  Boundary anchors (Diff_specific_n1 + Diff_shared_n2plus + Static): {n_anchor}")
print(f"    Diff_any (collapsed from Diff_specific_n1 + Diff_shared_n2plus): {n_diff_any}")
print(f"      of which Diff_specific_n1: {(bin_df['group_primary'] == 'Diff_specific_n1').sum()}")
print(f"      of which Diff_shared_n2plus: {(bin_df['group_primary'] == 'Diff_shared_n2plus').sum()}")
print(f"    Static: {n_static}")

# ============================================================
# Step 3: Reconstruct TAD domains as intervals between adjacent anchors
# ============================================================
print(f"[{time.time()-t0:.1f}s] Reconstructing TAD domains ...")
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
            # v4: collapsed labels (Diff_any / Static) for downstream burden compatibility
            'left_group_primary': left_anchor['group_primary_collapsed'],
            'right_group_primary': right_anchor['group_primary_collapsed'],
            'domain_size_bp': domain_end - domain_start,
        })

dom_df = pd.DataFrame(domains)
dom_df['domain_id'] = ['DOM_' + str(i+1).zfill(7) for i in range(len(dom_df))]
print(f"  Reconstructed TAD domains: {len(dom_df)}")
print(f"  Mean domain size: {dom_df['domain_size_bp'].mean()/1000:.1f} kb")
print(f"  Median domain size: {dom_df['domain_size_bp'].median()/1000:.1f} kb")

# Bounded domain QC (one-flank-only count)
n_diff_diff = ((dom_df['left_group_primary'] == 'Diff_any') &
                (dom_df['right_group_primary'] == 'Diff_any')).sum()
n_static_static = ((dom_df['left_group_primary'] == 'Static') &
                    (dom_df['right_group_primary'] == 'Static')).sum()
n_mixed = len(dom_df) - n_diff_diff - n_static_static
print(f"  QC: Diff_any-Diff_any flanks: {n_diff_diff}, "
      f"Static-Static: {n_static_static}, mixed: {n_mixed}")

out_bed = OUT_DIR / 'tad_domains_v4.bed'
dom_df_bed = dom_df[['chrom', 'start', 'end', 'domain_id', 'left_bin_id',
                      'right_bin_id', 'left_group_primary', 'right_group_primary',
                      'domain_size_bp']]
dom_df_bed.to_csv(out_bed, sep='\t', index=False)
print(f"  Saved {out_bed}")

# ============================================================
# Step 4: Map each boundary bin to its flanking TAD domains
# ============================================================
print(f"[{time.time()-t0:.1f}s] Mapping bins to flanking domains ...")
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
        # v4: collapsed Diff_any / Static label for downstream burden script
        'group_primary': b['group_primary_collapsed'],
        'left_domain_id': left_dom_id if left_dom_id else '',
        'right_domain_id': right_dom_id if right_dom_id else '',
        'has_one_flank_only': int(not (left_dom_id and right_dom_id)),
    })

bin_dom_df = pd.DataFrame(bin_to_dom)
n_one_flank = bin_dom_df['has_one_flank_only'].sum()
print(f"  Mapped bins: {len(bin_dom_df)}")
print(f"  One-flank-only bins (chromosome ends): {n_one_flank}")
print(f"  Bins by group_primary: {bin_dom_df['group_primary'].value_counts().to_dict()}")

out_bin_dom = OUT_DIR / 'bin_to_flanking_domains_v4.tsv'
bin_dom_df.to_csv(out_bin_dom, sep='\t', index=False)
print(f"  Saved {out_bin_dom}")

print(f"\n[Done in {time.time()-t0:.1f}s]")
