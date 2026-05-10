#!/usr/bin/env python3
"""
14_compute_per_sample_burden_static_v1.py

処理内容 (v231 round, M3 対応 — Static-bounded null control):
- v9 (Diff_any-bounded) と完全に同じ B' framework / SV interval / exclusion ロジックを
  Static-bounded domain (M3 反論用 null control) に適用。
- 違いは:
   * BIN_DOM = bin_to_static_bounded_domains_v1.tsv (Static anchor のみ)
   * GENES_PER_DOM = genes_per_static_domain_v1.tsv.gz
   * group_primary フィルター: 'Static' のみ採用 (v9 では 'Diff_any')
   * 5 burden variant 名: static_bounded, prox250/500/1000kb, static_bounded_excl
   * 出力ファイル名: *_static_v1.tsv

入力:
- output_v1/bin_to_static_bounded_domains_v1.tsv (NEW v231)
- output_v1/genes_per_static_domain_v1.tsv.gz (NEW v231)
- WGS event-bin: 04_wgs_sv_boundary_overlap/output_v10/sample_boundary_event_overlap_v10.tsv.gz
  (列: sample_id, ..., sv_chr, sv_start0, sv_end, ..., bin_id, ...)
- WGS sample cov: 05_wgs_sample_burden/output_v3/sample_burden_L2_and_specificity_v3.tsv
- arrayCGH event-bin: 08_arraycgh_sample_burden/output_v22/sample_event_bin_overlap_v22.tsv.gz
- arrayCGH cov: 08_arraycgh_sample_burden/output_v22/sample_covariates_v22.tsv
- MSSNG event-bin: output_v1/mssng_event_bins_dumped_v1.tsv.gz
- MSSNG sample cov: output_v1/mssng_sample_covariates_dumped_v1.tsv.gz

出力 (3 cohort files; 各 file に 5 burden variant の N_* 列):
- output_v1/sample_constraint_burden_discovery_static_v1.tsv
- output_v1/sample_constraint_burden_arraycgh_static_v1.tsv
- output_v1/sample_constraint_burden_mssng_static_v1.tsv

実行時間記録あり

期待される結果 (M3 反論用):
- Static-bounded domain は Diff_any anchor を通らないので、
  effect が attenuate (OR ≈ 1, P > 0.1) すれば「Diff_any anchor の生物学的選択性」が
  enrichment signal の本質であることを示せる。
- 逆に Static-bounded で同等の OR が出てしまうと circular reasoning と判定される。
"""

import re
import time
import gzip
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
from bisect import bisect_left, bisect_right

t0 = time.time()

BASE_DIR = Path('/lustre12/home/kushima-pg/tad04212026')
OUT_DIR = BASE_DIR / '14_constraint_enrichment_v1/output_v1'
ANNOT_DIR = Path('/lustre12/home/kushima-pg/annotationInfo')
GENCODE_GTF = ANNOT_DIR / 'gencode.v46.annotation.gtf.gz'
GNOMAD_LOEUF = ANNOT_DIR / 'gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz'

BIN_DOM = OUT_DIR / 'bin_to_static_bounded_domains_v1.tsv'
GENES_PER_DOM = OUT_DIR / 'genes_per_static_domain_v1.tsv.gz'

WGS_EVENT_BIN = BASE_DIR / '04_wgs_sv_boundary_overlap/output_v10/sample_boundary_event_overlap_v10.tsv.gz'
WGS_COV = BASE_DIR / '05_wgs_sample_burden/output_v3/sample_burden_L2_and_specificity_v3.tsv'

ARRAYCGH_EVENT_BIN = BASE_DIR / '08_arraycgh_sample_burden/output_v22/sample_event_bin_overlap_v22.tsv.gz'
ARRAYCGH_COV = BASE_DIR / '08_arraycgh_sample_burden/output_v22/sample_covariates_v22.tsv'

MSSNG_EVENT_BIN = OUT_DIR / 'mssng_event_bins_dumped_v1.tsv.gz'
MSSNG_COV = OUT_DIR / 'mssng_sample_covariates_dumped_v1.tsv.gz'

COVERAGE_FRAC = 0.10

LOEUF_STRICT = 0.35
LOEUF_RELAXED = 0.6
PLI_CUTOFF = 0.9

NEIGHBORHOOD_KB = [250, 500, 1000]  # ±X kb

# v231 M3: anchor filter ('Static' のみ。Diff_any は除外)
ANCHOR_GROUP_PRIMARY = 'Static'

WGS_REQUIRED_COV_COLS = {'sample_id', 'Diagnosis'}
ACGH_REQUIRED_COV_COLS = {'sample_id', 'diagnosis'}
MSSNG_REQUIRED_COV_COLS = {'sample_id', 'Diagnosis', 'FAMILYID',
                            'Sex', 'log1p_total_del_bases', 'log1p_total_gene_DEL'}
MSSNG_HARD_FAIL_NA_COLS = ['Sex', 'log1p_total_del_bases', 'log1p_total_gene_DEL']


def safe_strip_ab_suffix(sid):
    s = str(sid)
    return s[:-2] if s.endswith('AB') else s


def validate_required_cols(df, required, name):
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(
            f"[{name}] missing required columns: {sorted(missing)}\n"
            f"  Available: {sorted(df.columns)}"
        )


def sanitize_label(s):
    s = str(s)
    s = re.sub(r'[^0-9A-Za-z_]', '_', s)
    if s and s[0].isdigit():
        s = '_' + s
    s = re.sub(r'_+', '_', s).strip('_')
    return s if s else 'missing'


# ============================================================
# Step 1: Load Static-bounded domain mapping (v231 M3)
# ============================================================
print(f"[{time.time()-t0:.1f}s] Loading bin-to-Static-bounded-domain mapping (static_v1) ...")
bin_dom = pd.read_csv(BIN_DOM, sep='\t')
genes_per_dom = pd.read_csv(GENES_PER_DOM, sep='\t', compression='gzip')

# domain_id -> [(gene_id, gene_symbol, is_constraint_known, ...)]
DOM_GENES = {}
for did, g in genes_per_dom.groupby('domain_id'):
    DOM_GENES[did] = list(zip(
        g['gene_id'].tolist(),
        g['gene_symbol'].tolist(),
        g['is_constraint_known'].astype(int).tolist(),
        g['is_constrained_strict'].astype(int).tolist(),
        g['is_unconstrained_strict_known'].astype(int).tolist(),
        g['is_constrained_relaxed'].astype(int).tolist(),
        g['is_unconstrained_relaxed_known'].astype(int).tolist(),
        g['is_pLI_known'].astype(int).tolist(),
        g['is_pLI_constrained'].astype(int).tolist(),
        g['is_pLI_unconstrained_known'].astype(int).tolist(),
    ))

# bin_id -> (group_primary, chrom, bin_start, bin_end, left_dom, right_dom)
BIN_INFO = {}
for _, b in bin_dom.iterrows():
    BIN_INFO[b['bin_id']] = (
        b['group_primary'], b['chrom'], int(b['start0']), int(b['end']),
        b['left_domain_id'] if isinstance(b['left_domain_id'], str) else '',
        b['right_domain_id'] if isinstance(b['right_domain_id'], str) else '',
    )


# ============================================================
# Step 2: Load gene table for ±X kb neighborhood + direct-overlap exclusion
# ============================================================
print(f"[{time.time()-t0:.1f}s] Loading gnomAD constraint ...")
gnomad_df = pd.read_csv(GNOMAD_LOEUF, sep='\t', compression='gzip',
                         usecols=['gene', 'oe_lof_upper', 'pLI'])
gnomad_df = gnomad_df.dropna(subset=['gene'])
gnomad_df = gnomad_df.sort_values('oe_lof_upper').drop_duplicates(subset='gene', keep='first')
gnomad_df = gnomad_df.rename(columns={'oe_lof_upper': 'LOEUF'})

print(f"[{time.time()-t0:.1f}s] Parsing GENCODE v46 protein-coding genes ...")
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
            'gene_chr': chrom,
            'gene_start': start,
            'gene_end': end,
            'gene_symbol': attrs.get('gene_name', ''),
            'gene_id': attrs.get('gene_id', ''),
        })
gene_df = pd.DataFrame(genes)
gene_df = gene_df.merge(gnomad_df[['gene', 'LOEUF', 'pLI']],
                          left_on='gene_symbol', right_on='gene', how='left').drop(columns=['gene'])

# Constraint flags (consistent with annotate_static_v1)
gene_df['is_constraint_known'] = gene_df['LOEUF'].notna().astype(int)
gene_df['is_constrained_strict'] = ((gene_df['LOEUF'] < LOEUF_STRICT) & gene_df['is_constraint_known'].astype(bool)).fillna(False).astype(int)
gene_df['is_unconstrained_strict_known'] = ((gene_df['LOEUF'] >= LOEUF_STRICT) & gene_df['is_constraint_known'].astype(bool)).fillna(False).astype(int)
gene_df['is_constrained_relaxed'] = ((gene_df['LOEUF'] < LOEUF_RELAXED) & gene_df['is_constraint_known'].astype(bool)).fillna(False).astype(int)
gene_df['is_unconstrained_relaxed_known'] = ((gene_df['LOEUF'] >= LOEUF_RELAXED) & gene_df['is_constraint_known'].astype(bool)).fillna(False).astype(int)
pli_known = gene_df['pLI'].notna()
gene_df['is_pLI_known'] = pli_known.astype(int)
gene_df['is_pLI_constrained'] = ((gene_df['pLI'] > PLI_CUTOFF) & pli_known).fillna(False).astype(int)
gene_df['is_pLI_unconstrained_known'] = ((gene_df['pLI'] <= PLI_CUTOFF) & pli_known).fillna(False).astype(int)
gene_df['gene_mid'] = (gene_df['gene_start'] + gene_df['gene_end']) // 2

# Per-chrom sorted gene midpoints (for binary search ±X kb neighborhood)
GENE_BY_CHR = {}
for chrom, sub in gene_df.groupby('gene_chr'):
    sub = sub.sort_values('gene_mid').reset_index(drop=True)
    GENE_BY_CHR[chrom] = {
        'gene_mid': sub['gene_mid'].values,
        'records': sub[['gene_id', 'gene_symbol', 'is_constraint_known',
                         'is_constrained_strict', 'is_unconstrained_strict_known',
                         'is_constrained_relaxed', 'is_unconstrained_relaxed_known',
                         'is_pLI_known', 'is_pLI_constrained',
                         'is_pLI_unconstrained_known',
                         'gene_start', 'gene_end']].to_dict('records'),
    }
print(f"  Indexed {len(gene_df)} protein-coding genes by chromosome.")

# gene_id -> (chrom, gene_start, gene_end) lookup (for exact direct-overlap exclusion)
GENE_ID_TO_INTERVAL = {}
for _, g in gene_df.iterrows():
    GENE_ID_TO_INTERVAL[g['gene_id']] = (g['gene_chr'], int(g['gene_start']),
                                          int(g['gene_end']))


def parse_arraycgh_event_id(eid):
    """e.g., 'A134:chr2:106439107-107874624:DEL' -> ('chr2', 106439107, 107874624)."""
    if not isinstance(eid, str):
        return None, None, None
    parts = eid.split(':')
    if len(parts) < 4:
        return None, None, None
    chrom = parts[1]
    try:
        start_str, end_str = parts[2].split('-')
        return chrom, int(start_str), int(end_str)
    except (ValueError, IndexError):
        return None, None, None


def gene_count_record(gene_records):
    """Aggregate gene flags into per-sample record fields."""
    n_total = len(gene_records)
    n_known = sum(g['is_constraint_known'] for g in gene_records)
    n_unknown = n_total - n_known
    return {
        'N_total_genes': n_total,
        'N_constraint_known': n_known,
        'N_constraint_unknown': n_unknown,
        'N_constr_strict': sum(g['is_constrained_strict'] for g in gene_records),
        'N_unconstr_strict_known': sum(g['is_unconstrained_strict_known'] for g in gene_records),
        'N_constr_relaxed': sum(g['is_constrained_relaxed'] for g in gene_records),
        'N_unconstr_relaxed_known': sum(g['is_unconstrained_relaxed_known'] for g in gene_records),
        'N_pLI_known': sum(g['is_pLI_known'] for g in gene_records),
        'N_pLI_constr': sum(g['is_pLI_constrained'] for g in gene_records),
        'N_pLI_unconstr_known': sum(g['is_pLI_unconstrained_known'] for g in gene_records),
    }


def aggregate_burden_per_sample(df, sample_col, bin_col, sv_chr_col, sv_start_col, sv_end_col,
                                  label):
    """For each sample, compute 5 burden variants over Static disrupted bins (M3 null).

    Variants:
      (1) static_bounded:        Static-bounded domain (PRIMARY null control)
      (2) prox250kb / prox500kb / prox1000kb: ±X kb gene neighborhood
      (3) static_bounded_excl:   PRIMARY null - exact direct-overlap exclusion
    """
    print(f"[{time.time()-t0:.1f}s] Aggregating 5-variant burden for {label} "
          f"({len(df)} records) ...")
    if any(c is None for c in [sv_chr_col, sv_start_col, sv_end_col]):
        raise RuntimeError(
            f"static_v1 requires SV interval columns (sv_chr/sv_start/sv_end). "
            f"None passed for {label}."
        )

    rows = []
    for sample, sub in df.groupby(sample_col):
        # Set of disrupted bins + per-bin SV intervals (Static only; M3 null)
        disrupted_bin_ids = []
        sv_intervals_by_chrom = defaultdict(list)  # for exact direct-overlap exclusion
        seen_bin_set = set()
        for _, r in sub.iterrows():
            bid = r[bin_col]
            info = BIN_INFO.get(bid)
            if info is None:
                continue
            gp, b_chrom, b_start, b_end, ldid, rdid = info
            if gp != ANCHOR_GROUP_PRIMARY:  # 'Static'
                continue
            # collect bin (dedup)
            if bid not in seen_bin_set:
                seen_bin_set.add(bid)
                disrupted_bin_ids.append((bid, b_chrom, b_start, b_end, ldid, rdid))
            # collect SV interval (per-record; same (sample, bin) may have multi SV)
            try:
                sv_chr = r[sv_chr_col]
                sv_s = int(r[sv_start_col])
                sv_e = int(r[sv_end_col])
                if sv_chr and sv_s >= 0 and sv_e > sv_s:
                    sv_intervals_by_chrom[sv_chr].append((sv_s, sv_e))
            except (ValueError, TypeError, KeyError):
                continue

        # ----- Variant 1: static_bounded (PRIMARY null) -----
        seen_genes_dom = {}
        seen_dom_ids = set()
        for bid, b_chrom, b_start, b_end, ldid, rdid in disrupted_bin_ids:
            for did in (ldid, rdid):
                if did and did not in seen_dom_ids:
                    seen_dom_ids.add(did)
                    for entry in DOM_GENES.get(did, []):
                        gid = entry[0]
                        if gid not in seen_genes_dom:
                            seen_genes_dom[gid] = {
                                'gene_id': gid,
                                'gene_symbol': entry[1],
                                'is_constraint_known': entry[2],
                                'is_constrained_strict': entry[3],
                                'is_unconstrained_strict_known': entry[4],
                                'is_constrained_relaxed': entry[5],
                                'is_unconstrained_relaxed_known': entry[6],
                                'is_pLI_known': entry[7],
                                'is_pLI_constrained': entry[8],
                                'is_pLI_unconstrained_known': entry[9],
                            }

        # ----- Variant 5: static_bounded_excl (exact direct-overlap exclusion) -----
        # gene_start < sv_end AND gene_end > sv_start
        excluded_gene_ids = set()
        for gid in seen_genes_dom:
            interval = GENE_ID_TO_INTERVAL.get(gid)
            if interval is None:
                continue
            g_chrom, g_start, g_end = interval
            for sv_s, sv_e in sv_intervals_by_chrom.get(g_chrom, []):
                if g_start < sv_e and g_end > sv_s:
                    excluded_gene_ids.add(gid)
                    break

        seen_genes_dom_excl = {gid: rec for gid, rec in seen_genes_dom.items()
                                if gid not in excluded_gene_ids}

        # ----- Variants 2/3/4: prox250/500/1000 kb -----
        prox_genes = {kb: {} for kb in NEIGHBORHOOD_KB}
        for bid, b_chrom, b_start, b_end, ldid, rdid in disrupted_bin_ids:
            chrom_data = GENE_BY_CHR.get(b_chrom)
            if chrom_data is None:
                continue
            bin_mid = (b_start + b_end) / 2
            mids = chrom_data['gene_mid']
            recs = chrom_data['records']
            for kb in NEIGHBORHOOD_KB:
                X = kb * 1000
                lo = bisect_left(mids, bin_mid - X)
                hi = bisect_right(mids, bin_mid + X)
                for i in range(lo, hi):
                    g = recs[i]
                    gid = g['gene_id']
                    if gid not in prox_genes[kb]:
                        prox_genes[kb][gid] = g

        # Build record
        record = {'sample_id': sample}
        for k, v in gene_count_record(list(seen_genes_dom.values())).items():
            record[f'{k}_static_bounded'] = v
        for k, v in gene_count_record(list(seen_genes_dom_excl.values())).items():
            record[f'{k}_static_bounded_excl'] = v
        record['n_excluded_direct_overlap'] = len(excluded_gene_ids)
        for kb in NEIGHBORHOOD_KB:
            label_suf = f'prox{kb}kb'
            for k, v in gene_count_record(list(prox_genes[kb].values())).items():
                record[f'{k}_{label_suf}'] = v

        rows.append(record)

    out_df = pd.DataFrame(rows)
    for c in out_df.columns:
        if c.startswith('N_') or c == 'n_excluded_direct_overlap':
            out_df[c] = out_df[c].fillna(0).astype(int)
    print(f"  Samples processed: {len(out_df)}")
    return out_df


def left_merge_with_zero_fill(cov_df, exposure_df, name):
    out = cov_df.merge(exposure_df, on='sample_id', how='left', suffixes=('', '_exp'))
    n_cols = [c for c in out.columns if c.startswith('N_')]
    for c in n_cols:
        out[c] = out[c].fillna(0).astype(int)
    if 'n_excluded_direct_overlap' in out.columns:
        out['n_excluded_direct_overlap'] = out['n_excluded_direct_overlap'].fillna(0).astype(int)
    return out


# ============================================================
# Step 3: Discovery WGS
# ============================================================
print(f"\n[{time.time()-t0:.1f}s] === Discovery WGS (Static M3 null) ===")
wgs_df = pd.read_csv(WGS_EVENT_BIN, sep='\t', compression='gzip')
print(f"  WGS event-bin records: {len(wgs_df)}")
wgs_filt = wgs_df[(wgs_df['sv_type_norm'] == 'DEL') &
                   (wgs_df['overlap_frac_boundary'] >= COVERAGE_FRAC)].copy()
print(f"  DEL >={COVERAGE_FRAC} overlap: {len(wgs_filt)}")

exposure_wgs_raw = aggregate_burden_per_sample(
    wgs_filt, sample_col='sample_id', bin_col='bin_id',
    sv_chr_col='sv_chr', sv_start_col='sv_start0', sv_end_col='sv_end',
    label='discovery_wgs_static')

wgs_cov = pd.read_csv(WGS_COV, sep='\t')
validate_required_cols(wgs_cov, WGS_REQUIRED_COV_COLS, 'WGS_COV')
print(f"  WGS covariate samples: {len(wgs_cov)}")

miss_in_cov = set(exposure_wgs_raw['sample_id']) - set(wgs_cov['sample_id'])
if miss_in_cov:
    raise RuntimeError(
        f"WGS event-bin samples missing from cov: {len(miss_in_cov)} "
        f"examples={list(sorted(miss_in_cov))[:5]}"
    )

exposure_wgs = left_merge_with_zero_fill(wgs_cov, exposure_wgs_raw, 'WGS')
print(f"  WGS exposure rows: {len(exposure_wgs)}")
print(f"    ASD: {(exposure_wgs['Diagnosis']=='ASD').sum()}, "
      f"SZ: {(exposure_wgs['Diagnosis']=='SZ').sum()}, "
      f"Healthy: {(exposure_wgs['Diagnosis']=='Healthy').sum()}")

out_wgs = OUT_DIR / 'sample_constraint_burden_discovery_static_v1.tsv'
exposure_wgs.to_csv(out_wgs, sep='\t', index=False)
print(f"  Saved {out_wgs}")
# Sanity per variant
for var in ['static_bounded', 'prox250kb', 'prox500kb', 'prox1000kb', 'static_bounded_excl']:
    col = f'N_constr_strict_{var}'
    if col in exposure_wgs.columns:
        nz = (exposure_wgs[col] > 0).sum()
        print(f"  WGS non-zero {col}: {nz}")


# ============================================================
# Step 4: arrayCGH
# ============================================================
print(f"\n[{time.time()-t0:.1f}s] === arrayCGH (Static M3 null) ===")
acgh_df = pd.read_csv(ARRAYCGH_EVENT_BIN, sep='\t', compression='gzip')
acgh_filt = acgh_df[(acgh_df['sv_type'] == 'DEL') & (acgh_df['pattern'] == 'A')].copy()
print(f"  arrayCGH DEL Pattern A: {len(acgh_filt)}")

# parse event_id to extract SV interval (sv_chr, sv_start, sv_end)
parsed = acgh_filt['event_id'].apply(parse_arraycgh_event_id)
acgh_filt['sv_chr'] = parsed.apply(lambda x: x[0])
acgh_filt['sv_start'] = parsed.apply(lambda x: x[1])
acgh_filt['sv_end'] = parsed.apply(lambda x: x[2])
n_unparseable = acgh_filt['sv_chr'].isna().sum()
if n_unparseable > 0:
    print(f"  WARN: {n_unparseable} arrayCGH event_id failed to parse "
          f"(will lack SV interval, treated as no exclusion)")

exposure_acgh_raw = aggregate_burden_per_sample(
    acgh_filt, sample_col='sample_id', bin_col='bin_id',
    sv_chr_col='sv_chr', sv_start_col='sv_start', sv_end_col='sv_end',
    label='arraycgh_static')

acgh_cov = pd.read_csv(ARRAYCGH_COV, sep='\t')
validate_required_cols(acgh_cov, ACGH_REQUIRED_COV_COLS, 'ARRAYCGH_COV')

miss_in_cov_a = set(exposure_acgh_raw['sample_id']) - set(acgh_cov['sample_id'])
if miss_in_cov_a:
    raise RuntimeError(f"arrayCGH event-bin samples missing from cov: {len(miss_in_cov_a)}")

exposure_acgh = left_merge_with_zero_fill(acgh_cov, exposure_acgh_raw, 'arrayCGH')
if 'diagnosis' in exposure_acgh.columns and 'Diagnosis' not in exposure_acgh.columns:
    exposure_acgh = exposure_acgh.rename(columns={'diagnosis': 'Diagnosis'})

print(f"  arrayCGH exposure rows: {len(exposure_acgh)}")
print(f"    ASD: {(exposure_acgh['Diagnosis']=='ASD').sum()}, "
      f"CONT: {(exposure_acgh['Diagnosis']=='CONT').sum()}")

out_acgh = OUT_DIR / 'sample_constraint_burden_arraycgh_static_v1.tsv'
exposure_acgh.to_csv(out_acgh, sep='\t', index=False)
print(f"  Saved {out_acgh}")
for var in ['static_bounded', 'prox250kb', 'prox500kb', 'prox1000kb', 'static_bounded_excl']:
    col = f'N_constr_strict_{var}'
    if col in exposure_acgh.columns:
        nz = (exposure_acgh[col] > 0).sum()
        print(f"  arrayCGH non-zero {col}: {nz}")


# ============================================================
# Step 5: MSSNG
# ============================================================
print(f"\n[{time.time()-t0:.1f}s] === MSSNG (Static M3 null) ===")
if not MSSNG_EVENT_BIN.exists():
    raise RuntimeError(f"MSSNG event-bin not found: {MSSNG_EVENT_BIN}")
if not MSSNG_COV.exists():
    raise RuntimeError(f"MSSNG cov not found: {MSSNG_COV}")

mssng_df = pd.read_csv(MSSNG_EVENT_BIN, sep='\t', compression='gzip')
print(f"  MSSNG event-bin records: {len(mssng_df)}")
mssng_filt = mssng_df[(mssng_df['sv_type'] == 'DEL') & (mssng_df['pattern'] == 'A')].copy()

# require SV interval columns from v19 patch v4
required_mssng_cols = {'sv_chr', 'sv_start', 'sv_end'}
missing_sv_cols = required_mssng_cols - set(mssng_filt.columns)
if missing_sv_cols:
    raise RuntimeError(
        f"MSSNG event-bin missing SV interval columns: {sorted(missing_sv_cols)}. "
        f"Run v19 with patch v4 (PATCH_v18_to_v19_v4.md) to regenerate dump.\n"
        f"Available columns: {list(mssng_filt.columns)}"
    )

exposure_mssng_raw = aggregate_burden_per_sample(
    mssng_filt, sample_col='sample_id', bin_col='bin_id',
    sv_chr_col='sv_chr', sv_start_col='sv_start', sv_end_col='sv_end',
    label='mssng_static')

mssng_cov = pd.read_csv(MSSNG_COV, sep='\t', compression='gzip')
validate_required_cols(mssng_cov, MSSNG_REQUIRED_COV_COLS, 'MSSNG_COV')
print(f"  MSSNG covariate dump samples: {len(mssng_cov)}")

miss_in_cov_m = set(exposure_mssng_raw['sample_id']) - set(mssng_cov['sample_id'])
if miss_in_cov_m:
    raise RuntimeError(f"MSSNG event-bin samples missing from cov: {len(miss_in_cov_m)}")

# Platform / ancestry を sanitize 列名で dummy 化
dummy_added_cols = []
for cat_col, prefix in [('Platform', 'platform'), ('ancestry', 'ancestry')]:
    if cat_col in mssng_cov.columns:
        cat_series = mssng_cov[cat_col].fillna('missing').astype(str).replace('', 'missing')
        cat_series_clean = cat_series.apply(sanitize_label)
        d = pd.get_dummies(cat_series_clean, prefix=prefix, drop_first=True, dtype=int)
        d.columns = [sanitize_label(c) for c in d.columns]
        keep = [c for c in d.columns if d[c].sum() > 0]
        d = d[keep]
        mssng_cov = pd.concat([mssng_cov.drop(columns=[cat_col]), d], axis=1)
        dummy_added_cols.extend(d.columns.tolist())
        print(f"    {cat_col} → dummies: {keep}")
print(f"  Total dummies created: {dummy_added_cols}")

exposure_mssng = left_merge_with_zero_fill(mssng_cov, exposure_mssng_raw, 'MSSNG')

# Hard fail on required cov NaN
for c in MSSNG_HARD_FAIL_NA_COLS:
    n_na = exposure_mssng[c].isna().sum()
    print(f"    {c}: NaN count = {n_na} / {len(exposure_mssng)}")
    if n_na > 0:
        raise RuntimeError(f"MSSNG required cov '{c}' has {n_na} NaN values")

n_na_fid = exposure_mssng['FAMILYID'].isna().sum()
if n_na_fid > 0:
    exposure_mssng.loc[exposure_mssng['FAMILYID'].isna(), 'FAMILYID'] = (
        exposure_mssng.loc[exposure_mssng['FAMILYID'].isna(), 'sample_id']
        .apply(lambda s: f'_solo_{s}')
    )
    print(f"    FAMILYID: {n_na_fid} NaN filled with _solo_<sample_id>")

print(f"  MSSNG exposure rows: {len(exposure_mssng)}")
print(f"    ASD: {(exposure_mssng['Diagnosis']=='ASD').sum()}, "
      f"Sibling: {(exposure_mssng['Diagnosis']=='Sibling').sum()}")

out_mssng = OUT_DIR / 'sample_constraint_burden_mssng_static_v1.tsv'
exposure_mssng.to_csv(out_mssng, sep='\t', index=False)
print(f"  Saved {out_mssng}")
for var in ['static_bounded', 'prox250kb', 'prox500kb', 'prox1000kb', 'static_bounded_excl']:
    col = f'N_constr_strict_{var}'
    if col in exposure_mssng.columns:
        nz = (exposure_mssng[col] > 0).sum()
        print(f"  MSSNG non-zero {col}: {nz}")

print(f"\n[Done in {time.time()-t0:.1f}s]")
