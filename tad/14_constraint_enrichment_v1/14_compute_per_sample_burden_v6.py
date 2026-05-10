#!/usr/bin/env python3
"""
14_compute_per_sample_burden_v6.py

処理内容 (v5 → v6 変更点; v5 round 実行で判明した schema mismatch fix):
- BIN_DOM input path を bin_to_flanking_domains_v3.tsv → _v4.tsv に変更
  (v4 round で TAD anchor を Diff_specific_n1 + Diff_shared_n2plus + Static
  に修正し、group_primary 列を Diff_any/Static に二値化済み)
- 出力ファイル名 suffix を _v5 → _v6 に変更
- 既存修正の継承 (v5 から):
  - MSSNG required covariate NaN hard fail (Priority 1-2)
  - dummy column 名 sanitize (Priority 1-3 prep)
  - anti-join check (event-bin samples ⊆ cov samples)
  - WGS/arrayCGH/MSSNG 全 sample に covariate を left merge
  - MSSNG Platform/ancestry を pd.get_dummies で dummy 化
  - LOEUF unknown を unconstrained に混入させない
  - safe_strip_ab_suffix(), dedup by gene_id

入力:
- output_v1/bin_to_flanking_domains_v4.tsv (v4, schema 修正後)
- output_v1/genes_per_domain_v3.tsv.gz (v3, LOEUF 分離 logic は v2 不変、入力 TAD bed のみ v4 化)
- WGS event-bin: 04_wgs_sv_boundary_overlap/output_v10/sample_boundary_event_overlap_v10.tsv.gz
- WGS sample cov: 05_wgs_sample_burden/output_v3/sample_burden_L2_and_specificity_v3.tsv
- arrayCGH event-bin + cov: 08_arraycgh_sample_burden/output_v22/...
- MSSNG event-bin: output_v1/mssng_event_bins_dumped_v1.tsv.gz (v19 patch v3 出力)
- MSSNG sample cov: output_v1/mssng_sample_covariates_dumped_v1.tsv.gz (v19 patch v3 出力)

出力:
- output_v1/sample_constrained_burden_discovery_v6.tsv
- output_v1/sample_constrained_burden_arraycgh_v6.tsv
- output_v1/sample_constrained_burden_mssng_v6.tsv (Platform/ancestry sanitized dummies 含む)

実行時間記録あり
"""

import re
import time
import sys
import pandas as pd
import numpy as np
from pathlib import Path

t0 = time.time()

BASE_DIR = Path('/lustre12/home/kushima-pg/tad04212026')
OUT_DIR = BASE_DIR / '14_constraint_enrichment_v1/output_v1'

BIN_DOM = OUT_DIR / 'bin_to_flanking_domains_v4.tsv'
GENES_PER_DOM = OUT_DIR / 'genes_per_domain_v3.tsv.gz'

WGS_EVENT_BIN = BASE_DIR / '04_wgs_sv_boundary_overlap/output_v10/sample_boundary_event_overlap_v10.tsv.gz'
WGS_COV = BASE_DIR / '05_wgs_sample_burden/output_v3/sample_burden_L2_and_specificity_v3.tsv'

ARRAYCGH_EVENT_BIN = BASE_DIR / '08_arraycgh_sample_burden/output_v22/sample_event_bin_overlap_v22.tsv.gz'
ARRAYCGH_COV = BASE_DIR / '08_arraycgh_sample_burden/output_v22/sample_covariates_v22.tsv'

MSSNG_EVENT_BIN = OUT_DIR / 'mssng_event_bins_dumped_v1.tsv.gz'
MSSNG_COV = OUT_DIR / 'mssng_sample_covariates_dumped_v1.tsv.gz'

COVERAGE_FRAC = 0.10

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
    """v5 Priority 1-3 prep: dummy column 名を formula-safe にする.
    英数字・アンダースコア以外を _ に置換し、先頭が数字なら _ を prefix.
    例: '1.0' → '_1_0', 'GSA-MD' → 'GSA_MD'.
    """
    s = str(s)
    s = re.sub(r'[^0-9A-Za-z_]', '_', s)
    if s and s[0].isdigit():
        s = '_' + s
    s = re.sub(r'_+', '_', s).strip('_')
    return s if s else 'missing'


# ============================================================
# Load TAD domain mapping
# ============================================================
print(f"[{time.time()-t0:.1f}s] Loading bin-to-domain mapping ...")
bin_dom = pd.read_csv(BIN_DOM, sep='\t')
genes_per_dom = pd.read_csv(GENES_PER_DOM, sep='\t', compression='gzip')

dom_genes = {}
for did, g in genes_per_dom.groupby('domain_id'):
    dom_genes[did] = list(zip(
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

bin_dom_map = dict(zip(
    bin_dom['bin_id'],
    zip(bin_dom['group_primary'], bin_dom['left_domain_id'], bin_dom['right_domain_id']),
))


def aggregate_exposure(df, sample_col, bin_col, label):
    print(f"[{time.time()-t0:.1f}s] Aggregating exposure for {label} ({len(df)} records) ...")
    sample_domain_set = {}
    for sample, sub in df.groupby(sample_col):
        gp_doms = {'Diff_any': set(), 'Static': set()}
        for bin_id in sub[bin_col].unique():
            mapping = bin_dom_map.get(bin_id)
            if mapping is None:
                continue
            group_primary, left_dom, right_dom = mapping
            if group_primary not in gp_doms:
                continue
            if left_dom and isinstance(left_dom, str):
                gp_doms[group_primary].add(left_dom)
            if right_dom and isinstance(right_dom, str):
                gp_doms[group_primary].add(right_dom)
        sample_domain_set[sample] = gp_doms

    rows = []
    for sample, gp_doms in sample_domain_set.items():
        record = {'sample_id': sample}
        for gp_label, doms in gp_doms.items():
            seen = {}
            for did in doms:
                for entry in dom_genes.get(did, []):
                    gid = entry[0]
                    if gid not in seen:
                        seen[gid] = entry
            n_total = len(seen)
            n_known = sum(e[2] for e in seen.values())
            n_unknown = n_total - n_known
            n_constr_strict = sum(e[3] for e in seen.values())
            n_unconstr_strict_known = sum(e[4] for e in seen.values())
            n_constr_relaxed = sum(e[5] for e in seen.values())
            n_unconstr_relaxed_known = sum(e[6] for e in seen.values())
            n_pli_known = sum(e[7] for e in seen.values())
            n_pli_constr = sum(e[8] for e in seen.values())
            n_pli_unconstr_known = sum(e[9] for e in seen.values())
            tag = gp_label.lower()
            record[f'N_total_genes_{tag}'] = n_total
            record[f'N_constraint_known_{tag}'] = n_known
            record[f'N_constraint_unknown_{tag}'] = n_unknown
            record[f'N_constr_strict_{tag}'] = n_constr_strict
            record[f'N_unconstr_strict_known_{tag}'] = n_unconstr_strict_known
            record[f'N_constr_relaxed_{tag}'] = n_constr_relaxed
            record[f'N_unconstr_relaxed_known_{tag}'] = n_unconstr_relaxed_known
            record[f'N_pLI_known_{tag}'] = n_pli_known
            record[f'N_pLI_constr_{tag}'] = n_pli_constr
            record[f'N_pLI_unconstr_known_{tag}'] = n_pli_unconstr_known
        rows.append(record)
    out_df = pd.DataFrame(rows)
    for col in out_df.columns:
        if col.startswith('N_'):
            out_df[col] = out_df[col].fillna(0).astype(int)
    print(f"  Samples processed: {len(out_df)}")
    return out_df


def left_merge_with_zero_fill(cov_df, exposure_df, name):
    """Left merge cov base with exposure; fill missing N_* with 0."""
    out = cov_df.merge(exposure_df, on='sample_id', how='left', suffixes=('', '_exp'))
    n_cols = [c for c in out.columns
              if c.startswith('N_') and (c.endswith('_diff_any') or c.endswith('_static'))]
    for c in n_cols:
        out[c] = out[c].fillna(0).astype(int)
    return out


# ============================================================
# Step 1: Discovery WGS
# ============================================================
print(f"\n[{time.time()-t0:.1f}s] === Discovery WGS (Pattern A) ===")
wgs_df = pd.read_csv(WGS_EVENT_BIN, sep='\t', compression='gzip')
print(f"  WGS event-bin records: {len(wgs_df)}")
wgs_filt = wgs_df[(wgs_df['sv_type_norm'] == 'DEL') &
                   (wgs_df['overlap_frac_boundary'] >= COVERAGE_FRAC)].copy()
print(f"  DEL >={COVERAGE_FRAC} overlap: {len(wgs_filt)}")

exposure_wgs_raw = aggregate_exposure(wgs_filt, 'sample_id', 'bin_id', 'discovery_wgs')

wgs_cov = pd.read_csv(WGS_COV, sep='\t')
validate_required_cols(wgs_cov, WGS_REQUIRED_COV_COLS, 'WGS_COV')
print(f"  WGS covariate samples: {len(wgs_cov)}")

# v5 Priority 1: anti-join check (event-bin samples must be subset of cov samples)
miss_in_cov = set(exposure_wgs_raw['sample_id']) - set(wgs_cov['sample_id'])
if miss_in_cov:
    raise RuntimeError(
        f"WGS event-bin samples missing from cov denominator: {len(miss_in_cov)} "
        f"examples={list(sorted(miss_in_cov))[:5]}"
    )

exposure_wgs = left_merge_with_zero_fill(wgs_cov, exposure_wgs_raw, 'WGS')
print(f"  WGS exposure rows after left merge: {len(exposure_wgs)}")
print(f"    ASD: {(exposure_wgs['Diagnosis']=='ASD').sum()}, "
      f"SZ: {(exposure_wgs['Diagnosis']=='SZ').sum()}, "
      f"Healthy: {(exposure_wgs['Diagnosis']=='Healthy').sum()}")

out_wgs = OUT_DIR / 'sample_constrained_burden_discovery_v6.tsv'
exposure_wgs.to_csv(out_wgs, sep='\t', index=False)
print(f"  Saved {out_wgs}")


# ============================================================
# Step 2: arrayCGH
# ============================================================
print(f"\n[{time.time()-t0:.1f}s] === arrayCGH (Pattern A) ===")
acgh_df = pd.read_csv(ARRAYCGH_EVENT_BIN, sep='\t', compression='gzip')
acgh_filt = acgh_df[(acgh_df['sv_type'] == 'DEL') & (acgh_df['pattern'] == 'A')].copy()
print(f"  arrayCGH DEL Pattern A: {len(acgh_filt)}")

exposure_acgh_raw = aggregate_exposure(acgh_filt, 'sample_id', 'bin_id', 'arraycgh')

acgh_cov = pd.read_csv(ARRAYCGH_COV, sep='\t')
validate_required_cols(acgh_cov, ACGH_REQUIRED_COV_COLS, 'ARRAYCGH_COV')

# v5 Priority 1: anti-join check
miss_in_cov_a = set(exposure_acgh_raw['sample_id']) - set(acgh_cov['sample_id'])
if miss_in_cov_a:
    raise RuntimeError(
        f"arrayCGH event-bin samples missing from cov denominator: {len(miss_in_cov_a)}"
    )

exposure_acgh = left_merge_with_zero_fill(acgh_cov, exposure_acgh_raw, 'arrayCGH')
if 'diagnosis' in exposure_acgh.columns and 'Diagnosis' not in exposure_acgh.columns:
    exposure_acgh = exposure_acgh.rename(columns={'diagnosis': 'Diagnosis'})

print(f"  arrayCGH exposure rows: {len(exposure_acgh)}")
print(f"    ASD: {(exposure_acgh['Diagnosis']=='ASD').sum()}, "
      f"CONT: {(exposure_acgh['Diagnosis']=='CONT').sum()}")

out_acgh = OUT_DIR / 'sample_constrained_burden_arraycgh_v6.tsv'
exposure_acgh.to_csv(out_acgh, sep='\t', index=False)
print(f"  Saved {out_acgh}")


# ============================================================
# Step 3: MSSNG (v5 hard fail on missing required cov)
# ============================================================
print(f"\n[{time.time()-t0:.1f}s] === MSSNG (v5 hard fail strict) ===")
if not MSSNG_EVENT_BIN.exists():
    raise RuntimeError(f"MSSNG event-bin not found: {MSSNG_EVENT_BIN}")
if not MSSNG_COV.exists():
    raise RuntimeError(f"MSSNG sample covariate dump not found: {MSSNG_COV}")

mssng_df = pd.read_csv(MSSNG_EVENT_BIN, sep='\t', compression='gzip')
print(f"  MSSNG event-bin records: {len(mssng_df)}")
mssng_filt = mssng_df[(mssng_df['sv_type'] == 'DEL') & (mssng_df['pattern'] == 'A')].copy()
exposure_mssng_raw = aggregate_exposure(mssng_filt, 'sample_id', 'bin_id', 'mssng')

mssng_cov = pd.read_csv(MSSNG_COV, sep='\t', compression='gzip')
validate_required_cols(mssng_cov, MSSNG_REQUIRED_COV_COLS, 'MSSNG_COV')
print(f"  MSSNG covariate dump samples: {len(mssng_cov)}")

# v5 Priority 1: anti-join check
miss_in_cov_m = set(exposure_mssng_raw['sample_id']) - set(mssng_cov['sample_id'])
if miss_in_cov_m:
    raise RuntimeError(
        f"MSSNG event-bin samples missing from cov denominator: {len(miss_in_cov_m)} "
        f"examples={list(sorted(miss_in_cov_m))[:5]}\n"
        f"  Both files should originate from the same v19 patch run; mismatch indicates "
        f"a denominator integrity bug."
    )

# v5 Block 2: Platform / ancestry を sanitize 列名で dummy 化
dummy_added_cols = []
for cat_col, prefix in [('Platform', 'platform'), ('ancestry', 'ancestry')]:
    if cat_col in mssng_cov.columns:
        cat_series = mssng_cov[cat_col].fillna('missing').astype(str).replace('', 'missing')
        # v5 Priority 1-3: sanitize each unique value first
        cat_series_clean = cat_series.apply(sanitize_label)
        d = pd.get_dummies(cat_series_clean, prefix=prefix, drop_first=True, dtype=int)
        # Strip remaining unsafe chars from column names just in case
        d.columns = [sanitize_label(c) for c in d.columns]
        keep = [c for c in d.columns if d[c].sum() > 0]
        d = d[keep]
        mssng_cov = pd.concat([mssng_cov.drop(columns=[cat_col]), d], axis=1)
        dummy_added_cols.extend(d.columns.tolist())
        print(f"    {cat_col} → dummies (drop_first, sum>0 retained, sanitized): {keep}")
print(f"  Total dummies created: {dummy_added_cols}")

# v5 Block 1 (continued): cov base に left merge
exposure_mssng = left_merge_with_zero_fill(mssng_cov, exposure_mssng_raw, 'MSSNG')

# v5 Priority 1-2: required covariate hard fail on NaN (Sex, log1p_total_del_bases,
# log1p_total_gene_DEL). FAMILYID は _solo_ fallback OK.
for c in MSSNG_HARD_FAIL_NA_COLS:
    n_na = exposure_mssng[c].isna().sum()
    print(f"    {c}: NaN count = {n_na} / {len(exposure_mssng)}")
    if n_na > 0:
        bad = exposure_mssng.loc[exposure_mssng[c].isna(), 'sample_id'].head(10).tolist()
        raise RuntimeError(
            f"MSSNG required covariate '{c}' has {n_na} NaN values. "
            f"Examples: {bad}\n"
            f"  Required-non-null cols: {MSSNG_HARD_FAIL_NA_COLS}\n"
            f"  Re-check v19 patch v3 Sex / DEL-burden mapping logic."
        )

# FAMILYID NaN は _solo_ で埋める (cov dump 側で _solo_ になっているはずだが念のため)
n_na_fid = exposure_mssng['FAMILYID'].isna().sum()
if n_na_fid > 0:
    exposure_mssng.loc[exposure_mssng['FAMILYID'].isna(), 'FAMILYID'] = (
        exposure_mssng.loc[exposure_mssng['FAMILYID'].isna(), 'sample_id']
        .apply(lambda s: f'_solo_{s}')
    )
    print(f"    FAMILYID: {n_na_fid} NaN filled with _solo_<sample_id>")

print(f"  MSSNG exposure rows after left-merge with cov: {len(exposure_mssng)}")
print(f"    ASD: {(exposure_mssng['Diagnosis']=='ASD').sum()}, "
      f"Sibling: {(exposure_mssng['Diagnosis']=='Sibling').sum()}")

# Sanity: non-zero exposed sample に covariate が付いているか
ne_col = 'N_constr_strict_diff_any'
if ne_col in exposure_mssng.columns:
    exposed = exposure_mssng[exposure_mssng[ne_col] > 0]
    print(f"  Sanity check: non-zero {ne_col} rows = {len(exposed)}")
    if len(exposed) > 0:
        for c in MSSNG_HARD_FAIL_NA_COLS + ['FAMILYID']:
            n_na_e = exposed[c].isna().sum() if c in exposed.columns else 'N/A'
            print(f"    of which {c} NaN: {n_na_e}")

out_mssng = OUT_DIR / 'sample_constrained_burden_mssng_v6.tsv'
exposure_mssng.to_csv(out_mssng, sep='\t', index=False)
print(f"  Saved {out_mssng}")
print(f"  Columns ({len(exposure_mssng.columns)}): {list(exposure_mssng.columns)[:30]}{' ...' if len(exposure_mssng.columns) > 30 else ''}")

print(f"\n[Done in {time.time()-t0:.1f}s]")
