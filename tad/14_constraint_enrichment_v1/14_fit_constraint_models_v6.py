#!/usr/bin/env python3
"""
14_fit_constraint_models_v6.py

処理内容 (v5 → v6 変更点):
- 入出力 path のみ v5 → v6 に変更 (v4 round で TAD anchor schema 修正済みの
  burden_v6 出力を読む)。
- 統計設計・require_fit_ok・GEE matrix API・arrayCGH stratified primary・
  contrast meta・Holm 補正等は v5 から不変。

v5 で導入された fix (継承; ChatGPT v4 review safety patch):
- Priority 1-1: primary model failure / meta skip を **hard fail**.
  require_fit_ok() で fit_status / non-finite beta/se / non-positive SE をチェック。
  Test 1 / Test 3 arrayCGH / Test 3 MSSNG が ok でなければ RuntimeError.
  Primary IVW meta が作れない場合も RuntimeError (silent skip 廃止).
- Priority 1-3: GEE を formula API (smf.gee) → **matrix API (sm.GEE)** に変更.
  これにより dummy column 名に `.`, `-`, `/`, 数値 prefix が含まれても安全。
  さらに analytic subset 後の zero-variance covariate を自動 drop (perfect collinearity 回避).
- Priority 1-4: arrayCGH **platform-stratified IVW を primary** に変更
  (manuscript main framework との整合性; ChatGPT 推奨).
  pooled platform-adjusted logistic は sensitivity として保持.
  Test 4 external IVW meta は arrayCGH platform-stratified IVW estimate + MSSNG GEE で構成.
- Priority 1-5: contrast (β_constr - β_unconstr) の **2-cohort external IVW meta** を
  supporting row として出力. arrayCGH stratified-IVW の contrast は Agilent + NimbleGen 内
  で IVW し、それと MSSNG contrast を IVW meta する.
- Priority 2: outcome filter を明示 (assert_cols), MSSNG GEE では Sex を numeric にキャスト.
- 既存修正の継承 (v4 から):
  - scipy.stats を使用 (sm.distributions ではない)
  - Logistic 側 contrast (β_constr - β_unconstr)

入力:
- output_v1/sample_constrained_burden_discovery_v6.tsv (covariate joined)
- output_v1/sample_constrained_burden_arraycgh_v6.tsv (covariate joined)
- output_v1/sample_constrained_burden_mssng_v6.tsv  (covariate joined, sanitized dummies)

出力:
- output_v1/constraint_enrichment_results_v6.tsv

実行時間記録あり
"""

import time
import math
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

t0 = time.time()

BASE_DIR = Path('/lustre12/home/kushima-pg/tad04212026')
OUT_DIR = BASE_DIR / '14_constraint_enrichment_v1/output_v1'

EXPOSURE_DISC = OUT_DIR / 'sample_constrained_burden_discovery_v6.tsv'
EXPOSURE_ACGH = OUT_DIR / 'sample_constrained_burden_arraycgh_v6.tsv'
EXPOSURE_MSSNG = OUT_DIR / 'sample_constrained_burden_mssng_v6.tsv'

WGS_COVS = ['Sex_numeric', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10',
             'log1p_total_del_bases', 'log1p_total_gene_DEL']
ACGH_COVS_POOLED = ['sex', 'platform_nimblegen', 'log1p_total_del_bases_A', 'log1p_total_gene_DEL_A']
ACGH_COVS_STRAT = ['sex', 'log1p_total_del_bases_A', 'log1p_total_gene_DEL_A']
MSSNG_REQUIRED_COVS = ['Sex', 'log1p_total_del_bases', 'log1p_total_gene_DEL', 'FAMILYID']


def assert_cols(df, required, name):
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(
            f"[{name}] required columns missing: {missing}\n"
            f"  Available: {list(df.columns)}"
        )


def require_fit_ok(res, name, allowed=('ok', 'ok_gee')):
    """v5 Priority 1-1: primary model failure を hard fail."""
    status = str(res.get('fit_status', ''))
    ok = any(status == a or status.startswith(a) for a in allowed)
    if not ok:
        raise RuntimeError(f"{name} fit failed: {status}; full result={res}")
    for k in ('beta', 'se'):
        v = res.get(k)
        if v is None or (isinstance(v, float) and not np.isfinite(v)):
            raise RuntimeError(f"{name} has non-finite {k}: {v}")
    if 'se' in res and res['se'] is not None and res['se'] <= 0:
        raise RuntimeError(f"{name} has non-positive SE: {res['se']}")


def drop_zero_variance(sub, cols):
    """v5 Priority 1-3: analytic subset 後の zero-variance covariate を drop."""
    keep = []
    dropped = []
    for c in cols:
        if c not in sub.columns:
            dropped.append((c, 'absent'))
            continue
        if sub[c].nunique(dropna=True) > 1:
            keep.append(c)
        else:
            dropped.append((c, f'zero_variance(unique={sub[c].nunique()})'))
    if dropped:
        print(f"    Dropped covariates due to absent/zero-variance: {dropped}")
    return keep


def fit_logistic(df, exposure_constr, exposure_unconstr,
                  extra_static_constr=None, extra_static_unconstr=None,
                  covariates=None, label='', outcome='ASD_vs_Control'):
    cov = covariates or []
    cols_needed = ['case', exposure_constr, exposure_unconstr] + cov
    if extra_static_constr:
        cols_needed += [extra_static_constr]
    if extra_static_unconstr:
        cols_needed += [extra_static_unconstr]
    sub = df[[c for c in cols_needed if c in df.columns]].dropna().copy()
    if sub.empty or len(sub) < 50:
        return {'label': label, 'n_complete': len(sub), 'fit_status': 'too_few_samples'}

    # v5 Priority 1-3: drop zero-variance covariates after analytic filtering
    cov_kept = drop_zero_variance(sub, cov)
    static_kept = []
    for sc in [extra_static_constr, extra_static_unconstr]:
        if sc and sc in sub.columns and sub[sc].nunique(dropna=True) > 1:
            static_kept.append(sc)

    used_cols = [exposure_constr, exposure_unconstr] + cov_kept + static_kept
    X = sm.add_constant(sub[used_cols].astype(float), has_constant='add')
    y = sub['case'].astype(int)
    try:
        res = sm.Logit(y, X).fit(disp=False, method='newton', maxiter=200)
    except Exception as e:
        return {'label': label, 'n_complete': len(sub), 'fit_status': f'fit_error:{e}'}

    if exposure_constr not in res.params.index:
        return {'label': label, 'n_complete': len(sub), 'fit_status': 'param_missing'}

    beta = float(res.params[exposure_constr])
    se = float(res.bse[exposure_constr])
    p_two = float(res.pvalues[exposure_constr])
    z = beta / se
    p_one = float(stats.norm.sf(z))

    out = {
        'label': label, 'outcome': outcome, 'exposure': exposure_constr,
        'n_complete': len(sub),
        'beta': beta, 'se': se,
        'OR': math.exp(beta),
        'OR_lo95': math.exp(beta - 1.96 * se),
        'OR_hi95': math.exp(beta + 1.96 * se),
        'p_two_sided': p_two, 'p_one_sided': p_one,
        'fit_status': 'ok',
    }

    if exposure_unconstr in res.params.index:
        beta_un = float(res.params[exposure_unconstr])
        cov_mat = res.cov_params()
        var_d = (cov_mat.loc[exposure_constr, exposure_constr]
                 + cov_mat.loc[exposure_unconstr, exposure_unconstr]
                 - 2 * cov_mat.loc[exposure_constr, exposure_unconstr])
        delta = beta - beta_un
        if var_d > 0:
            se_d = math.sqrt(var_d)
            z_d = delta / se_d
            out['contrast_constr_minus_unconstr_delta'] = delta
            out['contrast_constr_minus_unconstr_se'] = se_d
            out['contrast_z'] = z_d
            out['contrast_p_two_sided'] = float(2 * stats.norm.sf(abs(z_d)))
            out['contrast_p_one_sided'] = float(stats.norm.sf(z_d))

    if extra_static_constr and extra_static_constr in res.params.index:
        cov_mat = res.cov_params()
        var_d = (cov_mat.loc[exposure_constr, exposure_constr]
                 + cov_mat.loc[extra_static_constr, extra_static_constr]
                 - 2 * cov_mat.loc[exposure_constr, extra_static_constr])
        delta = beta - float(res.params[extra_static_constr])
        if var_d > 0:
            se_d = math.sqrt(var_d)
            out['spec_diff_minus_static_delta'] = delta
            out['spec_se'] = se_d
            out['spec_z'] = delta / se_d
            out['spec_p_two_sided'] = float(2 * stats.norm.sf(abs(delta / se_d)))
            out['spec_p_one_sided'] = float(stats.norm.sf(delta / se_d))
    return out


def fit_gee_matrix(df, exposure_constr, exposure_unconstr,
                    family_id_col='FAMILYID', covariates=None,
                    label='', outcome='ASD_vs_Sibling'):
    """v5 Priority 1-3: matrix API (sm.GEE) で formula parser に依存しない実装.
    dummy column 名が `.`, `-`, 数値 prefix を含んでも安全."""
    cov = covariates or []
    cols_needed = ['case', exposure_constr, exposure_unconstr, family_id_col] + cov
    sub = df[[c for c in cols_needed if c in df.columns]].dropna().copy()
    if sub.empty or len(sub) < 50:
        return {'label': label, 'n_complete': len(sub), 'fit_status': 'too_few_samples'}

    # v5 Priority 1-3: drop zero-variance covariates
    cov_kept = drop_zero_variance(sub, cov)

    used_cols = [exposure_constr, exposure_unconstr] + cov_kept
    # numeric cast (Sex が string で入る可能性に備える)
    X = sub[used_cols].apply(pd.to_numeric, errors='coerce')
    if X.isna().any().any():
        bad_cols = X.columns[X.isna().any()].tolist()
        return {'label': label, 'n_complete': len(sub),
                 'fit_status': f'numeric_cast_error:{bad_cols}'}
    X = sm.add_constant(X.astype(float), has_constant='add')
    y = sub['case'].astype(int)
    groups = sub[family_id_col].astype(str)

    try:
        fam = sm.families.Binomial()
        ind = sm.cov_struct.Independence()
        res = sm.GEE(y, X, groups=groups, family=fam, cov_struct=ind).fit()
    except Exception as e:
        return {'label': label, 'n_complete': len(sub), 'fit_status': f'fit_error:{e}'}

    if exposure_constr not in res.params.index:
        return {'label': label, 'n_complete': len(sub), 'fit_status': 'param_missing'}

    beta = float(res.params[exposure_constr])
    se = float(res.bse[exposure_constr])
    p_two = float(res.pvalues.get(exposure_constr, float('nan')))
    z = beta / se if se else float('nan')
    p_one = float(stats.norm.sf(z)) if not np.isnan(z) else float('nan')

    out = {
        'label': label, 'outcome': outcome, 'exposure': exposure_constr,
        'n_complete': len(sub),
        'beta': beta, 'se': se,
        'OR': math.exp(beta) if np.isfinite(beta) else float('nan'),
        'OR_lo95': math.exp(beta - 1.96 * se) if np.isfinite(se) else float('nan'),
        'OR_hi95': math.exp(beta + 1.96 * se) if np.isfinite(se) else float('nan'),
        'p_two_sided': p_two, 'p_one_sided': p_one,
        'fit_status': 'ok_gee',
    }

    if exposure_unconstr in res.params.index:
        try:
            cov_mat = res.cov_params()
            var_d = (cov_mat.loc[exposure_constr, exposure_constr]
                     + cov_mat.loc[exposure_unconstr, exposure_unconstr]
                     - 2 * cov_mat.loc[exposure_constr, exposure_unconstr])
            delta = beta - float(res.params[exposure_unconstr])
            if var_d > 0:
                se_d = math.sqrt(var_d)
                out['contrast_constr_minus_unconstr_delta'] = delta
                out['contrast_constr_minus_unconstr_se'] = se_d
                out['contrast_z'] = delta / se_d
                out['contrast_p_two_sided'] = float(2 * stats.norm.sf(abs(delta / se_d)))
                out['contrast_p_one_sided'] = float(stats.norm.sf(delta / se_d))
        except Exception as e:
            out['contrast_status'] = f'contrast_error:{e}'
    return out


def ivw_meta(beta_a, se_a, beta_b, se_b):
    w_a = 1.0 / (se_a ** 2)
    w_b = 1.0 / (se_b ** 2)
    beta_meta = (beta_a * w_a + beta_b * w_b) / (w_a + w_b)
    se_meta = math.sqrt(1.0 / (w_a + w_b))
    z = beta_meta / se_meta
    p_two = float(2 * stats.norm.sf(abs(z)))
    p_one = float(stats.norm.sf(z))
    q = w_a * (beta_a - beta_meta) ** 2 + w_b * (beta_b - beta_meta) ** 2
    p_q = float(stats.chi2.sf(q, df=1))
    i2 = max(0.0, (q - 1) / q) if q > 0 else 0.0
    return {'beta_meta': beta_meta, 'se_meta': se_meta,
             'OR_meta': math.exp(beta_meta) if np.isfinite(beta_meta) else float('nan'),
             'OR_lo95_meta': math.exp(beta_meta - 1.96 * se_meta) if np.isfinite(se_meta) else float('nan'),
             'OR_hi95_meta': math.exp(beta_meta + 1.96 * se_meta) if np.isfinite(se_meta) else float('nan'),
             'p_meta_two_sided': p_two, 'p_meta_one_sided': p_one,
             'cochran_Q': q, 'p_cochran_Q': p_q, 'I2': i2}


# ============================================================
# Discovery WGS
# ============================================================
print(f"[{time.time()-t0:.1f}s] === Discovery WGS ===")
df_disc = pd.read_csv(EXPOSURE_DISC, sep='\t')
assert_cols(df_disc, ['sample_id', 'Diagnosis', 'N_constr_strict_diff_any'] + WGS_COVS, 'WGS')
df_disc_asd = df_disc[df_disc['Diagnosis'].isin(['ASD', 'Healthy'])].copy()
df_disc_asd['case'] = (df_disc_asd['Diagnosis'] == 'ASD').astype(int)
print(f"  WGS analytic: ASD={df_disc_asd['case'].sum()}, "
      f"Healthy={(df_disc_asd['case']==0).sum()}")

results = []
r1 = fit_logistic(df_disc_asd, 'N_constr_strict_diff_any',
                   'N_unconstr_strict_known_diff_any',
                   covariates=WGS_COVS, label='Test1_Discovery_LOEUF0.35_PRIMARY')
results.append(r1)
require_fit_ok(r1, 'Test1 Discovery WGS', allowed=('ok',))
print(f"  Test1: ok | β={r1['beta']:.4g} SE={r1['se']:.4g} P_two={r1['p_two_sided']:.4g} "
      f"P_contrast={r1.get('contrast_p_two_sided','N/A')}")

r1b = fit_logistic(df_disc_asd, 'N_constr_relaxed_diff_any',
                    'N_unconstr_relaxed_known_diff_any',
                    covariates=WGS_COVS, label='Test1b_Discovery_LOEUF0.6_sensitivity')
results.append(r1b)

r1c = fit_logistic(df_disc_asd, 'N_pLI_constr_diff_any',
                    'N_pLI_unconstr_known_diff_any',
                    covariates=WGS_COVS, label='Test1c_Discovery_pLI0.9_sensitivity')
results.append(r1c)

r2 = fit_logistic(df_disc_asd, 'N_constr_strict_diff_any',
                   'N_unconstr_strict_known_diff_any',
                   extra_static_constr='N_constr_strict_static',
                   extra_static_unconstr='N_unconstr_strict_known_static',
                   covariates=WGS_COVS, label='Test2_Discovery_Specificity_LOEUF0.35')
results.append(r2)


# ============================================================
# arrayCGH (v5 Priority 1-4: stratified IVW = primary)
# ============================================================
print(f"\n[{time.time()-t0:.1f}s] === arrayCGH (v5: stratified IVW PRIMARY) ===")
df_acgh = pd.read_csv(EXPOSURE_ACGH, sep='\t')
assert_cols(df_acgh, ['sample_id', 'Diagnosis', 'N_constr_strict_diff_any',
                       'platform_nimblegen'] + ACGH_COVS_STRAT, 'arrayCGH')
df_acgh_asd = df_acgh[df_acgh['Diagnosis'].isin(['ASD', 'CONT'])].copy()
df_acgh_asd['case'] = (df_acgh_asd['Diagnosis'] == 'ASD').astype(int)
print(f"  arrayCGH analytic: ASD={df_acgh_asd['case'].sum()}, "
      f"CONT={(df_acgh_asd['case']==0).sum()}")

# v5 Primary: platform-stratified Agilent + NimbleGen → within-array IVW
strat_results = []
for plat_label, plat_filter in [
    ('Agilent', df_acgh_asd['platform_nimblegen'] == 0),
    ('NimbleGen', df_acgh_asd['platform_nimblegen'] == 1),
]:
    sub = df_acgh_asd[plat_filter].copy()
    print(f"  arrayCGH stratum {plat_label}: n={len(sub)}")
    if len(sub) < 50:
        raise RuntimeError(f"arrayCGH {plat_label} too few samples ({len(sub)})")
    r_strat = fit_logistic(sub, 'N_constr_strict_diff_any',
                            'N_unconstr_strict_known_diff_any',
                            covariates=ACGH_COVS_STRAT,
                            label=f'arrayCGH_stratified_{plat_label}_LOEUF0.35')
    strat_results.append((plat_label, r_strat))
    results.append(r_strat)
    require_fit_ok(r_strat, f'arrayCGH stratum {plat_label}', allowed=('ok',))
    print(f"    {plat_label}: β={r_strat['beta']:.4g} SE={r_strat['se']:.4g} "
          f"P_two={r_strat['p_two_sided']:.4g}")

# v5 Primary arrayCGH: within-array IVW (Agilent + NimbleGen)
r_a = strat_results[0][1]
r_n = strat_results[1][1]
acgh_iv = ivw_meta(r_a['beta'], r_a['se'], r_n['beta'], r_n['se'])
r_acgh_primary = {
    'label': 'Test3_arrayCGH_within-array_IVW_LOEUF0.35_PRIMARY',
    'fit_status': 'ok',
    'beta': acgh_iv['beta_meta'], 'se': acgh_iv['se_meta'],
    'OR': acgh_iv['OR_meta'],
    'OR_lo95': acgh_iv['OR_lo95_meta'], 'OR_hi95': acgh_iv['OR_hi95_meta'],
    'p_two_sided': acgh_iv['p_meta_two_sided'],
    'p_one_sided': acgh_iv['p_meta_one_sided'],
    'cochran_Q': acgh_iv['cochran_Q'], 'p_cochran_Q': acgh_iv['p_cochran_Q'],
    'I2': acgh_iv['I2'],
}
results.append(r_acgh_primary)
require_fit_ok(r_acgh_primary, 'Test3 arrayCGH within-array IVW', allowed=('ok',))
print(f"  arrayCGH within-array IVW (PRIMARY): β={r_acgh_primary['beta']:.4g} "
      f"P_two={r_acgh_primary['p_two_sided']:.4g} I²={r_acgh_primary['I2']:.3g}")

# Pooled platform-adjusted as sensitivity
r_acgh_pooled = fit_logistic(df_acgh_asd, 'N_constr_strict_diff_any',
                              'N_unconstr_strict_known_diff_any',
                              covariates=ACGH_COVS_POOLED,
                              label='arrayCGH_pooled_LOEUF0.35_sensitivity')
results.append(r_acgh_pooled)
print(f"  arrayCGH pooled sensitivity: {r_acgh_pooled.get('fit_status')} "
      f"β={r_acgh_pooled.get('beta','N/A')}")


# ============================================================
# MSSNG GEE (v5: matrix API)
# ============================================================
print(f"\n[{time.time()-t0:.1f}s] === MSSNG (v5: matrix API GEE) ===")
exp_mssng = pd.read_csv(EXPOSURE_MSSNG, sep='\t')
assert_cols(exp_mssng, ['sample_id', 'Diagnosis', 'N_constr_strict_diff_any']
                       + MSSNG_REQUIRED_COVS, 'MSSNG')
df_mssng_asd = exp_mssng[exp_mssng['Diagnosis'].isin(['ASD', 'Sibling'])].copy()
df_mssng_asd['case'] = (df_mssng_asd['Diagnosis'] == 'ASD').astype(int)
print(f"  MSSNG analytic: ASD={df_mssng_asd['case'].sum()}, "
      f"Sibling={(df_mssng_asd['case']==0).sum()}")

# Available covariates (Sex + DEL burden + sanitized dummies)
avail = [c for c in MSSNG_REQUIRED_COVS if c != 'FAMILYID']
for c in df_mssng_asd.columns:
    if c.startswith('ancestry_') or c.startswith('platform_'):
        avail.append(c)
print(f"  MSSNG covariates available: {avail}")

r_mssng = fit_gee_matrix(df_mssng_asd, 'N_constr_strict_diff_any',
                          'N_unconstr_strict_known_diff_any',
                          family_id_col='FAMILYID', covariates=avail,
                          label='Test3_MSSNG_LOEUF0.35_GEE_PRIMARY')
results.append(r_mssng)
require_fit_ok(r_mssng, 'Test3 MSSNG GEE', allowed=('ok_gee',))
print(f"  MSSNG GEE: ok | β={r_mssng['beta']:.4g} SE={r_mssng['se']:.4g} "
      f"P_two={r_mssng['p_two_sided']:.4g} "
      f"P_contrast={r_mssng.get('contrast_p_two_sided','N/A')}")


# ============================================================
# Test 4: 2-cohort External IVW Meta (PRIMARY)
# arrayCGH within-array IVW estimate + MSSNG GEE estimate
# ============================================================
print(f"\n[{time.time()-t0:.1f}s] === Test 4: 2-cohort external IVW meta (PRIMARY) ===")
meta = ivw_meta(r_acgh_primary['beta'], r_acgh_primary['se'],
                 r_mssng['beta'], r_mssng['se'])
meta_record = {'label': 'Test4_2cohort_external_IVW_meta_LOEUF0.35_PRIMARY'}
meta_record.update(meta)
meta_record['fit_status'] = 'ok'
results.append(meta_record)
require_fit_ok(meta_record, 'Test4 External IVW Meta', allowed=('ok',))
print(f"  Test4 External IVW Meta: β={meta['beta_meta']:.4g} "
      f"P_two={meta['p_meta_two_sided']:.4g} I²={meta['I2']:.3g}")


# ============================================================
# Holm correction across pre-specified primary tests (Test 1 + Test 4)
# ============================================================
print(f"\n[{time.time()-t0:.1f}s] === Holm correction (Test 1 + Test 4) ===")
p1 = r1['p_one_sided']
p4 = meta['p_meta_one_sided']
_, p_holm, _, _ = multipletests([p1, p4], method='holm')
results.append({'label': 'Holm_PrimaryFamily_Test1_Discovery',
                 'p_one_sided': p1, 'p_holm_adjusted': float(p_holm[0])})
results.append({'label': 'Holm_PrimaryFamily_Test4_Meta',
                 'p_one_sided': p4, 'p_holm_adjusted': float(p_holm[1])})
print(f"  Holm Test1: {p_holm[0]:.4g}, Holm Test4: {p_holm[1]:.4g}")


# ============================================================
# v5 Priority 1-5: Contrast (β_constr - β_unconstr) external IVW meta
# arrayCGH stratified contrast (Agilent + NimbleGen IVW) + MSSNG contrast
# Reported as supporting analysis to the cover letter "across cohorts" claim
# ============================================================
print(f"\n[{time.time()-t0:.1f}s] === Contrast external IVW meta (supporting) ===")

def get_contrast(res):
    d = res.get('contrast_constr_minus_unconstr_delta')
    s = res.get('contrast_constr_minus_unconstr_se')
    if d is None or s is None or not np.isfinite(d) or not np.isfinite(s) or s <= 0:
        return None
    return d, s

# arrayCGH within-array contrast IVW (Agilent + NimbleGen)
c_a = get_contrast(r_a)
c_n = get_contrast(r_n)
if c_a and c_n:
    acgh_contrast_iv = ivw_meta(c_a[0], c_a[1], c_n[0], c_n[1])
    rec = {'label': 'arrayCGH_within-array_contrast_IVW_supporting'}
    rec.update(acgh_contrast_iv)
    results.append(rec)
    print(f"  arrayCGH within-array contrast IVW: Δ={acgh_contrast_iv['beta_meta']:.4g} "
          f"P_two={acgh_contrast_iv['p_meta_two_sided']:.4g}")
    # 2-cohort external contrast IVW (arrayCGH IVW + MSSNG)
    c_m = get_contrast(r_mssng)
    if c_m:
        ext_contrast_iv = ivw_meta(acgh_contrast_iv['beta_meta'], acgh_contrast_iv['se_meta'],
                                     c_m[0], c_m[1])
        rec2 = {'label': 'Test4_contrast_external_IVW_meta_supporting'}
        rec2.update(ext_contrast_iv)
        results.append(rec2)
        print(f"  Test4 contrast external IVW meta: Δ={ext_contrast_iv['beta_meta']:.4g} "
              f"P_two={ext_contrast_iv['p_meta_two_sided']:.4g} I²={ext_contrast_iv['I2']:.3g}")
    else:
        print(f"  MSSNG contrast unavailable; external contrast meta skipped.")
else:
    print(f"  arrayCGH stratified contrasts unavailable; external contrast meta skipped.")


# ============================================================
# Save results
# ============================================================
out_results = OUT_DIR / 'constraint_enrichment_results_v6.tsv'
res_df = pd.DataFrame(results)
res_df.to_csv(out_results, sep='\t', index=False)
print(f"\n  Saved {out_results}")
print(f"  Total rows: {len(res_df)}")
print(f"\n[Done in {time.time()-t0:.1f}s]")
