#!/usr/bin/env python3
"""
14_fit_constraint_models_v8.py

処理内容 (v7 → v8 変更点; ChatGPT v7 review Block 1, 2 fix):
- Block 1 fix: require_fit_ok() を IVW meta result の beta_meta/se_meta key
  にも対応させる (v7 では beta/se key のみ要求 → primary Test 4 meta が
  必ず hard fail で落ちる bug)。さらに np.isfinite() で NumPy scalar/int も
  安全に判定。
- Block 2 fix: primary variant で arrayCGH within-array IVW 不成立、
  MSSNG GEE 不成立、Test 4 external IVW meta 不成立を **明示的に hard fail**
  させる。v7 では r_acgh_iv=None / r_meta=None で silent skip → Holm が
  skip されるだけだったが、primary が欠落したまま結果ファイルが作成される
  事故を防ぐ。
- ChatGPT 推奨に従い meta record に beta/se aliases (beta_meta → beta) を追加
  し、Holm に必要な p_meta_one_sided も alias として p_one_sided に複写。

v7 → v8 変更点 (実装):
1. require_fit_ok() を堅くし、beta/beta_meta どちらでも検出可能に
2. run_3cohort_meta() 内で primary に対する 3 つの hard fail check を追加:
   (a) r_acgh_iv is None
   (b) r_mssng fit_status not ok
   (c) r_meta is None
3. r_meta record に beta/se/OR/p_one_sided/p_two_sided alias を追加

v5 から継承 (ChatGPT v4 review safety patch):
- GEE matrix API (formula parser 非依存)
- arrayCGH within-array IVW を Test 3 primary
- contrast (β_constr - β_unconstr) IVW meta supporting
- zero-variance covariate auto drop
- Holm 補正: pre-specified primary 2-test family のみ

入力:
- output_v1/sample_constraint_burden_discovery_v7.tsv (burden v7 出力)
- output_v1/sample_constraint_burden_arraycgh_v7.tsv
- output_v1/sample_constraint_burden_mssng_v7.tsv

出力:
- output_v1/constraint_enrichment_results_v8.tsv

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

EXPOSURE_DISC = OUT_DIR / 'sample_constraint_burden_discovery_v7.tsv'
EXPOSURE_ACGH = OUT_DIR / 'sample_constraint_burden_arraycgh_v7.tsv'
EXPOSURE_MSSNG = OUT_DIR / 'sample_constraint_burden_mssng_v7.tsv'

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


def require_fit_ok(res, name, allowed=('ok', 'ok_gee'), strict=True):
    """v8 Block 1 fix: cohort-level result (beta/se) と meta-level result
    (beta_meta/se_meta) の両方で正しく動くよう拡張。NumPy scalar/int も安全判定。
    """
    status = str(res.get('fit_status', ''))
    ok = any(status == a or status.startswith(a) for a in allowed)
    if not ok:
        msg = f"{name} fit failed: {status}; full result={res}"
        if strict:
            raise RuntimeError(msg)
        print(f"  WARN: {msg}")
        return False

    # v8: cohort-level なら beta/se、meta-level なら beta_meta/se_meta を見る
    beta = res.get('beta', res.get('beta_meta'))
    se = res.get('se', res.get('se_meta'))

    for k_label, v in [('beta(or beta_meta)', beta), ('se(or se_meta)', se)]:
        if v is None:
            msg = f"{name} has missing {k_label}: {v}"
            if strict:
                raise RuntimeError(msg)
            print(f"  WARN: {msg}")
            return False
        try:
            if not np.isfinite(float(v)):
                msg = f"{name} has non-finite {k_label}: {v}"
                if strict:
                    raise RuntimeError(msg)
                print(f"  WARN: {msg}")
                return False
        except (TypeError, ValueError):
            msg = f"{name} has non-numeric {k_label}: {v}"
            if strict:
                raise RuntimeError(msg)
            print(f"  WARN: {msg}")
            return False

    if float(se) <= 0:
        msg = f"{name} has non-positive SE: {se}"
        if strict:
            raise RuntimeError(msg)
        print(f"  WARN: {msg}")
        return False
    return True


def drop_zero_variance(sub, cols):
    keep = []
    for c in cols:
        if c not in sub.columns:
            continue
        if sub[c].nunique(dropna=True) > 1:
            keep.append(c)
    return keep


def fit_logistic(df, exposure_constr, exposure_unconstr,
                  covariates=None, label='', outcome='ASD_vs_Control'):
    cov = covariates or []
    cols_needed = ['case', exposure_constr, exposure_unconstr] + cov
    sub = df[[c for c in cols_needed if c in df.columns]].dropna().copy()
    if sub.empty or len(sub) < 50:
        return {'label': label, 'n_complete': len(sub), 'fit_status': 'too_few_samples'}
    cov_kept = drop_zero_variance(sub, cov)
    used_cols = [exposure_constr, exposure_unconstr] + cov_kept
    if sub[exposure_constr].nunique() < 2:
        return {'label': label, 'n_complete': len(sub),
                 'fit_status': f'zero_variance_exposure:{exposure_constr}'}
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
            out['contrast_constr_minus_unconstr_delta'] = delta
            out['contrast_constr_minus_unconstr_se'] = se_d
            out['contrast_z'] = delta / se_d
            out['contrast_p_two_sided'] = float(2 * stats.norm.sf(abs(delta / se_d)))
            out['contrast_p_one_sided'] = float(stats.norm.sf(delta / se_d))
    return out


def fit_gee_matrix(df, exposure_constr, exposure_unconstr,
                    family_id_col='FAMILYID', covariates=None,
                    label='', outcome='ASD_vs_Sibling'):
    cov = covariates or []
    cols_needed = ['case', exposure_constr, exposure_unconstr, family_id_col] + cov
    sub = df[[c for c in cols_needed if c in df.columns]].dropna().copy()
    if sub.empty or len(sub) < 50:
        return {'label': label, 'n_complete': len(sub), 'fit_status': 'too_few_samples'}
    if sub[exposure_constr].nunique() < 2:
        return {'label': label, 'n_complete': len(sub),
                 'fit_status': f'zero_variance_exposure:{exposure_constr}'}
    cov_kept = drop_zero_variance(sub, cov)
    used_cols = [exposure_constr, exposure_unconstr] + cov_kept
    X = sub[used_cols].apply(pd.to_numeric, errors='coerce')
    if X.isna().any().any():
        bad = X.columns[X.isna().any()].tolist()
        return {'label': label, 'n_complete': len(sub),
                 'fit_status': f'numeric_cast_error:{bad}'}
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
    if not (np.isfinite(beta_a) and np.isfinite(se_a) and np.isfinite(beta_b) and np.isfinite(se_b)):
        return {'fit_status': 'meta_input_nonfinite'}
    if se_a <= 0 or se_b <= 0:
        return {'fit_status': 'meta_input_se_nonpositive'}
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
             'cochran_Q': q, 'p_cochran_Q': p_q, 'I2': i2,
             'fit_status': 'ok'}


# ============================================================
# Load exposures
# ============================================================
print(f"[{time.time()-t0:.1f}s] Loading v7 exposure files ...")
df_disc = pd.read_csv(EXPOSURE_DISC, sep='\t')
df_acgh = pd.read_csv(EXPOSURE_ACGH, sep='\t')
df_mssng = pd.read_csv(EXPOSURE_MSSNG, sep='\t')

assert_cols(df_disc, ['sample_id', 'Diagnosis', 'N_constr_strict_diffany_bounded'] + WGS_COVS, 'WGS')
assert_cols(df_acgh, ['sample_id', 'Diagnosis', 'N_constr_strict_diffany_bounded',
                      'platform_nimblegen'] + ACGH_COVS_STRAT, 'arrayCGH')
assert_cols(df_mssng, ['sample_id', 'Diagnosis', 'N_constr_strict_diffany_bounded']
                       + MSSNG_REQUIRED_COVS, 'MSSNG')

df_disc_asd = df_disc[df_disc['Diagnosis'].isin(['ASD', 'Healthy'])].copy()
df_disc_asd['case'] = (df_disc_asd['Diagnosis'] == 'ASD').astype(int)

df_acgh_asd = df_acgh[df_acgh['Diagnosis'].isin(['ASD', 'CONT'])].copy()
df_acgh_asd['case'] = (df_acgh_asd['Diagnosis'] == 'ASD').astype(int)

df_mssng_asd = df_mssng[df_mssng['Diagnosis'].isin(['ASD', 'Sibling'])].copy()
df_mssng_asd['case'] = (df_mssng_asd['Diagnosis'] == 'ASD').astype(int)

print(f"  WGS analytic: ASD={df_disc_asd['case'].sum()}, "
      f"Healthy={(df_disc_asd['case']==0).sum()}")
print(f"  arrayCGH analytic: ASD={df_acgh_asd['case'].sum()}, "
      f"CONT={(df_acgh_asd['case']==0).sum()}")
print(f"  MSSNG analytic: ASD={df_mssng_asd['case'].sum()}, "
      f"Sibling={(df_mssng_asd['case']==0).sum()}")

# MSSNG covariates list (Sex + DEL burden + dummies)
mssng_cov_list = [c for c in MSSNG_REQUIRED_COVS if c != 'FAMILYID']
for c in df_mssng_asd.columns:
    if c.startswith('ancestry_') or c.startswith('platform_'):
        mssng_cov_list.append(c)
print(f"  MSSNG covariates: {mssng_cov_list}")


# ============================================================
# Burden variant configurations
# ============================================================
# (variant_suffix, exposure_constr_col, exposure_unconstr_col, label_prefix, primary_flag)
VARIANTS = [
    ('diffany_bounded', 'N_constr_strict_diffany_bounded', 'N_unconstr_strict_known_diffany_bounded',
     'DiffAnyBounded_LOEUF0.35', True),  # PRIMARY
    ('prox500kb', 'N_constr_strict_prox500kb', 'N_unconstr_strict_known_prox500kb',
     'Prox500kb_LOEUF0.35', False),  # primary supplement
    ('prox250kb', 'N_constr_strict_prox250kb', 'N_unconstr_strict_known_prox250kb',
     'Prox250kb_LOEUF0.35', False),
    ('prox1000kb', 'N_constr_strict_prox1000kb', 'N_unconstr_strict_known_prox1000kb',
     'Prox1Mb_LOEUF0.35', False),
    ('diffany_bounded_excl', 'N_constr_strict_diffany_bounded_excl',
     'N_unconstr_strict_known_diffany_bounded_excl',
     'DiffAnyBounded_DirectExcl_LOEUF0.35', False),
]

# Sensitivity cutoffs (only on primary diffany_bounded variant)
CUTOFF_SENSITIVITIES = [
    ('relaxed', 'N_constr_relaxed_diffany_bounded', 'N_unconstr_relaxed_known_diffany_bounded',
     'DiffAnyBounded_LOEUF0.6'),
    ('pli', 'N_pLI_constr_diffany_bounded', 'N_pLI_unconstr_known_diffany_bounded',
     'DiffAnyBounded_pLI0.9'),
]


def run_3cohort_meta(label_prefix, exposure_constr, exposure_unconstr, results,
                      strict_primary=False):
    """Run WGS Test 1 + arrayCGH stratified IVW + MSSNG GEE + 2-cohort external IVW Meta.
    Append all results to `results` list. Returns (r1, r_acgh_iv, r_mssng, r_meta).
    """
    print(f"\n[{time.time()-t0:.1f}s] === {label_prefix} ===")

    # Test 1: WGS Discovery
    r1 = fit_logistic(df_disc_asd, exposure_constr, exposure_unconstr,
                       covariates=WGS_COVS, label=f'Test1_WGS_{label_prefix}',
                       outcome='ASD_vs_Healthy')
    results.append(r1)
    require_fit_ok(r1, f'Test1 WGS {label_prefix}', allowed=('ok',), strict=strict_primary)
    print(f"  WGS Test1: status={r1.get('fit_status')} β={r1.get('beta','N/A')} "
          f"P_two={r1.get('p_two_sided','N/A')}")

    # Test 3 arrayCGH: stratified IVW (Agilent + NimbleGen)
    strat_results = []
    for plat_label, plat_filter in [
        ('Agilent', df_acgh_asd['platform_nimblegen'] == 0),
        ('NimbleGen', df_acgh_asd['platform_nimblegen'] == 1),
    ]:
        sub = df_acgh_asd[plat_filter].copy()
        if len(sub) < 50:
            print(f"  arrayCGH stratum {plat_label} too few samples ({len(sub)})")
            strat_results.append((plat_label, {'fit_status': 'too_few_samples'}))
            continue
        r_strat = fit_logistic(sub, exposure_constr, exposure_unconstr,
                                covariates=ACGH_COVS_STRAT,
                                label=f'arrayCGH_stratified_{plat_label}_{label_prefix}',
                                outcome='ASD_vs_CONT')
        strat_results.append((plat_label, r_strat))
        results.append(r_strat)
        print(f"  arrayCGH {plat_label}: {r_strat.get('fit_status')} "
              f"β={r_strat.get('beta','N/A')}")

    # arrayCGH within-array IVW (if both ok)
    r_acgh_iv = None
    if (len(strat_results) == 2
            and all(r.get('fit_status') == 'ok' for _, r in strat_results)):
        r_a = strat_results[0][1]
        r_n = strat_results[1][1]
        acgh_iv = ivw_meta(r_a['beta'], r_a['se'], r_n['beta'], r_n['se'])
        if acgh_iv.get('fit_status') == 'ok':
            r_acgh_iv = {
                'label': f'Test3_arrayCGH_within-array_IVW_{label_prefix}',
                'fit_status': 'ok',
                'outcome': 'ASD_vs_CONT',
                'beta': acgh_iv['beta_meta'], 'se': acgh_iv['se_meta'],
                'OR': acgh_iv['OR_meta'],
                'OR_lo95': acgh_iv['OR_lo95_meta'], 'OR_hi95': acgh_iv['OR_hi95_meta'],
                'p_two_sided': acgh_iv['p_meta_two_sided'],
                'p_one_sided': acgh_iv['p_meta_one_sided'],
                'cochran_Q': acgh_iv['cochran_Q'],
                'p_cochran_Q': acgh_iv['p_cochran_Q'],
                'I2': acgh_iv['I2'],
            }
            results.append(r_acgh_iv)
            require_fit_ok(r_acgh_iv, f'Test3 arrayCGH IVW {label_prefix}',
                            allowed=('ok',), strict=strict_primary)
            print(f"  arrayCGH within-array IVW: β={r_acgh_iv['beta']:.4g} "
                  f"P_two={r_acgh_iv['p_two_sided']:.4g}")

    # v8 Block 2 fix: primary 用 hard fail (a) arrayCGH IVW 不成立
    if strict_primary and r_acgh_iv is None:
        raise RuntimeError(
            f"Primary arrayCGH within-array IVW was not available for {label_prefix}. "
            f"Stratum results: {strat_results}"
        )

    # MSSNG GEE
    r_mssng = fit_gee_matrix(df_mssng_asd, exposure_constr, exposure_unconstr,
                              family_id_col='FAMILYID', covariates=mssng_cov_list,
                              label=f'Test3_MSSNG_GEE_{label_prefix}',
                              outcome='ASD_vs_Sibling')
    results.append(r_mssng)
    require_fit_ok(r_mssng, f'Test3 MSSNG GEE {label_prefix}',
                    allowed=('ok_gee',), strict=strict_primary)
    print(f"  MSSNG GEE: {r_mssng.get('fit_status')} β={r_mssng.get('beta','N/A')} "
          f"P_two={r_mssng.get('p_two_sided','N/A')}")

    # v8 Block 2 fix: primary 用 hard fail (b) MSSNG GEE 不成立
    if strict_primary and not str(r_mssng.get('fit_status', '')).startswith('ok'):
        raise RuntimeError(
            f"Primary MSSNG GEE failed for {label_prefix}: {r_mssng}"
        )

    # Test 4: 2-cohort external IVW meta
    r_meta = None
    if (r_acgh_iv and r_acgh_iv.get('fit_status') == 'ok'
            and r_mssng.get('fit_status', '').startswith('ok')):
        meta = ivw_meta(r_acgh_iv['beta'], r_acgh_iv['se'],
                         r_mssng['beta'], r_mssng['se'])
        if meta.get('fit_status') == 'ok':
            # v8 Block 1 fix: meta record に beta/se/OR/p alias を追加し
            # require_fit_ok() が beta/se で検出可能にする
            r_meta = {'label': f'Test4_2cohort_external_IVW_meta_{label_prefix}'}
            r_meta.update(meta)
            r_meta['beta'] = meta['beta_meta']
            r_meta['se'] = meta['se_meta']
            r_meta['OR'] = meta['OR_meta']
            r_meta['OR_lo95'] = meta['OR_lo95_meta']
            r_meta['OR_hi95'] = meta['OR_hi95_meta']
            r_meta['p_two_sided'] = meta['p_meta_two_sided']
            r_meta['p_one_sided'] = meta['p_meta_one_sided']
            results.append(r_meta)
            require_fit_ok(r_meta, f'Test4 IVW Meta {label_prefix}',
                            allowed=('ok',), strict=strict_primary)
            print(f"  Test4 External IVW Meta: β={meta['beta_meta']:.4g} "
                  f"P_two={meta['p_meta_two_sided']:.4g} I²={meta['I2']:.3g}")

    # v8 Block 2 fix: primary 用 hard fail (c) Test 4 meta 不成立
    if strict_primary and r_meta is None:
        raise RuntimeError(
            f"Primary Test4 external IVW meta was not generated for {label_prefix}."
        )
    return r1, r_acgh_iv, r_mssng, r_meta


# ============================================================
# Run all variants
# ============================================================
results = []
primary_r1 = None
primary_meta = None

for variant_suf, exp_c, exp_u, label_prefix, is_primary in VARIANTS:
    r1, r_acgh_iv, r_mssng, r_meta = run_3cohort_meta(
        label_prefix, exp_c, exp_u, results, strict_primary=is_primary)
    if is_primary:
        primary_r1 = r1
        primary_meta = r_meta

# Cutoff sensitivities (LOEUF<0.6, pLI>0.9) on primary variant
for suf, exp_c, exp_u, label_prefix in CUTOFF_SENSITIVITIES:
    run_3cohort_meta(label_prefix, exp_c, exp_u, results, strict_primary=False)


# ============================================================
# Holm correction across pre-specified primary 2-test family
# (Test 1 WGS DiffAnyBounded LOEUF<0.35) + (Test 4 IVW Meta DiffAnyBounded LOEUF<0.35)
# ============================================================
print(f"\n[{time.time()-t0:.1f}s] === Holm correction (pre-specified 2-test primary family) ===")
if primary_r1 and primary_meta:
    p1 = primary_r1.get('p_one_sided', float('nan'))
    p4 = primary_meta.get('p_meta_one_sided', float('nan'))
    if not (np.isnan(p1) or np.isnan(p4)):
        _, p_holm, _, _ = multipletests([p1, p4], method='holm')
        results.append({'label': 'Holm_PrimaryFamily_Test1_WGS_DiffAnyBounded_LOEUF0.35',
                         'p_one_sided': p1, 'p_holm_adjusted': float(p_holm[0])})
        results.append({'label': 'Holm_PrimaryFamily_Test4_Meta_DiffAnyBounded_LOEUF0.35',
                         'p_one_sided': p4, 'p_holm_adjusted': float(p_holm[1])})
        print(f"  Holm Test1: {p_holm[0]:.4g}, Holm Test4: {p_holm[1]:.4g}")
    else:
        print(f"  Holm skipped: NaN p-values (p1={p1}, p4={p4})")
else:
    print(f"  Holm skipped: primary results not available")


# ============================================================
# Contrast IVW external meta (β_constr - β_unconstr) for PRIMARY variant
# ============================================================
print(f"\n[{time.time()-t0:.1f}s] === Contrast IVW external meta (PRIMARY DiffAnyBounded) ===")

def get_contrast(res):
    if res is None:
        return None
    d = res.get('contrast_constr_minus_unconstr_delta')
    s = res.get('contrast_constr_minus_unconstr_se')
    if d is None or s is None or not np.isfinite(d) or not np.isfinite(s) or s <= 0:
        return None
    return d, s

# Find arrayCGH stratified contrast (Agilent + NimbleGen) and MSSNG contrast for primary variant
primary_label_prefix = 'DiffAnyBounded_LOEUF0.35'
r_a_strat = next((r for r in results
                   if r.get('label') == f'arrayCGH_stratified_Agilent_{primary_label_prefix}'), None)
r_n_strat = next((r for r in results
                   if r.get('label') == f'arrayCGH_stratified_NimbleGen_{primary_label_prefix}'), None)
r_mssng_primary = next((r for r in results
                         if r.get('label') == f'Test3_MSSNG_GEE_{primary_label_prefix}'), None)

c_a = get_contrast(r_a_strat)
c_n = get_contrast(r_n_strat)
c_m = get_contrast(r_mssng_primary)
if c_a and c_n:
    acgh_contrast_iv = ivw_meta(c_a[0], c_a[1], c_n[0], c_n[1])
    if acgh_contrast_iv.get('fit_status') == 'ok':
        rec = {'label': 'arrayCGH_within-array_contrast_IVW_DiffAnyBounded_supporting'}
        rec.update(acgh_contrast_iv)
        results.append(rec)
        print(f"  arrayCGH within-array contrast IVW: Δ={acgh_contrast_iv['beta_meta']:.4g} "
              f"P_two={acgh_contrast_iv['p_meta_two_sided']:.4g}")
        if c_m:
            ext_contrast = ivw_meta(acgh_contrast_iv['beta_meta'], acgh_contrast_iv['se_meta'],
                                      c_m[0], c_m[1])
            if ext_contrast.get('fit_status') == 'ok':
                rec2 = {'label': 'Test4_contrast_external_IVW_meta_DiffAnyBounded_supporting'}
                rec2.update(ext_contrast)
                results.append(rec2)
                print(f"  Test4 contrast external IVW: Δ={ext_contrast['beta_meta']:.4g} "
                      f"P_two={ext_contrast['p_meta_two_sided']:.4g}")
        else:
            print(f"  MSSNG contrast unavailable; external contrast meta skipped.")
    else:
        print(f"  arrayCGH contrast IVW skipped: {acgh_contrast_iv.get('fit_status')}")
else:
    print(f"  arrayCGH stratified contrasts unavailable; contrast meta skipped.")


# ============================================================
# Save results
# ============================================================
out_results = OUT_DIR / 'constraint_enrichment_results_v8.tsv'
res_df = pd.DataFrame(results)
res_df.to_csv(out_results, sep='\t', index=False)
print(f"\n  Saved {out_results}")
print(f"  Total rows: {len(res_df)}")
print(f"\n[Done in {time.time()-t0:.1f}s]")
