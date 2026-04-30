#!/usr/bin/env python3
"""
compute_l2_jaccard_v2.py - L2 boundary class annotation-level redundancy analysis (unified)

処理内容 (v244 round Phase B; ChatGPT review v159 指摘 全反映 + 1 script 化):
- bin_l2_annotation_v2.tsv.gz から L2 differential boundary bin (Diff_any union; n_L2_diff_support>0)
  を抽出。これは primary boundary-burden analysis で使用された post-filtered Diff_any
  union (4,980 bins) と同一の inclusion/exclusion mask を継承する。
- 10 classes 間の per-bin Jaccard overlap matrix (10x10) を計算
- Binary class-membership vectors の Pearson 相関行列 (10x10) を計算
- 相関行列の固有値から Meff (effective number of independent tests) を 3 法で算出:
    - Galwey (2009)
    - Li-Ji (2005, simplified)
    - Cheverud-Nyholt (1999)
- 統計サマリ: off-diagonal Jaccard (min/median/mean/max + max-pair),
  multi-class bin count, mean class membership
- Heatmap figure 生成 (annotation 圧縮 3 行 + ASCII symbols + 余白拡張)
- 全 outputs を一つの output ディレクトリに保存
- Reproducibility assertions: matrix shape, symmetry, max pair, multi-class count

Inputs (CLI args):
  --input PATH        bin_l2_annotation_v2.tsv.gz (REQUIRED)
  --outdir PATH       output directory (REQUIRED; auto-created)
  --out-prefix STR    figure file prefix (default: Supplementary_Fig_L2_Jaccard)

Outputs (saved to --outdir):
  l2_jaccard_matrix.tsv          pairwise Jaccard (10x10)
  l2_correlation_matrix.tsv      Pearson correlation between binary class vectors (10x10)
  l2_meff_summary.tsv            Meff (3 methods) + eigenvalues
  l2_jaccard_stats.tsv           off-diagonal stats + multi-class count
  l2_summary.json                machine-readable full summary
  <out-prefix>.svg / .pdf        heatmap figure

注意 (ChatGPT review 反映):
- Meff は Pearson correlation matrix of 10 binary class-membership vectors の固有値から
  計算する (Jaccard matrix そのものを固有値解析しない)
- 表現は "limited annotation-level redundancy" (annotation level) であり、
  association test 全体の独立性を主張しない
- Figure 内 symbols は ASCII (>=, <->)
- Source data reference: bin_l2_annotation_v2.tsv (Supp Methods §12.1 + §12.6)

Example:
  python compute_l2_jaccard_v2.py \
      --input  /home/kushima-pg/tad04292026/02_bin_l2_annotation/output_v2/bin_l2_annotation_v2.tsv.gz \
      --outdir /home/kushima-pg/tad04292026/15_l2_jaccard_v1/output \
      --out-prefix Supplementary_Fig_L2_Jaccard
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


CLASSES_DEFAULT = [
    'HPC_Exc-CA', 'HPC_Exc-DG', 'HPC_Exc-ENT',
    'HPC_Inh-CGE', 'HPC_Inh-MGE',
    'PFC_Astro',
    'PFC_Exc-DL', 'PFC_Exc-UL',
    'PFC_Inh-CGE', 'PFC_Inh-MGE',
]


def load_membership(input_path, classes):
    """
    Load Diff bins membership matrix from bin_l2_annotation_v2.tsv.gz.
    Diff bins = bins with n_L2_diff_support > 0 (post-filtered Diff_any union
    used in primary boundary-burden analysis; expected n=4,980).
    Returns (M, bin_ids) where M is (n_diff_bins, n_classes) binary matrix.
    """
    df = pd.read_csv(input_path, sep='\t', compression='infer')
    diff_mask = df['n_L2_diff_support'] > 0
    df_diff = df.loc[diff_mask].reset_index(drop=True)

    mem_cols = [f'membership_{c}' for c in classes]
    missing = [c for c in mem_cols if c not in df_diff.columns]
    if missing:
        raise ValueError(f'Missing membership columns: {missing}')

    M = df_diff[mem_cols].astype(int).values
    bin_ids = df_diff['bin_id'].tolist()
    return M, bin_ids


def compute_jaccard_matrix(M):
    """Compute pairwise Jaccard overlap |A∩B|/|A∪B| between class membership."""
    n = M.shape[1]
    J = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            a = M[:, i].astype(bool)
            b = M[:, j].astype(bool)
            inter = (a & b).sum()
            union = (a | b).sum()
            J[i, j] = inter / union if union > 0 else 0.0
    return J


def compute_meff(M):
    """
    Compute effective number of independent tests from Pearson correlation matrix
    of binary class-membership vectors.

    Returns dict with three Meff estimates:
      - Galwey (2009):           (sum sqrt(pos_eig))^2 / sum(pos_eig)
      - Li-Ji (2005, simplified): sum_i [I(λ_i ≥ 1) + (λ_i - floor(λ_i))]
      - Cheverud-Nyholt (1999):  1 + (M-1)*(1 - var(λ)/M)

    Note: Meff is computed from the Pearson/phi correlation matrix of binary
    class-membership vectors across the post-filtered Diff_any bin universe.
    """
    C = np.corrcoef(M.T)
    eigvals = np.linalg.eigvalsh(C)
    eigvals = np.sort(eigvals)[::-1]  # descending
    n = len(eigvals)

    pos_eig = np.maximum(eigvals, 0)

    # Galwey (2009)
    sqrt_sum = np.sqrt(pos_eig).sum()
    sum_pos = pos_eig.sum()
    Meff_galwey = (sqrt_sum ** 2) / sum_pos if sum_pos > 0 else float('nan')

    # Li-Ji (2005) simplified
    Meff_liji = float(sum(int(e >= 1) + (e - np.floor(e)) for e in eigvals if e > 0))

    # Cheverud-Nyholt (1999)
    var_eig = float(np.var(eigvals))
    Meff_chev = 1.0 + (n - 1) * (1.0 - var_eig / n)

    return {
        'eigenvalues_correlation': eigvals.tolist(),
        'M_classes': int(n),
        'Meff_Galwey_2009': float(Meff_galwey),
        'Meff_LiJi_2005_simplified': float(Meff_liji),
        'Meff_CheverudNyholt_1999': float(Meff_chev),
        'correlation_matrix_source': 'Pearson correlation of 10 binary class-membership vectors across the post-filtered Diff_any bin universe',
    }


def compute_stats(M, J, classes):
    """Compute off-diagonal Jaccard stats + multi-class membership counts."""
    n = J.shape[0]
    upper = []
    max_pair = (None, None, -1.0)
    for i in range(n):
        for j in range(i + 1, n):
            v = J[i, j]
            upper.append(v)
            if v > max_pair[2]:
                max_pair = (classes[i], classes[j], v)
    upper = np.array(upper)

    multi_class = int((M.sum(axis=1) >= 2).sum())
    n_diff = int(M.shape[0])
    mean_class_count = float(M.sum(axis=1).mean())

    return {
        'n_diff_bins': n_diff,
        'n_pairs': len(upper),
        'jaccard_min': float(upper.min()),
        'jaccard_median': float(np.median(upper)),
        'jaccard_mean': float(upper.mean()),
        'jaccard_max': float(upper.max()),
        'max_pair_class_a': max_pair[0],
        'max_pair_class_b': max_pair[1],
        'multi_class_bins_ge2': multi_class,
        'mean_class_membership_per_bin': mean_class_count,
    }


def make_heatmap(J, classes, stats, meff, out_svg, out_pdf):
    """Generate heatmap figure (improved layout per ChatGPT review)."""
    n = len(classes)
    disp = [c.replace('_', ' ') for c in classes]

    fig = plt.figure(figsize=(8.5, 9.5))

    # Heatmap (raised to give more annotation space below)
    ax = fig.add_axes([0.20, 0.42, 0.60, 0.54])

    J_plot = J.copy()
    for i in range(n):
        J_plot[i, i] = np.nan

    cmap = plt.cm.viridis
    im = ax.imshow(J_plot, cmap=cmap, aspect='equal', vmin=0.0, vmax=0.30,
                   interpolation='nearest')

    ax.set_xticks(range(n))
    ax.set_xticklabels(disp, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(n))
    ax.set_yticklabels(disp, fontsize=9)

    # Cell annotations
    for i in range(n):
        for j in range(n):
            if i == j:
                ax.text(j, i, '-', ha='center', va='center',
                        fontsize=8, color='#888888')
            else:
                val = J[i, j]
                color = 'white' if val > 0.18 else 'black'
                ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                        fontsize=7.5, color=color)

    # Colorbar
    cbar_ax = fig.add_axes([0.82, 0.42, 0.025, 0.54])
    cb = fig.colorbar(im, cax=cbar_ax, ticks=[0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
    cb.ax.set_ylabel('Jaccard overlap (per-bin)', fontsize=9, rotation=270, labelpad=12)
    cb.ax.tick_params(labelsize=8)

    # Annotation panel (3 lines, ASCII symbols, generous spacing below x-axis labels)
    ann_ax = fig.add_axes([0.10, 0.04, 0.85, 0.18])
    ann_ax.set_axis_off()
    pct = 100 * stats['multi_class_bins_ge2'] / stats['n_diff_bins']
    max_pair_disp_a = stats['max_pair_class_a'].replace('_', ' ')
    max_pair_disp_b = stats['max_pair_class_b'].replace('_', ' ')
    ann_lines = [
        f'Off-diagonal Jaccard (n={stats["n_pairs"]} pairs): median = {stats["jaccard_median"]:.3f}, mean = {stats["jaccard_mean"]:.3f}, max = {stats["jaccard_max"]:.3f} ({max_pair_disp_a} <-> {max_pair_disp_b}).',
        f'Effective number of independent tests: {meff["Meff_Galwey_2009"]:.2f} (Galwey 2009), {meff["Meff_LiJi_2005_simplified"]:.2f} (Li-Ji 2005), {meff["Meff_CheverudNyholt_1999"]:.2f} (Cheverud-Nyholt 1999).',
        f'Multi-class diff bins (>=2 classes): {stats["multi_class_bins_ge2"]:,} of {stats["n_diff_bins"]:,} = {pct:.1f}%.',
    ]
    y0 = 0.92
    for i, ln in enumerate(ann_lines):
        ann_ax.text(0.0, y0 - i * 0.32, ln, ha='left', va='top',
                    fontsize=9, family='sans-serif', wrap=True)

    plt.savefig(out_svg, format='svg', bbox_inches='tight')
    plt.savefig(out_pdf, format='pdf', bbox_inches='tight')
    plt.close(fig)


def run_assertions(M, J, stats, meff, classes):
    """Reproducibility checks (assert-based)."""
    # Matrix shape
    assert M.shape[1] == 10, f'Expected 10 classes, got {M.shape[1]}'
    assert J.shape == (10, 10), f'Jaccard matrix shape != 10x10 ({J.shape})'

    # Symmetry
    assert np.allclose(J, J.T, atol=1e-10), 'Jaccard matrix is not symmetric'

    # Diagonal
    assert np.allclose(np.diag(J), 1.0, atol=1e-10), \
        'Diagonal of Jaccard matrix should be 1.0'

    # Max pair (expected: PFC Exc-DL <-> PFC Exc-UL by prior NIG computation)
    expected_max_pair = {'PFC_Exc-DL', 'PFC_Exc-UL'}
    actual_max_pair = {stats['max_pair_class_a'], stats['max_pair_class_b']}
    assert expected_max_pair == actual_max_pair, \
        f'Expected max pair {expected_max_pair}, got {actual_max_pair}'

    # Stats sanity
    assert stats['n_diff_bins'] > 0
    assert 0.0 <= stats['jaccard_min'] <= stats['jaccard_max'] <= 1.0
    assert stats['jaccard_median'] <= stats['jaccard_mean'] + 0.10  # mean ~= median + slight skew

    # Meff range
    assert meff['Meff_Galwey_2009'] > 0
    assert meff['Meff_Galwey_2009'] <= meff['M_classes']

    print('  All reproducibility assertions PASSED.')


def main():
    parser = argparse.ArgumentParser(
        description='L2 boundary class annotation-level redundancy analysis')
    parser.add_argument('--input', required=True,
                        help='Path to bin_l2_annotation_v2.tsv.gz')
    parser.add_argument('--outdir', required=True,
                        help='Output directory (auto-created)')
    parser.add_argument('--out-prefix', default='Supplementary_Fig_L2_Jaccard',
                        help='Output figure file prefix')
    args = parser.parse_args()

    t0 = time.time()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    classes = CLASSES_DEFAULT

    # Load
    print(f'Loading membership data from {args.input} ...')
    M, bin_ids = load_membership(args.input, classes)
    print(f'  Loaded {M.shape[0]} Diff bins x {M.shape[1]} classes')

    # Compute
    print('Computing Jaccard matrix ...')
    J = compute_jaccard_matrix(M)

    print('Computing Pearson correlation + Meff (3 methods) ...')
    meff = compute_meff(M)

    print('Computing summary stats ...')
    stats = compute_stats(M, J, classes)

    # Reproducibility checks
    print('Running reproducibility assertions ...')
    run_assertions(M, J, stats, meff, classes)

    # Print summary
    print()
    print('=== Summary ===')
    print(f'Diff bins (post-filtered Diff_any union): {stats["n_diff_bins"]:,}')
    print(f'Off-diagonal Jaccard (n={stats["n_pairs"]}): '
          f'min={stats["jaccard_min"]:.3f}, '
          f'median={stats["jaccard_median"]:.3f}, '
          f'mean={stats["jaccard_mean"]:.3f}, '
          f'max={stats["jaccard_max"]:.3f}')
    print(f'Max pair: {stats["max_pair_class_a"]} <-> {stats["max_pair_class_b"]}')
    print(f'Multi-class bins (>=2): {stats["multi_class_bins_ge2"]:,}/{stats["n_diff_bins"]:,} '
          f'= {100*stats["multi_class_bins_ge2"]/stats["n_diff_bins"]:.1f}%')
    print(f'Mean class membership per diff bin: {stats["mean_class_membership_per_bin"]:.3f}')
    print(f'Meff (Galwey 2009):         {meff["Meff_Galwey_2009"]:.3f}')
    print(f'Meff (Li-Ji 2005):          {meff["Meff_LiJi_2005_simplified"]:.3f}')
    print(f'Meff (Cheverud-Nyholt 1999): {meff["Meff_CheverudNyholt_1999"]:.3f}')
    print()

    # Save outputs
    print('Saving outputs ...')

    # Jaccard matrix
    pd.DataFrame(J, index=classes, columns=classes).to_csv(
        outdir / 'l2_jaccard_matrix.tsv', sep='\t', float_format='%.6f')

    # Correlation matrix
    C = np.corrcoef(M.T)
    pd.DataFrame(C, index=classes, columns=classes).to_csv(
        outdir / 'l2_correlation_matrix.tsv', sep='\t', float_format='%.6f')

    # Meff summary
    pd.DataFrame([
        {'method': 'Galwey_2009', 'Meff': meff['Meff_Galwey_2009']},
        {'method': 'LiJi_2005_simplified', 'Meff': meff['Meff_LiJi_2005_simplified']},
        {'method': 'CheverudNyholt_1999', 'Meff': meff['Meff_CheverudNyholt_1999']},
    ]).to_csv(outdir / 'l2_meff_summary.tsv', sep='\t', index=False, float_format='%.6f')

    # Stats
    pd.DataFrame([stats]).to_csv(outdir / 'l2_jaccard_stats.tsv',
                                 sep='\t', index=False, float_format='%.6f')

    # Full summary JSON
    summary = {
        'input_file': str(args.input),
        'classes': classes,
        'stats': stats,
        'meff': meff,
        'note': 'Meff was computed from the Pearson correlation matrix of 10 binary class-membership vectors across the post-filtered Diff_any bin universe used in the primary boundary-burden analysis. Jaccard overlap and Meff describe annotation-level redundancy among the 10 L2 boundary classes; they do not, on their own, establish full independence of the underlying association tests, which depend additionally on per-individual deletion burden, SV length, gene density, and covariates.',
    }
    with open(outdir / 'l2_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    # Figure
    out_svg = outdir / f'{args.out_prefix}.svg'
    out_pdf = outdir / f'{args.out_prefix}.pdf'
    make_heatmap(J, classes, stats, meff, str(out_svg), str(out_pdf))
    print(f'  Saved {out_svg}')
    print(f'  Saved {out_pdf}')

    print()
    print(f'Elapsed: {time.time() - t0:.2f}s')


if __name__ == '__main__':
    main()
