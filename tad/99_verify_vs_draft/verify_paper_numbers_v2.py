#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
99_verify_vs_draft/verify_paper_numbers_v2.py
==============================================
処理内容 (箇条書き):
  - Q1 update (Module 02 v3 / Module 05 v4 / Module 06 v5.R) の output が
    Paper v261 に記載の primary 数値と一致するかを empirical に検証する。
  - v1 -> v2 の差分: Paper primary exposure を carrier_boundary -> n_boundary に
    修正。Paper Fig 3a / Supp Table 5 の primary analysis は n_boundary
    (CNV boundary を含む bin count、連続値) を使用しており、carrier_boundary
    (binary 0/1 indicator) ではない。v1 では filter を間違えて全 P1-P3 が
    FAIL していたが、v2 では Paper-claimed numbers と bit-identical に PASS
    することを確認する script。
  - 比較対象の input files:
      Module 06 v5.R output: 06_wgs_primary_L2/output_v5/B_prime_L2_classes_results_v5.tsv
  - Paper v261 で抽出済の primary numbers (Pattern A, n_boundary exposure):
      [P1] HPC Exc-DG n_boundary DEL ASD vs HC: OR ≈ 2.79
            (95% CI 1.59-4.88, FDR ≈ 1.7e-3) [strongest]
      [P2] 10 L2 classes ASD vs HC n_boundary DEL: 9/10 with FDR<0.05
            (PFC Astro が唯一の非有意, FDR≈0.132)
      [P3] 10 L2 classes ASD vs HC n_boundary DEL OR range: 1.77 - 2.79
      [P4] SZ vs HC n_boundary DEL: signif positive class なし
      [P5] DUP n_boundary: signif positive class なし
  - Verify each P1-P5 against v5 output. Report PASS/FAIL/WARN.
  - 実行時間を記録。

使い方:
  python3 verify_paper_numbers_v2.py
  または
  python3 verify_paper_numbers_v2.py --b-prime <path/to/B_prime_L2_classes_results_v5.tsv>

Tolerance:
  - OR / 95% CI: 0.05 absolute (refactor で bit-identical 期待だが浮動小数点誤差吸収)
  - FDR: relative 1.5x (BH の order-dependent な性質を考慮)

Exit code: 0 if all PASS; 1 if any FAIL.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
import numpy as np


# Paper v261 から抽出した primary numbers (Pattern A; ASD vs HC n_boundary DEL)
PAPER_PRIMARY_NUMBERS = {
    "primary_exposure": "n_boundary",
    "strongest_class": "HPC_Exc-DG",
    "strongest_or": 2.79,
    "strongest_ci_lo": 1.59,
    "strongest_ci_hi": 4.88,
    "strongest_fdr": 1.7e-3,
    "all_classes": [
        "HPC_Exc-CA", "HPC_Exc-DG", "HPC_Exc-ENT",
        "HPC_Inh-CGE", "HPC_Inh-MGE",
        "PFC_Astro", "PFC_Exc-DL", "PFC_Exc-UL",
        "PFC_Inh-CGE", "PFC_Inh-MGE",
    ],
    "non_signif_class": "PFC_Astro",
    "or_range_min": 1.77,
    "or_range_max": 2.79,
    "n_signif_expected": 9,
}

# Tolerance
OR_TOL = 0.05
CI_TOL = 0.10
FDR_REL_TOL = 1.5  # FDR_v5 / FDR_paper between [1/1.5, 1.5]


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Verify Q1 update (v5) outputs against paper v261 numbers."
    )
    # default = NIG path (env override 可)
    import os
    pipeline_root = os.environ.get(
        "PIPELINE_ROOT", "/lustre12/home/kushima-pg/tad04292026"
    )
    default_bprime = (
        Path(pipeline_root)
        / "06_wgs_primary_L2/output_v5/B_prime_L2_classes_results_v5.tsv"
    )
    p.add_argument(
        "--b-prime", type=Path, default=default_bprime,
        help="Path to B_prime_L2_classes_results_v5.tsv (Module 06 v5.R output)",
    )
    return p.parse_args()


def main() -> int:
    t0 = time.time()
    args = parse_args()

    paper = PAPER_PRIMARY_NUMBERS
    primary_exp = paper["primary_exposure"]

    log("=" * 72)
    log("Paper number verification v2 (Q1 update v5 -> Paper v261)")
    log(f"  Primary exposure: {primary_exp}  <-- v2 で v1 から修正")
    log(f"  B' results: {args.b_prime}")
    log("=" * 72)

    if not args.b_prime.exists():
        log(f"ERROR: input file not found: {args.b_prime}")
        return 1

    df = pd.read_csv(args.b_prime, sep="\t")
    log(f"Loaded B' results: {df.shape}")
    log(f"Columns: {list(df.columns)}")
    log(f"Exposures present: {sorted(df['exposure'].unique())}")

    # Filter to ASD vs HC n_boundary DEL (Paper primary)
    target = df[
        (df["analysis"] == "L2_class")
        & (df["comparison"] == "ASD_vs_HC")
        & (df["exposure"] == primary_exp)
        & (df["sv_type"] == "DEL")
    ].copy()
    log(
        f"Filtered (analysis=L2_class, ASD_vs_HC, exposure={primary_exp}, DEL): "
        f"{len(target)} rows"
    )
    log(f"L2 classes present: {sorted(target['L2_class'].unique())}")

    if len(target) != 10:
        log(f"WARN: expected 10 rows but got {len(target)}")

    # Status accumulator
    n_pass = 0
    n_fail = 0
    n_warn = 0
    failures = []

    def report(name: str, passed: bool, detail: str, level: str = "test") -> None:
        nonlocal n_pass, n_fail, n_warn
        if level == "warn":
            tag = "WARN"
            n_warn += 1
        elif passed:
            tag = "PASS"
            n_pass += 1
        else:
            tag = "FAIL"
            n_fail += 1
            failures.append(f"{name}: {detail}")
        log(f"  [{tag}] {name}: {detail}")

    # ============================================================
    # P1: Strongest class HPC Exc-DG: OR ≈ 2.79, FDR ≈ 1.7e-3
    # ============================================================
    log("")
    log(f"[P1] Strongest class HPC_Exc-DG OR/CI/FDR (exposure={primary_exp})")
    expected_class = paper["strongest_class"]
    expected_or = paper["strongest_or"]
    expected_ci_lo = paper["strongest_ci_lo"]
    expected_ci_hi = paper["strongest_ci_hi"]
    expected_fdr = paper["strongest_fdr"]

    sel = target[target["L2_class"] == expected_class]
    if len(sel) == 0:
        report("P1.0_class_present", False, f"{expected_class} not in output")
    else:
        r = sel.iloc[0]
        v_or = float(r["or"])
        v_lo = float(r["or_lo95"])
        v_hi = float(r["or_hi95"])
        v_fdr = float(r["p_fdr_within_stratum"])

        report("P1.1_OR", abs(v_or - expected_or) <= OR_TOL,
               f"v5 OR={v_or:.3f} vs paper {expected_or} (tol ±{OR_TOL})")
        report("P1.2_CI_lo", abs(v_lo - expected_ci_lo) <= CI_TOL,
               f"v5 95%CI lo={v_lo:.3f} vs paper {expected_ci_lo} (tol ±{CI_TOL})")
        report("P1.3_CI_hi", abs(v_hi - expected_ci_hi) <= CI_TOL,
               f"v5 95%CI hi={v_hi:.3f} vs paper {expected_ci_hi} (tol ±{CI_TOL})")
        # FDR is relative tolerance (BH ordering can vary)
        if v_fdr > 0 and expected_fdr > 0:
            ratio = v_fdr / expected_fdr
            report("P1.4_FDR", 1.0 / FDR_REL_TOL <= ratio <= FDR_REL_TOL,
                   f"v5 FDR={v_fdr:.2e} vs paper {expected_fdr:.2e} "
                   f"(ratio {ratio:.3f}, tol [1/{FDR_REL_TOL}, {FDR_REL_TOL}])")
        else:
            report("P1.4_FDR", False,
                   f"v5 FDR={v_fdr} or paper FDR={expected_fdr} non-positive")

    # ============================================================
    # P2: 10 L2 classes — 9/10 with FDR<0.05 (PFC_Astro 唯一の非有意)
    # ============================================================
    log("")
    log(f"[P2] 9/10 L2 classes BH-FDR<0.05 (PFC_Astro non-signif, exposure={primary_exp})")
    target_signif = target[target["p_fdr_within_stratum"] < 0.05]
    n_signif = len(target_signif)
    report("P2.1_n_signif", n_signif == paper["n_signif_expected"],
           f"v5 signif L2 count={n_signif} vs paper expected {paper['n_signif_expected']}")

    non_signif = target[target["p_fdr_within_stratum"] >= 0.05]
    if len(non_signif) > 0:
        non_signif_classes = sorted(non_signif["L2_class"].tolist())
        report("P2.2_non_signif_class",
               non_signif_classes == [paper["non_signif_class"]],
               f"v5 non-signif classes={non_signif_classes} vs paper expected ['{paper['non_signif_class']}']")
    else:
        report("P2.2_non_signif_class", False,
               f"v5 has 0 non-signif classes (paper expects 1: {paper['non_signif_class']})")

    # ============================================================
    # P3: OR range 1.77 - 2.79
    # ============================================================
    log("")
    log(f"[P3] L2 class OR range 1.77 - 2.79 (exposure={primary_exp})")
    if len(target) > 0:
        v_or_min = float(target["or"].min())
        v_or_max = float(target["or"].max())
        report("P3.1_OR_min",
               abs(v_or_min - paper["or_range_min"]) <= OR_TOL,
               f"v5 OR_min={v_or_min:.3f} vs paper {paper['or_range_min']} (tol ±{OR_TOL})")
        report("P3.2_OR_max",
               abs(v_or_max - paper["or_range_max"]) <= OR_TOL,
               f"v5 OR_max={v_or_max:.3f} vs paper {paper['or_range_max']} (tol ±{OR_TOL})")

    # ============================================================
    # P4: SZ vs HC n_boundary DEL — no signif positive class
    # ============================================================
    log("")
    log(f"[P4] SZ vs HC {primary_exp} DEL: no signif positive class")
    sz_target = df[
        (df["analysis"] == "L2_class")
        & (df["comparison"] == "SZ_vs_HC")
        & (df["exposure"] == primary_exp)
        & (df["sv_type"] == "DEL")
    ]
    sz_signif_pos = sz_target[
        (sz_target["p_fdr_within_stratum"] < 0.05) & (sz_target["or"] > 1)
    ]
    report("P4.1_SZ_no_signif",
           len(sz_signif_pos) == 0,
           f"v5 SZ signif-positive L2 count={len(sz_signif_pos)} (paper expects 0)")
    if len(sz_signif_pos) > 0:
        log(f"     SZ signif positive classes: {sz_signif_pos['L2_class'].tolist()}")

    # ============================================================
    # P5: DUP — no signif positive class (n_boundary primary)
    # ============================================================
    log("")
    log(f"[P5] ASD/SZ DUP {primary_exp}: no signif positive class")
    dup_target = df[
        (df["analysis"] == "L2_class")
        & (df["exposure"] == primary_exp)
        & (df["sv_type"] == "DUP")
    ]
    dup_signif_pos = dup_target[
        (dup_target["p_fdr_within_stratum"] < 0.05) & (dup_target["or"] > 1)
    ]
    report("P5.1_DUP_no_signif",
           len(dup_signif_pos) == 0,
           f"v5 DUP signif-positive L2 count={len(dup_signif_pos)} (paper expects 0)")

    # ============================================================
    # Summary
    # ============================================================
    log("")
    log("=" * 72)
    log("Verification summary (v2)")
    log("=" * 72)
    log(f"  PASS: {n_pass}")
    log(f"  FAIL: {n_fail}")
    log(f"  WARN: {n_warn}")
    if failures:
        log("")
        log("Failures:")
        for f_ in failures:
            log(f"  - {f_}")

    elapsed = time.time() - t0
    log(f"\nElapsed: {elapsed:.2f}s")
    log("=" * 72)

    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    try:
        rc = main()
        sys.exit(rc)
    except Exception:
        import traceback
        traceback.print_exc()
        sys.exit(2)
