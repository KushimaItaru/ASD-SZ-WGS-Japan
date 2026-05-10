#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ファイル名: tad_dynamic_boundary_specificity_test_v7.py
# 処理内容:
#   - Heffel の L2 diffbound BED から 10 developmental differential boundary classes を読み込む
#   - raw h5ad から boundary value == 2 かつ support>=2 の static boundary pool を構築し、全 diffbound union を除外する
#   - 各 bin に support count / mean boundary probability / gene density / segdup% を付与する
#   - sample-level の DEL burden について、dynamic bins と matched static bins を比較する
#   - 主解析は sample-level B' logistic regression（ASD vs Healthy）
#   - matched static null は MATCH_LEVEL に応じたマッチング変数で resampling する
#   - class別 + 10-class union の empirical P, null CI, QC を出力する
#   - event-level sensitivity（observed + matched static の OR）を出力する
#   - 実行時間をログ出力する
#
# v2 修正点:
#   1. OUTDIR NameError を修正（_OUTDIR を使用）
#   2. sample_to_idx を analysis_df の行順に合わせて修正
#   3. gene_density / segdup% 計算前に chromosome 内で bin を start 順にソート
#   4. 二重 resampling を解消（matched_sum / matched_exposed_sum を1回目のループで累積）
#   5. burden table の反復読み込みを解消
#   6. observed event-level OR を追加
#   7. matched_exposed_mean を matched_exposed_prob に修正
#
# v3 修正点:
#   8. event_df / event_to_idx の順序ずれを修正（unique_events順で event_df を構築）
#   9. compute_observed_event_or_from_dynamic_bins の overlap table 再読み込みを解消
#      （pre-loaded overlap data を引数で受け取る方式に変更）
#   10. SBATCH ヘッダ追加、N_RESAMPLES デフォルトを 100 に変更（試走用）
#
# v4 修正点:
#   11. MATCH_LEVEL 環境変数を追加（1/2/3 の 3 段階マッチング感度ラダー）
#       - Level 1 (loosest): chromosome only → 距離重みは gene_density + segdup_pct
#       - Level 2 (intermediate): chromosome + gene_bin + segdup_bin → 距離重みは gene_density + segdup_pct
#       - Level 3 (strictest, v3 と同等): chromosome + support_bin + meanprob_bin + gene_bin + segdup_bin
#         → 距離重みは support_count + mean_prob + gene_density + segdup_pct
#   12. OUTDIR / 出力ファイル名に _L{MATCH_LEVEL} サフィックスを追加
#   13. build_matching_bins / compute_candidate_pool_for_dynamic_bin を MATCH_LEVEL 対応に修正
#   14. 距離重み付け変数を MATCH_LEVEL に応じて動的に選択
#
# v5 修正点:
#   15. static pool 構築で (X > 0) → (X == 2) に変更
#       - h5ad の X matrix は TopDom 3値分類: 0=domain, 1=gap, 2=boundary
#       - (X > 0) だと gap (value=1) も含まれていた → boundary (value=2) のみに限定
#       - 原稿 Methods 記載「scHiCluster/TopDom boundary value of 2」と整合
#   16. dynamic_union の observed exposure を overlap table から 10-class union bins ベースで直接計算
#       - burden table の n_boundary_dev_union_DEL（11クラスunion）は使用しない → 10/11クラス不整合を解消
#       - sig8 union は DISCOVERY_SIGNIFICANT_CLASSES として定義（将来の sensitivity 用）
#   17. デフォルトを MATCH_LEVEL=2, N_RESAMPLES=1000 に変更（引数なしで本番実行可能に）
#   18. 出力ディレクトリ名を v5 に更新
#
# v6 修正点 (tad04212026 pipeline 移行):
#   19. common.paths_v1 から中央管理パスを取得（_BASEDIR / _L2_DIFFBOUND_DIR / _RAW_H5AD /
#       _OVERLAP_TABLE / _BURDEN_TABLE / _GTF_FILE / _SEGDUP_BED / _OUTDIR をすべて paths_v1
#       の定数に置換）
#   20. burden table の読み込み直後に列名 normalize を追加:
#       `n_boundary_HPC_Exc-DG_DEL` (Step 5 CamelCase-dash) →
#       `n_boundary_hpc_exc_dg_DEL` (snake_case, PRIMARY_BOUNDARY_CLASSES と整合)
#       対象プレフィクス: n_boundary_, n_events_, carrier_boundary_
#   21. 出力ファイル名サフィックスを v5 → v6 に更新（_L{MATCH_LEVEL} は維持）
#   22. 出力先は paths_v1.OUT_07_MATCHED_STATIC に統一
#
# v7 修正点 (WGS top-1% CNV count sample-level QC 伝播; path override only):
#   23. _OVERLAP_TABLE: paths_v1.F_04_EVENT_OVERLAP (output_v9) ->
#       hardcoded /lustre12/.../04_wgs_sv_boundary_overlap/output_v10/
#       sample_boundary_event_overlap_v10.tsv.gz
#       (Step 04 v10 で sample-level top-1% CNV count QC を適用済み)
#   24. _BURDEN_TABLE: paths_v1.F_05_SAMPLE_BURDEN_L2 (output_v2) ->
#       hardcoded /lustre12/.../05_wgs_sample_burden/output_v3/
#       sample_burden_L2_and_specificity_v3.tsv
#       (Step 05 v3 で post-QC データに基づき再集計済み)
#   25. paths_v1.py 自体は変更しない (中央管理ファイルへの侵襲を避ける)。
#       Step 04/05 の paths_v1 定数は将来的に bump する余地を残す。
#   26. 出力ファイル名サフィックス: _v6_L{MATCH_LEVEL} -> _v7_L{MATCH_LEVEL}
#   27. 出力先: output_tad_dynamic_boundary_specificity_v7_L{MATCH_LEVEL} に変更。
#   28. 解析ロジック (matched-static / dynamic boundary / Logit / event-level OR /
#       flanking control / class-wise + 10-class union, MATCH_LEVEL ladder) は
#       v6 と完全に同一。変更は path / 出力名のみ。
#   29. arrayCGH / MSSNG / GENCODE / SEGDUP / Heffel boundary などの非 WGS 入力は
#       v6 と完全同一 (paths_v1 経由のまま)。
#
# 実行例（デフォルト: Level 2, 1000 resamples）:
#   python tad_dynamic_boundary_specificity_test_v7.py
#
# 実行例（Level 3, 10000 resamples）:
#   MATCH_LEVEL=3 N_RESAMPLES=10000 python tad_dynamic_boundary_specificity_test_v7.py
#
# 実行例（Level 1, 探索）:
#   MATCH_LEVEL=1 N_RESAMPLES=1000 python tad_dynamic_boundary_specificity_test_v7.py
#
# 環境変数:
#   MATCH_LEVEL=2     (1=chrom only, 2=chrom+gene+segdup, 3=all 5 vars)  ← v5: デフォルト変更
#   TARGET_CLASS=      (空=全クラス実行, 指定=1クラスのみ。並列実行用)
#   OUTDIR=/path/to/output_dir
#   N_RESAMPLES=1000  ← v5: デフォルト変更
#   RANDOM_SEED=20260402
#   RUN_EVENT_LEVEL=1
#   RUN_FLANKING_CONTROL=0
#
# --- SBATCH (sbatch wrapper 用) ---
# #SBATCH -p ncbn-cpu
# #SBATCH --account=ncbn-cpu
# #SBATCH --cpus-per-task=8
# #SBATCH --mem=64G
# #SBATCH --time=24:00:00
# #SBATCH --output=specificity_test_v7_L%e_%j.log

import os
import re
import sys
import json
import time
import gzip
import math
import warnings
import traceback
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import statsmodels.api as sm

# v6: tad04212026 pipeline - centralized paths via common.paths_v1
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths_v1 import (
    HEFFEL_L2_DIFF_DOMAIN_BOUNDARIES,
    HEFFEL_DOMAIN_BOUNDARIES,
    F_04_EVENT_OVERLAP,
    F_05_SAMPLE_BURDEN_L2,
    GENCODE_GTF,
    SEGDUP_BED,
    OUT_07_MATCHED_STATIC,
    ensure_output_dirs,
)
# v6 (2026-04-21): common.naming_v1 で列名正規化（SV suffix 保持バグ修正）。
from common.naming_v1 import normalize_l2_burden_columns
# ensure_output_dirs() は main() から呼ぶ（import 時の副作用を避ける）。

warnings.filterwarnings("ignore")

# ============================================================
# CONFIG
# ============================================================

_L2_DIFFBOUND_DIR = str(HEFFEL_L2_DIFF_DOMAIN_BOUNDARIES)
_RAW_H5AD = str(HEFFEL_DOMAIN_BOUNDARIES / "BrainDev_raw.boundary.h5ad")
_OVERLAP_TABLE = "/lustre12/home/kushima-pg/tad04212026/04_wgs_sv_boundary_overlap/output_v10/sample_boundary_event_overlap_v10.tsv.gz"  # v7: top-1% CNV count QC 適用後 (paths_v1 を override)
_BURDEN_TABLE = "/lustre12/home/kushima-pg/tad04212026/05_wgs_sample_burden/output_v3/sample_burden_L2_and_specificity_v3.tsv"  # v7: top-1% CNV count QC 適用後 (paths_v1 を override)
_GTF_FILE = str(GENCODE_GTF)
_SEGDUP_BED = str(SEGDUP_BED)

# ---- v4: MATCH_LEVEL ----
MATCH_LEVEL = int(os.environ.get("MATCH_LEVEL", 2))  # v5: デフォルト 3→2
assert MATCH_LEVEL in (1, 2, 3), f"MATCH_LEVEL must be 1, 2, or 3 (got {MATCH_LEVEL})"

# v6: default OUTDIR = paths_v1.OUT_07_MATCHED_STATIC / match_level suffix
_OUTDIR = os.environ.get(
    "OUTDIR",
    str(OUT_07_MATCHED_STATIC / f"output_tad_dynamic_boundary_specificity_v7_L{MATCH_LEVEL}")
)

# 主解析
ANALYSIS_DIAGNOSIS_CASE = "ASD"
ANALYSIS_DIAGNOSIS_CONTROL = "Healthy"
SV_TYPE = "DEL"

# 10クラス（HPC_Astroは除外）
PRIMARY_BOUNDARY_CLASSES = [
    "hpc_exc_ca", "hpc_exc_dg", "hpc_exc_ent",
    "hpc_inh_cge", "hpc_inh_mge",
    "pfc_astro", "pfc_exc_dl", "pfc_exc_ul",
    "pfc_inh_cge", "pfc_inh_mge",
]

# v5: discovery で有意だった8クラス（dynamic_union の定義に使用）
# hpc_exc_ca と pfc_astro は discovery で non-significant → union から除外
DISCOVERY_SIGNIFICANT_CLASSES = [
    "hpc_exc_dg", "hpc_exc_ent",
    "hpc_inh_cge", "hpc_inh_mge",
    "pfc_exc_dl", "pfc_exc_ul",
    "pfc_inh_cge", "pfc_inh_mge",
]

# static pool = boundary value==2 かつ support>=2 かつ all diffbound union 除外
MIN_STATIC_SUPPORT = 2

# matching / resampling
N_RESAMPLES = int(os.environ.get("N_RESAMPLES", 1000))  # v5: デフォルト 100→1000
RANDOM_SEED = int(os.environ.get("RANDOM_SEED", 20260402))
N_PROB_BINS = int(os.environ.get("N_PROB_BINS", 5))
N_GENE_BINS = int(os.environ.get("N_GENE_BINS", 5))
N_SEGDUP_BINS = int(os.environ.get("N_SEGDUP_BINS", 3))
MIN_EXACT_POOL = int(os.environ.get("MIN_EXACT_POOL", 10))

# optional analyses
RUN_EVENT_LEVEL = os.environ.get("RUN_EVENT_LEVEL", "1") == "1"
RUN_FLANKING_CONTROL = os.environ.get("RUN_FLANKING_CONTROL", "0") == "1"

# ---- v4.1: TARGET_CLASS (並列実行用) ----
# 未設定 or 空文字 → 全クラスを実行（従来通り）
# 設定時 → 指定クラスのみ実行。出力ファイルにクラス名を付加
TARGET_CLASS = os.environ.get("TARGET_CLASS", "").strip()

# ---- v4: MATCH_LEVEL 別マッチング変数定義 ----
# Level 1: chromosome only → 距離重みは gene_density, segdup_pct
# Level 2: chromosome + gene_bin + segdup_bin → 距離重みは gene_density, segdup_pct
# Level 3: chromosome + support_bin + meanprob_bin + gene_bin + segdup_bin → 距離重みは全4変数
MATCH_LEVEL_CONFIG = {
    1: {
        "exact_key_vars": [],  # chrom is always included
        "distance_vars": ["gene_density", "segdup_pct"],
        "label": "chrom_only",
        "description": "Chromosome only (loosest matching)",
    },
    2: {
        "exact_key_vars": ["gene_bin", "segdup_bin"],
        "distance_vars": ["gene_density", "segdup_pct"],
        "label": "chrom_gene_segdup",
        "description": "Chromosome + gene density + segdup% (intermediate)",
    },
    3: {
        "exact_key_vars": ["support_bin", "meanprob_bin", "gene_bin", "segdup_bin"],
        "distance_vars": ["support_count", "mean_prob", "gene_density", "segdup_pct"],
        "label": "all_5_vars",
        "description": "Chromosome + support + mean_prob + gene density + segdup% (strictest)",
    },
}

_ML_CFG = MATCH_LEVEL_CONFIG[MATCH_LEVEL]

# output file names (v4: include L{MATCH_LEVEL}, and optionally TARGET_CLASS)
_CLASS_SUFFIX = f"_{TARGET_CLASS}" if TARGET_CLASS else ""
_MAIN_TSV = f"specificity_test_results_main_v7_L{MATCH_LEVEL}{_CLASS_SUFFIX}.tsv"
_PER_SAMPLE_TSV = f"specificity_test_per_sample_v7_L{MATCH_LEVEL}{_CLASS_SUFFIX}.tsv"
_FIGURE_SUMMARY_TSV = f"specificity_test_figure_summary_v7_L{MATCH_LEVEL}{_CLASS_SUFFIX}.tsv"
_QC_JSON = f"specificity_test_qc_v7_L{MATCH_LEVEL}{_CLASS_SUFFIX}.json"
_EVENT_TSV = f"specificity_test_event_level_v7_L{MATCH_LEVEL}{_CLASS_SUFFIX}.tsv"

# ============================================================
# LOGGING
# ============================================================

_T0 = time.time()


def log(msg: str) -> None:
    elapsed = time.time() - _T0
    print(f"[{elapsed:8.1f}s] {msg}", flush=True)


# ============================================================
# UTIL
# ============================================================

def norm_chrom(c) -> str:
    s = str(c).replace("chr", "")
    if s.isdigit():
        return f"chr{int(s)}"
    if s.upper() in ("X", "Y", "M", "MT"):
        return f"chr{s.upper()}"
    return f"chr{s}"


def detect_column(df: pd.DataFrame, candidates: List[str],
                  required: bool = True) -> Optional[str]:
    col_map = {str(c).lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in col_map:
            return col_map[cand.lower()]
    if required:
        raise KeyError(
            f"Missing column. candidates={candidates}, "
            f"actual={list(df.columns)}"
        )
    return None


def assert_unique_keys(df: pd.DataFrame, key_cols: List[str],
                       label: str, n_show: int = 10) -> None:
    dup_mask = df.duplicated(subset=key_cols, keep=False)
    if dup_mask.any():
        n_dup = int(dup_mask.sum())
        ex = df.loc[dup_mask, key_cols].head(n_show).to_string(index=False)
        raise ValueError(f"[{label}] duplicate rows found: n={n_dup}\n{ex}")


def make_bin_id(chrom: pd.Series, start0: pd.Series,
                end: pd.Series) -> pd.Series:
    return chrom.astype(str) + ":" + start0.astype(str) + "-" + end.astype(str)


def safe_qcut_edges(x: pd.Series, q: int) -> np.ndarray:
    vals = pd.to_numeric(x, errors="coerce")
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return np.array([0.0, 1.0])
    probs = np.linspace(0.0, 1.0, q + 1)
    edges = np.quantile(vals, probs)
    edges = np.unique(edges)
    if len(edges) < 2:
        v = float(vals.iloc[0] if hasattr(vals, "iloc") else vals[0])
        return np.array([v - 1e-9, v + 1e-9])
    return edges


def make_quantile_bins_from_edges(x: pd.Series,
                                  edges: np.ndarray) -> pd.Series:
    vals = pd.to_numeric(x, errors="coerce")
    idx = np.digitize(vals, edges[1:-1], right=False)
    return pd.Series(idx, index=x.index)


def merge_intervals(
    intervals: List[Tuple[int, int]]
) -> List[Tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(int(s), int(e)) for s, e in merged]


def odds_ratio_2x2(a: int, b: int, c: int, d: int) -> float:
    """Haldane-Anscombe correction."""
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))


# ============================================================
# LOAD DIFFBOUND BINS
# ============================================================

def infer_lineage_key_from_filename(path: Path) -> str:
    name = path.name.replace("_diffbound.bed.gz", "")
    return name.lower().replace("-", "_")


def load_diffbound_bins(
    diff_dir: str,
) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Returns (class_bins, analysis_union, all_diff_union, sig8_union)."""
    log("Loading differential boundary bins ...")
    diff_files = sorted(Path(diff_dir).glob("*_diffbound.bed.gz"))
    if len(diff_files) == 0:
        raise FileNotFoundError(f"No diffbound files found in: {diff_dir}")

    class_bins: Dict[str, pd.DataFrame] = {}
    all_bins_list: List[pd.DataFrame] = []
    analysis_bins_list: List[pd.DataFrame] = []
    sig8_bins_list: List[pd.DataFrame] = []  # v5: sig8 union 用

    for f in diff_files:
        lineage_key = infer_lineage_key_from_filename(f)
        df = pd.read_csv(
            f, sep="\t", compression="gzip",
            usecols=["chrom", "start", "end", "chi2_sc"]
        )
        df["chrom"] = df["chrom"].map(norm_chrom)
        df["start0"] = pd.to_numeric(df["start"], errors="coerce").astype(int)
        df["end"] = pd.to_numeric(df["end"], errors="coerce").astype(int)
        df["bin_id"] = make_bin_id(df["chrom"], df["start0"], df["end"])
        df["lineage_key"] = lineage_key
        df = df[
            ["chrom", "start0", "end", "bin_id", "chi2_sc", "lineage_key"]
        ].drop_duplicates("bin_id").copy()

        all_bins_list.append(df.copy())
        if lineage_key in PRIMARY_BOUNDARY_CLASSES:
            class_bins[lineage_key] = df.copy()
            analysis_bins_list.append(df.copy())
        # v5: sig8 union 構築
        if lineage_key in DISCOVERY_SIGNIFICANT_CLASSES:
            sig8_bins_list.append(df.copy())

        log(f"  {lineage_key}: {df['bin_id'].nunique()} bins")

    all_diff_union = pd.concat(all_bins_list, axis=0, ignore_index=True)
    all_diff_union = all_diff_union[
        ["chrom", "start0", "end", "bin_id"]
    ].drop_duplicates("bin_id").copy()

    analysis_union = pd.concat(analysis_bins_list, axis=0, ignore_index=True)
    analysis_union = analysis_union[
        ["chrom", "start0", "end", "bin_id"]
    ].drop_duplicates("bin_id").copy()

    # v5: sig8 union（discovery有意8クラスの union）
    sig8_union = pd.concat(sig8_bins_list, axis=0, ignore_index=True)
    sig8_union = sig8_union[
        ["chrom", "start0", "end", "bin_id"]
    ].drop_duplicates("bin_id").copy()

    log(f"  all diffbound union bins (11 files): "
        f"{all_diff_union['bin_id'].nunique()}")
    log(f"  analysis diffbound union bins (10 classes): "
        f"{analysis_union['bin_id'].nunique()}")
    log(f"  sig8 diffbound union bins (8 classes): "
        f"{sig8_union['bin_id'].nunique()}")

    return class_bins, analysis_union, all_diff_union, sig8_union


# ============================================================
# LOAD H5AD BOUNDARY METRICS
# ============================================================

def load_h5ad_support_metrics(raw_h5ad: str) -> pd.DataFrame:
    log("Loading raw h5ad support metrics ...")
    import anndata as ad
    from scipy import sparse as sp

    adata = ad.read_h5ad(raw_h5ad)
    var = adata.var.reset_index(drop=True).copy()
    for col in ["chrom", "start", "end"]:
        if col not in var.columns:
            raise ValueError(
                f"raw h5ad var is missing required column: {col}"
            )

    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    else:
        X = X.tocsr()

    # v5: TopDom 3値分類 (0=domain, 1=gap, 2=boundary)
    # boundary (value==2) のみをカウント（v4 では X>0 で gap も含んでいた）
    support_count = np.asarray(
        (X == 2).sum(axis=0)
    ).ravel().astype(np.int32)
    mean_prob = np.asarray(X.mean(axis=0)).ravel().astype(np.float64)
    max_prob = np.asarray(
        X.max(axis=0).toarray()
    ).ravel().astype(np.float64)

    out = pd.DataFrame({
        "chrom": var["chrom"].map(norm_chrom).astype(str),
        "start0": pd.to_numeric(var["start"], errors="coerce").astype(int),
        "end": pd.to_numeric(var["end"], errors="coerce").astype(int),
        "support_count": support_count,
        "mean_prob": mean_prob,
        "max_prob": max_prob,
    })
    out["bin_id"] = make_bin_id(out["chrom"], out["start0"], out["end"])
    out = out[
        ["chrom", "start0", "end", "bin_id",
         "support_count", "mean_prob", "max_prob"]
    ].copy()
    assert_unique_keys(out, ["bin_id"], "h5ad_bin_metrics")
    log(f"  h5ad bins: {len(out)}")
    return out


# ============================================================
# GENE DENSITY / SEGDUP%
# ============================================================

def build_gene_index(
    gtf_path: str,
) -> Dict[str, List[Tuple[int, int]]]:
    log("Building gene interval index ...")
    gene_intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "gene":
                continue
            attrs = fields[8]
            if ('gene_type "protein_coding"' not in attrs
                    and 'gene_biotype "protein_coding"' not in attrs):
                continue
            chrom = norm_chrom(fields[0])
            try:
                s1 = int(fields[3])
                e1 = int(fields[4])
            except Exception:
                continue
            if e1 < s1:
                s1, e1 = e1, s1
            s0 = max(0, s1 - 1)  # GTF 1-based closed -> 0-based half-open
            if e1 <= s0:
                continue
            gene_intervals[chrom].append((s0, e1))

    for chrom in list(gene_intervals.keys()):
        gene_intervals[chrom].sort()

    n_total = sum(len(v) for v in gene_intervals.values())
    log(f"  protein-coding genes indexed: {n_total}")
    return gene_intervals


def annotate_gene_density(
    bin_df: pd.DataFrame,
    gene_index: Dict[str, List[Tuple[int, int]]],
) -> pd.DataFrame:
    log("Annotating gene density per bin ...")
    out = bin_df.copy()
    gene_density = np.zeros(len(out), dtype=np.int32)

    chrom_to_rows: Dict[str, List[Tuple[int, int, int]]] = defaultdict(list)
    for i, (chrom, s, e) in enumerate(
        zip(out["chrom"], out["start0"], out["end"])
    ):
        chrom_to_rows[chrom].append((i, int(s), int(e)))

    for chrom, rows in chrom_to_rows.items():
        rows.sort(key=lambda x: x[1])  # v2 fix: start 順にソート
        genes = gene_index.get(chrom, [])
        if len(genes) == 0:
            continue
        gi = 0
        for i, s, e in rows:
            while gi < len(genes) and genes[gi][1] <= s:
                gi += 1
            gj = gi
            c = 0
            while gj < len(genes) and genes[gj][0] < e:
                if genes[gj][1] > s:
                    c += 1
                gj += 1
            gene_density[i] = c

    out["gene_density"] = gene_density
    return out


def load_segdup_merged(
    segdup_bed: str,
) -> Dict[str, List[Tuple[int, int]]]:
    log("Loading and merging segdup BED ...")
    intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    with open(segdup_bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                chrom = norm_chrom(parts[0])
                s = int(parts[1])
                e = int(parts[2])
            except Exception:
                continue
            if e > s:
                intervals[chrom].append((s, e))

    merged = {}
    for chrom, ints in intervals.items():
        merged[chrom] = merge_intervals(ints)
    return merged


def compute_segdup_pct_for_bins(
    bin_df: pd.DataFrame,
    segdup_merged: Dict[str, List[Tuple[int, int]]],
) -> pd.DataFrame:
    log("Annotating segdup% per bin ...")
    out = bin_df.copy()
    segdup_pct = np.zeros(len(out), dtype=np.float64)

    chrom_to_rows: Dict[str, List[Tuple[int, int, int]]] = defaultdict(list)
    for i, (chrom, s, e) in enumerate(
        zip(out["chrom"], out["start0"], out["end"])
    ):
        chrom_to_rows[chrom].append((i, int(s), int(e)))

    for chrom, rows in chrom_to_rows.items():
        rows.sort(key=lambda x: x[1])  # v2 fix: start 順にソート
        ints = segdup_merged.get(chrom, [])
        if len(ints) == 0:
            continue
        j = 0
        for i, s, e in rows:
            bin_len = e - s
            if bin_len <= 0:
                continue
            while j < len(ints) and ints[j][1] <= s:
                j += 1
            k = j
            cov = 0
            while k < len(ints) and ints[k][0] < e:
                ov = max(0, min(e, ints[k][1]) - max(s, ints[k][0]))
                cov += ov
                k += 1
            segdup_pct[i] = 100.0 * cov / bin_len

    out["segdup_pct"] = segdup_pct
    return out


# ============================================================
# BUILD STATIC POOL AND BIN METRICS
# ============================================================

def build_static_pool(
    h5ad_metrics: pd.DataFrame,
    all_diff_union: pd.DataFrame,
) -> pd.DataFrame:
    log("Building static boundary pool ...")
    diff_ids = set(all_diff_union["bin_id"].tolist())
    static_df = h5ad_metrics.loc[
        (h5ad_metrics["support_count"] >= MIN_STATIC_SUPPORT)
        & (~h5ad_metrics["bin_id"].isin(diff_ids))
    ].copy()
    static_df["is_static_pool"] = 1
    log(f"  static pool bins: {len(static_df)}")
    return static_df


def annotate_bin_metrics_for_dynamic_and_static(
    class_bins: Dict[str, pd.DataFrame],
    analysis_union: pd.DataFrame,
    static_df: pd.DataFrame,
    h5ad_metrics: pd.DataFrame,
    gene_index: Dict[str, List[Tuple[int, int]]],
    segdup_merged: Dict[str, List[Tuple[int, int]]],
) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame]:
    """v5: dynamic_union = 10-class union (primary)。observed exposure は overlap table から直接計算。"""
    log("Annotating dynamic/static bins with metrics ...")

    dynamic_sets: Dict[str, pd.DataFrame] = {}
    for cls, df in class_bins.items():
        tmp = df.merge(
            h5ad_metrics[["bin_id", "support_count", "mean_prob", "max_prob"]],
            on="bin_id", how="left", validate="1:1",
        )
        if tmp["support_count"].isna().any():
            raise ValueError(
                f"Missing h5ad metrics for dynamic class {cls}"
            )
        dynamic_sets[cls] = tmp.copy()

    # v5: dynamic_union = 10-class union (primary)
    # observed exposure は overlap table から直接計算（burden table の 11-class 列は不使用）
    union_df = analysis_union.merge(
        h5ad_metrics[["bin_id", "support_count", "mean_prob", "max_prob"]],
        on="bin_id", how="left", validate="1:1",
    )
    if union_df["support_count"].isna().any():
        raise ValueError("Missing h5ad metrics for analysis union bins")
    dynamic_sets["dynamic_union"] = union_df.copy()

    all_dynamic = pd.concat(
        [df.assign(set_name=cls) for cls, df in dynamic_sets.items()],
        axis=0, ignore_index=True,
    )
    all_dynamic = annotate_gene_density(all_dynamic, gene_index)
    all_dynamic = compute_segdup_pct_for_bins(all_dynamic, segdup_merged)

    updated_dynamic_sets: Dict[str, pd.DataFrame] = {}
    for cls in dynamic_sets.keys():
        updated_dynamic_sets[cls] = (
            all_dynamic.loc[all_dynamic["set_name"] == cls]
            .drop(columns=["set_name"])
            .copy()
        )

    static_annot = annotate_gene_density(static_df.copy(), gene_index)
    static_annot = compute_segdup_pct_for_bins(static_annot, segdup_merged)

    return (
        updated_dynamic_sets,
        updated_dynamic_sets["dynamic_union"],
        static_annot,
    )


# ============================================================
# BURDEN TABLE / OBSERVED MODEL
# ============================================================

def load_burden_analysis_frame(
    burden_path: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (analysis_df, full_burden_df)."""
    log("Loading burden table ...")
    df = pd.read_csv(burden_path, sep="\t", low_memory=False)
    # v6 (2026-04-21): use common.naming_v1.normalize_l2_burden_columns which
    #   (a) lowercases the L2-class body only (preserves SV suffix _DEL/_DUP uppercase),
    #   (b) skips group-scheme columns (the body contains "__"),
    #   (c) skips unknown L2 classes — so PRIMARY_BOUNDARY_CLASSES keys match
    #       get_observed_exposure_column() which asks for `n_boundary_{class}_DEL`.
    df, _rename_map = normalize_l2_burden_columns(df)
    if _rename_map:
        log(f"  column normalize: renamed {len(_rename_map)} L2 burden columns "
            f"(e.g. '{next(iter(_rename_map))}' -> '{_rename_map[next(iter(_rename_map))]}')")

    dx_col = detect_column(df, ["Diagnosis"])
    sex_col = detect_column(df, ["Sex_numeric"])
    total_bases_col = detect_column(df, ["log1p_total_del_bases"])
    total_gene_col = detect_column(df, ["log1p_total_gene_DEL"])

    pc_cols = sorted(
        [c for c in df.columns if re.fullmatch(r"PC\d+", str(c))],
        key=lambda x: int(str(x).replace("PC", "")),
    )[:10]
    if len(pc_cols) < 10:
        raise ValueError(
            f"Need PC1-PC10 in burden table, found only {pc_cols}"
        )

    needed = (
        ["sample_id", dx_col, sex_col, total_bases_col, total_gene_col]
        + pc_cols
    )
    out = df[needed].copy()
    out = out.rename(columns={
        dx_col: "Diagnosis",
        sex_col: "Sex_numeric",
        total_bases_col: "log1p_total_del_bases",
        total_gene_col: "log1p_total_gene_DEL",
    })

    for c in (
        ["Sex_numeric", "log1p_total_del_bases", "log1p_total_gene_DEL"]
        + pc_cols
    ):
        out[c] = pd.to_numeric(out[c], errors="coerce")

    out = out.loc[
        out["Diagnosis"].isin(
            [ANALYSIS_DIAGNOSIS_CASE, ANALYSIS_DIAGNOSIS_CONTROL]
        )
    ].copy()
    out["is_case"] = (
        out["Diagnosis"] == ANALYSIS_DIAGNOSIS_CASE
    ).astype(int)
    out = out.dropna(
        subset=(
            ["Sex_numeric", "log1p_total_del_bases", "log1p_total_gene_DEL"]
            + pc_cols
        )
    ).copy()

    assert_unique_keys(out, ["sample_id"], "analysis_burden_df")
    n_case = int(out["is_case"].sum())
    n_ctrl = int((1 - out["is_case"]).sum())
    log(f"  analysis samples: {len(out)} "
        f"({ANALYSIS_DIAGNOSIS_CASE}={n_case}, "
        f"{ANALYSIS_DIAGNOSIS_CONTROL}={n_ctrl})")
    return out, df


def get_observed_exposure_column(class_name: str) -> Optional[str]:
    # v5: dynamic_union は overlap table から直接計算するため None を返す
    if class_name == "dynamic_union":
        return None
    return f"n_boundary_{class_name}_DEL"


def fit_bprime_logit(
    analysis_df: pd.DataFrame,
    exposure: np.ndarray,
) -> Dict[str, float]:
    """Model B': Logit(is_case) ~ exposure + Sex + log1p_bases + log1p_gene + PC1-10."""
    tmp = analysis_df.copy()
    tmp["exposure"] = pd.to_numeric(
        pd.Series(exposure, index=tmp.index), errors="coerce"
    )
    tmp = tmp.dropna(subset=["exposure"]).copy()

    y = tmp["is_case"].astype(int).to_numpy()
    X_cols = (
        ["exposure", "Sex_numeric",
         "log1p_total_del_bases", "log1p_total_gene_DEL"]
        + [c for c in tmp.columns if re.fullmatch(r"PC\d+", str(c))]
    )
    X = tmp[X_cols].astype(float)
    X = sm.add_constant(X, has_constant="add")

    res = {
        "n_complete": int(len(tmp)),
        "n_case": int(tmp["is_case"].sum()),
        "n_control": int(len(tmp) - tmp["is_case"].sum()),
        "n_exposed_case": int(
            np.sum((tmp["exposure"] > 0) & (tmp["is_case"] == 1))
        ),
        "n_exposed_control": int(
            np.sum((tmp["exposure"] > 0) & (tmp["is_case"] == 0))
        ),
        "case_mean_burden": float(
            tmp.loc[tmp["is_case"] == 1, "exposure"].mean()
        ),
        "control_mean_burden": float(
            tmp.loc[tmp["is_case"] == 0, "exposure"].mean()
        ),
        "beta": np.nan,
        "se": np.nan,
        "or": np.nan,
        "ci_lo": np.nan,
        "ci_hi": np.nan,
        "p_value": np.nan,
        "converged": False,
        "note": "",
    }

    if tmp["exposure"].nunique() < 2 or tmp["is_case"].nunique() < 2:
        res["note"] = "No variation in exposure/outcome"
        return res

    for method in ["newton", "lbfgs", "bfgs"]:
        try:
            fit = sm.Logit(y, X).fit(disp=0, maxiter=200, method=method)
            if np.isfinite(fit.params["exposure"]):
                beta = float(fit.params["exposure"])
                se = float(fit.bse["exposure"])
                p = float(fit.pvalues["exposure"])
                res["beta"] = beta
                res["se"] = se
                res["or"] = float(np.exp(beta))
                res["ci_lo"] = float(np.exp(beta - 1.96 * se))
                res["ci_hi"] = float(np.exp(beta + 1.96 * se))
                res["p_value"] = p
                res["converged"] = bool(
                    fit.mle_retvals.get("converged", False)
                )
                return res
        except Exception:
            continue

    res["note"] = "All optimizers failed"
    return res


# ============================================================
# OVERLAP TABLE / STATIC & EVENT MAPS
# ============================================================

def load_overlap_and_build_maps(
    overlap_path: str,
    analysis_samples: List[str],
) -> Tuple[
    Dict[str, np.ndarray],   # static_bin_to_sample_idxs
    Dict[str, np.ndarray],   # static_bin_to_event_idxs
    pd.DataFrame,            # event_df (aligned with event_to_idx)
    pd.DataFrame,            # ov_filtered (pre-loaded, reusable)
]:
    """
    Load overlap table once. Build:
    - static bin -> sample index map
    - static bin -> event index map
    - event_df aligned with event_to_idx ordering (v3 fix)
    - filtered overlap table for reuse by observed event OR
    """
    log("Loading overlap table and building maps ...")
    ov = pd.read_csv(
        overlap_path, sep="\t", compression="gzip", low_memory=False,
    )

    sample_col = detect_column(ov, ["sample_id"])
    event_col = detect_column(ov, ["event_id"])
    svtype_col = detect_column(ov, ["sv_type_norm"])
    bin_col = detect_column(ov, ["bin_id"])
    dx_col = detect_column(ov, ["Diagnosis"])
    group_col = detect_column(ov, ["group_key"], required=False)
    overlaps_diff_col = detect_column(ov, ["overlaps_diffbound"])

    # DEL + analysis samples のみに絞り込み
    ov_filtered = ov.loc[
        (ov[svtype_col].astype(str) == SV_TYPE)
        & (ov[sample_col].isin(analysis_samples))
        & (ov[dx_col].isin(
            [ANALYSIS_DIAGNOSIS_CASE, ANALYSIS_DIAGNOSIS_CONTROL]
        ))
    ].copy()

    # static overlap records
    if (group_col is not None
            and "static_raw_support2"
            in ov_filtered[group_col].astype(str).unique()):
        ov_static = ov_filtered.loc[
            ov_filtered[group_col].astype(str) == "static_raw_support2"
        ].copy()
    else:
        ov_static = ov_filtered.loc[
            ov_filtered[overlaps_diff_col] == 0
        ].copy()

    # v2 fix: analysis_df の行順に合わせる（sorted なし）
    sample_to_idx = {sid: i for i, sid in enumerate(analysis_samples)}

    # event index: sorted unique events
    unique_events = sorted(
        ov_filtered[[event_col]].drop_duplicates()[event_col].tolist()
    )
    event_to_idx = {eid: i for i, eid in enumerate(unique_events)}

    # v3 fix: event_df を unique_events と同じ順序で構築
    event_map = (
        ov_filtered[[event_col, dx_col]]
        .drop_duplicates(subset=[event_col])
        .rename(columns={event_col: "event_id", dx_col: "Diagnosis"})
    )
    event_df = (
        pd.DataFrame({"event_id": unique_events})
        .merge(event_map, on="event_id", how="left", validate="1:1")
    )
    event_df["is_case"] = (
        event_df["Diagnosis"] == ANALYSIS_DIAGNOSIS_CASE
    ).astype(int)

    # static bin -> sample/event indices
    static_bin_to_sample_idxs: Dict[str, np.ndarray] = {}
    static_bin_to_event_idxs: Dict[str, np.ndarray] = {}

    for bid, grp in ov_static.groupby(bin_col):
        sidx = np.array(
            sorted({
                sample_to_idx[s]
                for s in grp[sample_col].tolist()
                if s in sample_to_idx
            }),
            dtype=np.int32,
        )
        eidx = np.array(
            sorted({
                event_to_idx[e]
                for e in grp[event_col].tolist()
                if e in event_to_idx
            }),
            dtype=np.int32,
        )
        static_bin_to_sample_idxs[str(bid)] = sidx
        static_bin_to_event_idxs[str(bid)] = eidx

    log(f"  static bins with >=1 DEL overlap: "
        f"{len(static_bin_to_sample_idxs)}")
    log(f"  unique DEL events in analysis set: {len(event_df)}")
    log(f"  overlap table filtered rows (reusable): "
        f"{len(ov_filtered)}")

    return (
        static_bin_to_sample_idxs,
        static_bin_to_event_idxs,
        event_df,
        ov_filtered,
    )


def compute_observed_event_or(
    dynamic_bin_ids: set,
    ov_filtered: pd.DataFrame,
    event_df: pd.DataFrame,
) -> Dict[str, float]:
    """
    v3: pre-loaded ov_filtered を使い、overlap table の再読み込みを回避。
    event_df は load_overlap_and_build_maps で構築済みの全 event 一覧。
    """
    bin_col = detect_column(ov_filtered, ["bin_id"])
    event_col = detect_column(ov_filtered, ["event_id"])

    exposed_events = set(
        ov_filtered.loc[
            ov_filtered[bin_col].isin(dynamic_bin_ids),
            event_col,
        ].tolist()
    )

    exposed = event_df["event_id"].isin(exposed_events).astype(int).to_numpy()
    case = event_df["is_case"].to_numpy()

    a = int(np.sum((case == 1) & (exposed == 1)))
    b = int(np.sum((case == 1) & (exposed == 0)))
    c = int(np.sum((case == 0) & (exposed == 1)))
    d = int(np.sum((case == 0) & (exposed == 0)))

    return {
        "observed_n_exposed_case_events": a,
        "observed_n_exposed_ctrl_events": c,
        "observed_event_or": float(odds_ratio_2x2(a, b, c, d)),
    }


# ============================================================
# MATCHING PREP (v4: MATCH_LEVEL 対応)
# ============================================================

def build_matching_bins(
    dynamic_df: pd.DataFrame,
    static_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    """
    v4: MATCH_LEVEL に応じて exact_key に含める変数を変更。
    - Level 1: chrom のみ（exact_key_vars = []）
    - Level 2: chrom + gene_bin + segdup_bin
    - Level 3: chrom + support_bin + meanprob_bin + gene_bin + segdup_bin
    quantile bin は全レベルで常に計算する（距離重みで使う可能性があるため）。
    """
    log(f"Preparing matching bins (MATCH_LEVEL={MATCH_LEVEL}: "
        f"{_ML_CFG['description']}) ...")

    combined = pd.concat([
        dynamic_df[["support_count", "mean_prob", "gene_density", "segdup_pct"]],
        static_df[["support_count", "mean_prob", "gene_density", "segdup_pct"]],
    ], axis=0, ignore_index=True)

    support_edges = safe_qcut_edges(
        combined["support_count"], max(3, N_PROB_BINS)
    )
    meanprob_edges = safe_qcut_edges(combined["mean_prob"], N_PROB_BINS)
    gene_edges = safe_qcut_edges(combined["gene_density"], N_GENE_BINS)
    segdup_edges = safe_qcut_edges(combined["segdup_pct"], N_SEGDUP_BINS)

    def annotate_bins(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        out["support_bin"] = make_quantile_bins_from_edges(
            out["support_count"], support_edges
        ).astype(int)
        out["meanprob_bin"] = make_quantile_bins_from_edges(
            out["mean_prob"], meanprob_edges
        ).astype(int)
        out["gene_bin"] = make_quantile_bins_from_edges(
            out["gene_density"], gene_edges
        ).astype(int)
        out["segdup_bin"] = make_quantile_bins_from_edges(
            out["segdup_pct"], segdup_edges
        ).astype(int)
        return out

    dyn = annotate_bins(dynamic_df)
    sta = annotate_bins(static_df)

    # v4: exact_key の構築を MATCH_LEVEL に応じて変更
    exact_key_vars = _ML_CFG["exact_key_vars"]  # chrom は常に先頭に付加

    static_lookup: Dict[tuple, List[int]] = defaultdict(list)
    static_chr_lookup: Dict[str, List[int]] = defaultdict(list)

    for idx, row in sta.reset_index(drop=True).iterrows():
        chrom = row["chrom"]
        exact_key = tuple([chrom] + [int(row[v]) for v in exact_key_vars])
        static_lookup[exact_key].append(idx)
        static_chr_lookup[chrom].append(idx)

    return (
        dyn.reset_index(drop=True),
        sta.reset_index(drop=True),
        {
            "support_edges": support_edges,
            "meanprob_edges": meanprob_edges,
            "gene_edges": gene_edges,
            "segdup_edges": segdup_edges,
            "static_lookup": static_lookup,
            "static_chr_lookup": static_chr_lookup,
            "exact_key_vars": exact_key_vars,
        },
    )


def compute_candidate_pool_for_dynamic_bin(
    dyn_row: pd.Series,
    static_df: pd.DataFrame,
    lookup_bundle: Dict[str, object],
) -> Tuple[np.ndarray, np.ndarray, str]:
    """
    v4: MATCH_LEVEL に応じた progressive relaxation と距離重み付け。

    Level 1:
      - Start: chrom_only → 距離重みは gene_density + segdup_pct
      - No relaxation needed (already at loosest)

    Level 2:
      - Start: chrom + gene_bin + segdup_bin (exact)
      - Relax 1: chrom + gene_bin
      - Relax 2: chrom_only
      - 距離重みは gene_density + segdup_pct

    Level 3:
      - Start: chrom + support_bin + meanprob_bin + gene_bin + segdup_bin (exact)
      - Relax 1: chrom + support_bin + gene_bin
      - Relax 2: chrom + support_bin
      - Relax 3: chrom_only
      - 距離重みは support_count + mean_prob + gene_density + segdup_pct
    """
    chrom = dyn_row["chrom"]

    static_lookup = lookup_bundle["static_lookup"]
    static_chr_lookup = lookup_bundle["static_chr_lookup"]
    exact_key_vars = lookup_bundle["exact_key_vars"]

    # --- v4: MATCH_LEVEL 別 progressive relaxation ---
    if MATCH_LEVEL == 1:
        # Level 1: chromosome only — no exact_key lookup, go straight to chrom
        cand = list(static_chr_lookup.get(chrom, []))
        match_level = "chrom_only"

    elif MATCH_LEVEL == 2:
        # Level 2: chrom + gene_bin + segdup_bin
        gene_bin = int(dyn_row["gene_bin"])
        segdup_bin = int(dyn_row["segdup_bin"])

        # exact: chrom + gene_bin + segdup_bin
        exact_key = (chrom, gene_bin, segdup_bin)
        cand = list(static_lookup.get(exact_key, []))
        match_level = "exact_L2"

        if len(cand) < MIN_EXACT_POOL:
            # relax 1: chrom + gene_bin
            mask = (
                (static_df["chrom"] == chrom)
                & (static_df["gene_bin"] == gene_bin)
            )
            cand = static_df.index[mask].to_numpy().tolist()
            match_level = "chrom_gene"

        if len(cand) < MIN_EXACT_POOL:
            # relax 2: chrom only
            cand = list(static_chr_lookup.get(chrom, []))
            match_level = "chrom_only"

    else:
        # Level 3: chrom + support_bin + meanprob_bin + gene_bin + segdup_bin (v3 同等)
        support_bin = int(dyn_row["support_bin"])
        meanprob_bin = int(dyn_row["meanprob_bin"])
        gene_bin = int(dyn_row["gene_bin"])
        segdup_bin = int(dyn_row["segdup_bin"])

        exact_key = (chrom, support_bin, meanprob_bin, gene_bin, segdup_bin)
        cand = list(static_lookup.get(exact_key, []))
        match_level = "exact_L3"

        if len(cand) < MIN_EXACT_POOL:
            mask = (
                (static_df["chrom"] == chrom)
                & (static_df["support_bin"] == support_bin)
                & (static_df["gene_bin"] == gene_bin)
            )
            cand = static_df.index[mask].to_numpy().tolist()
            match_level = "chrom_support_gene"

        if len(cand) < MIN_EXACT_POOL:
            mask = (
                (static_df["chrom"] == chrom)
                & (static_df["support_bin"] == support_bin)
            )
            cand = static_df.index[mask].to_numpy().tolist()
            match_level = "chrom_support"

        if len(cand) < MIN_EXACT_POOL:
            cand = list(static_chr_lookup.get(chrom, []))
            match_level = "chrom_only"

    if len(cand) == 0:
        raise ValueError(
            f"No static candidates for dynamic bin: "
            f"{dyn_row['bin_id']} on {chrom}"
        )

    # --- v4: 距離重み付けは MATCH_LEVEL に応じた変数のみ ---
    dist_vars = _ML_CFG["distance_vars"]

    cand_df = static_df.loc[cand, dist_vars].copy()
    ref = dyn_row[dist_vars].astype(float).to_numpy()
    X = cand_df.astype(float).to_numpy()

    med = np.nanmedian(X, axis=0)
    mad = np.nanmedian(np.abs(X - med), axis=0)
    mad[mad == 0] = 1.0

    ref_z = (ref - med) / mad
    X_z = (X - med) / mad
    d2 = np.sum((X_z - ref_z) ** 2, axis=1)
    w = np.exp(-0.5 * d2)
    if np.all(w <= 0) or not np.isfinite(w).any():
        w = np.ones(len(cand), dtype=float)
    w = w / w.sum()

    return np.asarray(cand, dtype=np.int32), w.astype(np.float64), match_level


def prepare_candidate_map_for_class(
    dynamic_df: pd.DataFrame,
    static_df: pd.DataFrame,
    lookup_bundle: Dict[str, object],
) -> Tuple[List[np.ndarray], List[np.ndarray], List[str]]:
    log(f"Preparing candidate map for set with "
        f"{len(dynamic_df)} dynamic bins ...")
    cand_idx_list: List[np.ndarray] = []
    cand_w_list: List[np.ndarray] = []
    match_levels: List[str] = []

    for _, row in dynamic_df.iterrows():
        cand_idx, cand_w, match_level = compute_candidate_pool_for_dynamic_bin(
            row, static_df, lookup_bundle,
        )
        cand_idx_list.append(cand_idx)
        cand_w_list.append(cand_w)
        match_levels.append(match_level)

    return cand_idx_list, cand_w_list, match_levels


# ============================================================
# RESAMPLING / MAIN ANALYSIS
# ============================================================

def draw_matched_static_set(
    rng: np.random.Generator,
    cand_idx_list: List[np.ndarray],
    cand_w_list: List[np.ndarray],
) -> np.ndarray:
    chosen = []
    for cand_idx, cand_w in zip(cand_idx_list, cand_w_list):
        pick = rng.choice(cand_idx, size=1, replace=True, p=cand_w)[0]
        chosen.append(int(pick))
    # unique bins for burden counting
    chosen = np.array(sorted(set(chosen)), dtype=np.int32)
    return chosen


def compute_sample_burden_from_static_selection(
    selected_static_idx: np.ndarray,
    static_df: pd.DataFrame,
    static_bin_to_sample_idxs: Dict[str, np.ndarray],
    n_samples: int,
) -> np.ndarray:
    burden = np.zeros(n_samples, dtype=np.int32)
    if len(selected_static_idx) == 0:
        return burden

    selected_bin_ids = static_df.iloc[selected_static_idx]["bin_id"].tolist()
    for bid in selected_bin_ids:
        sidx = static_bin_to_sample_idxs.get(str(bid))
        if sidx is not None and len(sidx) > 0:
            np.add.at(burden, sidx, 1)
    return burden


def compute_event_overlap_or_from_static_selection(
    selected_static_idx: np.ndarray,
    static_df: pd.DataFrame,
    static_bin_to_event_idxs: Dict[str, np.ndarray],
    event_df: pd.DataFrame,
) -> Dict[str, float]:
    exposed = np.zeros(len(event_df), dtype=np.int8)
    selected_bin_ids = static_df.iloc[selected_static_idx]["bin_id"].tolist()
    for bid in selected_bin_ids:
        eidx = static_bin_to_event_idxs.get(str(bid))
        if eidx is not None and len(eidx) > 0:
            exposed[eidx] = 1

    case = event_df["is_case"].to_numpy()

    a = int(np.sum((case == 1) & (exposed == 1)))
    b = int(np.sum((case == 1) & (exposed == 0)))
    c = int(np.sum((case == 0) & (exposed == 1)))
    d = int(np.sum((case == 0) & (exposed == 0)))

    return {
        "n_exposed_case_events": a,
        "n_exposed_ctrl_events": c,
        "event_or": float(odds_ratio_2x2(a, b, c, d)),
    }


def empirical_p_upper_tail(
    observed: float,
    null_values: np.ndarray,
) -> float:
    null_values = np.asarray(null_values, dtype=float)
    ok = np.isfinite(null_values)
    if ok.sum() == 0 or not np.isfinite(observed):
        return np.nan
    return float((np.sum(null_values[ok] >= observed) + 1) / (ok.sum() + 1))


def run_class_specificity_test(
    class_name: str,
    dynamic_df: pd.DataFrame,
    static_df: pd.DataFrame,
    lookup_bundle: Dict[str, object],
    analysis_df: pd.DataFrame,
    observed_exposure: np.ndarray,
    static_bin_to_sample_idxs: Dict[str, np.ndarray],
    static_bin_to_event_idxs: Dict[str, np.ndarray],
    event_df: pd.DataFrame,
    rng_seed: int,
    observed_event_stats: Optional[Dict[str, float]] = None,
) -> Tuple[Dict[str, object], pd.DataFrame]:

    rng = np.random.default_rng(rng_seed)

    # observed model
    observed_fit = fit_bprime_logit(analysis_df, observed_exposure)

    # candidate map
    cand_idx_list, cand_w_list, match_levels = (
        prepare_candidate_map_for_class(dynamic_df, static_df, lookup_bundle)
    )

    # resampling arrays
    n_samples = len(analysis_df)
    null_or = np.full(N_RESAMPLES, np.nan, dtype=np.float64)
    null_beta = np.full(N_RESAMPLES, np.nan, dtype=np.float64)
    null_case_mean = np.full(N_RESAMPLES, np.nan, dtype=np.float64)
    null_ctrl_mean = np.full(N_RESAMPLES, np.nan, dtype=np.float64)
    null_unique_bins = np.full(N_RESAMPLES, np.nan, dtype=np.float64)

    # v2 fix: accumulate per-sample burden in single pass
    matched_sum = np.zeros(n_samples, dtype=np.float64)
    matched_exposed_sum = np.zeros(n_samples, dtype=np.float64)

    if RUN_EVENT_LEVEL:
        null_event_or = np.full(N_RESAMPLES, np.nan, dtype=np.float64)
    else:
        null_event_or = None

    log(f"  Starting {N_RESAMPLES} resamples ...")
    t_resample_start = time.time()

    for r in range(N_RESAMPLES):
        selected_static_idx = draw_matched_static_set(
            rng, cand_idx_list, cand_w_list,
        )
        matched_burden = compute_sample_burden_from_static_selection(
            selected_static_idx,
            static_df,
            static_bin_to_sample_idxs,
            n_samples,
        )

        # accumulate for per-sample output (v2 fix)
        matched_sum += matched_burden
        matched_exposed_sum += (matched_burden > 0).astype(float)

        fit = fit_bprime_logit(analysis_df, matched_burden)

        null_or[r] = fit["or"]
        null_beta[r] = fit["beta"]
        null_case_mean[r] = fit["case_mean_burden"]
        null_ctrl_mean[r] = fit["control_mean_burden"]
        null_unique_bins[r] = len(selected_static_idx)

        if RUN_EVENT_LEVEL:
            ev = compute_event_overlap_or_from_static_selection(
                selected_static_idx,
                static_df,
                static_bin_to_event_idxs,
                event_df,
            )
            null_event_or[r] = ev["event_or"]

        if (r + 1) % max(1, N_RESAMPLES // 10) == 0:
            elapsed = time.time() - t_resample_start
            log(f"    resample {r+1}/{N_RESAMPLES} "
                f"({elapsed:.1f}s elapsed)")

    log(f"  Resampling done in {time.time() - t_resample_start:.1f}s")

    # empirical P
    emp_p_or = empirical_p_upper_tail(observed_fit["or"], null_or)
    emp_p_beta = empirical_p_upper_tail(observed_fit["beta"], null_beta)

    def _pctl(arr):
        ok = np.isfinite(arr)
        if ok.any():
            return np.nanpercentile(arr[ok], [2.5, 50, 97.5])
        return [np.nan, np.nan, np.nan]

    or_ci = _pctl(null_or)
    beta_ci = _pctl(null_beta)
    case_mean_ci = _pctl(null_case_mean)
    ctrl_mean_ci = _pctl(null_ctrl_mean)
    unique_bin_ci = _pctl(null_unique_bins)

    match_level_counts = pd.Series(match_levels).value_counts().to_dict()

    result = {
        "class": class_name,
        "match_level_setting": MATCH_LEVEL,  # v4: 記録用
        "n_dynamic_bins": int(dynamic_df["bin_id"].nunique()),
        "n_static_pool_bins": int(static_df["bin_id"].nunique()),
        "observed_n_complete": observed_fit["n_complete"],
        "observed_n_case": observed_fit["n_case"],
        "observed_n_control": observed_fit["n_control"],
        "observed_n_exposed_case": observed_fit["n_exposed_case"],
        "observed_n_exposed_control": observed_fit["n_exposed_control"],
        "observed_case_mean_burden": observed_fit["case_mean_burden"],
        "observed_control_mean_burden": observed_fit["control_mean_burden"],
        "observed_beta": observed_fit["beta"],
        "observed_se": observed_fit["se"],
        "observed_or": observed_fit["or"],
        "observed_ci_lo": observed_fit["ci_lo"],
        "observed_ci_hi": observed_fit["ci_hi"],
        "observed_p_value": observed_fit["p_value"],
        "matched_static_or_mean": float(np.nanmean(null_or)),
        "matched_static_or_median": float(np.nanmedian(null_or)),
        "matched_static_or_ci_lo": float(or_ci[0]),
        "matched_static_or_ci_hi": float(or_ci[2]),
        "matched_static_beta_mean": float(np.nanmean(null_beta)),
        "matched_static_beta_ci_lo": float(beta_ci[0]),
        "matched_static_beta_ci_hi": float(beta_ci[2]),
        "matched_static_case_mean_burden_mean": float(
            np.nanmean(null_case_mean)
        ),
        "matched_static_case_mean_burden_ci_lo": float(case_mean_ci[0]),
        "matched_static_case_mean_burden_ci_hi": float(case_mean_ci[2]),
        "matched_static_control_mean_burden_mean": float(
            np.nanmean(null_ctrl_mean)
        ),
        "matched_static_control_mean_burden_ci_lo": float(ctrl_mean_ci[0]),
        "matched_static_control_mean_burden_ci_hi": float(ctrl_mean_ci[2]),
        "matched_static_unique_bins_mean": float(
            np.nanmean(null_unique_bins)
        ),
        "matched_static_unique_bins_ci_lo": float(unique_bin_ci[0]),
        "matched_static_unique_bins_ci_hi": float(unique_bin_ci[2]),
        "effect_size_diff_or": (
            float(observed_fit["or"] - np.nanmean(null_or))
            if np.isfinite(observed_fit["or"]) else np.nan
        ),
        "effect_size_diff_beta": (
            float(observed_fit["beta"] - np.nanmean(null_beta))
            if np.isfinite(observed_fit["beta"]) else np.nan
        ),
        "empirical_p_or": emp_p_or,
        "empirical_p_beta": emp_p_beta,
        "match_level_counts": match_level_counts,
        "n_resamples": N_RESAMPLES,
    }

    # per-sample output (v2 fix: single pass)
    matched_mean = matched_sum / N_RESAMPLES
    matched_exposed_prob = matched_exposed_sum / N_RESAMPLES

    # event-level
    if RUN_EVENT_LEVEL and observed_event_stats is not None:
        event_or_ci = _pctl(null_event_or)
        result.update({
            "observed_event_or": observed_event_stats.get(
                "observed_event_or", np.nan
            ),
            "observed_n_exposed_case_events": observed_event_stats.get(
                "observed_n_exposed_case_events", np.nan
            ),
            "observed_n_exposed_ctrl_events": observed_event_stats.get(
                "observed_n_exposed_ctrl_events", np.nan
            ),
            "matched_static_event_or_mean": float(
                np.nanmean(null_event_or)
            ),
            "matched_static_event_or_ci_lo": float(event_or_ci[0]),
            "matched_static_event_or_ci_hi": float(event_or_ci[2]),
            "empirical_p_event_or": empirical_p_upper_tail(
                observed_event_stats.get("observed_event_or", np.nan),
                null_event_or,
            ),
        })

    tmp_ps = analysis_df[["sample_id", "Diagnosis"]].copy()
    tmp_ps["class"] = class_name
    tmp_ps["observed_dynamic_burden"] = observed_exposure
    tmp_ps["matched_static_burden_mean"] = matched_mean
    tmp_ps["dynamic_minus_static_mean"] = (
        tmp_ps["observed_dynamic_burden"]
        - tmp_ps["matched_static_burden_mean"]
    )
    tmp_ps["observed_exposed"] = (
        tmp_ps["observed_dynamic_burden"] > 0
    ).astype(int)
    tmp_ps["matched_exposed_prob"] = matched_exposed_prob  # v2 fix

    return result, tmp_ps


# ============================================================
# FLANKING CONTROL (OPTIONAL)
# ============================================================

def build_flanking_control_bins(
    dynamic_union_df: pd.DataFrame,
    h5ad_metrics: pd.DataFrame,
    all_diff_union: pd.DataFrame,
) -> Dict[str, pd.DataFrame]:
    log("Building flanking control bins ...")
    diff_ids = set(all_diff_union["bin_id"].tolist())
    h5ad_ids = set(h5ad_metrics["bin_id"].tolist())
    bin_size = 25000

    out: Dict[str, pd.DataFrame] = {}
    for flank in [1, 2]:
        rows = []
        for _, row in dynamic_union_df.iterrows():
            chrom = row["chrom"]
            s = int(row["start0"])
            e = int(row["end"])
            for shift in [-flank, flank]:
                ns = s + shift * bin_size
                ne = e + shift * bin_size
                if ns < 0 or ne <= ns:
                    continue
                bid = f"{chrom}:{ns}-{ne}"
                if bid in diff_ids:
                    continue
                if bid not in h5ad_ids:
                    continue
                rows.append((chrom, ns, ne, bid))
        df = pd.DataFrame(
            rows, columns=["chrom", "start0", "end", "bin_id"]
        ).drop_duplicates("bin_id")
        out[f"flank_{flank}"] = df
        log(f"  flank_{flank}: {len(df)} bins")
    return out


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    total_start = time.time()
    ensure_output_dirs()  # v6 (2026-04-21): moved from module top-level.
    outdir = Path(_OUTDIR)
    outdir.mkdir(parents=True, exist_ok=True)

    qc = {
        "script": "tad_dynamic_boundary_specificity_test_v7.py",
        "version": "v6",
        "match_level": MATCH_LEVEL,
        "match_level_config": {
            "label": _ML_CFG["label"],
            "description": _ML_CFG["description"],
            "exact_key_vars": _ML_CFG["exact_key_vars"],
            "distance_vars": _ML_CFG["distance_vars"],
        },
        "inputs": {
            "diffbound_dir": _L2_DIFFBOUND_DIR,
            "raw_h5ad": _RAW_H5AD,
            "overlap_table": _OVERLAP_TABLE,
            "burden_table": _BURDEN_TABLE,
            "gtf": _GTF_FILE,
            "segdup_bed": _SEGDUP_BED,
        },
        "analysis": {
            "case": ANALYSIS_DIAGNOSIS_CASE,
            "control": ANALYSIS_DIAGNOSIS_CONTROL,
            "sv_type": SV_TYPE,
            "n_resamples": N_RESAMPLES,
            "random_seed": RANDOM_SEED,
            "run_event_level": RUN_EVENT_LEVEL,
            "run_flanking_control": RUN_FLANKING_CONTROL,
            "matching_variables": (
                ["chromosome"] + _ML_CFG["exact_key_vars"]
            ),
            "distance_weighting_variables": _ML_CFG["distance_vars"],
            "matching_strategy": (
                f"MATCH_LEVEL={MATCH_LEVEL}: "
                f"{_ML_CFG['description']}. "
                "Progressive chromosome-constrained relaxation "
                "+ Gaussian kernel weighted sampling."
            ),
        },
        "class_qc": {},
    }

    log("=" * 80)
    log("tad_dynamic_boundary_specificity_test_v7.py")
    log(f"  MATCH_LEVEL = {MATCH_LEVEL} ({_ML_CFG['description']})")
    log(f"  Exact key vars: chromosome + {_ML_CFG['exact_key_vars']}")
    log(f"  Distance weighting vars: {_ML_CFG['distance_vars']}")
    log("  Main analysis: sample-level B' logistic regression")
    log("  Null: matched static boundary resampling")
    log(f"  Target: {ANALYSIS_DIAGNOSIS_CASE} vs "
        f"{ANALYSIS_DIAGNOSIS_CONTROL}, {SV_TYPE}")
    log(f"  N_RESAMPLES = {N_RESAMPLES}")
    log("=" * 80)

    # 1. load diffbound bins (v5: sig8_union も返す)
    class_bins, analysis_union, all_diff_union, sig8_union = load_diffbound_bins(
        _L2_DIFFBOUND_DIR
    )

    # 2. load h5ad metrics and build static pool
    h5ad_metrics = load_h5ad_support_metrics(_RAW_H5AD)
    static_df = build_static_pool(h5ad_metrics, all_diff_union)

    # 3. annotate dynamic/static bins (v5: 10-class union = primary)
    gene_index = build_gene_index(_GTF_FILE)
    segdup_merged = load_segdup_merged(_SEGDUP_BED)
    dynamic_sets, dynamic_union_annot, static_annot = (
        annotate_bin_metrics_for_dynamic_and_static(
            class_bins, analysis_union, static_df,
            h5ad_metrics, gene_index, segdup_merged,
        )
    )

    if RUN_FLANKING_CONTROL:
        flanking_sets = build_flanking_control_bins(
            dynamic_union_annot, h5ad_metrics, all_diff_union,
        )
        for fk, fdf in flanking_sets.items():
            fdf = fdf.merge(
                h5ad_metrics[
                    ["bin_id", "support_count", "mean_prob", "max_prob"]
                ],
                on="bin_id", how="left", validate="1:1",
            )
            fdf = annotate_gene_density(fdf, gene_index)
            fdf = compute_segdup_pct_for_bins(fdf, segdup_merged)
            dynamic_sets[fk] = fdf.copy()

    # 4. burden df and full burden table (v2: 1回だけ読む)
    analysis_df, burden_full = load_burden_analysis_frame(_BURDEN_TABLE)
    analysis_samples = analysis_df["sample_id"].tolist()

    # 5. overlap table -> static/event maps + pre-loaded filtered table (v3)
    (
        static_bin_to_sample_idxs,
        static_bin_to_event_idxs,
        event_df,
        ov_filtered,
    ) = load_overlap_and_build_maps(_OVERLAP_TABLE, analysis_samples)

    # 6. run class-wise specificity tests
    rng_master = np.random.default_rng(RANDOM_SEED)
    main_results: List[Dict] = []
    per_sample_all: List[pd.DataFrame] = []
    event_level_rows: List[Dict] = []

    full_class_order = PRIMARY_BOUNDARY_CLASSES + ["dynamic_union"]
    if RUN_FLANKING_CONTROL:
        full_class_order += ["flank_1", "flank_2"]

    # v4.1: TARGET_CLASS 指定時は1クラスのみ実行
    if TARGET_CLASS:
        if TARGET_CLASS not in full_class_order:
            raise ValueError(
                f"TARGET_CLASS='{TARGET_CLASS}' is not in class_order: "
                f"{full_class_order}"
            )
        # rng_master を TARGET_CLASS の位置まで送る（再現性保証）
        target_idx = full_class_order.index(TARGET_CLASS)
        for _ in range(target_idx):
            rng_master.integers(0, 2**31 - 1)  # skip seeds for prior classes
        class_order = [TARGET_CLASS]
        log(f"  TARGET_CLASS={TARGET_CLASS} (index {target_idx}/{len(full_class_order)})")
    else:
        class_order = full_class_order

    for class_name in class_order:
        log("-" * 80)
        log(f"Running specificity test for: {class_name}")

        if class_name == "dynamic_union":
            dyn_df = dynamic_union_annot.copy()
        else:
            dyn_df = dynamic_sets[class_name].copy()

        dyn_df, static_match_df, lookup_bundle = build_matching_bins(
            dyn_df, static_annot.copy(),
        )

        # observed exposure
        # v5: dynamic_union は overlap table から 10-class union bins で直接計算
        #     （burden table の n_boundary_dev_union_DEL は 11クラスunion なので使わない）
        obs_col = (
            get_observed_exposure_column(class_name)
            if class_name not in {"flank_1", "flank_2"}
            else None
        )

        if obs_col is not None and obs_col in burden_full.columns:
            tmp = (
                burden_full[["sample_id", obs_col]]
                .copy()
                .rename(columns={obs_col: "observed_exposure"})
            )
            obs_df = analysis_df[["sample_id"]].merge(
                tmp, on="sample_id", how="left",
            )
            observed_exposure = (
                pd.to_numeric(obs_df["observed_exposure"], errors="coerce")
                .fillna(0)
                .to_numpy()
            )
        else:
            # dynamic_union / flanking controls: compute from ov_filtered (v5/v3)
            bin_col = detect_column(ov_filtered, ["bin_id"])
            sample_col = detect_column(ov_filtered, ["sample_id"])
            dyn_bin_set = set(dyn_df["bin_id"].tolist())
            ov_sub = ov_filtered.loc[
                ov_filtered[bin_col].isin(dyn_bin_set)
            ].copy()
            tmp_cnt = (
                ov_sub.groupby(sample_col)[bin_col]
                .nunique()
                .reset_index()
                .rename(columns={
                    sample_col: "sample_id",
                    bin_col: "observed_exposure",
                })
            )
            obs_df = analysis_df[["sample_id"]].merge(
                tmp_cnt, on="sample_id", how="left",
            )
            observed_exposure = (
                pd.to_numeric(obs_df["observed_exposure"], errors="coerce")
                .fillna(0)
                .to_numpy()
            )

        # observed event OR (v3: pre-loaded ov_filtered を使用)
        observed_event_stats = None
        if RUN_EVENT_LEVEL:
            dyn_bin_set = set(dyn_df["bin_id"].tolist())
            observed_event_stats = compute_observed_event_or(
                dyn_bin_set, ov_filtered, event_df,
            )

        class_seed = int(rng_master.integers(0, 2**31 - 1))

        result, per_sample_df = run_class_specificity_test(
            class_name=class_name,
            dynamic_df=dyn_df,
            static_df=static_match_df,
            lookup_bundle=lookup_bundle,
            analysis_df=analysis_df,
            observed_exposure=observed_exposure,
            static_bin_to_sample_idxs=static_bin_to_sample_idxs,
            static_bin_to_event_idxs=static_bin_to_event_idxs,
            event_df=event_df,
            rng_seed=class_seed,
            observed_event_stats=observed_event_stats,
        )

        main_results.append(result)
        per_sample_all.append(per_sample_df)

        qc["class_qc"][class_name] = {
            "n_dynamic_bins": int(result["n_dynamic_bins"]),
            "n_static_pool_bins": int(result["n_static_pool_bins"]),
            "observed_case_mean_burden": result["observed_case_mean_burden"],
            "observed_control_mean_burden": result[
                "observed_control_mean_burden"
            ],
            "observed_or": result["observed_or"],
            "empirical_p_or": result["empirical_p_or"],
            "match_level_counts": result["match_level_counts"],
        }

        if RUN_EVENT_LEVEL:
            event_level_rows.append({
                "class": class_name,
                "observed_event_or": result.get(
                    "observed_event_or", np.nan
                ),
                "observed_n_exposed_case_events": result.get(
                    "observed_n_exposed_case_events", np.nan
                ),
                "observed_n_exposed_ctrl_events": result.get(
                    "observed_n_exposed_ctrl_events", np.nan
                ),
                "matched_static_event_or_mean": result.get(
                    "matched_static_event_or_mean", np.nan
                ),
                "matched_static_event_or_ci_lo": result.get(
                    "matched_static_event_or_ci_lo", np.nan
                ),
                "matched_static_event_or_ci_hi": result.get(
                    "matched_static_event_or_ci_hi", np.nan
                ),
                "empirical_p_event_or": result.get(
                    "empirical_p_event_or", np.nan
                ),
            })

        log(f"  {class_name}: observed_OR={result['observed_or']:.3f}, "
            f"null_OR_median={result['matched_static_or_median']:.3f}, "
            f"null_OR_mean={result['matched_static_or_mean']:.3f}, "
            f"emp_P={result['empirical_p_or']:.4f}")

    # --- output ---
    main_df = pd.DataFrame(main_results)
    per_sample_df = pd.concat(per_sample_all, axis=0, ignore_index=True)

    # v5: figure_summary に median 列を追加（本文・figure caption が median を使用するため）
    figure_df = main_df[[
        "class", "match_level_setting",
        "observed_or", "observed_ci_lo", "observed_ci_hi",
        "matched_static_or_median",
        "matched_static_or_mean",
        "matched_static_or_ci_lo", "matched_static_or_ci_hi",
        "empirical_p_or",
        "observed_case_mean_burden", "observed_control_mean_burden",
        "matched_static_case_mean_burden_mean",
        "matched_static_control_mean_burden_mean",
    ]].copy()

    main_path = outdir / _MAIN_TSV
    per_sample_path = outdir / _PER_SAMPLE_TSV
    figure_path = outdir / _FIGURE_SUMMARY_TSV
    qc_path = outdir / _QC_JSON

    main_df.to_csv(
        main_path, sep="\t", index=False, float_format="%.6g",
    )
    per_sample_df.to_csv(
        per_sample_path, sep="\t", index=False, float_format="%.6g",
    )
    figure_df.to_csv(
        figure_path, sep="\t", index=False, float_format="%.6g",
    )

    if RUN_EVENT_LEVEL and len(event_level_rows) > 0:
        event_df_out = pd.DataFrame(event_level_rows)
        event_df_out.to_csv(
            outdir / _EVENT_TSV, sep="\t", index=False, float_format="%.6g",
        )

    qc["outputs"] = {
        "main_tsv": str(main_path),
        "per_sample_tsv": str(per_sample_path),
        "figure_summary_tsv": str(figure_path),
        "event_tsv": (
            str(outdir / _EVENT_TSV) if RUN_EVENT_LEVEL else None
        ),
    }
    qc["runtime_seconds"] = time.time() - total_start

    with open(qc_path, "w") as f:
        json.dump(qc, f, indent=2, default=str)

    log("=" * 80)
    log("DONE")
    log(f"  MATCH_LEVEL:       {MATCH_LEVEL} ({_ML_CFG['label']})")
    log(f"  Main results:      {main_path}")
    log(f"  Per-sample table:  {per_sample_path}")
    log(f"  Figure summary:    {figure_path}")
    if RUN_EVENT_LEVEL:
        log(f"  Event-level table: {outdir / _EVENT_TSV}")
    log(f"  QC JSON:           {qc_path}")
    log(f"  Total elapsed:     {time.time() - total_start:.1f}s")
    log("=" * 80)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log("FATAL ERROR")
        log(str(e))
        traceback.print_exc()
        sys.exit(1)
