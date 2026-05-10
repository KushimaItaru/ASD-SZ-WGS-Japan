#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ファイル名: 34_intersect_sv_with_heffel_master_v11.py
# - 処理内容:
#   - Pattern B (NAHR-included) / Pattern C (NAHR-excluded + |SV_length| < 1 MB)
#     の感度分析 overlap table を作成する (v10 = Pattern A の sibling)。
#   - v10 -> v11 の差分:
#       1. PATTERN 環境変数 (必須) で patB / patC を切替
#       2. 入力 SV table を Step 03 v9 の Pattern 別出力から読む:
#            patB: OUT_V9/wgs_rare_sv_events_v9_patB.tsv.gz
#            patC: OUT_V9/wgs_rare_sv_events_v9_patC.tsv.gz
#          (paths_v1.F_03_WGS_SV_EVENTS は Pattern A 用なので hardcode で上書き)
#       3. 出力ディレクトリ: output_v11_{PATTERN}/
#          出力ファイル suffix: _v10 -> _v11_{PATTERN}
#            sample_boundary_event_overlap_v11_patB.tsv.gz
#            sample_mechanistic_burden_v11_patB.tsv
#            sample_mechanistic_burden_summary_v11_patB.tsv
#            (patC 同様)
#       4. Pattern C の size cap は上流 (Step 03 v9) で既に適用済みのため、
#          本 script では追加フィルタ不要 (入力段階で既に <1MB のみ)。
#       5. 解析ロジック (sample-level top-1% QC / gene covariate / overlap /
#          3 系統 burden / matching_level stratification) は v10 と完全一致。
#   - 元 v10 の処理内容:
#   - rare SV events と Heffel boundary master v9 を overlap し、event-level overlap table を作成する
#   - 2026-04-22 update (v9 -> v10, sample-level top-1% CNV burden QC 実装):
#       Script 03 v8 (34a_extract_rare_sv_events_v8.py) は event-level でのみ
#       top-1% CNV burden sample の event を除外していた。v10 では sample universe
#       からも同じ top-1% sample を除外 (sinfo + sv + overlap_df + full_sample
#       全てに filter を適用) し、downstream の Script 05/06/11/12/10 が
#       完全に QC 済 sample universe を継承できるようにする。
#       MSSNG (Script 09 v18) および arrayCGH (Script 08 v22) には適用しない
#       (WGS discovery 限定の QC)。
#   - 2026-04-21 update (in-place, ChatGPT review blockers 2 & 3): (v9 由来)
#       Blocker 2: protein-coding gene overlap カウント (total_gene_DEL/DUP, log1p_total_gene_DEL/DUP)
#                  を sample_mechanistic_burden_v10.tsv にマージする
#       Blocker 3: bin overlap filter を `overlap_bp >= 1` から
#                  `overlap_frac_boundary >= 0.10` に変更
#   - master v9 の新列 (matching_level, group_key_L2_from_obs) を overlap_df に伝播する
#   - sample info + PCA から ASD/SZ/Healthy の full cohort sample table を構築する
#   - carrier sample の DEL/DUP 別 burden（total bases / total counts / total exon counts）を統合し、non-carrier は 0 埋めする
#   - sample-level mechanistic burden table を full cohort ベースで作成する
#   - excitatory / astro / inhibitory(CGE/MGE) の count と carrier flag を、SV type 別（*_DEL / *_DUP）に出力する
#   - 列名は固定列番号ではなく列名ベースで動的に選択する
#   - v10 からの主な変更点 (v9 -> v10):
#     (1) sample-level top-1% CNV burden QC の追加
#         (Script 03 v8 と同じ quantile threshold を CNV_SAMPLE_COUNTS から計算)
#         - high-burden sample 集合を sinfo, sv, overlap_df, full_sample 全てから除外
#         - 除外前後のサンプル数を Diagnosis 別に log 出力
#     (2) 出力 suffix を v9 -> v10 に更新
#         - sample_boundary_event_overlap_v9.tsv.gz -> _v10.tsv.gz
#         - sample_mechanistic_burden_v9.tsv        -> _v10.tsv
#         - sample_mechanistic_burden_summary_v9.tsv-> _v10.tsv
#     (3) default outdir を OUT_04_SV_BOUNDARY_OVERLAP (output_v9) から
#         output_v10 に明示変更 (paths_v1.py は不変)
#     (4) CLI 引数 --cnv-sample-counts, --high-burden-percentile を追加
#         (default: paths_v1.CNV_SAMPLE_COUNTS, 0.99)
#   - v9 からの主な変更点 (v9 由来):
#     (1) tad04212026/ パイプラインに移行、common/paths_v1.py を使用
#     (2) 出力 suffix を v8 → v9 に更新
#     (3) 引数なし実行可能（defaults を paths_v1 から提供）
#     (4) 2026-04-21 追加: gene covariates 組込 + overlap_frac_boundary 0.10 default
#   - v7 からの主な変更点:
#     (1) master v9 入力に対応（matching_level, group_key_L2_from_obs 列の存在確認と伝播）
#     (2) burden 集計を 3 系統出力:
#         - unsuffixed: row-level, matching_level ∈ {exact, normalized, l2_aggregated}（全 match、current primary output）
#         - _strict:    row-level, matching_level ∈ {exact, normalized} のみ（v7 相当、sensitivity）
#         - _l2unique:  L2-unique (bin_id × group_key_L2_from_obs), 全 match（L2 lineage 単位、最も生物学的に妥当）
#     (3) unique event 集計も strict 系列を追加
#     (4) 出力 suffix を _v7 → _v8 に統一
#     (5) summary に matching_level breakdown を追加
#   - 実行時間を記録する

import argparse
import gzip
import os
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

# ----- common paths -----
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths_v1 import (
    F_03_WGS_SV_EVENTS, F_01_BOUNDARY_MASTER,
    SAMPLE_INFO, PCA_EIGENVEC, GENCODE_GTF,
    OUT_04_SV_BOUNDARY_OVERLAP, ensure_output_dirs,
    CNV_SAMPLE_COUNTS, PIPELINE_ROOT,  # v10: top-1% QC 用
)

# v10: paths_v1.py の OUT_04 は output_v9 を指しているため、
# v10 用の output_v10 ディレクトリを明示的に構成する
OUT_04_V10 = PIPELINE_ROOT / "04_wgs_sv_boundary_overlap" / "output_v10"

# ============================================================
# v11: PATTERN 環境変数で Pattern B / Pattern C を切替
# ============================================================
_PATTERN = os.environ.get("PATTERN", "").strip()
if _PATTERN not in ("patB", "patC"):
    raise RuntimeError(
        f"v11 requires PATTERN env var = 'patB' or 'patC'. "
        f"Got: {_PATTERN!r}. Example: PATTERN=patB python3 34_...v11.py"
    )

# Step 03 v9 の出力 Pattern 別 rare SV event table
F_03_WGS_SV_EVENTS_V11 = (
    PIPELINE_ROOT / "03_wgs_sv_events" / "output_v9"
    / f"wgs_rare_sv_events_v9_{_PATTERN}.tsv.gz"
)

# v11 Pattern 別 output dir
OUT_04_V11 = PIPELINE_ROOT / "04_wgs_sv_boundary_overlap" / f"output_v11_{_PATTERN}"

# 出力 suffix
V11_SUFFIX = f"v11_{_PATTERN}"


def now_str() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    print(f"[{now_str()}] {msg}", file=sys.stderr, flush=True)


def detect_column(df: pd.DataFrame, candidates: List[str], required: bool = True) -> Optional[str]:
    cols_lower = {str(c).lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    if required:
        raise KeyError(f"候補列が見つかりません: {candidates}\n実際の列: {list(df.columns)}")
    return None


def read_table_auto(path: Path) -> pd.DataFrame:
    suffixes = "".join(path.suffixes).lower()
    compression = "gzip" if suffixes.endswith(".gz") else None

    try:
        return pd.read_csv(path, sep="\t", compression=compression, low_memory=False)
    except Exception:
        pass

    try:
        return pd.read_csv(path, sep=None, engine="python", compression=compression, low_memory=False)
    except Exception as e:
        raise RuntimeError(f"ファイルを読めませんでした: {path}\n{e}")


def normalize_chr(x: str) -> str:
    s = str(x)
    if s.startswith("chr"):
        return s
    return f"chr{s}"


def load_pca_file(path: Path, n_pcs: int = 10) -> pd.DataFrame:
    pca = pd.read_csv(path, sep=r"\s+", header=0, low_memory=False)
    if "#IID" in pca.columns:
        pca_id_col = "#IID"
    elif "IID" in pca.columns:
        pca_id_col = "IID"
    else:
        pca_id_col = pca.columns[1] if len(pca.columns) > 1 else pca.columns[0]

    pc_cols = [c for c in pca.columns if str(c).startswith("PC") or str(c).startswith("pc")][:n_pcs]
    out = pca[[pca_id_col] + pc_cols].copy()
    out = out.rename(columns={pca_id_col: "sample_id"})
    return out


def make_event_id(df: pd.DataFrame, sample_col: str, chr_col: str, start_col: str, end_col: str, svtype_col: str) -> pd.Series:
    return (
        df[sample_col].astype(str) + "|" +
        df[chr_col].astype(str) + ":" +
        df[start_col].astype(str) + "-" +
        df[end_col].astype(str) + "|" +
        df[svtype_col].astype(str)
    )


# =========================================================================
# v10: sample-level top-1% CNV burden QC helper
# =========================================================================

def load_high_burden_samples(
    cnv_counts_path: Path,
    percentile: float,
) -> tuple[set[str], float]:
    """
    Script 03 v8 (34a_extract_rare_sv_events_v8.py) と同じロジックで
    top-percentile 以上の CNV 数を持つ sample を抽出する。

    Returns:
        (high_burden_samples, thr)
    """
    cnt = read_table_auto(cnv_counts_path)
    cnt_sample_col = detect_column(cnt, ["sampleID", "sample_id", "Sample_ID", "IID"])
    cnt_count_col = detect_column(cnt, ["cnv_count", "count", "n_cnv"])
    cnt = cnt.rename(columns={cnt_sample_col: "sample_id", cnt_count_col: "cnv_count"})
    thr = float(cnt["cnv_count"].quantile(percentile))
    high_burden = set(cnt.loc[cnt["cnv_count"] >= thr, "sample_id"].astype(str))
    return high_burden, thr


# =========================================================================
# 2026-04-21: Blocker 2 - protein-coding gene overlap index
#   Copied from unified_heffel_tad_pipeline_v18.run_step3 /
#   build_gene_interval_index / get_overlapping_genes.
#   Important: this Script 4 normalizes SV chrom to "chr*" prefix form, so
#   the gene_index here also uses "chr*" keys (different from v18 which
#   stripped "chr"). SV and master master use "chr*" consistently.
# =========================================================================

def parse_gtf_attributes(attr_str: str) -> dict:
    attrs = {}
    for field in attr_str.strip().rstrip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        parts = field.split(" ", 1)
        if len(parts) == 2:
            attrs[parts[0].strip()] = parts[1].strip().strip('"')
    return attrs


def build_gene_interval_index(gtf_path: Path):
    """Build protein-coding gene interval index from GENCODE GTF.

    gene_key = 'gene_id|gene_name' to dedup per unique gene across different
    SVs from the same sample. Chromosome keys use the 'chr*' form consumed
    by the rest of this script (see normalize_chr).
    """
    gene_records: Dict[str, List[tuple]] = defaultdict(list)
    n_gene_rows = 0
    gtf_path = str(gtf_path)
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = parse_gtf_attributes(fields[8])
            gene_type = attrs.get("gene_type", attrs.get("gene_biotype", ""))
            if gene_type != "protein_coding":
                continue
            chrom = normalize_chr(fields[0])  # chr*
            try:
                s1, e1 = int(fields[3]), int(fields[4])
            except Exception:
                continue
            if e1 < s1:
                s1, e1 = e1, s1
            s0 = max(0, s1 - 1)
            if e1 <= s0:
                continue
            gene_id = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", gene_id if gene_id else "NA")
            gene_key = f"{gene_id}|{gene_name}" if gene_id else gene_name
            gene_records[chrom].append((s0, e1, gene_key))
            n_gene_rows += 1

    log(f"  Parsed {n_gene_rows} protein-coding genes from GTF")
    gene_index = {}
    for chrom, lst in gene_records.items():
        arr = np.array(lst, dtype=object)
        starts = arr[:, 0].astype(np.int64)
        ends = arr[:, 1].astype(np.int64)
        genes = arr[:, 2].astype(str)
        order = np.argsort(starts)
        gene_index[chrom] = {
            "start": starts[order],
            "end":   ends[order],
            "gene":  genes[order],
        }
    return gene_index, n_gene_rows


def get_overlapping_genes(chrom: str, sv_start: int, sv_end: int, gene_index) -> set:
    """Return set of protein-coding gene_keys overlapping the SV interval."""
    if chrom not in gene_index:
        return set()
    starts = gene_index[chrom]["start"]
    ends = gene_index[chrom]["end"]
    genes = gene_index[chrom]["gene"]
    # candidate window: genes with start < sv_end
    idx_hi = int(np.searchsorted(starts, sv_end, side="left"))
    if idx_hi <= 0:
        return set()
    mask = ends[:idx_hi] > sv_start
    if not np.any(mask):
        return set()
    return set(genes[:idx_hi][mask].tolist())


def compute_sample_gene_overlap_counts(
    sv_df: pd.DataFrame,
    gene_index,
    sample_col: str,
    chr_col: str,
    start_col: str,
    end_col: str,
    svtype_col: str,
) -> pd.DataFrame:
    """
    Return DataFrame[sample_id, total_gene_DEL, total_gene_DUP,
                     log1p_total_gene_DEL, log1p_total_gene_DUP].
    Unique gene counting dedups identical gene_keys hit by multiple SVs in
    the same sample × SV_type stratum.
    """
    sample_gene_sets: Dict[str, Dict[str, set]] = {"DEL": defaultdict(set), "DUP": defaultdict(set)}
    for row in sv_df.itertuples(index=False):
        svt = str(getattr(row, svtype_col)).upper().strip()
        if svt not in ("DEL", "DUP"):
            continue
        sid = getattr(row, sample_col)
        chrom = getattr(row, chr_col)
        s0 = int(getattr(row, start_col))
        e = int(getattr(row, end_col))
        overlapped = get_overlapping_genes(chrom, s0, e, gene_index)
        if overlapped:
            sample_gene_sets[svt][sid].update(overlapped)

    all_sids = set(sv_df[sample_col].astype(str).tolist())
    for svt in ("DEL", "DUP"):
        for sid in all_sids:
            sample_gene_sets[svt].setdefault(sid, set())

    rows = []
    for sid in all_sids:
        n_del = len(sample_gene_sets["DEL"][sid])
        n_dup = len(sample_gene_sets["DUP"][sid])
        rows.append({
            "sample_id": sid,
            "total_gene_DEL": int(n_del),
            "total_gene_DUP": int(n_dup),
        })
    df = pd.DataFrame(rows)
    df["log1p_total_gene_DEL"] = np.log1p(df["total_gene_DEL"])
    df["log1p_total_gene_DUP"] = np.log1p(df["total_gene_DUP"])

    for svt in ("DEL", "DUP"):
        col = f"total_gene_{svt}"
        nz = int((df[col] > 0).sum())
        log(f"  gene overlap covariate {col}: nonzero={nz}, "
            f"mean={df[col].mean():.2f}, max={int(df[col].max())}")
    return df


# ---------------------------------------------------------------------------
# 3 系統 burden 集計 (v8 で導入、現行 v10 でも踏襲)
# - unsuffixed: row-level, matching_level ∈ {exact, normalized, l2_aggregated}
# - _strict:    row-level, matching_level ∈ {exact, normalized}
# - _l2unique:  L2-unique (bin_id × group_key_L2_from_obs), all matches
# ---------------------------------------------------------------------------
STRICT_LEVELS = {"exact", "normalized"}


def _compute_row_level_metrics(sub: pd.DataFrame, sample_col: str, suffix: str) -> pd.DataFrame:
    """
    row-level burden 集計（v7 互換ロジック）。
    suffix は追加する列名接尾辞（"" または "_strict"）。
    """
    if sub.empty:
        return pd.DataFrame(columns=[sample_col])

    grp = sub.groupby(sample_col, as_index=False).agg(
        **{
            f"n_boundary_support2_all{suffix}":      ("bin_id", "size"),
            f"n_boundary_consensus_support2{suffix}": ("is_consensus_support2", "sum"),
            f"n_boundary_overlaps_diffbound{suffix}": ("overlaps_diffbound", "sum"),
            f"n_boundary_exc_rg_gain{suffix}":       ("is_exc_rg_gain", "sum"),
            f"n_boundary_astro_rg_loss{suffix}":     ("is_astro_rg_loss", "sum"),
            f"n_boundary_inh_cge_gain{suffix}":      ("is_inh_cge_gain", "sum"),
            f"n_boundary_inh_mge_gain{suffix}":      ("is_inh_mge_gain", "sum"),
        }
    )

    ev_all = (
        sub.groupby(sample_col)["event_id"].nunique()
        .rename(f"n_unique_events_overlap_any{suffix}").reset_index()
    )
    ev_exc = (
        sub.loc[sub["is_exc_rg_gain"] == 1]
        .groupby(sample_col)["event_id"].nunique()
        .rename(f"n_unique_events_overlap_exc_primary{suffix}").reset_index()
    )
    ev_inh_cge = (
        sub.loc[sub["is_inh_cge_gain"] == 1]
        .groupby(sample_col)["event_id"].nunique()
        .rename(f"n_unique_events_overlap_inh_cge_primary{suffix}").reset_index()
    )
    ev_inh_mge = (
        sub.loc[sub["is_inh_mge_gain"] == 1]
        .groupby(sample_col)["event_id"].nunique()
        .rename(f"n_unique_events_overlap_inh_mge_primary{suffix}").reset_index()
    )

    tmp = grp.merge(ev_all, on=sample_col, how="left") \
             .merge(ev_exc, on=sample_col, how="left") \
             .merge(ev_inh_cge, on=sample_col, how="left") \
             .merge(ev_inh_mge, on=sample_col, how="left")
    return tmp


def _compute_l2unique_metrics(sub: pd.DataFrame, sample_col: str) -> pd.DataFrame:
    """
    L2-unique burden 集計:
    (sample_id, bin_id, group_key_L2_from_obs) で unique 化してから flag sum を取る。
    同じ L2 lineage 内の複数 fine cluster による二重カウントを避ける。
    """
    if sub.empty:
        return pd.DataFrame(columns=[sample_col])

    dedup_keys = [sample_col, "bin_id", "group_key_L2_from_obs"]
    agg_map = {
        "is_consensus_support2": "max",
        "overlaps_diffbound":   "max",
        "is_exc_rg_gain":       "max",
        "is_astro_rg_loss":     "max",
        "is_inh_cge_gain":      "max",
        "is_inh_mge_gain":      "max",
        "event_id":             "first",
    }
    deduped = sub.groupby(dedup_keys, as_index=False).agg(agg_map)

    grp = deduped.groupby(sample_col, as_index=False).agg(
        n_boundary_support2_all_l2unique=("bin_id", "size"),
        n_boundary_consensus_support2_l2unique=("is_consensus_support2", "sum"),
        n_boundary_overlaps_diffbound_l2unique=("overlaps_diffbound", "sum"),
        n_boundary_exc_rg_gain_l2unique=("is_exc_rg_gain", "sum"),
        n_boundary_astro_rg_loss_l2unique=("is_astro_rg_loss", "sum"),
        n_boundary_inh_cge_gain_l2unique=("is_inh_cge_gain", "sum"),
        n_boundary_inh_mge_gain_l2unique=("is_inh_mge_gain", "sum"),
    )
    return grp


def aggregate_sample_burden_by_type(
    overlap_df: pd.DataFrame,
    sample_meta_df: pd.DataFrame,
    sample_col: str,
    svtype_col: str = "sv_type_norm",
) -> pd.DataFrame:
    """
    3 系統の burden を SV type × sample 単位で集計する (v8 で導入、現行 v10 でも踏襲)。
    """
    out = sample_meta_df.copy().set_index(sample_col, drop=False)
    sv_types = ["DEL", "DUP"]

    metric_templates_all = [
        "n_boundary_support2_all", "n_boundary_consensus_support2",
        "n_boundary_overlaps_diffbound", "n_boundary_exc_rg_gain",
        "n_boundary_astro_rg_loss", "n_boundary_inh_cge_gain", "n_boundary_inh_mge_gain",
        "n_unique_events_overlap_any", "n_unique_events_overlap_exc_primary",
        "n_unique_events_overlap_inh_cge_primary", "n_unique_events_overlap_inh_mge_primary",
        "carrier_boundary_any_overlap", "carrier_boundary_exc_rg_gain",
        "carrier_boundary_astro_rg_loss", "carrier_boundary_inh_cge_gain",
        "carrier_boundary_inh_mge_gain",
    ]
    metric_templates_strict = [m + "_strict" for m in metric_templates_all]
    metric_templates_l2unique = [
        "n_boundary_support2_all_l2unique", "n_boundary_consensus_support2_l2unique",
        "n_boundary_overlaps_diffbound_l2unique", "n_boundary_exc_rg_gain_l2unique",
        "n_boundary_astro_rg_loss_l2unique", "n_boundary_inh_cge_gain_l2unique",
        "n_boundary_inh_mge_gain_l2unique",
        "carrier_boundary_exc_rg_gain_l2unique", "carrier_boundary_astro_rg_loss_l2unique",
        "carrier_boundary_inh_cge_gain_l2unique", "carrier_boundary_inh_mge_gain_l2unique",
    ]
    all_metric_templates = metric_templates_all + metric_templates_strict + metric_templates_l2unique

    for svt in sv_types:
        for mt in all_metric_templates:
            out[f"{mt}_{svt}"] = 0

    if overlap_df.empty:
        return out.reset_index(drop=True)

    for svt in sv_types:
        sub_all = overlap_df.loc[overlap_df[svtype_col] == svt].copy()
        if sub_all.empty:
            continue

        tmp_all = _compute_row_level_metrics(sub_all, sample_col, suffix="")
        sub_strict = sub_all.loc[sub_all["matching_level"].isin(STRICT_LEVELS)].copy()
        tmp_strict = _compute_row_level_metrics(sub_strict, sample_col, suffix="_strict")
        tmp_l2u = _compute_l2unique_metrics(sub_all, sample_col)

        tmp = pd.DataFrame({sample_col: out.index})
        for t in [tmp_all, tmp_strict, tmp_l2u]:
            if not t.empty:
                tmp = tmp.merge(t, on=sample_col, how="left")

        numeric_cols = [c for c in tmp.columns if c != sample_col]
        for c in numeric_cols:
            tmp[c] = tmp[c].fillna(0).astype(int)

        tmp["carrier_boundary_any_overlap"] = (tmp["n_unique_events_overlap_any"] > 0).astype(int)
        tmp["carrier_boundary_exc_rg_gain"] = (tmp["n_boundary_exc_rg_gain"] > 0).astype(int)
        tmp["carrier_boundary_astro_rg_loss"] = (tmp["n_boundary_astro_rg_loss"] > 0).astype(int)
        tmp["carrier_boundary_inh_cge_gain"] = (tmp["n_boundary_inh_cge_gain"] > 0).astype(int)
        tmp["carrier_boundary_inh_mge_gain"] = (tmp["n_boundary_inh_mge_gain"] > 0).astype(int)

        tmp["carrier_boundary_any_overlap_strict"] = (tmp["n_unique_events_overlap_any_strict"] > 0).astype(int)
        tmp["carrier_boundary_exc_rg_gain_strict"] = (tmp["n_boundary_exc_rg_gain_strict"] > 0).astype(int)
        tmp["carrier_boundary_astro_rg_loss_strict"] = (tmp["n_boundary_astro_rg_loss_strict"] > 0).astype(int)
        tmp["carrier_boundary_inh_cge_gain_strict"] = (tmp["n_boundary_inh_cge_gain_strict"] > 0).astype(int)
        tmp["carrier_boundary_inh_mge_gain_strict"] = (tmp["n_boundary_inh_mge_gain_strict"] > 0).astype(int)

        tmp["carrier_boundary_exc_rg_gain_l2unique"] = (tmp["n_boundary_exc_rg_gain_l2unique"] > 0).astype(int)
        tmp["carrier_boundary_astro_rg_loss_l2unique"] = (tmp["n_boundary_astro_rg_loss_l2unique"] > 0).astype(int)
        tmp["carrier_boundary_inh_cge_gain_l2unique"] = (tmp["n_boundary_inh_cge_gain_l2unique"] > 0).astype(int)
        tmp["carrier_boundary_inh_mge_gain_l2unique"] = (tmp["n_boundary_inh_mge_gain_l2unique"] > 0).astype(int)

        tmp = tmp.set_index(sample_col, drop=True)

        for c in tmp.columns:
            new_col = f"{c}_{svt}"
            out.loc[tmp.index, new_col] = tmp[c].astype(int)

    return out.reset_index(drop=True)


def main():
    t0 = time.time()
    ensure_output_dirs()
    # v11: output_v11_{PATTERN} ディレクトリを確実に作成
    OUT_04_V11.mkdir(parents=True, exist_ok=True)

    parser = argparse.ArgumentParser(
        description=f"v11 ({_PATTERN}): Intersect rare DEL/DUP events with Heffel boundary "
                    f"master v9 and build full-cohort sample-level burden table with "
                    f"matching_level stratification. Input: Step 03 v9 {_PATTERN} output. "
                    f"Same analysis logic as v10 (Pattern A). "
                    f"PATTERN env var required ('patB' or 'patC')."
    )
    # v11: default input は Step 03 v9 の Pattern 別出力
    parser.add_argument("--sv", default=str(F_03_WGS_SV_EVENTS_V11))
    parser.add_argument("--master", default=str(F_01_BOUNDARY_MASTER))
    parser.add_argument("--sample-info", default=str(SAMPLE_INFO))
    parser.add_argument("--pca", default=str(PCA_EIGENVEC))
    parser.add_argument("--gtf", default=str(GENCODE_GTF),
                        help="GENCODE GTF for protein-coding gene interval index "
                             "(default: paths_v1.GENCODE_GTF)")
    # v10: top-1% QC 用引数 (v11 でも同一)
    parser.add_argument("--cnv-sample-counts", default=str(CNV_SAMPLE_COUNTS),
                        help="Per-sample CNV count file used to define top-percentile "
                             "high-burden samples (matches Script 03 v9 input).")
    parser.add_argument("--high-burden-percentile", type=float, default=0.99,
                        help="Percentile threshold for sample-level high-burden exclusion "
                             "(default 0.99; matches Script 03 v9 event-level QC).")
    # v11: デフォルト outdir を output_v11_{PATTERN} に変更 (paths_v1 は不変)
    parser.add_argument("--outdir", default=str(OUT_04_V11))
    parser.add_argument("--sample-col", default=None)
    parser.add_argument("--chr-col", default=None)
    parser.add_argument("--start-col", default=None)
    parser.add_argument("--end-col", default=None)
    parser.add_argument("--diag-col", default=None)
    parser.add_argument("--sex-col", default=None)
    parser.add_argument("--svtype-col", default=None)
    parser.add_argument("--keep-dx", default="ASD,SZ,Healthy")
    parser.add_argument("--min-overlap-bp", type=int, default=1,
                        help="Minimum overlap in bp (legacy absolute filter, default 1). "
                             "Kept for backward compatibility; combined with "
                             "--min-overlap-frac-boundary via AND.")
    parser.add_argument("--min-overlap-frac-boundary", type=float, default=0.10,
                        help="Minimum overlap_frac_boundary = overlap_bp / boundary_bin_width. "
                             "Default 0.10 (matches unified_heffel_tad_pipeline_v18.bin_overlap_thr).")
    args = parser.parse_args()

    log(f"Start 34_intersect_sv_with_heffel_master_v11.py [PATTERN={_PATTERN}]")
    log(f"  input SV: {args.sv}")
    log(f"  output dir: {args.outdir}")
    log(f"  output suffix: {V11_SUFFIX}")
    log(f"  overlap thresholds: --min-overlap-bp={args.min_overlap_bp}, "
        f"--min-overlap-frac-boundary={args.min_overlap_frac_boundary}")
    log(f"  sample-level QC: --cnv-sample-counts={args.cnv_sample_counts}, "
        f"--high-burden-percentile={args.high_burden_percentile}")

    sv_path = Path(args.sv)
    master_path = Path(args.master)
    sample_info_path = Path(args.sample_info)
    pca_path = Path(args.pca)
    gtf_path = Path(args.gtf)
    cnv_counts_path = Path(args.cnv_sample_counts)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for p in [sv_path, master_path, sample_info_path, pca_path, gtf_path, cnv_counts_path]:
        if not p.exists():
            raise FileNotFoundError(f"Not found: {p}")

    # =====================================================================
    # v10: Load high-burden sample set (top percentile by cnv_count)
    # =====================================================================
    log(f"[v10] Reading CNV sample counts: {cnv_counts_path}")
    high_burden_samples, thr = load_high_burden_samples(
        cnv_counts_path, args.high_burden_percentile
    )
    log(f"[v10] Sample-level top-{(1-args.high_burden_percentile)*100:.0f}% threshold: "
        f"cnv_count >= {thr:.0f}  -> {len(high_burden_samples)} samples flagged for exclusion")

    log(f"Reading SV table: {sv_path}")
    sv = read_table_auto(sv_path)

    sample_col = args.sample_col or detect_column(sv, ["sample_id", "Sample_ID", "sample", "IID", "id", "sampleid"])
    chr_col = args.chr_col or detect_column(sv, ["chrom", "chr", "chromosome", "CHROM", "SV_chrom"])
    start_col = args.start_col or detect_column(sv, ["start0", "start", "POS", "sv_start", "BEGIN"])
    end_col = args.end_col or detect_column(sv, ["end", "stop", "sv_end", "END", "SV_end"])
    diag_col = args.diag_col or detect_column(sv, ["Diagnosis", "diagnosis", "DX", "case_control", "phenotype"], required=False)
    sex_col = args.sex_col or detect_column(sv, ["Sex", "sex", "SEX", "gender"], required=False)
    svtype_col = args.svtype_col or detect_column(sv, ["sv_type_norm", "SV_type", "sv_type", "type"])

    log(f"Detected SV columns: sample={sample_col}, chr={chr_col}, start={start_col}, end={end_col}, svtype={svtype_col}, diagnosis={diag_col}, sex={sex_col}")

    sv = sv.copy()
    sv[chr_col] = sv[chr_col].map(normalize_chr)
    sv[start_col] = pd.to_numeric(sv[start_col], errors="coerce")
    sv[end_col] = pd.to_numeric(sv[end_col], errors="coerce")
    sv[svtype_col] = sv[svtype_col].astype(str).str.upper().str.strip()

    sv = sv.dropna(subset=[sample_col, chr_col, start_col, end_col]).copy()
    sv[start_col] = sv[start_col].astype(int)
    sv[end_col] = sv[end_col].astype(int)
    sv = sv.loc[sv[end_col] > sv[start_col]].copy()
    sv["event_id"] = make_event_id(sv, sample_col, chr_col, start_col, end_col, svtype_col)
    sv["sv_len"] = sv[end_col] - sv[start_col]

    # =====================================================================
    # v10: Apply sample-level high-burden exclusion to SV events
    # (Script 03 v8 はこれらの event を既に event-level で除外しているが、
    #  念のため sample-level でも filter を明示適用して冪等にする)
    # =====================================================================
    n_sv_before = sv.shape[0]
    n_sample_before = sv[sample_col].astype(str).nunique()
    sv = sv.loc[~sv[sample_col].astype(str).isin(high_burden_samples)].copy()
    n_sv_after = sv.shape[0]
    n_sample_after = sv[sample_col].astype(str).nunique()
    log(f"[v10] High-burden sample exclusion on SV events: "
        f"{n_sv_before} -> {n_sv_after} events "
        f"({n_sample_before} -> {n_sample_after} unique samples)")

    # ---------------------------------------------------------------------
    # 2026-04-21: protein-coding gene covariate computation (Blocker 2)
    # ---------------------------------------------------------------------
    log(f"Building protein-coding gene interval index from: {gtf_path}")
    gene_index, n_gene_rows = build_gene_interval_index(gtf_path)
    log(f"  Loaded {n_gene_rows} protein-coding genes into interval index "
        f"across {len(gene_index)} chromosomes.")

    log("Computing per-sample × SV type unique gene overlap counts ...")
    gene_cov = compute_sample_gene_overlap_counts(
        sv_df=sv,
        gene_index=gene_index,
        sample_col=sample_col,
        chr_col=chr_col,
        start_col=start_col,
        end_col=end_col,
        svtype_col=svtype_col,
    )
    gene_cov = gene_cov.rename(columns={sample_col: "sample_id"}) \
        if sample_col != "sample_id" else gene_cov
    # free memory
    del gene_index

    # Preserve existing per-sample burden columns coming from the SV event table.
    burden_cols_present = [c for c in [
        "total_del_bases", "total_dup_bases",
        "total_del_count", "total_dup_count",
        "total_exon_del_count", "total_exon_dup_count",
        "cnv_count"
    ] if c in sv.columns]

    agg_dict = {c: "max" for c in burden_cols_present}
    carrier_sample_burden = sv.groupby(sample_col, as_index=False).agg(agg_dict) if agg_dict else sv[[sample_col]].drop_duplicates().copy()
    carrier_sample_burden["carrier_any_rare_sv"] = 1

    log(f"Reading boundary master v9: {master_path}")
    master = pd.read_csv(master_path, sep="\t", compression="gzip", low_memory=False)
    master["chrom"] = master["chrom"].map(normalize_chr)
    master["start0"] = pd.to_numeric(master["start0"], errors="coerce").astype(int)
    master["end"] = pd.to_numeric(master["end"], errors="coerce").astype(int)

    v9_required_cols = ["matching_level", "group_key_L2_from_obs"]
    missing_v9_cols = [c for c in v9_required_cols if c not in master.columns]
    if missing_v9_cols:
        raise KeyError(
            f"master v9 に必要列がありません: {missing_v9_cols}\n"
            f"master_v7 など古いバージョンを指定していませんか? v9 を指定してください。"
        )
    log(f"master rows: {master.shape[0]}")
    log(f"master matching_level distribution:")
    for lvl, cnt in master["matching_level"].value_counts().items():
        log(f"    {lvl:<15s}: {cnt}")

    log(f"Reading sample info: {sample_info_path}")
    sinfo = read_table_auto(sample_info_path)
    sinfo_sample_col = detect_column(sinfo, ["SampleID", "sampleID", "sample_id", "IID", "#IID", "ID"])
    sinfo_dx_col = detect_column(sinfo, ["Diagnosis", "diagnosis", "Diagnosis_DSM5", "phenotype", "disease"])
    sinfo_sex_col = detect_column(sinfo, ["Sex", "sex", "gender"])

    sinfo = sinfo.rename(columns={
        sinfo_sample_col: "sample_id",
        sinfo_dx_col: "Diagnosis",
        sinfo_sex_col: "Sex"
    })
    keep_dx = {x.strip() for x in args.keep_dx.split(",") if x.strip() != ""}
    sinfo = sinfo.loc[sinfo["Diagnosis"].isin(keep_dx)].copy()
    sinfo["Sex_numeric"] = sinfo["Sex"].map({"M": 1, "Male": 1, "F": 0, "Female": 0})

    # =====================================================================
    # v10: Apply sample-level high-burden exclusion to sinfo (sample universe)
    # =====================================================================
    log(f"[v10] sinfo counts BEFORE high-burden exclusion (Diagnosis breakdown):")
    for dx_, n_ in sinfo["Diagnosis"].value_counts().items():
        log(f"    {dx_}: {n_}")
    n_sinfo_before = sinfo.shape[0]
    sinfo = sinfo.loc[~sinfo["sample_id"].astype(str).isin(high_burden_samples)].copy()
    n_sinfo_after = sinfo.shape[0]
    log(f"[v10] sinfo sample universe: {n_sinfo_before} -> {n_sinfo_after} "
        f"(removed {n_sinfo_before - n_sinfo_after} high-burden samples)")
    log(f"[v10] sinfo counts AFTER high-burden exclusion (Diagnosis breakdown):")
    for dx_, n_ in sinfo["Diagnosis"].value_counts().items():
        log(f"    {dx_}: {n_}")

    log(f"Reading PCA: {pca_path}")
    pca = load_pca_file(pca_path, n_pcs=10)

    full_sample = (
        sinfo[["sample_id", "Diagnosis", "Sex", "Sex_numeric"]]
        .drop_duplicates(subset=["sample_id"])
        .merge(pca, on="sample_id", how="left")
    )

    full_sample = full_sample.merge(
        carrier_sample_burden.rename(columns={sample_col: "sample_id"}),
        on="sample_id",
        how="left"
    )

    # 2026-04-21: merge gene covariates (total_gene_DEL/DUP, log1p_*)
    full_sample = full_sample.merge(gene_cov, on="sample_id", how="left")

    numeric_zero_cols = [
        "total_del_bases", "total_dup_bases",
        "total_del_count", "total_dup_count",
        "total_exon_del_count", "total_exon_dup_count",
        "cnv_count", "carrier_any_rare_sv",
        "total_gene_DEL", "total_gene_DUP",
    ]
    for c in numeric_zero_cols:
        if c not in full_sample.columns:
            full_sample[c] = 0
        full_sample[c] = full_sample[c].fillna(0)

    # log1p_total_gene_* : fillna(0) then recompute log1p from the raw counts
    # to avoid NaN propagation for samples with no SVs (non-carriers).
    for svt in ("DEL", "DUP"):
        raw_col = f"total_gene_{svt}"
        log_col = f"log1p_total_gene_{svt}"
        full_sample[log_col] = np.log1p(full_sample[raw_col].astype(float))

    for c in ["total_del_bases", "total_dup_bases", "cnv_count"]:
        if c in full_sample.columns:
            full_sample[c] = full_sample[c].astype(float)

    for c in ["total_del_count", "total_dup_count", "total_exon_del_count", "total_exon_dup_count",
              "carrier_any_rare_sv", "total_gene_DEL", "total_gene_DUP"]:
        if c in full_sample.columns:
            full_sample[c] = full_sample[c].astype(int)

    log(f"[v10] Full cohort samples retained (after high-burden excl): {full_sample.shape[0]}")
    log(f"[v10] Full cohort Diagnosis breakdown:")
    for dx_, n_ in full_sample["Diagnosis"].value_counts().items():
        log(f"    {dx_}: {n_}")
    log(f"Carrier samples in full cohort: {int(full_sample['carrier_any_rare_sv'].sum())}")
    log(f"Samples with >=1 protein-coding gene DEL overlap: "
        f"{int((full_sample['total_gene_DEL'] > 0).sum())}")
    log(f"Samples with >=1 protein-coding gene DUP overlap: "
        f"{int((full_sample['total_gene_DUP'] > 0).sum())}")

    log("Running chromosome-wise interval overlap")
    overlap_chunks = []

    common_chrs = sorted(set(sv[chr_col].unique()) & set(master["chrom"].unique()))
    log(f"Common chromosomes: {len(common_chrs)}")

    for chrom in common_chrs:
        sv_chr = sv.loc[sv[chr_col] == chrom].copy()
        master_chr = master.loc[master["chrom"] == chrom].copy()

        if sv_chr.empty or master_chr.empty:
            continue

        sv_chr = sv_chr.sort_values(start_col).reset_index(drop=True)
        master_chr = master_chr.sort_values("start0").reset_index(drop=True)

        m_starts = master_chr["start0"].to_numpy()

        for _, row in sv_chr.iterrows():
            s = int(row[start_col])
            e = int(row[end_col])

            right = np.searchsorted(m_starts, e, side="left")
            if right == 0:
                continue

            cand = master_chr.iloc[:right]
            cand = cand.loc[cand["end"] > s]
            if cand.empty:
                continue

            tmp = cand.copy()
            tmp["sample_id"] = row[sample_col]
            tmp["event_id"] = row["event_id"]
            tmp["sv_chr"] = chrom
            tmp["sv_start0"] = s
            tmp["sv_end"] = e
            tmp["sv_len"] = row["sv_len"]
            tmp["sv_type_norm"] = row[svtype_col]

            if diag_col is not None:
                tmp["Diagnosis"] = row[diag_col]
            if sex_col is not None:
                tmp["Sex"] = row[sex_col]

            tmp["overlap_bp"] = np.minimum(tmp["end"].to_numpy(), e) - np.maximum(tmp["start0"].to_numpy(), s)
            tmp["overlap_frac_boundary"] = tmp["overlap_bp"] / (tmp["end"] - tmp["start0"])
            tmp["overlap_frac_sv"] = tmp["overlap_bp"] / max(1, (e - s))

            # 2026-04-21: Blocker 3 — combined filter: overlap_bp >= min_overlap_bp
            #            AND overlap_frac_boundary >= min_overlap_frac_boundary.
            tmp = tmp.loc[
                (tmp["overlap_bp"] >= args.min_overlap_bp)
                & (tmp["overlap_frac_boundary"] >= args.min_overlap_frac_boundary)
            ].copy()
            if tmp.empty:
                continue

            overlap_cols = [
                "sample_id", "event_id", "sv_type_norm", "sv_chr", "sv_start0", "sv_end", "sv_len",
                "bin_id", "group_key", "group_key_norm", "group_key_L2_from_obs",
                "chrom", "start0", "end",
                "raw_value", "impute_value", "is_consensus_support2",
                "overlaps_diffbound", "matching_level",
                "group_prob_in_diff", "dev_direction",
                "trajectory_class", "anchor_type", "region", "lineage_key",
                "is_exc_rg_gain", "is_astro_rg_loss",
                "is_inh_cge_gain", "is_inh_mge_gain",
                "overlap_bp", "overlap_frac_boundary", "overlap_frac_sv"
            ]
            if diag_col is not None:
                overlap_cols.append("Diagnosis")
            if sex_col is not None:
                overlap_cols.append("Sex")

            overlap_chunks.append(tmp[overlap_cols])

    if len(overlap_chunks) == 0:
        log("No overlaps found")
        overlap_df = pd.DataFrame(columns=["sample_id", "event_id", "sv_type_norm", "bin_id", "matching_level", "group_key_L2_from_obs"])
    else:
        overlap_df = pd.concat(overlap_chunks, axis=0, ignore_index=True)

    # v10: safety — overlap_df からも high-burden sample を除外 (sv は既に除外済だが
    # 冪等性確保のため)
    if not overlap_df.empty:
        n_ovl_before = overlap_df.shape[0]
        overlap_df = overlap_df.loc[~overlap_df["sample_id"].astype(str).isin(high_burden_samples)].copy()
        n_ovl_after = overlap_df.shape[0]
        if n_ovl_before != n_ovl_after:
            log(f"[v10] overlap_df high-burden excl: {n_ovl_before} -> {n_ovl_after} rows")

    # v11: 出力 suffix を v11_{PATTERN} に
    overlap_path = outdir / f"sample_boundary_event_overlap_{V11_SUFFIX}.tsv.gz"
    log(f"Writing overlap table: {overlap_path}")
    overlap_df.to_csv(overlap_path, sep="\t", index=False, compression="gzip")

    log(f"overlap_df shape: {overlap_df.shape}")
    if not overlap_df.empty:
        log("overlap_df matching_level distribution:")
        for lvl, cnt in overlap_df["matching_level"].value_counts().items():
            log(f"    {lvl:<15s}: {cnt}")

    log("Building full-cohort sample-level burden table (3 系統: unsuffixed/_strict/_l2unique)")
    sample_burden = aggregate_sample_burden_by_type(overlap_df, full_sample, "sample_id", svtype_col="sv_type_norm")

    for c in ["total_del_bases", "total_dup_bases", "total_del_count", "total_dup_count", "total_exon_del_count", "total_exon_dup_count"]:
        if c not in sample_burden.columns:
            sample_burden[c] = 0

    sample_burden["log1p_total_del_bases"] = np.log1p(sample_burden["total_del_bases"])
    sample_burden["log1p_total_dup_bases"] = np.log1p(sample_burden["total_dup_bases"])
    sample_burden["log1p_total_del_count"] = np.log1p(sample_burden["total_del_count"])
    sample_burden["log1p_total_dup_count"] = np.log1p(sample_burden["total_dup_count"])
    sample_burden["log1p_total_exon_del_count"] = np.log1p(sample_burden["total_exon_del_count"])
    sample_burden["log1p_total_exon_dup_count"] = np.log1p(sample_burden["total_exon_dup_count"])
    # gene covariates' log1p cols already present from full_sample merge, but recompute
    # to ensure consistency if burden path overwrites them.
    for svt in ("DEL", "DUP"):
        raw_col = f"total_gene_{svt}"
        log_col = f"log1p_total_gene_{svt}"
        if raw_col in sample_burden.columns:
            sample_burden[log_col] = np.log1p(sample_burden[raw_col].astype(float))

    # v11: 出力 suffix を v11_{PATTERN} に
    burden_path = outdir / f"sample_mechanistic_burden_{V11_SUFFIX}.tsv"
    log(f"Writing sample burden table: {burden_path}")
    sample_burden.to_csv(burden_path, sep="\t", index=False)

    summary = {
        "n_full_cohort_samples": int(full_sample["sample_id"].nunique()),
        "n_carrier_samples_input": int(full_sample["carrier_any_rare_sv"].sum()),
        "n_sv_events_input": int(sv["event_id"].nunique()),
        "n_overlap_rows": int(overlap_df.shape[0]),
        "n_overlap_events": int(overlap_df["event_id"].nunique()) if not overlap_df.empty else 0,
        "n_samples_with_any_overlap": int(overlap_df["sample_id"].nunique()) if not overlap_df.empty else 0,
        # v10: QC 指標
        "v10_high_burden_threshold_cnv_count": thr,
        "v10_high_burden_percentile": args.high_burden_percentile,
        "v10_n_high_burden_samples_excluded": len(high_burden_samples),
        "v10_n_sinfo_samples_after_excl": int(n_sinfo_after),
        "v10_sinfo_samples_excluded": int(n_sinfo_before - n_sinfo_after),
        # gene covariate sanity (v18 anchor: total_gene_DEL nonzero~2322 mean~0.45, total_gene_DUP nonzero~3388 mean~1.12)
        "n_samples_total_gene_DEL_gt0": int((full_sample["total_gene_DEL"] > 0).sum()),
        "n_samples_total_gene_DUP_gt0": int((full_sample["total_gene_DUP"] > 0).sum()),
        "mean_total_gene_DEL": float(full_sample["total_gene_DEL"].mean()),
        "mean_total_gene_DUP": float(full_sample["total_gene_DUP"].mean()),
        "min_overlap_bp_used": int(args.min_overlap_bp),
        "min_overlap_frac_boundary_used": float(args.min_overlap_frac_boundary),
    }

    if not overlap_df.empty:
        for lvl in ["exact", "normalized", "l2_aggregated"]:
            sub_lvl = overlap_df.loc[overlap_df["matching_level"] == lvl]
            summary[f"n_overlap_rows_{lvl}"] = int(sub_lvl.shape[0])
            summary[f"n_overlap_events_{lvl}"] = int(sub_lvl["event_id"].nunique())
            summary[f"n_samples_with_overlap_{lvl}"] = int(sub_lvl["sample_id"].nunique())
    else:
        for lvl in ["exact", "normalized", "l2_aggregated"]:
            summary[f"n_overlap_rows_{lvl}"] = 0
            summary[f"n_overlap_events_{lvl}"] = 0
            summary[f"n_samples_with_overlap_{lvl}"] = 0

    if not overlap_df.empty:
        for svt in ["DEL", "DUP"]:
            sub = overlap_df.loc[overlap_df["sv_type_norm"] == svt]
            summary[f"n_overlap_events_{svt}"] = int(sub["event_id"].nunique())
            summary[f"n_samples_with_any_overlap_{svt}"] = int(sub["sample_id"].nunique())

        for svt in ["DEL", "DUP"]:
            summary[f"n_full_samples_exc_gain_carrier_{svt}"] = int(sample_burden[f"carrier_boundary_exc_rg_gain_{svt}"].sum())
            summary[f"n_full_samples_inh_cge_gain_carrier_{svt}"] = int(sample_burden[f"carrier_boundary_inh_cge_gain_{svt}"].sum())
            summary[f"n_full_samples_inh_mge_gain_carrier_{svt}"] = int(sample_burden[f"carrier_boundary_inh_mge_gain_{svt}"].sum())
            summary[f"n_full_samples_exc_gain_carrier_strict_{svt}"] = int(sample_burden[f"carrier_boundary_exc_rg_gain_strict_{svt}"].sum())
            summary[f"n_full_samples_inh_cge_gain_carrier_strict_{svt}"] = int(sample_burden[f"carrier_boundary_inh_cge_gain_strict_{svt}"].sum())
            summary[f"n_full_samples_inh_mge_gain_carrier_strict_{svt}"] = int(sample_burden[f"carrier_boundary_inh_mge_gain_strict_{svt}"].sum())
            summary[f"n_full_samples_exc_gain_carrier_l2unique_{svt}"] = int(sample_burden[f"carrier_boundary_exc_rg_gain_l2unique_{svt}"].sum())
            summary[f"n_full_samples_inh_cge_gain_carrier_l2unique_{svt}"] = int(sample_burden[f"carrier_boundary_inh_cge_gain_l2unique_{svt}"].sum())
            summary[f"n_full_samples_inh_mge_gain_carrier_l2unique_{svt}"] = int(sample_burden[f"carrier_boundary_inh_mge_gain_l2unique_{svt}"].sum())
    else:
        zero_keys = [
            "n_overlap_events_DEL", "n_overlap_events_DUP",
            "n_samples_with_any_overlap_DEL", "n_samples_with_any_overlap_DUP",
        ]
        for svt in ["DEL", "DUP"]:
            for tier in ["", "_strict", "_l2unique"]:
                zero_keys.extend([
                    f"n_full_samples_exc_gain_carrier{tier}_{svt}",
                    f"n_full_samples_inh_cge_gain_carrier{tier}_{svt}",
                    f"n_full_samples_inh_mge_gain_carrier{tier}_{svt}",
                ])
        for key in zero_keys:
            summary[key] = 0

    # v11: pattern タグを summary に追加
    summary["v11_pattern"] = _PATTERN
    summary_df = pd.DataFrame([summary])
    # v11: 出力 suffix を v11_{PATTERN} に
    summary_path = outdir / f"sample_mechanistic_burden_summary_{V11_SUFFIX}.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    elapsed = time.time() - t0
    log(f"Done. elapsed_sec={elapsed:.2f}")


if __name__ == "__main__":
    main()
