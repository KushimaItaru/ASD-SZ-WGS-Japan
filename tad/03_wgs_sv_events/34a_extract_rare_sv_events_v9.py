#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ファイル名: 34a_extract_rare_sv_events_v9.py
# - 処理内容:
#   - v8 (Pattern A = NAHR-excluded) を継承し、1 回の実行で Pattern A/B/C の
#     3 rare SV event table を生成する。感度分析 (unified_heffel_tad_pipeline_v20
#     の Pattern B/C 解析) を新パイプライン (tad04292026/) に移植したもの。
#
#   - 共通 QC (v8 と同一、事前に全 Pattern 共通):
#       1. Annotation_mode == full
#       2. DEL/DUP のみ
#       3. autosome (chr1-chr22)
#       4. |SV_length| >= 25 kb
#       5. event-level dedup (sample_id + chrom + start0 + end + sv_type)
#       6. ASD/SZ/Healthy のみ (multi-sample row check も同時)
#       7. high-burden sample 除外 (cnv_count >= 99 percentile) = 事前 QC layer 1
#       8. segdup >= 50% 除外
#       9. common CNV one-direction overlap >= 0.30 除外
#      10. exclusion BED 被覆率 >= 50% 除外
#      11. NAHR GD locus 注釈 (is_nahr_gd_cnv 列追加) — 除外はしない
#
#   - Pattern 別フィルタ (分岐):
#       Pattern A (NAHR-excluded, main analysis):
#           is_nahr_gd_cnv == 0
#           suffix = "_patA"
#       Pattern B (NAHR-included, sensitivity):
#           全イベントをそのまま保持
#           suffix = "_patB"
#       Pattern C (NAHR-excluded + |SV_length| < 1,000,000, sensitivity):
#           is_nahr_gd_cnv == 0 かつ sv_len_abs < 1_000_000
#           suffix = "_patC"
#
#   - Pattern ごとに独立に以下を計算:
#       - PCA merge (sample_id ベース)
#       - GENCODE v46 basic exon union から per-SV の ExonU_Count
#         (exon_index は全 Pattern 共通で 1 回だけ構築)
#       - per-sample burden (DEL/DUP 別 total bases / counts / exon counts)
#         これは Pattern ごとに event set が異なるため必ず再計算
#       - 出力 TSV + summary TSV
#
#   - 出力 (output_v9/ 配下):
#       wgs_rare_sv_events_v9_patA.tsv.gz
#       wgs_rare_sv_events_v9_patA_summary.tsv
#       wgs_rare_sv_events_v9_patB.tsv.gz
#       wgs_rare_sv_events_v9_patB_summary.tsv
#       wgs_rare_sv_events_v9_patC.tsv.gz
#       wgs_rare_sv_events_v9_patC_summary.tsv
#       wgs_rare_sv_events_v9_overall_qc.tsv  (全 Pattern の共通 QC step ログ)
#
#   - v8 -> v9 変更点:
#       1. 3 Pattern を 1 回の実行で生成 (I/O 効率化, AnnotSV 読み込み 1 回)
#       2. 出力 suffix: v8 -> v9 + _patA/_patB/_patC
#       3. 出力ディレクトリ: output_v8 -> output_v9 (hardcoded override,
#          paths_v1.py は変更しない)
#       4. Pattern C の size cap: sv_len_abs < 1_000_000 bp (1 MB strict)
#       5. 共通 QC ログを overall_qc.tsv に記録 (transparency 向上)
#       6. Pattern B の is_nahr_gd_cnv 列は保持するが除外しない
#
#   - v8 と同一の挙動 (互換性維持):
#       - NAHR GD locus 判定ロジック (15q11.2 は target cov >= 50%,
#         それ以外は reciprocal overlap >= 50%)
#       - sample-level top-1% exclusion (high-burden_percentile=0.99)
#       - segdup / common CNV / exclusion BED の閾値
#       - event-level dedup
#       - multi-sample row check
#       - GENCODE basic exon union
#       - PCA 10 列の merge
#
#   - Sample-level QC (追加):
#       本スクリプトは event-level QC の一部として top-1% cnv_count サンプルを
#       除外する (v7/v8 と同一挙動)。Step 04 v10 / Step 05 v3 の sample-level
#       top-1% QC (CNV count) との関係は下記:
#         - Step 03 v9 で top-1% サンプル全 event を本 script から除外
#         - Step 04 v10 の sample-level top-1% QC は冗長だが依然有効
#           (anchor は同じ cnv_count テーブル)
#
#   - 実行時間を記録する

import argparse
import gzip
import json
import re
import sys
import time
from bisect import bisect_left, bisect_right
from intervaltree import IntervalTree, Interval
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# ----- common paths -----
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths_v1 import (
    ANNOTSV_TABLE, SAMPLE_INFO, COMMON_CNV_TABLE, SEGDUP_BED, EXCLUSION_BED,
    PCA_EIGENVEC, CNV_SAMPLE_COUNTS, GENCODE_GTF, CURATED_GD_FILE,
    PIPELINE_ROOT,
)

# ==============================================================
# v9: hardcoded output directory override
#     paths_v1.py の OUT_03_WGS_SV_EVENTS (output_v8) は変更せず、
#     v9 専用 output_v9/ を直接指定する
# ==============================================================
_OUT_V9 = PIPELINE_ROOT / "03_wgs_sv_events" / "output_v9"

# v9 専用 Pattern 定義
PATTERNS = ["patA", "patB", "patC"]
PATTERN_DESC = {
    "patA": "NAHR-excluded (main analysis, v8 behavior)",
    "patB": "NAHR-included (sensitivity)",
    "patC": "NAHR-excluded + |SV_length| < 1,000,000 bp (sensitivity)",
}
PATC_SIZE_MAX = 1_000_000  # 1 MB strict, v20 patC と同一


# ==============================================================
# Utility (v8 と同一)
# ==============================================================
def now_str() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    print(f"[{now_str()}] {msg}", file=sys.stderr, flush=True)


def detect_column(df: pd.DataFrame, candidates: List[str], required: bool = True) -> Optional[str]:
    lower_map = {str(c).lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    if required:
        raise KeyError(
            f"候補列が見つかりません: {candidates}\n実際の列: {list(df.columns)}"
        )
    return None


def read_table_auto(path: Path) -> pd.DataFrame:
    suffix = "".join(path.suffixes).lower()
    compression = "gzip" if suffix.endswith(".gz") else None
    try:
        return pd.read_csv(path, sep="\t", compression=compression, low_memory=False)
    except Exception:
        pass
    try:
        return pd.read_csv(path, sep=None, engine="python", compression=compression)
    except Exception as e:
        raise RuntimeError(f"ファイルを読めませんでした: {path}\n{e}")


def norm_chrom(x) -> str:
    s = str(x).replace("chr", "")
    if s.isdigit():
        return f"chr{int(s)}"
    if s.upper() in ("X", "Y", "M", "MT"):
        return f"chr{s.upper()}"
    return f"chr{s}"


def is_autosome_chr(chrom: str) -> bool:
    m = re.fullmatch(r"chr([0-9]+)", str(chrom))
    if not m:
        return False
    n = int(m.group(1))
    return 1 <= n <= 22


def check_multi_sample_rows(df: pd.DataFrame, col_sample: str, n_show: int = 20) -> None:
    raw_ids = df[col_sample].astype(str)
    token_counts = raw_ids.str.split(r"[;,| ]+").apply(
        lambda tokens: len([t for t in tokens if t != ""]) if isinstance(tokens, list) else 1
    )
    multi_mask = token_counts >= 2
    n_multi = int(multi_mask.sum())
    if n_multi > 0:
        examples = df.loc[multi_mask, col_sample].head(n_show).tolist()
        raise ValueError(
            f"[multi-sample row check] {n_multi} rows contain multiple sample IDs in '{col_sample}'.\n"
            f"First {min(n_show, n_multi)} examples:\n"
            + "\n".join(f"  {ex}" for ex in examples)
            + "\n\nAnnotSV full mode should yield single-sample rows."
        )


def parse_first_sample_id(x: str) -> str:
    s = str(x).strip()
    if s == "" or s.lower() == "nan":
        return s
    parts = re.split(r"[;,| ]+", s)
    parts = [p for p in parts if p != ""]
    return parts[0] if len(parts) > 0 else s


def make_event_id(sample_id: str, chrom: str, start0: int, end: int, sv_type: str) -> str:
    return f"{sample_id}|{chrom}:{start0}-{end}|{sv_type}"


def one_dir_overlap(sv_s: int, sv_e: int, db_s: int, db_e: int) -> float:
    ov = max(0, min(sv_e, db_e) - max(sv_s, db_s))
    sv_len = sv_e - sv_s
    return ov / sv_len if sv_len > 0 else 0.0


def merge_intervals_sorted(starts, ends):
    if len(starts) == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    ms, me = [starts[0]], [ends[0]]
    for s, e in zip(starts[1:], ends[1:]):
        if s <= me[-1]:
            me[-1] = max(me[-1], e)
        else:
            ms.append(s)
            me.append(e)
    return np.array(ms, dtype=np.int64), np.array(me, dtype=np.int64)


def build_merged_interval_index(intervals_by_chrom):
    index = {}
    for chrom, lst in intervals_by_chrom.items():
        arr = np.asarray(lst, dtype=np.int64)
        if arr.size == 0:
            continue
        order = np.argsort(arr[:, 0])
        ms, me = merge_intervals_sorted(arr[order, 0], arr[order, 1])
        index[chrom] = {"start": ms, "end": me}
    return index


def compute_coverage_from_index(chrom: str, sv_start: int, sv_end: int, merged_index) -> float:
    sv_len = sv_end - sv_start
    if sv_len <= 0:
        return 0.0
    if chrom not in merged_index:
        return 0.0
    ms = merged_index[chrom]["start"]
    me = merged_index[chrom]["end"]
    if ms.size == 0:
        return 0.0
    idx_start = bisect_right(ms, sv_end) - 1
    idx_end = bisect_left(ms, sv_start)
    cov = 0
    for i in range(max(0, idx_end - 1), min(len(ms), idx_start + 2)):
        ov = max(0, min(sv_end, int(me[i])) - max(sv_start, int(ms[i])))
        cov += ov
    return min(100.0, 100.0 * cov / sv_len)


def load_bed3(path: Path) -> Dict[str, List[Tuple[int, int]]]:
    out = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = norm_chrom(parts[0])
            try:
                s = int(parts[1])
                e = int(parts[2])
            except ValueError:
                continue
            if e > s:
                out[chrom].append((s, e))
    return out


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


def build_exon_index_from_gtf(gtf_path: Path):
    exon_intervals = defaultdict(list)
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "exon":
                continue
            if 'tag "basic"' not in fields[8]:
                continue
            chrom = norm_chrom(fields[0])
            try:
                s1 = int(fields[3])
                e1 = int(fields[4])
            except ValueError:
                continue
            if e1 < s1:
                s1, e1 = e1, s1
            s0 = max(0, s1 - 1)
            if e1 > s0:
                exon_intervals[chrom].append((s0, e1))
    return build_merged_interval_index(exon_intervals)


def count_union_exons_overlapped(chrom: str, sv_start: int, sv_end: int, exon_index) -> int:
    if chrom not in exon_index:
        return 0
    es = exon_index[chrom]["start"]
    ee = exon_index[chrom]["end"]
    if es.size == 0:
        return 0
    return max(
        0,
        int(np.searchsorted(es, sv_end, side="left"))
        - int(np.searchsorted(ee, sv_start, side="right"))
    )


# ==============================================================
# NAHR Genomic Disorder CNV helpers (v7 と同一)
# ==============================================================
_RE_15Q11_2 = re.compile(r"15q11\.2", re.IGNORECASE)
_RE_Q13 = re.compile(r"q13", re.IGNORECASE)


def _is_15q11_2(gd_id: str) -> bool:
    s = str(gd_id).strip()
    return bool(_RE_15Q11_2.search(s) and not _RE_Q13.search(s))


def _as_bool_gd(x) -> bool:
    s = str(x).strip().upper()
    return s in {"TRUE", "T", "1", "YES", "Y"}


def _calc_ro(s1, e1, s2, e2) -> float:
    ov = max(0, min(e1, e2) - max(s1, s2))
    l1, l2 = e1 - s1, e2 - s2
    if l1 <= 0 or l2 <= 0 or ov <= 0:
        return 0.0
    return min(ov / l1, ov / l2)


def _calc_target_cov(cnv_s, cnv_e, locus_s, locus_e) -> float:
    ov = max(0, min(cnv_e, locus_e) - max(cnv_s, locus_s))
    locus_len = locus_e - locus_s
    return ov / locus_len if locus_len > 0 else 0.0


def load_nahr_gd_loci(path: str):
    df = pd.read_csv(path, sep="\t", dtype=str)
    col_map = {c.lower(): c for c in df.columns}
    nahr_col = col_map.get("nahr", None)
    if nahr_col is None:
        raise RuntimeError(f"Column 'nahr' not found in {path}")
    df["nahr_bool"] = df[nahr_col].apply(_as_bool_gd)
    df = df[df["nahr_bool"]].copy()

    chr_col = col_map.get("chr", None)
    start_col = col_map.get("start", None)
    end_col = col_map.get("end", None)
    cnv_col = col_map.get("cnv", None)
    gd_id_col = col_map.get("gd_id", None)

    trees = {}
    for _, r in df.iterrows():
        c = str(r[chr_col]).replace("chr", "")
        c_norm = norm_chrom("chr" + c)
        try:
            s, e = int(r[start_col]), int(r[end_col])
        except Exception:
            continue
        cnv_type = str(r[cnv_col]).upper().strip()
        gid = str(r[gd_id_col]).strip()
        trees.setdefault(cnv_type, {}).setdefault(c_norm, IntervalTree()).add(Interval(s, e, gid))
    n_loci = sum(sum(len(t) for t in chrom_trees.values()) for chrom_trees in trees.values())
    return trees, n_loci


def is_nahr_gd_cnv(chrom: str, start: int, end: int, svtype: str, gd_trees: dict) -> bool:
    svtype_upper = svtype.upper().strip()
    if svtype_upper not in gd_trees:
        return False
    if chrom not in gd_trees[svtype_upper]:
        return False
    for iv in gd_trees[svtype_upper][chrom].overlap(start, end):
        gid = str(iv.data)
        if _is_15q11_2(gid):
            if _calc_target_cov(start, end, iv.begin, iv.end) >= 0.5:
                return True
        else:
            if _calc_ro(start, end, iv.begin, iv.end) >= 0.5:
                return True
    return False


# ==============================================================
# Per-pattern output helpers
# ==============================================================
def compute_per_sample_burden(sv_df: pd.DataFrame, target_sv_types: List[str]) -> pd.DataFrame:
    """Pattern ごとに、指定された event set から per-sample DEL/DUP の
    total_bases / total_count / total_exon_count を算出して返す。"""
    per_sample_svt = (
        sv_df.groupby(["sample_id", "sv_type_norm"], as_index=False)
        .agg(
            total_bases=("sv_len_abs", "sum"),
            total_count=("event_id", "nunique"),
            total_exon_count=("ExonU_Count", "sum"),
        )
    )
    wide_base = per_sample_svt.pivot(
        index="sample_id", columns="sv_type_norm", values="total_bases"
    ).fillna(0)
    wide_count = per_sample_svt.pivot(
        index="sample_id", columns="sv_type_norm", values="total_count"
    ).fillna(0)
    wide_exon = per_sample_svt.pivot(
        index="sample_id", columns="sv_type_norm", values="total_exon_count"
    ).fillna(0)

    for svt in sorted(target_sv_types):
        if svt not in wide_base.columns:
            wide_base[svt] = 0
        if svt not in wide_count.columns:
            wide_count[svt] = 0
        if svt not in wide_exon.columns:
            wide_exon[svt] = 0

    wide_base = wide_base[[svt for svt in sorted(target_sv_types)]].copy()
    wide_count = wide_count[[svt for svt in sorted(target_sv_types)]].copy()
    wide_exon = wide_exon[[svt for svt in sorted(target_sv_types)]].copy()

    wide_base.columns = [f"total_{svt.lower()}_bases" for svt in wide_base.columns]
    wide_count.columns = [f"total_{svt.lower()}_count" for svt in wide_count.columns]
    wide_exon.columns = [f"total_exon_{svt.lower()}_count" for svt in wide_exon.columns]

    return pd.concat([wide_base, wide_count, wide_exon], axis=1).reset_index()


def write_pattern_outputs(
    sv_df: pd.DataFrame,
    pattern: str,
    outdir: Path,
    target_sv_types: List[str],
    keep_dx: set,
    n_input_rows: int,
    n_nahr_total_annotated: int,
    col_id: Optional[str],
    col_tad: Optional[str],
) -> dict:
    """Pattern ごとに event table / summary を書き出し、QC stat 辞書を返す。"""
    sub = sv_df.copy()
    suffix = pattern

    # final columns (v8 と同等 + pattern 列)
    final_cols = [
        "event_id", "sample_id", "Diagnosis", "Sex", "Sex_numeric",
        "sv_type_norm", "chrom", "start0", "end", "sv_len_abs",
        "segdup_pct", "exclusion_pct", "is_common_cnv_onedir30",
        "is_high_burden_sample", "is_nahr_gd_cnv", "ExonU_Count", "cnv_count",
        "AnnotSV_ID" if col_id is not None else None,
        "TAD_coordinate" if col_tad is not None else None,
    ] + [f"PC{i}" for i in range(1, 11) if f"PC{i}" in sub.columns]
    final_cols += [c for c in [
        "total_del_bases", "total_dup_bases",
        "total_del_count", "total_dup_count",
        "total_exon_del_count", "total_exon_dup_count"
    ] if c in sub.columns]
    final_cols = [c for c in final_cols if c is not None and c in sub.columns]
    out = sub[final_cols].copy()
    out.insert(0, "pattern", pattern)  # pattern 識別子を先頭列に追加

    out_path = outdir / f"wgs_rare_sv_events_v9_{suffix}.tsv.gz"
    summary_path = outdir / f"wgs_rare_sv_events_v9_{suffix}_summary.tsv"

    log(f"[{pattern}] Writing rare SV events: {out_path}")
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")

    # summary (v8 と同じ形式、pattern 列を追加)
    summary_rows = []
    for svt in sorted(target_sv_types):
        svsub = out[out["sv_type_norm"] == svt].copy()
        base_col = f"total_{svt.lower()}_bases"
        count_col = f"total_{svt.lower()}_count"
        exon_col = f"total_exon_{svt.lower()}_count"
        row_all = {
            "pattern": pattern,
            "SV_type": svt,
            "Diagnosis": "ALL",
            "n_input_rows_total": int(n_input_rows),
            "n_nahr_total_annotated": int(n_nahr_total_annotated),
            "n_final_events": int(svsub.shape[0]),
            "n_final_samples": int(svsub["sample_id"].nunique()) if svsub.shape[0] > 0 else 0,
            "median_sv_len": float(svsub["sv_len_abs"].median()) if svsub.shape[0] > 0 else np.nan,
            "mean_sv_len": float(svsub["sv_len_abs"].mean()) if svsub.shape[0] > 0 else np.nan,
            "max_sv_len": float(svsub["sv_len_abs"].max()) if svsub.shape[0] > 0 else np.nan,
            "median_total_bases_per_sample": float(svsub[base_col].dropna().median()) if base_col in svsub.columns and svsub.shape[0] > 0 else np.nan,
            "median_total_count_per_sample": float(svsub[count_col].dropna().median()) if count_col in svsub.columns and svsub.shape[0] > 0 else np.nan,
            "median_total_exon_count_per_sample": float(svsub[exon_col].dropna().median()) if exon_col in svsub.columns and svsub.shape[0] > 0 else np.nan,
        }
        summary_rows.append(row_all)
        for dx in sorted(keep_dx):
            sub_dx = svsub.loc[svsub["Diagnosis"] == dx]
            summary_rows.append({
                "pattern": pattern,
                "SV_type": svt,
                "Diagnosis": dx,
                "n_input_rows_total": int(n_input_rows),
                "n_final_events": int(sub_dx.shape[0]),
                "n_final_samples": int(sub_dx["sample_id"].nunique()) if sub_dx.shape[0] > 0 else 0,
                "median_sv_len": float(sub_dx["sv_len_abs"].median()) if sub_dx.shape[0] > 0 else np.nan,
                "mean_sv_len": float(sub_dx["sv_len_abs"].mean()) if sub_dx.shape[0] > 0 else np.nan,
                "max_sv_len": float(sub_dx["sv_len_abs"].max()) if sub_dx.shape[0] > 0 else np.nan,
                "median_total_bases_per_sample": float(sub_dx[base_col].dropna().median()) if base_col in sub_dx.columns and sub_dx.shape[0] > 0 else np.nan,
                "median_total_count_per_sample": float(sub_dx[count_col].dropna().median()) if count_col in sub_dx.columns and sub_dx.shape[0] > 0 else np.nan,
                "median_total_exon_count_per_sample": float(sub_dx[exon_col].dropna().median()) if exon_col in sub_dx.columns and sub_dx.shape[0] > 0 else np.nan,
            })

    summary = pd.DataFrame(summary_rows)
    log(f"[{pattern}] Writing summary: {summary_path}")
    summary.to_csv(summary_path, sep="\t", index=False)

    return {
        "pattern": pattern,
        "n_events": int(out.shape[0]),
        "n_samples": int(out["sample_id"].nunique()) if out.shape[0] > 0 else 0,
        "max_sv_len": float(out["sv_len_abs"].max()) if out.shape[0] > 0 else np.nan,
        "n_events_del": int((out["sv_type_norm"] == "DEL").sum()),
        "n_events_dup": int((out["sv_type_norm"] == "DUP").sum()),
    }


# ==============================================================
# Main
# ==============================================================
def main():
    parser = argparse.ArgumentParser(
        description="Extract rare DEL/DUP event tables for Pattern A/B/C (v9, sensitivity analysis)."
    )
    parser.add_argument("--annotsv", default=str(ANNOTSV_TABLE))
    parser.add_argument("--sample-info", default=str(SAMPLE_INFO))
    parser.add_argument("--common-cnv", default=str(COMMON_CNV_TABLE))
    parser.add_argument("--segdup-bed", default=str(SEGDUP_BED))
    parser.add_argument("--exclusion-bed", default=str(EXCLUSION_BED))
    parser.add_argument("--pca", default=str(PCA_EIGENVEC))
    parser.add_argument("--cnv-sample-counts", default=str(CNV_SAMPLE_COUNTS))
    parser.add_argument("--gencode-gtf", default=str(GENCODE_GTF))
    parser.add_argument("--curated-gd-file", default=str(CURATED_GD_FILE),
                        help="Path to curated_genomic_disorder_cnv_loci_v3.txt for NAHR annotation")
    parser.add_argument("--outdir", default=str(_OUT_V9),
                        help="v9: hardcoded override of paths_v1.OUT_03_WGS_SV_EVENTS")
    parser.add_argument("--sv-types", default="DEL,DUP",
                        help="comma-separated SV types to keep (default: DEL,DUP)")
    parser.add_argument("--min-sv-len", type=int, default=25000)
    parser.add_argument("--segdup-max-pct", type=float, default=50.0)
    parser.add_argument("--common-cnv-onedir-thr", type=float, default=0.30)
    parser.add_argument("--exclusion-overlap-thr", type=float, default=0.50)
    parser.add_argument("--high-burden-percentile", type=float, default=0.99)
    parser.add_argument("--patc-size-max", type=int, default=PATC_SIZE_MAX,
                        help="Pattern C upper size limit (|SV_length| < this)")
    args = parser.parse_args()

    t0 = time.time()
    log("Start 34a_extract_rare_sv_events_v9.py (Pattern A/B/C)")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    target_sv_types = [x.strip().upper() for x in args.sv_types.split(",") if x.strip() != ""]
    if len(target_sv_types) == 0:
        raise ValueError("No valid SV types specified in --sv-types")

    # ------------------------------------------------------------------
    # Load AnnotSV
    # ------------------------------------------------------------------
    qc_steps = []
    log(f"Reading AnnotSV: {args.annotsv}")
    sv = read_table_auto(Path(args.annotsv))
    col_id = detect_column(sv, ["AnnotSV_ID"], required=False)
    col_chr = detect_column(sv, ["SV_chrom", "chrom", "chr", "CHROM"])
    col_start = detect_column(sv, ["SV_start", "start", "POS"])
    col_end = detect_column(sv, ["SV_end", "end", "END"])
    col_len = detect_column(sv, ["SV_length", "length", "sv_len"])
    col_type = detect_column(sv, ["SV_type", "sv_type", "type"])
    col_sample = detect_column(sv, ["Samples_ID", "sampleID", "sample_id", "Sample_ID"])
    col_mode = detect_column(sv, ["Annotation_mode"], required=False)
    col_tad = detect_column(sv, ["TAD_coordinate"], required=False)

    sv = sv.copy()
    n_input_rows = sv.shape[0]
    qc_steps.append({"step": "load_annotsv", "n_rows": n_input_rows})

    if col_mode is not None:
        sv = sv.loc[sv[col_mode].astype(str).str.lower() == "full"].copy()
    log(f"After full mode: {sv.shape[0]}")
    qc_steps.append({"step": "full_mode", "n_rows": sv.shape[0]})

    check_multi_sample_rows(sv, col_sample, n_show=20)
    log("Multi-sample row check passed")

    sv["sample_id"] = sv[col_sample].map(parse_first_sample_id)
    sv["chrom"] = sv[col_chr].map(norm_chrom)
    sv["start0"] = np.minimum(
        pd.to_numeric(sv[col_start], errors="coerce"),
        pd.to_numeric(sv[col_end], errors="coerce")
    )
    sv["end"] = np.maximum(
        pd.to_numeric(sv[col_start], errors="coerce"),
        pd.to_numeric(sv[col_end], errors="coerce")
    )
    sv["sv_length_raw"] = pd.to_numeric(sv[col_len], errors="coerce")
    sv["sv_len_abs"] = sv["sv_length_raw"].abs()
    sv["sv_type_norm"] = sv[col_type].astype(str).str.upper().str.strip()

    sv = sv.dropna(subset=["sample_id", "chrom", "start0", "end", "sv_len_abs"]).copy()
    sv["start0"] = sv["start0"].astype(int)
    sv["end"] = sv["end"].astype(int)
    sv = sv.loc[sv["end"] > sv["start0"]].copy()

    sv = sv.loc[sv["sv_type_norm"].isin(target_sv_types)].copy()
    log(f"After SV_type filter ({','.join(sorted(target_sv_types))}): {sv.shape[0]}")
    qc_steps.append({"step": "sv_type_filter", "n_rows": sv.shape[0]})

    sv = sv.loc[sv["chrom"].map(is_autosome_chr)].copy()
    log(f"After autosome filter: {sv.shape[0]}")
    qc_steps.append({"step": "autosome_filter", "n_rows": sv.shape[0]})

    sv = sv.loc[sv["sv_len_abs"] >= args.min_sv_len].copy()
    log(f"After |SV_length| >= {args.min_sv_len}: {sv.shape[0]}")
    qc_steps.append({"step": f"min_sv_len_{args.min_sv_len}", "n_rows": sv.shape[0]})

    dedup_cols = ["sample_id", "chrom", "start0", "end", "sv_type_norm"]
    sv = sv.sort_values(dedup_cols).drop_duplicates(subset=dedup_cols, keep="first").copy()
    sv["event_id"] = [
        make_event_id(s, c, st, en, svt)
        for s, c, st, en, svt in zip(
            sv["sample_id"], sv["chrom"], sv["start0"], sv["end"], sv["sv_type_norm"]
        )
    ]
    log(f"After event dedup: {sv.shape[0]}")
    qc_steps.append({"step": "event_dedup", "n_rows": sv.shape[0]})

    # ------------------------------------------------------------------
    # Sample info
    # ------------------------------------------------------------------
    log(f"Reading sample info: {args.sample_info}")
    sinfo = read_table_auto(Path(args.sample_info))
    id_col = detect_column(sinfo, ["SampleID", "sampleID", "sample_id", "IID", "#IID", "ID"])
    dx_col = detect_column(sinfo, ["Diagnosis", "diagnosis", "Diagnosis_DSM5", "phenotype", "disease"])
    sex_col = detect_column(sinfo, ["Sex", "sex", "gender"])
    sinfo = sinfo.rename(columns={id_col: "sample_id", dx_col: "Diagnosis", sex_col: "Sex"})
    keep_dx = {"ASD", "SZ", "Healthy"}
    sinfo = sinfo.loc[sinfo["Diagnosis"].isin(keep_dx)].copy()
    sinfo["Sex_numeric"] = sinfo["Sex"].map({"M": 1, "Male": 1, "F": 0, "Female": 0})
    sv = sv.merge(
        sinfo[["sample_id", "Diagnosis", "Sex", "Sex_numeric"]].drop_duplicates(subset=["sample_id"]),
        on="sample_id",
        how="inner"
    )
    log(f"After sample info merge (ASD/SZ/Healthy only): {sv.shape[0]}")
    qc_steps.append({"step": "sample_info_merge", "n_rows": sv.shape[0]})

    # ------------------------------------------------------------------
    # High-burden sample exclusion (top-1% cnv_count)
    # ------------------------------------------------------------------
    log(f"Reading CNV sample counts: {args.cnv_sample_counts}")
    cnt = read_table_auto(Path(args.cnv_sample_counts))
    cnt_sample_col = detect_column(cnt, ["sampleID", "sample_id", "Sample_ID", "IID"])
    cnt_count_col = detect_column(cnt, ["cnv_count", "count", "n_cnv"])
    cnt = cnt.rename(columns={cnt_sample_col: "sample_id", cnt_count_col: "cnv_count"})
    thr_99 = cnt["cnv_count"].quantile(args.high_burden_percentile)
    high_burden_samples = set(cnt.loc[cnt["cnv_count"] >= thr_99, "sample_id"].astype(str))
    sv["is_high_burden_sample"] = sv["sample_id"].isin(high_burden_samples).astype(int)
    n_before = sv.shape[0]
    sv = sv.loc[sv["is_high_burden_sample"] == 0].copy()
    log(f"After high-burden exclusion (cnv_count >= {thr_99:.0f}): {n_before} -> {sv.shape[0]}")
    qc_steps.append({"step": "high_burden_excl", "n_rows": sv.shape[0], "thr": float(thr_99)})

    # ------------------------------------------------------------------
    # segdup filter
    # ------------------------------------------------------------------
    log(f"Reading segdup BED: {args.segdup_bed}")
    segdup_raw = load_bed3(Path(args.segdup_bed))
    segdup_index = build_merged_interval_index(segdup_raw)
    segdup_pcts = np.array([
        compute_coverage_from_index(c, s, e, segdup_index)
        for c, s, e in zip(sv["chrom"], sv["start0"], sv["end"])
    ])
    sv["segdup_pct"] = segdup_pcts
    n_before = sv.shape[0]
    sv = sv.loc[sv["segdup_pct"] < args.segdup_max_pct].copy()
    log(f"After segdup filter (< {args.segdup_max_pct}%): {n_before} -> {sv.shape[0]}")
    qc_steps.append({"step": "segdup_filter", "n_rows": sv.shape[0]})

    # ------------------------------------------------------------------
    # common CNV filter
    # ------------------------------------------------------------------
    log(f"Reading common CNV file: {args.common_cnv}")
    ccnv = read_table_auto(Path(args.common_cnv))
    cchr = detect_column(ccnv, ["chrom", "chr", "SV_chrom", "CHROM"])
    cstart = detect_column(ccnv, ["start", "start0", "SV_start", "POS", "min_start"])
    cend = detect_column(ccnv, ["end", "SV_end", "END", "max_end"])
    ctype = detect_column(ccnv, ["SV_type", "sv_type", "type"], required=False)
    ccnv["chrom_norm"] = ccnv[cchr].map(norm_chrom)
    ccnv["start0"] = pd.to_numeric(ccnv[cstart], errors="coerce")
    ccnv["end"] = pd.to_numeric(ccnv[cend], errors="coerce")
    ccnv = ccnv.dropna(subset=["chrom_norm", "start0", "end"]).copy()
    ccnv["start0"] = ccnv["start0"].astype(int)
    ccnv["end"] = ccnv["end"].astype(int)
    ccnv = ccnv.loc[ccnv["end"] > ccnv["start0"]].copy()
    if ctype is not None:
        ccnv["_type"] = ccnv[ctype].astype(str).str.upper().str.strip()
    else:
        ccnv["_type"] = "*"
    ccnv_dict = defaultdict(list)
    for row in ccnv[["chrom_norm", "_type", "start0", "end"]].itertuples(index=False, name=None):
        ccnv_dict[(row[0], row[1])].append((int(row[2]), int(row[3])))

    def is_common_cnv(chrom, start0, end, svtype):
        key_exact = (chrom, svtype)
        key_any = (chrom, "*")
        if key_exact in ccnv_dict:
            for ds, de in ccnv_dict[key_exact]:
                if one_dir_overlap(start0, end, ds, de) >= args.common_cnv_onedir_thr:
                    return True
        if key_any in ccnv_dict:
            for ds, de in ccnv_dict[key_any]:
                if one_dir_overlap(start0, end, ds, de) >= args.common_cnv_onedir_thr:
                    return True
        return False

    common_mask = np.array([
        is_common_cnv(c, s, e, svt)
        for c, s, e, svt in zip(sv["chrom"], sv["start0"], sv["end"], sv["sv_type_norm"])
    ])
    sv["is_common_cnv_onedir30"] = common_mask.astype(int)
    n_before = sv.shape[0]
    sv = sv.loc[sv["is_common_cnv_onedir30"] == 0].copy()
    log(f"After common CNV filter (one-dir >= {args.common_cnv_onedir_thr}): {n_before} -> {sv.shape[0]}")
    qc_steps.append({"step": "common_cnv_filter", "n_rows": sv.shape[0]})

    # ------------------------------------------------------------------
    # exclusion BED filter
    # ------------------------------------------------------------------
    log(f"Reading exclusion BED: {args.exclusion_bed}")
    excl_raw = load_bed3(Path(args.exclusion_bed))
    excl_index = build_merged_interval_index(excl_raw)
    excl_pcts = np.array([
        compute_coverage_from_index(c, s, e, excl_index)
        for c, s, e in zip(sv["chrom"], sv["start0"], sv["end"])
    ])
    sv["exclusion_pct"] = excl_pcts
    n_before = sv.shape[0]
    sv = sv.loc[sv["exclusion_pct"] < (args.exclusion_overlap_thr * 100.0)].copy()
    log(f"After exclusion BED filter (< {args.exclusion_overlap_thr*100:.0f}%): {n_before} -> {sv.shape[0]}")
    qc_steps.append({"step": "exclusion_bed_filter", "n_rows": sv.shape[0]})

    # ------------------------------------------------------------------
    # NAHR GD annotation (event-level, not filtered here)
    # ------------------------------------------------------------------
    log(f"Reading curated GD file: {args.curated_gd_file}")
    gd_trees, n_gd_loci = load_nahr_gd_loci(args.curated_gd_file)
    log(f"Loaded {n_gd_loci} NAHR GD loci")
    nahr_mask = np.array([
        is_nahr_gd_cnv(c, s, e, svt, gd_trees)
        for c, s, e, svt in zip(sv["chrom"], sv["start0"], sv["end"], sv["sv_type_norm"])
    ])
    sv["is_nahr_gd_cnv"] = nahr_mask.astype(int)
    n_nahr_total_annotated = int(nahr_mask.sum())
    log(f"NAHR GD CNV annotated (not filtered yet): {n_nahr_total_annotated}")
    qc_steps.append({"step": "nahr_annotation", "n_rows": sv.shape[0],
                     "n_nahr_annotated": n_nahr_total_annotated})

    # ------------------------------------------------------------------
    # PCA merge (全 Pattern 共通)
    # ------------------------------------------------------------------
    log(f"Reading PCA: {args.pca}")
    pca = load_pca_file(Path(args.pca), n_pcs=10)
    sv = sv.merge(pca, on="sample_id", how="left")
    log("Merged PCA")

    # ------------------------------------------------------------------
    # GENCODE exon union (全 Pattern 共通: per-SV 値)
    # ------------------------------------------------------------------
    log(f"Building exon index from GENCODE: {args.gencode_gtf}")
    exon_index = build_exon_index_from_gtf(Path(args.gencode_gtf))
    exon_counts = np.array([
        count_union_exons_overlapped(c, s, e, exon_index)
        for c, s, e in zip(sv["chrom"], sv["start0"], sv["end"])
    ], dtype=np.int32)
    sv["ExonU_Count"] = exon_counts
    log("Computed exon-overlap counts per SV")

    # AnnotSV_ID, TAD_coordinate rename (v8 と同一)
    if col_id is not None:
        sv = sv.rename(columns={col_id: "AnnotSV_ID"})
    if col_tad is not None:
        sv = sv.rename(columns={col_tad: "TAD_coordinate"})

    cnt_keep = cnt[["sample_id", "cnv_count"]].drop_duplicates(subset=["sample_id"]).copy()
    sv = sv.merge(cnt_keep, on="sample_id", how="left")

    # ------------------------------------------------------------------
    # Pattern branching
    # ------------------------------------------------------------------
    pattern_stats = {}
    for pattern in PATTERNS:
        if pattern == "patA":
            sv_pat = sv.loc[sv["is_nahr_gd_cnv"] == 0].copy()
        elif pattern == "patB":
            sv_pat = sv.copy()
        elif pattern == "patC":
            sv_pat = sv.loc[
                (sv["is_nahr_gd_cnv"] == 0) & (sv["sv_len_abs"] < args.patc_size_max)
            ].copy()
        else:
            raise RuntimeError(f"Unknown pattern: {pattern}")

        log(f"[{pattern}] events={sv_pat.shape[0]} samples={sv_pat['sample_id'].nunique()}"
            f" ({PATTERN_DESC[pattern]})")

        # per-sample burden (Pattern-specific)
        per_sample = compute_per_sample_burden(sv_pat, target_sv_types)
        sv_pat = sv_pat.merge(per_sample, on="sample_id", how="left")

        # 出力
        stat = write_pattern_outputs(
            sv_pat, pattern, outdir, target_sv_types, keep_dx,
            n_input_rows, n_nahr_total_annotated, col_id, col_tad,
        )
        pattern_stats[pattern] = stat

    # ------------------------------------------------------------------
    # Overall QC log + summary JSON
    # ------------------------------------------------------------------
    qc_df = pd.DataFrame(qc_steps)
    qc_path = outdir / "wgs_rare_sv_events_v9_overall_qc.tsv"
    qc_df.to_csv(qc_path, sep="\t", index=False)
    log(f"Wrote overall QC log: {qc_path}")

    # JSON summary for reproducibility
    summary_json = {
        "script": "34a_extract_rare_sv_events_v9.py",
        "timestamp": now_str(),
        "args": {k: str(v) for k, v in vars(args).items()},
        "n_input_rows_total": int(n_input_rows),
        "high_burden_threshold_99pctl_cnv_count": float(thr_99),
        "n_gd_loci_nahr": int(n_gd_loci),
        "n_nahr_total_annotated": int(n_nahr_total_annotated),
        "patterns": pattern_stats,
        "patc_size_max_bp": int(args.patc_size_max),
        "qc_steps": qc_steps,
    }
    json_path = outdir / "wgs_rare_sv_events_v9_summary.json"
    with open(json_path, "w") as fh:
        json.dump(summary_json, fh, indent=2, default=float)
    log(f"Wrote summary JSON: {json_path}")

    elapsed = time.time() - t0
    log(f"Done. elapsed_sec={elapsed:.2f}")


if __name__ == "__main__":
    main()
