#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ファイル名: 34a_extract_rare_sv_events_v8.py
# - 処理内容:
#   - AnnotSV の autosomal SV テーブルから、論文本番 TAD 解析（tad_label_permutation_v8_allmodels.sh）に整合する rare SV event table を作成する
#   - Annotation_mode == full、DEL/DUP のみ、|SV_length| >= 25kb、ASD/SZ/Healthy のみを抽出する
#   - high-burden sample を cnv_count の 0.99 quantile 以上で除外する
#   - segdup 被覆率 < 50%、common CNV one-direction overlap < 0.30、exclusion BED 被覆率 < 50% で QC を揃える
#   - common CNV テーブルの min_start / max_end 列にも対応する
#   - GENCODE v46 basic exon union を使って SV ごとの exon overlap count を計算する
#   - sample ごとの DEL/DUP 別 total bases / total counts / total exon counts を付与する
#   - DEL と DUP を1つの event table にまとめて出力し、後続の overlap / regression で SV_type ごとに分けて解析できる形にする
#   - 列名は固定列番号ではなく列名ベースで動的に選択する
#   - v8 変更点:
#     - tad04212026/ パイプラインに移行
#     - パスを common/paths_v1.py で一元管理
#     - argparse は維持（必要に応じて override 可能）だが全オプションに default を設定
#     - 引数なしで `python 34a_extract_rare_sv_events_v8.py` で実行可能
#     - 出力 suffix を v7 → v8 に更新
#   - v7 変更点:
#     - NAHR Genomic Disorder CNV除外フィルターを追加（curated_genomic_disorder_cnv_loci_v3.txt の NAHR=True ローカスとオーバーラップするCNVを除外）
#     - 15q11.2 はtarget coverage >= 50%、それ以外は reciprocal overlap >= 50% で判定
#     - intervaltree を使用した高速オーバーラップ判定
#   - v5 変更点:
#     - args.segdup-bed の構文バグを修正（args.segdup_bed に統一）
#     - multi-sample row（Samples_ID に複数サンプルが混入）を検出し、存在すれば ValueError で停止する
#     - summary TSV に diagnosis 別の event 数・unique sample 数を追加
#   - 実行時間を記録する

import argparse
import gzip
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
    OUT_03_WGS_SV_EVENTS, ensure_output_dirs,
)


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
    """
    v5 追加: Samples_ID 列に複数サンプルが混入している行を検出し、
    存在すれば ValueError で停止する。
    parse_first_sample_id() を安全に使うための事前チェック。
    """
    raw_ids = df[col_sample].astype(str)
    multi_mask = raw_ids.str.contains(r"[;,| ]{1,}", regex=True, na=False)
    # 空白1個はサンプルID内部にないと仮定し、split して token 数を確認
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
            + "\n\nAnnotSV full mode should yield single-sample rows. "
            "Please check the input or filtering logic."
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


# =====================================================================
# NAHR Genomic Disorder CNV exclusion helpers (v7)
# =====================================================================
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
    """NAHR=True のGDローカスを読み込み、IntervalTreeとして返す"""
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
    """CNVがNAHR GDローカスにマッチするかを判定"""
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



def main():
    ensure_output_dirs()
    parser = argparse.ArgumentParser(
        description="Extract rare DEL/DUP event table aligned to original TAD QC pipeline."
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
                        help="Path to curated_genomic_disorder_cnv_loci_v3.txt for NAHR exclusion")
    parser.add_argument("--outdir", default=str(OUT_03_WGS_SV_EVENTS))
    parser.add_argument("--sv-types", default="DEL,DUP",
                        help="comma-separated SV types to keep (default: DEL,DUP)")
    parser.add_argument("--min-sv-len", type=int, default=25000)
    parser.add_argument("--segdup-max-pct", type=float, default=50.0)
    parser.add_argument("--common-cnv-onedir-thr", type=float, default=0.30)
    parser.add_argument("--exclusion-overlap-thr", type=float, default=0.50)
    parser.add_argument("--high-burden-percentile", type=float, default=0.99)
    args = parser.parse_args()

    t0 = time.time()
    log("Start 34a_extract_rare_sv_events_v8.py (NAHR GD excluded)")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    target_sv_types = {x.strip().upper() for x in args.sv_types.split(",") if x.strip() != ""}
    if len(target_sv_types) == 0:
        raise ValueError("No valid SV types specified in --sv-types")

    # ------------------------------------------------------------------
    # Load AnnotSV
    # ------------------------------------------------------------------
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

    if col_mode is not None:
        sv = sv.loc[sv[col_mode].astype(str).str.lower() == "full"].copy()
    log(f"After full mode: {sv.shape[0]}")

    # v5: multi-sample row check (after full mode filter, before parse_first_sample_id)
    check_multi_sample_rows(sv, col_sample, n_show=20)
    log("Multi-sample row check passed: all rows contain single sample ID")

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

    sv = sv.loc[sv["chrom"].map(is_autosome_chr)].copy()
    log(f"After autosome filter: {sv.shape[0]}")

    sv = sv.loc[sv["sv_len_abs"] >= args.min_sv_len].copy()
    log(f"After |SV_length| >= {args.min_sv_len}: {sv.shape[0]}")

    dedup_cols = ["sample_id", "chrom", "start0", "end", "sv_type_norm"]
    sv = sv.sort_values(dedup_cols).drop_duplicates(subset=dedup_cols, keep="first").copy()
    sv["event_id"] = [
        make_event_id(s, c, st, en, svt)
        for s, c, st, en, svt in zip(
            sv["sample_id"], sv["chrom"], sv["start0"], sv["end"], sv["sv_type_norm"]
        )
    ]
    log(f"After event dedup: {sv.shape[0]}")
    log("SV counts after dedup by type:")
    print(sv["sv_type_norm"].value_counts().sort_index(), file=sys.stderr)

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
    log("SV counts by diagnosis and type:")
    print(sv.groupby(["Diagnosis", "sv_type_norm"]).size(), file=sys.stderr)

    # ------------------------------------------------------------------
    # High-burden exclusion
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
    log("SV counts by type after high-burden exclusion:")
    print(sv["sv_type_norm"].value_counts().sort_index(), file=sys.stderr)

    # ------------------------------------------------------------------
    # segdup (v5: fixed args.segdup_bed)
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
    log("SV counts by type after segdup filter:")
    print(sv["sv_type_norm"].value_counts().sort_index(), file=sys.stderr)

    # ------------------------------------------------------------------
    # common CNV
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
    log("SV counts by type after common CNV filter:")
    print(sv["sv_type_norm"].value_counts().sort_index(), file=sys.stderr)

    # ------------------------------------------------------------------
    # exclusion BED
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
    log("SV counts by type after exclusion BED filter:")
    print(sv["sv_type_norm"].value_counts().sort_index(), file=sys.stderr)


    # ------------------------------------------------------------------
    # NAHR Genomic Disorder CNV exclusion (v7)
    # ------------------------------------------------------------------
    log(f"Reading curated GD file: {args.curated_gd_file}")
    gd_trees, n_gd_loci = load_nahr_gd_loci(args.curated_gd_file)
    log(f"Loaded {n_gd_loci} NAHR GD loci")
    nahr_mask = np.array([
        is_nahr_gd_cnv(c, s, e, svt, gd_trees)
        for c, s, e, svt in zip(sv["chrom"], sv["start0"], sv["end"], sv["sv_type_norm"])
    ])
    sv["is_nahr_gd_cnv"] = nahr_mask.astype(int)
    n_before_nahr = sv.shape[0]
    n_nahr_excluded = int(nahr_mask.sum())
    sv = sv.loc[sv["is_nahr_gd_cnv"] == 0].copy()
    log(f"After NAHR GD exclusion: {n_before_nahr} -> {sv.shape[0]} (excluded {n_nahr_excluded})")
    log("SV counts by type after NAHR GD exclusion:")
    print(sv["sv_type_norm"].value_counts().sort_index(), file=sys.stderr)

    # ------------------------------------------------------------------
    # PCA
    # ------------------------------------------------------------------
    log(f"Reading PCA: {args.pca}")
    pca = load_pca_file(Path(args.pca), n_pcs=10)
    sv = sv.merge(pca, on="sample_id", how="left")
    log("Merged PCA")

    # ------------------------------------------------------------------
    # GENCODE exon union
    # ------------------------------------------------------------------
    log(f"Building exon index from GENCODE: {args.gencode_gtf}")
    exon_index = build_exon_index_from_gtf(Path(args.gencode_gtf))
    exon_counts = np.array([
        count_union_exons_overlapped(c, s, e, exon_index)
        for c, s, e in zip(sv["chrom"], sv["start0"], sv["end"])
    ], dtype=np.int32)
    sv["ExonU_Count"] = exon_counts
    log("Computed exon-overlap counts per SV")

    # ------------------------------------------------------------------
    # Per-sample burden by SV type
    # ------------------------------------------------------------------
    per_sample_svt = (
        sv.groupby(["sample_id", "sv_type_norm"], as_index=False)
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

    per_sample = pd.concat([wide_base, wide_count, wide_exon], axis=1).reset_index()
    sv = sv.merge(per_sample, on="sample_id", how="left")

    cnt_keep = cnt[["sample_id", "cnv_count"]].drop_duplicates(subset=["sample_id"]).copy()
    sv = sv.merge(cnt_keep, on="sample_id", how="left")

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------
    final_cols = [
        "event_id", "sample_id", "Diagnosis", "Sex", "Sex_numeric",
        "sv_type_norm", "chrom", "start0", "end", "sv_len_abs",
        "segdup_pct", "exclusion_pct", "is_common_cnv_onedir30",
        "is_high_burden_sample", "is_nahr_gd_cnv", "ExonU_Count", "cnv_count",
        "AnnotSV_ID" if col_id is not None else None,
        "TAD_coordinate" if col_tad is not None else None,
    ] + [f"PC{i}" for i in range(1, 11) if f"PC{i}" in sv.columns]
    final_cols += [c for c in [
        "total_del_bases", "total_dup_bases",
        "total_del_count", "total_dup_count",
        "total_exon_del_count", "total_exon_dup_count"
    ] if c in sv.columns]

    if col_id is not None:
        sv = sv.rename(columns={col_id: "AnnotSV_ID"})
    if col_tad is not None:
        sv = sv.rename(columns={col_tad: "TAD_coordinate"})

    final_cols = [c for c in final_cols if c is not None and c in sv.columns]
    out = sv[final_cols].copy()

    out_path = outdir / "wgs_rare_sv_events_v8.tsv.gz"
    summary_path = outdir / "wgs_rare_sv_events_v8_summary.tsv"

    log(f"Writing rare SV events: {out_path}")
    out.to_csv(out_path, sep="\t", index=False, compression="gzip")

    # Summary (v5: enhanced with diagnosis-level breakdown)
    summary_rows = []
    for svt in sorted(target_sv_types):
        sub = out[out["sv_type_norm"] == svt].copy()
        base_col = f"total_{svt.lower()}_bases"
        count_col = f"total_{svt.lower()}_count"
        exon_col = f"total_exon_{svt.lower()}_count"
        summary_rows.append({
            "SV_type": svt,
            "Diagnosis": "ALL",
            "n_input_rows_total": int(n_input_rows),
            "n_nahr_excluded": int(n_nahr_excluded),
            "n_final_events": int(sub.shape[0]),
            "n_final_samples": int(sub["sample_id"].nunique()) if sub.shape[0] > 0 else 0,
            "median_sv_len": float(sub["sv_len_abs"].median()) if sub.shape[0] > 0 else np.nan,
            "mean_sv_len": float(sub["sv_len_abs"].mean()) if sub.shape[0] > 0 else np.nan,
            "median_total_bases_per_sample": float(sub[base_col].dropna().median()) if base_col in sub.columns and sub.shape[0] > 0 else np.nan,
            "median_total_count_per_sample": float(sub[count_col].dropna().median()) if count_col in sub.columns and sub.shape[0] > 0 else np.nan,
            "median_total_exon_count_per_sample": float(sub[exon_col].dropna().median()) if exon_col in sub.columns and sub.shape[0] > 0 else np.nan,
        })
        # v5: per-diagnosis breakdown
        for dx in sorted(keep_dx):
            sub_dx = sub.loc[sub["Diagnosis"] == dx]
            summary_rows.append({
                "SV_type": svt,
                "Diagnosis": dx,
                "n_input_rows_total": int(n_input_rows),
                "n_final_events": int(sub_dx.shape[0]),
                "n_final_samples": int(sub_dx["sample_id"].nunique()) if sub_dx.shape[0] > 0 else 0,
                "median_sv_len": float(sub_dx["sv_len_abs"].median()) if sub_dx.shape[0] > 0 else np.nan,
                "mean_sv_len": float(sub_dx["sv_len_abs"].mean()) if sub_dx.shape[0] > 0 else np.nan,
                "median_total_bases_per_sample": float(sub_dx[base_col].dropna().median()) if base_col in sub_dx.columns and sub_dx.shape[0] > 0 else np.nan,
                "median_total_count_per_sample": float(sub_dx[count_col].dropna().median()) if count_col in sub_dx.columns and sub_dx.shape[0] > 0 else np.nan,
                "median_total_exon_count_per_sample": float(sub_dx[exon_col].dropna().median()) if exon_col in sub_dx.columns and sub_dx.shape[0] > 0 else np.nan,
            })

    summary = pd.DataFrame(summary_rows)
    log(f"Writing summary: {summary_path}")
    summary.to_csv(summary_path, sep="\t", index=False)

    elapsed = time.time() - t0
    log(f"Done. elapsed_sec={elapsed:.2f}")


if __name__ == "__main__":
    main()
