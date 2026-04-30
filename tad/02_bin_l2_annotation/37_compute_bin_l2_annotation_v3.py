#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
37_compute_bin_l2_annotation_v3.py
=====================================
処理内容:
  - 01_v9 master (heffel_boundary_master_v9.tsv.gz) の各 bin に対し、
    10 L2 diffbound classes (HPC_Astro 除外) の membership を付与
  - n_L2_diff_support (その bin が何個の L2 lineage で diff boundary として
    現れるか; 0-10) を計算
  - best_matching_level (exact > normalized > l2_aggregated > none) を bin 毎に抽出
  - group_primary を割り当て (S1 = paper の primary grouping):
      Diff_specific_n1 (n=1) / Diff_shared_n2plus (n>=2) / Static (diff=0)
      この union (Diff_specific_n1 ∪ Diff_shared_n2plus) が paper の Diff_any
  - 論文 Fig. 3a (10 L2 classes) と Diff_any-vs-Static specificity 解析 用の
    bin annotation table を出力
  - 実行時間を記録
入力:
  --master   : 01_v9 master tsv.gz (heffel_boundary_master_v9.tsv.gz)
  --lineage-xlsx : lineage_stage_clusters_v3.xlsx
  --outdir   : 出力ディレクトリ
出力:
  bin_l2_annotation_v3.tsv.gz (bin-level annotation, 1 row per bin)
  bin_l2_annotation_summary_v3.tsv (QC カウント)
入力/出力パス:
  入力  : F_01_BOUNDARY_MASTER (common.paths_v1)
          LINEAGE_XLSX         (common.paths_v1)
  出力先: OUT_02_BIN_L2_ANNOT
  v3 の変更点 (v2 → v3):
    - sensitivity grouping S2/S3/S4/S5 (group_s2, group_s3, group_s4, group_s5)
      を削除。paper では使用していないため (Diff_any union のみ paper integral)。
    - 関数 lbl_s2 / lbl_s3 / lbl_s4 / lbl_s5 を削除、lbl_primary のみ残す。
    - assign_group_columns() で group_primary 列のみ生成。
    - make_summary() / main() の group_cols list を group_primary のみに更新。
    - 出力 suffix を v2 → v3。
    - その他 (chunked aggregation, L2 class membership, n_L2_diff_support
      計算, matching_level rank) は v2 から完全不変。group_primary の数値結果は
      v2 と bit-identical になる想定。
"""
from __future__ import annotations
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


# =========================================================
# L2 classes (10, HPC_Astro 除外)
# =========================================================
L2_CLASSES_ALL = [
    "HPC_Astro", "HPC_Exc-CA", "HPC_Exc-DG", "HPC_Exc-ENT",
    "HPC_Inh-CGE", "HPC_Inh-MGE",
    "PFC_Astro", "PFC_Exc-DL", "PFC_Exc-UL",
    "PFC_Inh-CGE", "PFC_Inh-MGE",
]
L2_CLASSES_EXCLUDED = {"HPC_Astro"}   # v200 Methods: <500 bins で除外
L2_CLASSES = [c for c in L2_CLASSES_ALL if c not in L2_CLASSES_EXCLUDED]
assert len(L2_CLASSES) == 10, "Expected 10 L2 classes after HPC_Astro exclusion"

MATCHING_LEVEL_RANK = {
    "exact":          3,
    "normalized":     2,
    "l2_aggregated":  1,
    "none":           0,
    "":               0,
    None:             0,
}


def load_l2_classes_from_xlsx(xlsx_path: Path) -> list[str]:
    """lineage_stage_clusters_v3.xlsx の Summary シートから lineage リストを読む。
       期待: 11 個を返し、本スクリプト側で HPC_Astro を除外する。
    """
    try:
        import openpyxl  # noqa: F401
    except ImportError:
        log("openpyxl が無いため hard-coded L2 リストを使用します。")
        return list(L2_CLASSES_ALL)

    try:
        df = pd.read_excel(xlsx_path, sheet_name="Summary", engine="openpyxl")
    except Exception as e:
        log(f"xlsx 読み込み失敗 ({e})。hard-coded L2 リストで続行します。")
        return list(L2_CLASSES_ALL)

    if "Lineage (L2_diff)" not in df.columns:
        log(f"'Lineage (L2_diff)' 列が見つからない。hard-coded を使用。")
        return list(L2_CLASSES_ALL)

    xlsx_classes = (
        df["Lineage (L2_diff)"]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )
    xlsx_classes = [c for c in xlsx_classes if c in set(L2_CLASSES_ALL)]
    if set(xlsx_classes) != set(L2_CLASSES_ALL):
        log(
            f"xlsx 由来 lineage と hard-coded に差: "
            f"xlsx={sorted(set(xlsx_classes))}, hard={sorted(L2_CLASSES_ALL)}"
        )
    return xlsx_classes if xlsx_classes else list(L2_CLASSES_ALL)


# ----- common paths -----
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths_v1 import (
    F_01_BOUNDARY_MASTER,
    LINEAGE_XLSX,
    OUT_02_BIN_L2_ANNOT,
    ensure_output_dirs,
)


def parse_args() -> argparse.Namespace:
    """All arguments default to common/paths_v1 locations; argparse is retained
    only as an override hook so the script can still be run with overrides."""
    p = argparse.ArgumentParser(description="Build bin-level L2 class annotation from 01_v9 master.")
    p.add_argument("--master", default=F_01_BOUNDARY_MASTER, type=Path,
                   help="Path to heffel_boundary_master_v9.tsv.gz (33_v9 output); "
                        "defaults to common.paths_v1.F_01_BOUNDARY_MASTER")
    p.add_argument("--lineage-xlsx", default=LINEAGE_XLSX, type=Path,
                   help="Path to lineage_stage_clusters_v3.xlsx; "
                        "defaults to common.paths_v1.LINEAGE_XLSX")
    p.add_argument("--outdir", default=OUT_02_BIN_L2_ANNOT, type=Path,
                   help="Output directory (will be created); "
                        "defaults to common.paths_v1.OUT_02_BIN_L2_ANNOT")
    p.add_argument("--chunksize", type=int, default=500_000,
                   help="pandas read_csv chunksize for streaming aggregation.")
    return p.parse_args()


def aggregate_master_streaming(master_path: Path, chunksize: int) -> pd.DataFrame:
    """
    01_v9 master を chunk 読みして bin 単位に集約する。
    必要列: bin_id, chrom, start0, end, overlaps_diffbound, matching_level, lineage_key,
           raw_value, impute_value
    返り値 columns:
      bin_id, chrom, start0, end,
      overlaps_diffbound_any (0/1),
      raw_value_max, impute_value_max,
      best_matching_level,
      membership_<L2_CLASS> (0/1) for each of 10 L2_CLASSES,
      n_L2_diff_support (int 0-10)
    """
    use_cols = [
        "bin_id", "chrom", "start0", "end",
        "overlaps_diffbound", "matching_level", "lineage_key",
        "raw_value", "impute_value",
    ]
    dtype = {
        "bin_id": "string",
        "chrom": "string",
        "start0": "int64",
        "end": "int64",
        "overlaps_diffbound": "Int8",
        "matching_level": "string",
        "lineage_key": "string",
        "raw_value": "Int16",
        "impute_value": "Float32",
    }

    # bin_id -> coords
    coords: dict[str, tuple[str, int, int]] = {}
    # bin_id -> int (any diffbound row)
    diff_any: dict[str, int] = {}
    # bin_id -> set(lineage)
    bin_lineages: dict[str, set[str]] = {}
    # bin_id -> best matching rank + label
    best_ml_rank: dict[str, int] = {}
    best_ml_label: dict[str, str] = {}
    # bin_id -> max raw_value / impute_value
    raw_max: dict[str, int] = {}
    imp_max: dict[str, float] = {}

    total_rows = 0
    l2_class_set = set(L2_CLASSES)

    log(f"Reading master (chunksize={chunksize:,}) ...")
    for i, chunk in enumerate(
        pd.read_csv(
            master_path, sep="\t", compression="gzip",
            usecols=use_cols, dtype=dtype, chunksize=chunksize,
            na_values=["", "NA", "nan"],
        ),
        start=1,
    ):
        total_rows += len(chunk)

        # coords (first seen wins)
        unseen = chunk.loc[~chunk["bin_id"].isin(coords), ["bin_id", "chrom", "start0", "end"]].drop_duplicates("bin_id")
        for _, r in unseen.iterrows():
            coords[r["bin_id"]] = (r["chrom"], int(r["start0"]), int(r["end"]))

        # overlaps_diffbound (max)
        diff_chunk = chunk.groupby("bin_id")["overlaps_diffbound"].max()
        for bin_id, val in diff_chunk.items():
            if pd.isna(val):
                continue
            v = int(val)
            if v > diff_any.get(bin_id, 0):
                diff_any[bin_id] = v

        # lineage membership (diffbound rows only, excluded classes skipped)
        diff_rows = chunk[
            (chunk["overlaps_diffbound"] == 1)
            & (chunk["lineage_key"].notna())
            & (chunk["lineage_key"].isin(l2_class_set))
        ][["bin_id", "lineage_key"]]
        for bin_id, lin in zip(diff_rows["bin_id"].tolist(), diff_rows["lineage_key"].tolist()):
            s = bin_lineages.get(bin_id)
            if s is None:
                s = set()
                bin_lineages[bin_id] = s
            s.add(lin)

        # best matching_level
        chunk_ml = chunk[["bin_id", "matching_level"]].copy()
        chunk_ml["rank"] = chunk_ml["matching_level"].map(MATCHING_LEVEL_RANK).fillna(0).astype(int)
        idx = chunk_ml.groupby("bin_id")["rank"].idxmax()
        best_rows = chunk_ml.loc[idx]
        for bin_id, rank, label in zip(best_rows["bin_id"].tolist(),
                                       best_rows["rank"].tolist(),
                                       best_rows["matching_level"].tolist()):
            if rank > best_ml_rank.get(bin_id, -1):
                best_ml_rank[bin_id] = int(rank)
                best_ml_label[bin_id] = "none" if (label is None or label == "" or pd.isna(label)) else str(label)

        # raw / impute max
        rmax = chunk.groupby("bin_id")["raw_value"].max()
        for bin_id, val in rmax.items():
            if pd.isna(val):
                continue
            v = int(val)
            if v > raw_max.get(bin_id, -1):
                raw_max[bin_id] = v
        imax = chunk.groupby("bin_id")["impute_value"].max()
        for bin_id, val in imax.items():
            if pd.isna(val):
                continue
            v = float(val)
            if v > imp_max.get(bin_id, float("-inf")):
                imp_max[bin_id] = v

        log(f"  chunk {i}: rows={len(chunk):,} (cumulative {total_rows:,}); "
            f"bins seen so far={len(coords):,}")

    log(f"Total rows read: {total_rows:,}; unique bins: {len(coords):,}")

    # Build DataFrame
    bin_ids = sorted(coords.keys())
    records = []
    for b in bin_ids:
        c, s, e = coords[b]
        lineages = bin_lineages.get(b, set())
        row = {
            "bin_id": b,
            "chrom": c,
            "start0": s,
            "end": e,
            "overlaps_diffbound_any": diff_any.get(b, 0),
            "raw_value_max": raw_max.get(b, 0),
            "impute_value_max": imp_max.get(b, float("nan")),
            "best_matching_level": best_ml_label.get(b, "none"),
            "n_L2_diff_support": len(lineages),
        }
        for cls in L2_CLASSES:
            row[f"membership_{cls}"] = 1 if cls in lineages else 0
        records.append(row)
    return pd.DataFrame.from_records(records)


def assign_group_columns(anno: pd.DataFrame) -> pd.DataFrame:
    """Assign primary group label (S1 = paper's primary grouping).

    v3 で S2-S5 sensitivity grouping を削除済 (paper 不使用)。
    group_primary の union (Diff_specific_n1 ∪ Diff_shared_n2plus) が paper の
    Diff_any boundary set 定義 (Supp Methods Para 251 で規定)。
    """
    n = anno["n_L2_diff_support"].astype(int)
    diff = anno["overlaps_diffbound_any"].astype(int)

    def lbl_primary(di: int, ni: int) -> str:
        if di == 0:
            return "Static"
        if ni == 1:
            return "Diff_specific_n1"
        if ni >= 2:
            return "Diff_shared_n2plus"
        return "Static"

    anno["group_primary"] = [lbl_primary(di, ni) for di, ni in zip(diff, n)]
    return anno


def make_summary(anno: pd.DataFrame) -> pd.DataFrame:
    """QC: 各グループ定義でのカウント集計。"""
    rows = []
    rows.append(("total_bins", len(anno)))
    rows.append(("diff_bins", int((anno["overlaps_diffbound_any"] == 1).sum())))
    rows.append(("static_bins", int((anno["overlaps_diffbound_any"] == 0).sum())))
    for n in range(11):
        rows.append((f"n_L2_diff_support_eq_{n}",
                     int((anno["n_L2_diff_support"] == n).sum())))
    for cls in L2_CLASSES:
        rows.append((f"membership_{cls}",
                     int(anno[f"membership_{cls}"].sum())))
    # v3: group_primary のみ (S2-S5 削除済)
    for col in ["group_primary"]:
        vc = anno[col].value_counts()
        for k, v in vc.items():
            rows.append((f"{col}__{k}", int(v)))
    # matching_level QC
    vc_ml = anno["best_matching_level"].value_counts()
    for k, v in vc_ml.items():
        rows.append((f"best_matching_level__{k}", int(v)))
    return pd.DataFrame(rows, columns=["metric", "count"])


def main() -> None:
    t0 = time.time()
    ensure_output_dirs()

    args = parse_args()
    master_path  = args.master
    lineage_xlsx = args.lineage_xlsx
    outdir       = args.outdir
    chunksize    = args.chunksize

    log("=" * 60)
    log("Start 37_compute_bin_l2_annotation_v3.py (S2-S5 sensitivity grouping removed)")
    log(f"  master: {master_path}")
    log(f"  lineage_xlsx: {lineage_xlsx}")
    log(f"  outdir: {outdir}")

    xlsx_classes = load_l2_classes_from_xlsx(lineage_xlsx)
    log(f"L2 classes from xlsx (pre-filter): {xlsx_classes}")
    log(f"L2 classes used (n={len(L2_CLASSES)}, HPC_Astro excluded): {L2_CLASSES}")

    anno = aggregate_master_streaming(master_path, chunksize=chunksize)
    log(f"Aggregated bin table: {anno.shape}")

    anno = assign_group_columns(anno)

    # reorder columns (v3: group_primary のみ)
    id_cols = ["bin_id", "chrom", "start0", "end"]
    flag_cols = ["overlaps_diffbound_any", "raw_value_max", "impute_value_max",
                 "best_matching_level", "n_L2_diff_support"]
    membership_cols = [f"membership_{c}" for c in L2_CLASSES]
    group_cols = ["group_primary"]
    anno = anno[id_cols + flag_cols + membership_cols + group_cols]

    out_anno = outdir / "bin_l2_annotation_v3.tsv.gz"
    log(f"Writing annotation: {out_anno}")
    anno.to_csv(out_anno, sep="\t", index=False, compression="gzip")

    summary = make_summary(anno)
    out_sum = outdir / "bin_l2_annotation_summary_v3.tsv"
    log(f"Writing summary: {out_sum}")
    summary.to_csv(out_sum, sep="\t", index=False)

    # preview log
    log("Group primary counts:")
    for k, v in anno["group_primary"].value_counts().items():
        log(f"  {k}: {v:,}")
    log(f"Diff bins (overlaps_diffbound_any=1): {(anno['overlaps_diffbound_any']==1).sum():,}")
    log(f"  specific (n=1): {((anno['n_L2_diff_support']==1) & (anno['overlaps_diffbound_any']==1)).sum():,}")
    log(f"  shared (n>=2): {((anno['n_L2_diff_support']>=2) & (anno['overlaps_diffbound_any']==1)).sum():,}")
    log(f"Static bins: {(anno['overlaps_diffbound_any']==0).sum():,}")

    elapsed = time.time() - t0
    log(f"Done. Elapsed: {elapsed:.2f}s = {elapsed/60:.2f}min")
    log("=" * 60)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        import traceback
        traceback.print_exc()
        sys.exit(1)
