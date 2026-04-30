#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
38_compute_sample_burden_L2_and_specificity_v4.py
=====================================================
処理内容:
  - [v4 変更点] S2-S5 sensitivity grouping を削除 (paper 不使用):
    * GROUP_COLUMNS から group_s2 / group_s3 / group_s4 / group_s5 を削除、
      group_primary のみ残す。
    * 入力 bin annotation を 02_v2 → 02_v3 (S2-S5 列削除済) に切替。
      paths_v1.py 更新で F_02_BIN_L2_ANNOTATION が v3 を指すため自動的に切替。
    * 出力 dir / suffix を v3 → v4 に更新。
    * 統計・QC ロジック (DEL/DUP, n_boundary/n_events/carrier_boundary,
      sample-level top-1% QC 適用済 input) は v3 と完全一致。
    * group_primary に該当する数値結果は v3 と bit-identical になる想定。
  - [v3 変更点 (継承)] 入力を Script 04 v10 (sample-level top-1% CNV count QC
    適用済み) のファイルに切り替え:
    * event overlap : output_v10/sample_boundary_event_overlap_v10.tsv.gz
    * covariate burden: output_v10/sample_mechanistic_burden_v10.tsv
  - 04_v10 event overlap (sample_boundary_event_overlap_v10.tsv.gz) と
    02_v3 の bin-level annotation (bin_l2_annotation_v3.tsv.gz) を join し、
    サンプル x 10 L2 classes x (DEL / DUP) の bin-count / event-count /
    carrier burden を計算する + group_primary (S1) の Diff_specific_n1 /
    Diff_shared_n2plus / Static を計算する
  - 04_v10 の既存 burden (sample_mechanistic_burden_v10.tsv) から共変量
    (Diagnosis, Sex, Sex_numeric, PC1-10, total_bases/count, total_gene, log1p_*)
    を継承して最終テーブルに結合
  - 10 L2 classes (HPC_Astro 除外):
      HPC_Exc-CA, HPC_Exc-DG, HPC_Exc-ENT, HPC_Inh-CGE, HPC_Inh-MGE,
      PFC_Astro, PFC_Exc-DL, PFC_Exc-UL, PFC_Inh-CGE, PFC_Inh-MGE
  - group_primary (paper の primary specificity grouping):
      Diff_specific_n1 / Diff_shared_n2plus / Static
      union (Diff_specific_n1 ∪ Diff_shared_n2plus) が paper の Diff_any
  - Exposure の 3 種類を出力:
      bin-count   : n_boundary_<label>_<svtype>
      event-count : n_events_<label>_<svtype>
      carrier     : carrier_boundary_<label>_<svtype>
  - 実行時間記録
入力:
  --event-overlap   : sample_boundary_event_overlap_v10.tsv.gz (04_v10 output)
  --bin-annotation  : bin_l2_annotation_v3.tsv.gz (02_v3 output, group_primary only)
  --covariate-burden: sample_mechanistic_burden_v10.tsv (04_v10 output, 共変量継承用)
  --outdir          : 出力ディレクトリ (output_v4/)
出力:
  sample_burden_L2_and_specificity_v4.tsv
  sample_burden_L2_and_specificity_summary_v4.tsv (QC)
"""
from __future__ import annotations
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ----- common paths -----
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths_v1 import (
    F_02_BIN_L2_ANNOTATION,
    PIPELINE_ROOT,
    ensure_output_dirs,
)

# =========================================================
# v4: input は 04_v10 (sample-level top-1% QC 適用済), output は v4 dir
# =========================================================
OUT_04_V10 = PIPELINE_ROOT / "04_wgs_sv_boundary_overlap" / "output_v10"
F_04_V10_EVENT_OVERLAP = OUT_04_V10 / "sample_boundary_event_overlap_v10.tsv.gz"
F_04_V10_SAMPLE_BURDEN_COV = OUT_04_V10 / "sample_mechanistic_burden_v10.tsv"

# v4 output dir
OUT_05_V4 = PIPELINE_ROOT / "05_wgs_sample_burden" / "output_v4"


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


# =========================================================
# 共通: L2 classes と group 定義
# v4 変更: GROUP_COLUMNS から S2-S5 削除 (paper 不使用)
# =========================================================
L2_CLASSES = [
    "HPC_Exc-CA", "HPC_Exc-DG", "HPC_Exc-ENT",
    "HPC_Inh-CGE", "HPC_Inh-MGE",
    "PFC_Astro", "PFC_Exc-DL", "PFC_Exc-UL",
    "PFC_Inh-CGE", "PFC_Inh-MGE",
]
GROUP_COLUMNS = ["group_primary"]  # v4: S2-S5 削除済 (paper 不使用)
SV_TYPES = ["DEL", "DUP"]

# 共変量として継承する列名 (04_v10 burden 由来)
COVARIATE_COLS = [
    "sample_id",
    "Diagnosis", "Sex", "Sex_numeric",
    "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
    "total_del_bases", "total_dup_bases",
    "total_del_count", "total_dup_count",
    "total_exon_del_count", "total_exon_dup_count",
    "cnv_count", "carrier_any_rare_sv",
    "log1p_total_del_bases", "log1p_total_dup_bases",
    "log1p_total_del_count", "log1p_total_dup_count",
    "log1p_total_exon_del_count", "log1p_total_exon_dup_count",
    "total_gene_DEL", "log1p_total_gene_DEL",
    "total_gene_DUP", "log1p_total_gene_DUP",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute sample x (10 L2 classes + group_primary) x SV type burden "
                    "(v4: input=04_v10 + 02_v3; S2-S5 sensitivity grouping removed)."
    )
    p.add_argument("--event-overlap", default=F_04_V10_EVENT_OVERLAP, type=Path,
                   help="Path to sample_boundary_event_overlap_v10.tsv.gz (04_v10).")
    p.add_argument("--bin-annotation", default=F_02_BIN_L2_ANNOTATION, type=Path,
                   help="Path to bin_l2_annotation_v3.tsv.gz (02_v3, group_primary only).")
    p.add_argument("--covariate-burden", default=F_04_V10_SAMPLE_BURDEN_COV, type=Path,
                   help="Path to sample_mechanistic_burden_v10.tsv (04_v10, for covariates).")
    p.add_argument("--outdir", default=OUT_05_V4, type=Path,
                   help="Output directory (output_v4/, will be created).")
    p.add_argument("--chunksize", type=int, default=500_000,
                   help="Chunk size for event overlap reading (default: 500000).")
    return p.parse_args()


def load_bin_annotation(path: Path) -> pd.DataFrame:
    log(f"Loading bin annotation: {path}")
    anno = pd.read_csv(path, sep="\t", compression="gzip",
                       dtype={"bin_id": "string",
                              "chrom": "string",
                              "best_matching_level": "string"})
    expected = ["bin_id", "overlaps_diffbound_any", "n_L2_diff_support"] \
               + [f"membership_{c}" for c in L2_CLASSES] + GROUP_COLUMNS
    missing = [c for c in expected if c not in anno.columns]
    if missing:
        raise RuntimeError(f"bin_annotation missing columns: {missing}")
    log(f"  annotation shape: {anno.shape}")
    return anno


# Gene-burden covariates are produced by Step 4 (34_intersect_sv_with_heffel_master_v10.py)
# and consumed as hard requirements by Step 6 (39_fit_B_prime_L2_and_specificity_v5.R).
# If they are missing here, the R script will hard-fail; fail fast instead of warning.
REQUIRED_COVARIATE_COLS = [
    "sample_id", "Diagnosis", "Sex_numeric",
    "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
    "log1p_total_del_bases", "log1p_total_dup_bases",
    "log1p_total_gene_DEL",  "log1p_total_gene_DUP",
]


def load_covariate_burden(path: Path) -> pd.DataFrame:
    log(f"Loading covariate burden: {path}")
    df = pd.read_csv(path, sep="\t")
    missing_required = [c for c in REQUIRED_COVARIATE_COLS if c not in df.columns]
    if missing_required:
        raise ValueError(
            "Required covariates missing from Step 4 burden table "
            f"({path.name}): {missing_required}. "
            "Step 6 R fit (B\' logistic regression) will hard-fail without these. "
            "Re-run Step 4 (34_intersect_sv_with_heffel_master_v10.py) to generate "
            "log1p_total_gene_DEL / log1p_total_gene_DUP."
        )
    missing_optional = [c for c in COVARIATE_COLS if c not in df.columns]
    if missing_optional:
        log(f"  [INFO] optional covariate columns missing (will be dropped): {missing_optional}")
    keep = [c for c in COVARIATE_COLS if c in df.columns]
    df = df[keep].drop_duplicates("sample_id").copy()
    log(f"  covariate shape: {df.shape}")
    return df


def stream_event_overlap(
    path: Path,
    anno_small: pd.DataFrame,
    chunksize: int,
) -> pd.DataFrame:
    """event overlap を chunk 読みして bin annotation と join (v3 と同一)."""
    use_cols = ["sample_id", "event_id", "sv_type_norm", "bin_id"]
    dtype = {
        "sample_id": "string",
        "event_id": "string",
        "sv_type_norm": "string",
        "bin_id": "string",
    }
    needed_anno_cols = (
        ["bin_id"]
        + [f"membership_{c}" for c in L2_CLASSES]
        + GROUP_COLUMNS
    )
    anno_small = anno_small[needed_anno_cols].copy()
    anno_small = anno_small.set_index("bin_id")

    pieces = []
    total_rows = 0
    log(f"Reading event overlap in chunks: {path}")
    for i, chunk in enumerate(
        pd.read_csv(
            path, sep="\t", compression="gzip",
            usecols=use_cols, dtype=dtype, chunksize=chunksize,
        ),
        start=1,
    ):
        total_rows += len(chunk)
        chunk = chunk[chunk["sv_type_norm"].isin(SV_TYPES)]
        merged = chunk.join(anno_small, on="bin_id", how="inner")
        if len(merged):
            pieces.append(merged)
        log(f"  chunk {i}: rows={len(chunk):,} after sv filter, "
            f"joined={len(merged):,} (cumulative {total_rows:,})")
    if not pieces:
        raise RuntimeError("No overlap rows after filtering and joining.")
    out = pd.concat(pieces, ignore_index=True)
    log(f"Joined long-format shape: {out.shape}")
    return out


def compute_burden(
    long_df: pd.DataFrame,
    samples: pd.Series,
) -> pd.DataFrame:
    """sample x (L2 class or group) x SV type burden 集計.

    v4 では GROUP_COLUMNS = ['group_primary'] のみ (S2-S5 削除済)。
    """
    all_samples = pd.Index(samples.unique(), name="sample_id")
    n_samples = len(all_samples)
    log(f"Computing burden for {n_samples:,} samples")

    out = pd.DataFrame(index=all_samples)

    for sv in SV_TYPES:
        sub = long_df[long_df["sv_type_norm"] == sv]
        log(f"  SV={sv}: rows={len(sub):,}")

        # ---- 10 L2 classes ----
        for c in L2_CLASSES:
            col_mem = f"membership_{c}"
            sub_c = sub[sub[col_mem] == 1]
            if len(sub_c) == 0:
                out[f"n_boundary_{c}_{sv}"] = 0
                out[f"n_events_{c}_{sv}"] = 0
                out[f"carrier_boundary_{c}_{sv}"] = 0
                continue
            nbin = sub_c.groupby("sample_id")["bin_id"].nunique()
            nevt = sub_c.groupby("sample_id")["event_id"].nunique()
            out[f"n_boundary_{c}_{sv}"] = nbin.reindex(all_samples, fill_value=0).astype(int)
            out[f"n_events_{c}_{sv}"] = nevt.reindex(all_samples, fill_value=0).astype(int)
            out[f"carrier_boundary_{c}_{sv}"] = (out[f"n_boundary_{c}_{sv}"] >= 1).astype(int)

        # ---- group_primary (S1) のみ (v4: S2-S5 削除済) ----
        for gcol in GROUP_COLUMNS:
            labels = long_df[gcol].dropna().unique().tolist()
            for lbl in labels:
                sub_l = sub[sub[gcol] == lbl]
                if len(sub_l) == 0:
                    out[f"n_boundary_{gcol}__{lbl}_{sv}"] = 0
                    out[f"n_events_{gcol}__{lbl}_{sv}"] = 0
                    out[f"carrier_boundary_{gcol}__{lbl}_{sv}"] = 0
                    continue
                nbin = sub_l.groupby("sample_id")["bin_id"].nunique()
                nevt = sub_l.groupby("sample_id")["event_id"].nunique()
                out[f"n_boundary_{gcol}__{lbl}_{sv}"] = nbin.reindex(all_samples, fill_value=0).astype(int)
                out[f"n_events_{gcol}__{lbl}_{sv}"] = nevt.reindex(all_samples, fill_value=0).astype(int)
                out[f"carrier_boundary_{gcol}__{lbl}_{sv}"] = (out[f"n_boundary_{gcol}__{lbl}_{sv}"] >= 1).astype(int)

    out = out.reset_index()
    log(f"Burden table shape: {out.shape}")
    return out


def make_summary(burden: pd.DataFrame) -> pd.DataFrame:
    rows = []
    rows.append(("n_samples", len(burden)))
    for sv in SV_TYPES:
        for c in L2_CLASSES:
            col = f"n_boundary_{c}_{sv}"
            if col in burden.columns:
                rows.append((f"carrier_{c}_{sv}", int(burden[f"carrier_boundary_{c}_{sv}"].sum())))
                rows.append((f"total_bin_overlaps_{c}_{sv}", int(burden[col].sum())))
        for gcol in GROUP_COLUMNS:
            for col in burden.columns:
                if col.startswith(f"carrier_boundary_{gcol}__") and col.endswith(f"_{sv}"):
                    rows.append((col, int(burden[col].sum())))
    return pd.DataFrame(rows, columns=["metric", "count"])


def main() -> None:
    ensure_output_dirs()
    args = parse_args()
    t0 = time.time()
    log("=" * 60)
    log("Start 38_compute_sample_burden_L2_and_specificity_v4.py "
        "(Pattern A primary; S2-S5 sensitivity grouping removed)")
    log(f"  event_overlap: {args.event_overlap}")
    log(f"  bin_annotation: {args.bin_annotation}")
    log(f"  covariate_burden: {args.covariate_burden}")
    log(f"  outdir: {args.outdir}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    # ---- 1. Load bin annotation (v3, group_primary only) ----
    anno = load_bin_annotation(args.bin_annotation)

    # ---- 2. Stream-join event overlap ----
    long_df = stream_event_overlap(
        args.event_overlap, anno, chunksize=args.chunksize
    )

    # ---- 3. Load covariate burden ----
    cov = load_covariate_burden(args.covariate_burden)
    samples = cov["sample_id"].astype("string")

    # ---- 4. Compute burden ----
    burden = compute_burden(long_df, samples)

    # ---- 5. Merge with covariates ----
    in_cov = set(cov["sample_id"].astype(str))
    in_bur = set(burden["sample_id"].astype(str))
    missing_from_cov = in_bur - in_cov
    if missing_from_cov:
        log(f"  [WARN] {len(missing_from_cov)} samples in burden NOT in covariate table; dropping.")
    final = cov.merge(burden, how="left", on="sample_id")
    num_cols = [c for c in final.columns if c.startswith(("n_boundary_", "n_events_", "carrier_boundary_"))]
    final[num_cols] = final[num_cols].fillna(0).astype(int)
    log(f"Final merged table shape: {final.shape}")

    # ---- 6. Write ----
    out_final = args.outdir / "sample_burden_L2_and_specificity_v4.tsv"
    log(f"Writing: {out_final}")
    final.to_csv(out_final, sep="\t", index=False)

    summary = make_summary(final)
    out_sum = args.outdir / "sample_burden_L2_and_specificity_summary_v4.tsv"
    log(f"Writing summary: {out_sum}")
    summary.to_csv(out_sum, sep="\t", index=False)

    # ---- 7. Log preview ----
    log("Preview of carrier counts (primary groups):")
    for sv in SV_TYPES:
        for label in ["Diff_specific_n1", "Diff_shared_n2plus", "Static"]:
            col = f"carrier_boundary_group_primary__{label}_{sv}"
            if col in final.columns:
                by_dx = final.groupby("Diagnosis")[col].sum()
                log(f"  {col}: {by_dx.to_dict()}")

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
