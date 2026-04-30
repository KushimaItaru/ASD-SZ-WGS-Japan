#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
38_compute_sample_burden_L2_and_specificity_v5.py
=====================================================
処理内容:
  - [v5 変更点] Option II cascade: v4 から robustness を強化
    * [FIX] v4 の dead code を除去: `missing_from_cov = in_bur - in_cov` は
      burden が cov の sample 一覧から作られるため常に空だった。v5 では
      「event overlap には存在するが cov に居ない sample」を stream 段で
      明示的に pre-check し WARN ログする (silent sample loss 防止)。
    * [ADD] 「cov には存在するが event overlap に 1 件も出現しない sample」
      の件数を INFO ログ + summary に記録 (ゼロ exposure だが解析に含まれる点を透明化)。
    * [ADD] 出力 TSV を sample_id で sort して deterministic 化。
    * [ADD] summary に QC counter (n_samples_cov, n_samples_burden_only,
      n_samples_cov_only, n_total_events_kept) を追加。
    * それ以外は v4 と完全一致 (統計・カラム生成ロジック不変)。
  - v4 から継承する機能:
    * PATTERN 環境変数 (必須: 'patB' または 'patC') で切替
    * 入力: Step 04 v11 (Pattern 別 overlap + burden) の出力を参照
        event overlap   : output_v11_{PATTERN}/sample_boundary_event_overlap_v11_{PATTERN}.tsv.gz
        covariate burden: output_v11_{PATTERN}/sample_mechanistic_burden_v11_{PATTERN}.tsv
    * 出力: output_v5_{PATTERN}/
        sample_burden_L2_and_specificity_v5_{PATTERN}.tsv
        sample_burden_L2_and_specificity_summary_v5_{PATTERN}.tsv
    * 04_v11 event overlap と 02_v2 bin-level annotation を join
    * サンプル x (10 L2 classes + 5 group definitions) x (DEL / DUP) の
      bin-count / event-count / carrier burden を計算
    * 04_v11 の既存 burden から共変量 (Diagnosis, Sex, Sex_numeric, PC1-10,
      total_bases/count, total_gene, log1p_*) を継承
    * 10 L2 classes (HPC_Astro 除外):
        HPC_Exc-CA, HPC_Exc-DG, HPC_Exc-ENT, HPC_Inh-CGE, HPC_Inh-MGE,
        PFC_Astro, PFC_Exc-DL, PFC_Exc-UL, PFC_Inh-CGE, PFC_Inh-MGE
    * 5 group definitions (group_primary, group_s2, group_s3, group_s4, group_s5)
    * Exposure 3 種類: n_boundary_, n_events_, carrier_boundary_
  - 実行時間記録
"""
from __future__ import annotations
import argparse
import os
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

# ============================================================
# v5: PATTERN 環境変数で Pattern B / Pattern C を切替 (v4 と同様)
# ============================================================
_PATTERN = os.environ.get("PATTERN", "").strip()
if _PATTERN not in ("patB", "patC"):
    raise RuntimeError(
        f"v5 requires PATTERN env var = 'patB' or 'patC'. "
        f"Got: {_PATTERN!r}. Example: PATTERN=patB python3 38_...v5.py"
    )

# Step 04 v11 の Pattern 別出力ディレクトリ + ファイル
OUT_04_V11 = PIPELINE_ROOT / "04_wgs_sv_boundary_overlap" / f"output_v11_{_PATTERN}"
F_04_V11_EVENT_OVERLAP = OUT_04_V11 / f"sample_boundary_event_overlap_v11_{_PATTERN}.tsv.gz"
F_04_V11_SAMPLE_BURDEN_COV = OUT_04_V11 / f"sample_mechanistic_burden_v11_{_PATTERN}.tsv"

# v5 output dir (Pattern 別)
OUT_05_V5 = PIPELINE_ROOT / "05_wgs_sample_burden" / f"output_v5_{_PATTERN}"

# 出力 suffix
V5_SUFFIX = f"v5_{_PATTERN}"


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


# =========================================================
# 共通: L2 classes と group 定義 (v4 と同一)
# =========================================================
L2_CLASSES = [
    "HPC_Exc-CA", "HPC_Exc-DG", "HPC_Exc-ENT",
    "HPC_Inh-CGE", "HPC_Inh-MGE",
    "PFC_Astro", "PFC_Exc-DL", "PFC_Exc-UL",
    "PFC_Inh-CGE", "PFC_Inh-MGE",
]
GROUP_COLUMNS = ["group_primary", "group_s2", "group_s3", "group_s4", "group_s5"]
SV_TYPES = ["DEL", "DUP"]

# 共変量として継承する列名（04_v11 burden 由来、v4 と同一）
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

# Step 6 R fit (B') が必要とする hard-required covariates (v4 と同一)
REQUIRED_COVARIATE_COLS = [
    "sample_id", "Diagnosis", "Sex_numeric",
    "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
    "log1p_total_del_bases", "log1p_total_dup_bases",
    "log1p_total_gene_DEL",  "log1p_total_gene_DUP",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=f"Compute sample x (L2 classes + specificity groups) x SV type burden "
                    f"(v5: input=04_v11 {_PATTERN}; PATTERN env var required; "
                    f"silent-sample-loss fix)."
    )
    p.add_argument("--event-overlap", default=F_04_V11_EVENT_OVERLAP, type=Path,
                   help=f"Path to sample_boundary_event_overlap_v11_{_PATTERN}.tsv.gz (04_v11).")
    p.add_argument("--bin-annotation", default=F_02_BIN_L2_ANNOTATION, type=Path,
                   help="Path to bin_l2_annotation_v2.tsv.gz (37).")
    p.add_argument("--covariate-burden", default=F_04_V11_SAMPLE_BURDEN_COV, type=Path,
                   help=f"Path to sample_mechanistic_burden_v11_{_PATTERN}.tsv (04_v11).")
    p.add_argument("--outdir", default=OUT_05_V5, type=Path,
                   help=f"Output directory (output_v5_{_PATTERN}/, will be created).")
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


def load_covariate_burden(path: Path) -> pd.DataFrame:
    log(f"Loading covariate burden: {path}")
    df = pd.read_csv(path, sep="\t")
    missing_required = [c for c in REQUIRED_COVARIATE_COLS if c not in df.columns]
    if missing_required:
        raise ValueError(
            "Required covariates missing from Step 4 v11 burden table "
            f"({path.name}): {missing_required}. "
            "Step 6 R fit (B' logistic regression) will hard-fail without these. "
            "Re-run Step 4 v11 (34_intersect_sv_with_heffel_master_v11.py)."
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
    cov_samples: set,
    chunksize: int,
) -> tuple:
    """event overlap を chunk 読みして bin annotation と join.

    v5 変更:
      - cov に居ない sample を chunk 単位で検出し WARN ログ
      - 戻り値に event overlap 由来の QC counter を追加
    """
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
    samples_not_in_cov: set = set()  # event overlap にあるが cov に居ない sample
    samples_seen_in_events: set = set()  # event overlap に 1 件以上出現した sample
    total_events_kept = 0

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

        # v5 FIX: chunk 内の sample を cov 集合と突合
        chunk_samples = set(chunk["sample_id"].dropna().astype(str).unique())
        samples_seen_in_events.update(chunk_samples)
        bad = chunk_samples - cov_samples
        if bad:
            samples_not_in_cov.update(bad)
            # chunk 内で cov に居ない sample の event は保持せず除外
            chunk = chunk[chunk["sample_id"].isin(cov_samples)]

        merged = chunk.join(anno_small, on="bin_id", how="inner")
        if len(merged):
            pieces.append(merged)
            total_events_kept += len(merged)
        log(f"  chunk {i}: rows={len(chunk):,} after sv filter + cov filter, "
            f"joined={len(merged):,} (cumulative {total_rows:,})")

    if not pieces:
        raise RuntimeError("No overlap rows after filtering and joining.")

    out = pd.concat(pieces, ignore_index=True)
    log(f"Joined long-format shape: {out.shape}")

    if samples_not_in_cov:
        log(f"  [WARN] {len(samples_not_in_cov)} samples in event_overlap NOT in "
            f"covariate table — their events were DROPPED "
            f"(first 5: {sorted(samples_not_in_cov)[:5]}).")

    qc = {
        "n_raw_rows": int(total_rows),
        "n_events_kept": int(total_events_kept),
        "n_samples_seen_in_events": int(len(samples_seen_in_events)),
        "n_samples_event_only_dropped": int(len(samples_not_in_cov)),
        "samples_event_only_dropped": sorted(samples_not_in_cov),
        "samples_seen_in_events": samples_seen_in_events,
    }
    return out, qc


def compute_burden(
    long_df: pd.DataFrame,
    samples: pd.Series,
) -> pd.DataFrame:
    """sample x (L2 class or group) x SV type burden 集計 (v4 と同一ロジック)."""
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

        # ---- 5 group definitions ----
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


def make_summary(
    burden: pd.DataFrame,
    qc_counters: dict,
) -> pd.DataFrame:
    rows = []
    # ---- v5 QC counters を先頭に配置 ----
    rows.append(("pattern", _PATTERN))
    rows.append(("n_samples", len(burden)))
    rows.append(("n_samples_cov_total", qc_counters.get("n_samples_cov_total", "")))
    rows.append(("n_samples_cov_zero_events", qc_counters.get("n_samples_cov_zero_events", "")))
    rows.append(("n_samples_event_only_dropped", qc_counters.get("n_samples_event_only_dropped", "")))
    rows.append(("n_events_kept", qc_counters.get("n_events_kept", "")))
    rows.append(("n_raw_rows_event_overlap", qc_counters.get("n_raw_rows", "")))

    # ---- v4 互換の carrier / total 集計 ----
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
    log(f"Start 38_compute_sample_burden_L2_and_specificity_v5.py [PATTERN={_PATTERN}]")
    log(f"  event_overlap: {args.event_overlap}")
    log(f"  bin_annotation: {args.bin_annotation}")
    log(f"  covariate_burden: {args.covariate_burden}")
    log(f"  outdir: {args.outdir}")

    args.outdir.mkdir(parents=True, exist_ok=True)

    # ---- 1. Load bin annotation ----
    anno = load_bin_annotation(args.bin_annotation)

    # ---- 2. Load covariate burden (先に load して sample 集合を取得) ----
    cov = load_covariate_burden(args.covariate_burden)
    cov_samples_set = set(cov["sample_id"].astype(str))
    log(f"  cov samples: n={len(cov_samples_set)}")

    # ---- 3. Stream-join event overlap (v5: cov sample 突合を内部で実施) ----
    long_df, stream_qc = stream_event_overlap(
        args.event_overlap, anno,
        cov_samples=cov_samples_set,
        chunksize=args.chunksize,
    )

    # ---- 4. Compute burden (cov の sample を基準に) ----
    samples = cov["sample_id"].astype("string")
    burden = compute_burden(long_df, samples)

    # ---- 5. Merge with covariates ----
    # v5 FIX: v4 の dead code (missing_from_cov) を削除。代わりに v5 版の
    # 明示的カウンタを記録する。burden は cov sample ベースで生成されるため
    # cov と burden の sample 集合は等しいことを assert する。
    in_bur = set(burden["sample_id"].astype(str))
    in_cov = cov_samples_set
    assert in_bur == in_cov, (
        f"Internal invariant violated: burden samples ({len(in_bur)}) != "
        f"covariate samples ({len(in_cov)}). "
        f"burden_only={len(in_bur - in_cov)}, cov_only={len(in_cov - in_bur)}"
    )

    # ---- v5 ADD: cov には居るが event には 1 件も出現しなかった sample をカウント ----
    samples_seen = stream_qc.get("samples_seen_in_events", set())
    samples_cov_zero_events = cov_samples_set - samples_seen
    if samples_cov_zero_events:
        log(f"  [INFO] {len(samples_cov_zero_events)} samples in cov had zero events "
            f"in event_overlap (burden=0 across all classes). "
            f"This is expected for some controls; included in analysis with zero exposure.")

    final = cov.merge(burden, how="left", on="sample_id")
    num_cols = [c for c in final.columns if c.startswith(("n_boundary_", "n_events_", "carrier_boundary_"))]
    final[num_cols] = final[num_cols].fillna(0).astype(int)

    # ---- v5 ADD: deterministic output のため sample_id で sort ----
    final = final.sort_values("sample_id", kind="mergesort").reset_index(drop=True)
    log(f"Final merged table shape: {final.shape}")

    # ---- 6. Write ----
    out_final = args.outdir / f"sample_burden_L2_and_specificity_{V5_SUFFIX}.tsv"
    log(f"Writing: {out_final}")
    final.to_csv(out_final, sep="\t", index=False)

    qc_counters = {
        "n_samples_cov_total": int(len(cov_samples_set)),
        "n_samples_cov_zero_events": int(len(samples_cov_zero_events)),
        "n_samples_event_only_dropped": int(stream_qc.get("n_samples_event_only_dropped", 0)),
        "n_events_kept": int(stream_qc.get("n_events_kept", 0)),
        "n_raw_rows": int(stream_qc.get("n_raw_rows", 0)),
    }
    summary = make_summary(final, qc_counters)
    out_sum = args.outdir / f"sample_burden_L2_and_specificity_summary_{V5_SUFFIX}.tsv"
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

    # ---- 8. v5 summary log ----
    log("-" * 60)
    log("[v5 QC summary]")
    for k, v in qc_counters.items():
        log(f"  {k}: {v}")
    log("-" * 60)

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
