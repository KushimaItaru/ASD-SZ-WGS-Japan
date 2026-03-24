#!/usr/bin/env python3
"""
ehdn/04_merge_ehdn_novel_norm_v2.py

処理内容:
- EHdn profile の locus.tsv を全サンプルから読み込み
- gene_regions_1kb_pad.bed で genic に限定
- depths_all.tsv による深度正規化（TARGET_DEPTH=40x）
- 出力:
  merged_results_novel/novel_loci_genic_norm.tsv
  skipped_samples_depth_missing.log
- 実行時間を記録

入力（環境変数 or 既定）:
- EHDN_OUT_DIR（ehdn_output）
- SAMPLE_LIST_DIR/ehdn_all_samples.tsv
- DEPTHS_ALL_TSV（depth/depths_all.tsv）
- GENE_REGIONS_BED（resources/gene_regions_1kb_pad.bed）
"""

from __future__ import annotations

import os
import sys
import subprocess
from datetime import datetime
from pathlib import Path
import pandas as pd


def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def read_depths(depth_tsv: Path) -> dict[str, float]:
    df = pd.read_csv(depth_tsv, sep="\t", dtype={"sample": str})
    if "sample" not in df.columns or "depth" not in df.columns:
        # fallback: 1列目sample 2列目depth
        df = pd.read_csv(depth_tsv, sep="\t", header=0)
        df.columns = [c.lower() for c in df.columns]
    df["sample"] = df["sample"].astype(str)
    df["depth"] = pd.to_numeric(df["depth"], errors="coerce")
    df = df.dropna(subset=["sample", "depth"])
    return dict(zip(df["sample"], df["depth"].astype(float)))


def main() -> None:
    start = datetime.now()
    print(f"[{ts()}] [INFO] Start 04_merge_ehdn_novel_norm_v2.py")

    project_root = Path(__file__).resolve().parents[1]
    ehdn_dir = Path(os.environ.get("EHDN_OUT_DIR", str(project_root / "ehdn_output")))
    sample_list = Path(os.environ.get("SAMPLE_LIST_DIR", str(project_root / "sample_lists"))) / "ehdn_all_samples.tsv"
    depth_tsv = Path(os.environ.get("DEPTHS_ALL_TSV", str(project_root / "depth" / "depths_all.tsv")))
    gene_bed = Path(os.environ.get("GENE_REGIONS_BED", str(project_root / "resources" / "gene_regions_1kb_pad.bed")))
    out_dir = Path(os.environ.get("MERGED_NOVEL_DIR", str(project_root / "merged_results_novel")))
    out_dir.mkdir(exist_ok=True)
    target_depth = float(os.environ.get("TARGET_DEPTH", "40.0"))

    if not sample_list.exists():
        print(f"[{ts()}] [ERROR] sample list not found: {sample_list}", file=sys.stderr)
        sys.exit(1)
    if not depth_tsv.exists():
        print(f"[{ts()}] [ERROR] depth file not found: {depth_tsv}", file=sys.stderr)
        sys.exit(1)
    if not gene_bed.exists():
        print(f"[{ts()}] [ERROR] gene bed not found: {gene_bed}", file=sys.stderr)
        sys.exit(1)
    if not ehdn_dir.exists():
        print(f"[{ts()}] [ERROR] ehdn output dir not found: {ehdn_dir}", file=sys.stderr)
        sys.exit(1)

    depth_map = read_depths(depth_tsv)
    df_samples = pd.read_csv(sample_list, sep="\t", dtype=str)
    if not {"SampleID", "Group"}.issubset(df_samples.columns):
        print(f"[{ts()}] [ERROR] Sample list must contain SampleID and Group columns.", file=sys.stderr)
        sys.exit(1)

    # Collect all loci coords for genic filter
    print(f"[{ts()}] [INFO] Collecting locus coordinates from EHdn outputs...")
    temp_bed = out_dir / "tmp_all_loci.bed"
    # *.locus.tsv: headerあり。1行目スキップして contig start end を抽出
    cmd = f"find {ehdn_dir} -name '*.locus.tsv' -print0 | xargs -0 awk 'FNR>1{{print $1,$2,$3}}' OFS='\\t' > {temp_bed}"
    subprocess.run(cmd, shell=True, check=True)

    filtered_bed = out_dir / "genic_loci.bed"
    # 重要: bedtools stdin は -a - を使う（v1の改善点）
    cmd2 = f"sort -k1,1 -k2,2n {temp_bed} | uniq | bedtools intersect -a - -b {gene_bed} -u > {filtered_bed}"
    subprocess.run(cmd2, shell=True, check=True)

    valid = set()
    with open(filtered_bed, "r") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                valid.add(f"{p[0]}:{p[1]}-{p[2]}")
    print(f"[{ts()}] [INFO] Genic loci (unique): {len(valid):,}")

    results = []
    skipped = []
    processed = 0

    for r in df_samples.itertuples(index=False):
        sid = str(getattr(r, "SampleID"))
        grp = str(getattr(r, "Group"))

        f_path = ehdn_dir / sid / f"{sid}.locus.tsv"
        if not f_path.exists():
            continue

        d = depth_map.get(sid, None)
        if d is None or d <= 0:
            skipped.append(sid)
            continue

        norm_factor = target_depth / float(d)

        try:
            df = pd.read_csv(f_path, sep="\t")
            # EHdnの基本列（想定: contig,start,end,motif,num_anc_irrs）
            if not {"contig", "start", "end", "motif", "num_anc_irrs"}.issubset(df.columns):
                continue

            df["key"] = df["contig"].astype(str) + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
            df = df[df["key"].isin(valid)].copy()
            if df.empty:
                continue

            df["SampleID"] = sid
            df["Group"] = grp
            df["num_anc_irrs_norm"] = df["num_anc_irrs"].astype(float) * norm_factor

            df.drop(columns=["key"], inplace=True)
            results.append(df)
        except Exception:
            continue

        processed += 1
        if processed % 500 == 0:
            print(f"[{ts()}] [INFO] Processed {processed} samples...", end="\r")

    if skipped:
        skip_f = out_dir / "skipped_samples_depth_missing.log"
        skip_f.write_text("\n".join(sorted(set(skipped))) + "\n")
        print(f"\n[{ts()}] [WARN] Skipped samples (missing depth): {len(set(skipped))} -> {skip_f}")

    if not results:
        print(f"[{ts()}] [ERROR] No merged results generated.", file=sys.stderr)
        sys.exit(2)

    df_out = pd.concat(results, ignore_index=True)
    out_tsv = out_dir / "novel_loci_genic_norm.tsv"
    df_out.to_csv(out_tsv, sep="\t", index=False, float_format="%.6f")
    print(f"[{ts()}] [INFO] Wrote {out_tsv} (rows={len(df_out):,})")

    # cleanup temp bed
    try:
        temp_bed.unlink(missing_ok=True)
    except Exception:
        pass

    end = datetime.now()
    print(f"[{ts()}] [DONE] Elapsed={(end-start).total_seconds():.1f}s")


if __name__ == "__main__":
    main()
