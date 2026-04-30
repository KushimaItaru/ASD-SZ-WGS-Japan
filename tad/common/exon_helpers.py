"""
common/exon_helpers.py
======================
共通 helper module for exon-exclusion sensitivity analyses (Module 12 / 13)。

Module 12 (WGS exon-exclusion) と Module 13 (MSSNG exon-exclusion) で重複していた
以下の関数・定数を集約:
  - L2_MEMBERSHIP_COLS: 10 L2 class の membership 列名
  - build_merged_exon_bed(): GENCODE GTF を読み protein-coding exon を chrom 単位で
                            merge 済み区間として返す
  - overlaps_any_exon(): bisect-based 半開区間 overlap 判定
  - load_bin_l2_annotation(): bin_l2_annotation tsv.gz を読み (exon_free_bins,
                              exon_overlap_bins, any_l2_bins) を返す

Module-specific な処理 (count_disrupted_bins / fit_logit / fit_gee / fit_firth /
main 等) は各 module 内に残す。

Usage in each module:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from common.exon_helpers import (
        L2_MEMBERSHIP_COLS,
        build_merged_exon_bed,
        overlaps_any_exon,
        load_bin_l2_annotation,
    )

History:
  2026-04-29: Initial extraction from Module 12 v3 / Module 13 v3 (Q2 round).
              本 helper への切り出しにより、Module 12 v4 / Module 13 v4 で
              GENCODE GTF parsing + bin × exon overlap 計算のコード重複 (~150 行 × 2)
              を削減。helper 化前の数値結果と bit-identical (関数シグネチャ + 内部
              ロジック完全保持)。
"""

from __future__ import annotations

import bisect
import gzip
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple


# =========================================================
# Constants
# =========================================================
# 10 L2 cell-type-resolved boundary classes (HPC_Astro 除外、Methods §9 規定)
L2_MEMBERSHIP_COLS = [
    "membership_HPC_Exc-CA",
    "membership_HPC_Exc-DG",
    "membership_HPC_Exc-ENT",
    "membership_HPC_Inh-CGE",
    "membership_HPC_Inh-MGE",
    "membership_PFC_Astro",
    "membership_PFC_Exc-DL",
    "membership_PFC_Exc-UL",
    "membership_PFC_Inh-CGE",
    "membership_PFC_Inh-MGE",
]


# =========================================================
# Logging utility (shared between modules)
# Format matches Module 12 v3 / Module 13 v3 original style:
#   [YYYY-MM-DD HH:MM:SS] [   elapsed_s] msg
# Each module imports this via `from common.exon_helpers import log` so all
# log lines (helper-internal + module-specific) share a single timestamp/elapsed
# format. Each script invocation gets its own _T0 (initialized at import time).
# =========================================================
_T0 = time.time()


def log(msg: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    elapsed = time.time() - _T0
    print(f"[{stamp}] [{elapsed:7.1f}s] {msg}", flush=True)


# =========================================================
# Step 1: Build merged protein-coding exon intervals from GENCODE GTF
# =========================================================
def build_merged_exon_bed(gtf_path: Path) -> Dict[str, List[Tuple[int, int]]]:
    """GENCODE GTF を読み、protein_coding exon を chrom 単位で merge 済み区間として返す。

    Returns:
        Dict mapping chromosome name -> list of (start_0based, end_exclusive) tuples,
        sorted by start. Adjacent/overlapping exons are merged.
    """
    log(f"Reading GENCODE GTF: {gtf_path}")
    exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    n_exon = 0
    with gzip.open(str(gtf_path), "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            if fields[2] != "exon":
                continue
            if 'gene_type "protein_coding"' not in fields[8]:
                continue
            chrom = fields[0]
            s = int(fields[3]) - 1  # 1-based inclusive -> 0-based inclusive
            e = int(fields[4])
            exons[chrom].append((s, e))
            n_exon += 1
    log(f"  Raw protein-coding exon intervals: {n_exon}")

    merged: Dict[str, List[Tuple[int, int]]] = {}
    total = 0
    for chrom in sorted(exons.keys()):
        intervals = sorted(exons[chrom])
        m: List[Tuple[int, int]] = [intervals[0]]
        for s, e in intervals[1:]:
            if s <= m[-1][1]:
                m[-1] = (m[-1][0], max(m[-1][1], e))
            else:
                m.append((s, e))
        merged[chrom] = m
        total += len(m)
    log(f"  Merged exon intervals: {total}")
    return merged


# =========================================================
# Step 1b: bisect-based exon overlap test
# =========================================================
def overlaps_any_exon(
    chrom: str,
    start: int,
    end: int,
    merged_exons: Dict[str, List[Tuple[int, int]]],
    exon_starts: Dict[str, List[int]],
) -> bool:
    """半開区間 [start, end) が merged exon と重複するかを bisect で判定。

    `exon_starts` is a per-chrom list of exon start coordinates (precomputed
    from `merged_exons`); this allows O(log n) bisect to find the candidate
    region. Returns True at first overlap found, False if no overlap.
    """
    if chrom not in exon_starts:
        return False
    starts = exon_starts[chrom]
    intervals = merged_exons[chrom]
    idx = bisect.bisect_right(starts, start) - 1
    for i in range(max(0, idx), len(intervals)):
        es, ee = intervals[i]
        if es >= end:
            break
        if ee > start:
            return True
    return False


# =========================================================
# Step 2: Load bin_l2_annotation -> exon-free / exon-overlap / any-L2 sets
# =========================================================
def load_bin_l2_annotation(
    path: Path,
    merged_exons: Dict[str, List[Tuple[int, int]]],
    exon_starts: Dict[str, List[int]],
) -> Tuple[Set[str], Set[str], Set[str]]:
    """bin_l2_annotation tsv.gz を読み、bin を 3 つの set に分類。

    Args:
        path: Path to bin_l2_annotation_v3.tsv.gz (or v2 backward compat).
              Must contain columns: bin_id, chrom, start0, end, and 10
              membership_<L2_class> columns.
        merged_exons: Output of build_merged_exon_bed().
        exon_starts: Per-chrom list of exon start positions (for bisect).

    Returns:
        Tuple of three sets:
          - exon_free_bins: bins NOT overlapping any protein-coding exon
          - exon_overlap_bins: bins overlapping at least one protein-coding exon
          - any_l2_bins: bins with at least one membership_<L2_class> > 0
                         (10 L2 classes pool; HPC_Astro excluded)
    """
    log(f"Reading bin L2 annotation: {path}")
    exon_free_bins: Set[str] = set()
    exon_overlap_bins: Set[str] = set()
    any_l2_bins: Set[str] = set()

    with gzip.open(str(path), "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        req_cols = ["bin_id", "chrom", "start0", "end"] + L2_MEMBERSHIP_COLS
        missing = [c for c in req_cols if c not in col]
        if missing:
            raise RuntimeError(
                f"bin_l2_annotation missing columns: {missing}"
            )
        idx_bin = col["bin_id"]
        idx_chr = col["chrom"]
        idx_s = col["start0"]
        idx_e = col["end"]
        idx_mem = [col[c] for c in L2_MEMBERSHIP_COLS]

        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < max(idx_mem) + 1:
                continue
            bin_id = fields[idx_bin]
            is_l2 = False
            for mi in idx_mem:
                try:
                    v = float(fields[mi]) if fields[mi] else 0.0
                except ValueError:
                    v = 0.0
                if v > 0:
                    is_l2 = True
                    break
            if is_l2:
                any_l2_bins.add(bin_id)

            chrom = fields[idx_chr]
            try:
                bs = int(fields[idx_s])
                be = int(fields[idx_e])
            except ValueError:
                continue
            if overlaps_any_exon(chrom, bs, be, merged_exons, exon_starts):
                exon_overlap_bins.add(bin_id)
            else:
                exon_free_bins.add(bin_id)

    total = len(exon_free_bins) + len(exon_overlap_bins)
    log(
        f"  Total annotated bins: {total}  "
        f"exon-overlap={len(exon_overlap_bins)} "
        f"({100*len(exon_overlap_bins)/max(total,1):.1f}%)  "
        f"exon-free={len(exon_free_bins)} "
        f"({100*len(exon_free_bins)/max(total,1):.1f}%)"
    )
    log(f"  Any-L2 pool (10 classes, hpc_astro excluded): {len(any_l2_bins)}")
    return exon_free_bins, exon_overlap_bins, any_l2_bins
