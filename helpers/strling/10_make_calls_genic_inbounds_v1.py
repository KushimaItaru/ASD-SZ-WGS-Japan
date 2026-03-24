#!/usr/bin/env python3
# ============================================================================
# NOTE: This script contains hardcoded default paths specific to the original
# analysis environment (NIG supercomputer). To run in a different environment,
# update the paths below or set the corresponding environment variables
# (e.g. SAMPLE_INFO, PCA_EIGENVEC, CRAM_BASE_DIR1, CRAM_BASE_DIR2).
# ============================================================================
# 10_make_calls_genic_inbounds_v1.py
# - 処理内容:
#   - STRling の joint-bounds（genic + repeatunit 3–8bp）から IN_BOUNDS の locus key を作成
#   - calls_genic/ 配下の *-genotype.txt を、IN_BOUNDS のみ残すようにフィルタ
#   - 出力先 calls_genic_inbounds/ に *-genotype.txt を作成（ヘッダ保持）
#   - 各ファイルの入力/出力行数、IN/OUT 内訳を summary TSV に記録
#   - 実行時間（秒）を記録
#
# 使い方:
#   Run via the top-level wrapper: strling/02_strling_casecontrol_call_and_outliers.sh
#
# パスを明示する場合:
#   python3 strling/10_make_calls_genic_inbounds_v1.py \
#     --bounds /path/to/joint-bounds.genic_1kbpad.len3_8.txt \
#     --calls-dir /path/to/calls_genic \
#     --out-dir /path/to/calls_genic_inbounds

from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def eprint(msg: str) -> None:
    print(msg, file=sys.stderr)


def norm_colname(x: str) -> str:
    x = x.strip()
    if x.startswith("#"):
        x = x[1:]
    return x.strip().lower()


def detect_indices(header_fields: List[str]) -> Tuple[int, int, int, int]:
    """
    chrom/left/right/repeatunit の列位置をヘッダ名から推定。
    見つからなければ (0,1,2,3) を返す。
    """
    norm = [norm_colname(x) for x in header_fields]
    # STRling bounds: #chrom left right repeat ...
    # STRling genotype: #chrom left right repeatunit ...
    chrom_i = None
    left_i = None
    right_i = None
    rep_i = None

    for i, c in enumerate(norm):
        if c in {"chrom", "chr", "contig"} and chrom_i is None:
            chrom_i = i
        if c in {"left", "start"} and left_i is None:
            left_i = i
        if c in {"right", "end", "stop"} and right_i is None:
            right_i = i
        if c in {"repeatunit", "repeat_unit", "repeat"} and rep_i is None:
            rep_i = i

    if chrom_i is None or left_i is None or right_i is None or rep_i is None:
        return (0, 1, 2, 3)
    return (chrom_i, left_i, right_i, rep_i)


def make_key(chrom: str, left: str, right: str, rep: str) -> str:
    # key を 4列で一意化（タブ区切り）
    return f"{chrom}\t{left}\t{right}\t{rep}"


def load_inbounds_keys(bounds_path: Path) -> Tuple[set[str], List[str]]:
    """
    boundsファイル（ヘッダ付き）から IN_BOUNDS の key 集合を作る。
    戻り値: (keys_set, header_fields)
    """
    if not bounds_path.exists():
        raise FileNotFoundError(f"Bounds file not found: {bounds_path}")
    if bounds_path.stat().st_size == 0:
        raise ValueError(f"Bounds file is empty: {bounds_path}")

    with bounds_path.open("r") as f:
        header = f.readline().rstrip("\n")
        if not header:
            raise ValueError(f"Bounds header is empty: {bounds_path}")
        header_fields = header.split("\t")
        ci, li, ri, mi = detect_indices(header_fields)

        keys: set[str] = set()
        n = 0
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) <= max(ci, li, ri, mi):
                continue
            k = make_key(parts[ci], parts[li], parts[ri], parts[mi])
            keys.add(k)
            n += 1

    return keys, header_fields


@dataclass
class FileSummary:
    sample_id: str
    infile: str
    outfile: str
    input_data_lines: int
    output_data_lines: int
    in_bounds_lines: int
    out_of_bounds_lines: int
    status: str
    note: str


def sample_id_from_filename(p: Path) -> str:
    name = p.name
    if name.endswith("-genotype.txt"):
        return name[:-len("-genotype.txt")]
    # fallback
    return name.split("-", 1)[0]


def filter_one_genotype(
    infile: Path,
    outfile: Path,
    keys: set[str],
    overwrite: bool,
    skip_if_exists: bool,
) -> FileSummary:
    sid = sample_id_from_filename(infile)

    if not infile.exists():
        return FileSummary(sid, str(infile), str(outfile), 0, 0, 0, 0, "MISSING_INFILE", "infile not found")

    if infile.stat().st_size == 0:
        return FileSummary(sid, str(infile), str(outfile), 0, 0, 0, 0, "EMPTY_INFILE", "infile is empty")

    if outfile.exists() and skip_if_exists and outfile.stat().st_size > 0:
        # 既に作成済みならスキップ
        return FileSummary(sid, str(infile), str(outfile), 0, 0, 0, 0, "SKIP_EXISTS", "outfile exists")

    if outfile.exists() and (not overwrite) and (not skip_if_exists):
        return FileSummary(sid, str(infile), str(outfile), 0, 0, 0, 0, "SKIP_NO_OVERWRITE", "overwrite disabled")

    # read header
    try:
        with infile.open("r") as fin:
            header = fin.readline().rstrip("\n")
            if not header:
                return FileSummary(sid, str(infile), str(outfile), 0, 0, 0, 0, "BAD_HEADER", "empty header")
            hfields = header.split("\t")
            ci, li, ri, mi = detect_indices(hfields)

            # write
            outfile.parent.mkdir(parents=True, exist_ok=True)
            with outfile.open("w") as fout:
                fout.write(header + "\n")

                in_bounds = 0
                out_bounds = 0
                in_lines = 0
                out_lines = 0

                for line in fin:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) <= max(ci, li, ri, mi):
                        # 列崩れは out_of_bounds 扱い（安全側）
                        out_bounds += 1
                        in_lines += 1
                        continue

                    k = make_key(parts[ci], parts[li], parts[ri], parts[mi])
                    in_lines += 1
                    if k in keys:
                        fout.write(line + "\n")
                        out_lines += 1
                        in_bounds += 1
                    else:
                        out_bounds += 1

        return FileSummary(
            sample_id=sid,
            infile=str(infile),
            outfile=str(outfile),
            input_data_lines=in_lines,
            output_data_lines=out_lines,
            in_bounds_lines=in_bounds,
            out_of_bounds_lines=out_bounds,
            status="OK",
            note="",
        )

    except Exception as e:
        return FileSummary(sid, str(infile), str(outfile), 0, 0, 0, 0, "FAIL", f"{type(e).__name__}: {e}")


def main() -> None:
    ap = argparse.ArgumentParser()
    # デフォルトはあなたのプロジェクト構造に合わせる
    project_root = Path(__file__).resolve().parents[2]
    out_root = project_root / "strling_output_genomewide"

    ap.add_argument(
        "--bounds",
        default=str(out_root / "str-results" / "joint-bounds.genic_1kbpad.len3_8.txt"),
        help="IN_BOUNDS を定義する joint-bounds（genic+len3_8）",
    )
    ap.add_argument(
        "--calls-dir",
        default=str(out_root / "calls_genic"),
        help="入力 genotype があるディレクトリ（calls_genic）",
    )
    ap.add_argument(
        "--out-dir",
        default=str(out_root / "calls_genic_inbounds"),
        help="出力先ディレクトリ（calls_genic_inbounds）",
    )
    ap.add_argument(
        "--pattern",
        default="*-genotype.txt",
        help="入力ファイルパターン（デフォルト: *-genotype.txt）",
    )
    ap.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="既存出力を上書きする（デフォルト: しない）",
    )
    ap.add_argument(
        "--skip-if-exists",
        action="store_true",
        default=True,
        help="出力が既に存在し非空ならスキップ（デフォルト: 有効）",
    )
    ap.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="デバッグ用: 処理する最大ファイル数（0=全て）",
    )
    args = ap.parse_args()

    t0 = time.time()
    print(f"[{ts()}] [INFO] Start {Path(__file__).name}")
    print(f"[{ts()}] [INFO] bounds={args.bounds}")
    print(f"[{ts()}] [INFO] calls_dir={args.calls_dir}")
    print(f"[{ts()}] [INFO] out_dir={args.out_dir}")
    print(f"[{ts()}] [INFO] pattern={args.pattern}")
    print(f"[{ts()}] [INFO] overwrite={args.overwrite} skip_if_exists={args.skip_if_exists}")

    bounds_path = Path(args.bounds)
    calls_dir = Path(args.calls_dir)
    out_dir = Path(args.out_dir)

    if not calls_dir.exists():
        raise SystemExit(f"[ERROR] calls_dir not found: {calls_dir}")

    # load keys
    keys, _ = load_inbounds_keys(bounds_path)
    print(f"[{ts()}] [INFO] Loaded IN_BOUNDS keys: {len(keys):,}")

    files = sorted(calls_dir.glob(args.pattern))
    if not files:
        raise SystemExit(f"[ERROR] No genotype files matched: {calls_dir}/{args.pattern}")

    if args.max_files and args.max_files > 0:
        files = files[: args.max_files]
        print(f"[{ts()}] [INFO] max_files applied: {len(files)} files")

    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "calls_genic_inbounds.filter_summary.tsv"
    log_path = out_dir / f"calls_genic_inbounds.filter_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"

    summaries: List[FileSummary] = []
    n_ok = n_fail = n_skip = 0

    for i, infile in enumerate(files, start=1):
        sid = sample_id_from_filename(infile)
        outfile = out_dir / infile.name

        s = filter_one_genotype(
            infile=infile,
            outfile=outfile,
            keys=keys,
            overwrite=args.overwrite,
            skip_if_exists=args.skip_if_exists,
        )
        summaries.append(s)

        if s.status == "OK":
            n_ok += 1
        elif s.status.startswith("SKIP"):
            n_skip += 1
        else:
            n_fail += 1

        if i % 500 == 0 or i == 1 or i == len(files):
            print(f"[{ts()}] [INFO] Progress {i}/{len(files)} | OK={n_ok} SKIP={n_skip} FAIL={n_fail}", end="\r")

    print()
    # write summary TSV
    with summary_path.open("w") as w:
        w.write("\t".join([
            "SampleID", "infile", "outfile",
            "input_data_lines", "output_data_lines",
            "in_bounds_lines", "out_of_bounds_lines",
            "status", "note"
        ]) + "\n")
        for s in summaries:
            w.write("\t".join([
                s.sample_id, s.infile, s.outfile,
                str(s.input_data_lines), str(s.output_data_lines),
                str(s.in_bounds_lines), str(s.out_of_bounds_lines),
                s.status, s.note.replace("\t", " ")
            ]) + "\n")

    # write log
    elapsed = time.time() - t0
    with log_path.open("w") as w:
        w.write(f"Started: {ts()}\n")
        w.write(f"bounds: {bounds_path}\n")
        w.write(f"calls_dir: {calls_dir}\n")
        w.write(f"out_dir: {out_dir}\n")
        w.write(f"pattern: {args.pattern}\n")
        w.write(f"IN_BOUNDS keys: {len(keys)}\n")
        w.write(f"files_total: {len(files)}\n")
        w.write(f"OK: {n_ok}\n")
        w.write(f"SKIP: {n_skip}\n")
        w.write(f"FAIL: {n_fail}\n")
        w.write(f"Elapsed_sec: {elapsed:.1f}\n")

    print(f"[{ts()}] [DONE] Wrote summary: {summary_path}")
    print(f"[{ts()}] [DONE] Wrote log    : {log_path}")
    print(f"[{ts()}] [DONE] OK={n_ok} SKIP={n_skip} FAIL={n_fail} Elapsed={elapsed:.1f}s")


if __name__ == "__main__":
    main()
