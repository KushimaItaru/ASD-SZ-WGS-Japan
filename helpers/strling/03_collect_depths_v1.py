#!/usr/bin/env python3
"""
strling/03_collect_depths_v1.py

Description:
- Merge depth/depth_parts/*.tsv into depths_all.tsv
- Output format: sample \t depth (with header)
- Record execution time

Usage:
  Run via the top-level wrapper: ehdn/02_ehdn_depth.sh
"""

from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path
import pandas as pd


def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def main() -> None:
    start = datetime.now()
    print(f"[{ts()}] [INFO] Start 03_collect_depths_v1.py")

    project_root = Path(__file__).resolve().parents[2]
    parts_dir = Path(os.environ.get("DEPTH_PARTS_DIR", str(project_root / "depth" / "depth_parts")))
    out_tsv = Path(os.environ.get("DEPTHS_ALL_TSV", str(project_root / "depth" / "depths_all.tsv")))

    if not parts_dir.exists():
        raise SystemExit(f"[ERROR] depth parts dir not found: {parts_dir}")

    files = sorted(parts_dir.glob("*.tsv"))
    if not files:
        raise SystemExit(f"[ERROR] No depth part files: {parts_dir}/*.tsv")

    rows = []
    for f in files:
        try:
            line = f.read_text().strip().split("\t")
            if len(line) >= 2:
                sid = str(line[0])
                depth = float(line[1])
                rows.append((sid, depth))
        except Exception:
            continue

    df = pd.DataFrame(rows, columns=["sample", "depth"])
    df = df.drop_duplicates(subset=["sample"])
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[{ts()}] [INFO] Wrote {out_tsv} (n={len(df)})")

    end = datetime.now()
    print(f"[{ts()}] [DONE] Elapsed={(end-start).total_seconds():.1f}s")


if __name__ == "__main__":
    main()
