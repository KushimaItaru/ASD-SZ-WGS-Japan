#!/usr/bin/env python3
"""
strling/03_collect_depths_v1.py

処理内容:
- depth/depth_parts/*.tsv を結合して depths_all.tsv を作成
- 出力形式: sample \t depth（ヘッダ付き）
- 実行時間を記録

使い方:
  cd ~/str_12282025
  python3 strling/03_collect_depths_v1.py
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

    project_root = Path(__file__).resolve().parents[1]
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
