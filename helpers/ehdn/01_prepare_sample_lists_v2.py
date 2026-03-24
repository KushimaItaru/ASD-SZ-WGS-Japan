#!/usr/bin/env python3
# ============================================================================
# NOTE: This script contains hardcoded default paths specific to the original
# analysis environment (NIG supercomputer). To run in a different environment,
# update the paths below or set the corresponding environment variables
# (e.g. SAMPLE_INFO, PCA_EIGENVEC, CRAM_BASE_DIR1, CRAM_BASE_DIR2).
# ============================================================================
"""
ehdn/01_prepare_sample_lists_v2.py（v2）

処理内容:
- SampleInfo（TSV）から列名ベースで SampleID / Diagnosis / Father / Mother を自動検出
- CRAMパスを2つのベースディレクトリから探索し、存在するサンプルのみ採用
- 用途別にサンプルリストを分離して出力：
  1) EHdn実行・深度計算・マージ対象（ケース + Healthy + family_member）
     -> sample_lists/ehdn_all_samples.tsv
  2) case-control burden解析対象（ケース + Diagnosis=="Healthy" のみ）
     -> sample_lists/casecontrol_samples.tsv
  3) de novo/trio解析用（proband と父母のリンク、および必要個体のリスト）
     -> sample_lists/trio_links.tsv, sample_lists/trio_samples.tsv
- 実行時間を記録

注意:
- case-controlのcontrolsは Diagnosis=="Healthy" のみ（family_memberは除外）
- EHdnプロファイル作成では family_member も含める（de novo用に必要）

使い方:
  cd ~/str_12282025
  python3 ehdn/01_prepare_sample_lists_v2.py
"""

from __future__ import annotations

import os
import sys
from datetime import datetime
from pathlib import Path
import pandas as pd


def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def find_col(df: pd.DataFrame, candidates: list[str], required: bool = True) -> str | None:
    cols = list(df.columns)
    lower = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand in cols:
            return cand
        if cand.lower() in lower:
            return lower[cand.lower()]
    if required:
        raise KeyError(f"Required column not found. candidates={candidates}, columns={cols}")
    return None


def parse_parent_columns(df: pd.DataFrame) -> tuple[pd.Series, pd.Series] | tuple[None, None]:
    """
    Father/Mother の列名揺れに対応：
    - Father, Mother
    - FatherID, MotherID
    - Father:Mother（結合列）
    """
    cols = set(df.columns)

    if "Father" in cols and "Mother" in cols:
        return df["Father"].astype(str), df["Mother"].astype(str)

    if "FatherID" in cols and "MotherID" in cols:
        return df["FatherID"].astype(str), df["MotherID"].astype(str)

    if "Father:Mother" in cols:
        def split_parents(val: str) -> tuple[str, str]:
            s = str(val)
            if ":" in s:
                a, b = s.split(":", 1)
                return a, b
            return "NA", "NA"

        tmp = df["Father:Mother"].astype(str).apply(split_parents)
        father = tmp.apply(lambda x: x[0])
        mother = tmp.apply(lambda x: x[1])
        return father, mother

    return None, None


def find_cram(sample_id: str, base1: str, base2: str) -> str | None:
    p1 = Path(base1) / sample_id / f"{sample_id}.cram"
    p2 = Path(base2) / sample_id / f"{sample_id}.cram"
    if p1.exists():
        return str(p1)
    if p2.exists():
        return str(p2)
    return None


def main() -> None:
    start = datetime.now()
    print(f"[{ts()}] [INFO] Start 01_prepare_sample_lists_v2.py")

    project_root = Path(__file__).resolve().parents[1]

    sample_info = os.environ.get(
        "SAMPLE_INFO",
        "/lustre12/home/kushima-pg/sampleInfo/GRIFIN_srWGS_SampleInfo_11242025.txt"  # CONFIGURE
    )
    base1 = os.environ.get("CRAM_BASE_DIR1", "/lustre12/home/grifinpd-pg/analysis/parabricks")  # CONFIGURE
    base2 = os.environ.get("CRAM_BASE_DIR2", "/lustre12/home/ncbn-share-pg/control_genome/pb3.1.0/results")  # CONFIGURE
    out_dir = project_root / "sample_lists"
    out_dir.mkdir(exist_ok=True)

    if not Path(sample_info).exists():
        print(f"[{ts()}] [ERROR] SampleInfo not found: {sample_info}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(sample_info, sep="\t", dtype=str)
    col_sid = find_col(df, ["SampleID", "sample_id", "sample"])
    col_dx = find_col(df, ["Diagnosis", "diagnosis", "DX", "dx"], required=True)

    df[col_sid] = df[col_sid].astype(str)
    df[col_dx] = df[col_dx].astype(str)

    # Parent columns（なければ trio 出力はスキップ）
    father_col, mother_col = parse_parent_columns(df)

    # Groupの基本分類（Diagnosisの値をそのまま保持しつつ、用途別に使い分け）
    # - casecontrol controls: dx=="Healthy" のみ
    # - family_member は casecontrol から除外、ただし EHdn 実行対象としては含める
    df["Diagnosis_norm"] = df[col_dx].astype(str).str.strip()

    # CRAM探索（存在するもののみ）
    print(f"[{ts()}] [INFO] Finding CRAM paths for {len(df)} rows in SampleInfo...")
    df["CRAM_Path"] = df[col_sid].apply(lambda s: find_cram(str(s), base1, base2))
    df["CRAM_Exists"] = df["CRAM_Path"].notna()

    n_missing = int((~df["CRAM_Exists"]).sum())
    if n_missing > 0:
        print(f"[{ts()}] [WARN] Missing CRAM for {n_missing} samples. They will be excluded from all downstream lists.")

    df = df[df["CRAM_Exists"]].copy()

    # ---------- EHdn実行対象（プロファイル生成用） ----------
    # ASD / SZ / Healthy / family_member を含める
    ehdn_keep = df["Diagnosis_norm"].isin(["ASD", "SZ", "Healthy", "family_member"])
    df_ehdn = df[ehdn_keep].copy()
    df_ehdn_out = df_ehdn[[col_sid, "Diagnosis_norm", "CRAM_Path"]].rename(columns={col_sid: "SampleID"})
    df_ehdn_out = df_ehdn_out.rename(columns={"Diagnosis_norm": "Group"})  # Group列にASD/SZ/Healthy/family_member
    df_ehdn_out = df_ehdn_out.drop_duplicates(subset=["SampleID"])

    ehdn_file = out_dir / "ehdn_all_samples.tsv"
    df_ehdn_out.to_csv(ehdn_file, sep="\t", index=False)
    print(f"[{ts()}] [INFO] Wrote {ehdn_file} (n={len(df_ehdn_out)})")

    # ---------- case-control対象（burden解析用） ----------
    # cases: ASD/SZ
    # controls: Diagnosis=="Healthy" のみ（family_member除外）
    cc_keep = df["Diagnosis_norm"].isin(["ASD", "SZ", "Healthy"])
    df_cc = df[cc_keep].copy()
    df_cc_out = df_cc[[col_sid, "Diagnosis_norm", "CRAM_Path"]].rename(columns={col_sid: "SampleID"})
    df_cc_out = df_cc_out.rename(columns={"Diagnosis_norm": "Group"})  # Group列にASD/SZ/Healthyのみ
    df_cc_out = df_cc_out.drop_duplicates(subset=["SampleID"])

    cc_file = out_dir / "casecontrol_samples.tsv"
    df_cc_out.to_csv(cc_file, sep="\t", index=False)
    print(f"[{ts()}] [INFO] Wrote {cc_file} (n={len(df_cc_out)})")

    # ---------- trio（de novo）用 ----------
    # proband は ASD または SZ のみを対象（必要に応じて拡張可能）
    trio_links_file = out_dir / "trio_links.tsv"
    trio_samples_file = out_dir / "trio_samples.tsv"

    if father_col is None or mother_col is None:
        print(f"[{ts()}] [WARN] Parent columns not detected in SampleInfo. Skipping trio files.")
    else:
        df_tmp = df.copy()
        df_tmp["FatherID"] = father_col
        df_tmp["MotherID"] = mother_col

        probands = df_tmp[df_tmp["Diagnosis_norm"].isin(["ASD", "SZ"])].copy()
        probands = probands[[col_sid, "Diagnosis_norm", "FatherID", "MotherID"]].rename(columns={col_sid: "ProbandID"})
        # 欠損IDの正規化
        for c in ["FatherID", "MotherID"]:
            probands[c] = probands[c].astype(str).replace({"nan": "NA", "": "NA", ".": "NA"})

        # CRAM存在しているID集合
        available = set(df[col_sid].astype(str).tolist())

        def ok_id(x: str) -> bool:
            return (x is not None) and (str(x) not in ["NA", "nan", "", "."]) and (str(x) in available)

        probands["Father_hasCRAM"] = probands["FatherID"].apply(ok_id)
        probands["Mother_hasCRAM"] = probands["MotherID"].apply(ok_id)

        probands.to_csv(trio_links_file, sep="\t", index=False)
        print(f"[{ts()}] [INFO] Wrote {trio_links_file} (n_probands={len(probands)})")

        # trio_samples = proband + father + mother（CRAMありのみ）
        trio_ids = set()
        for r in probands.itertuples(index=False):
            trio_ids.add(str(r.ProbandID))
            if r.Father_hasCRAM:
                trio_ids.add(str(r.FatherID))
            if r.Mother_hasCRAM:
                trio_ids.add(str(r.MotherID))

        df_trio = df_ehdn_out[df_ehdn_out["SampleID"].isin(trio_ids)].copy()
        df_trio.to_csv(trio_samples_file, sep="\t", index=False)
        print(f"[{ts()}] [INFO] Wrote {trio_samples_file} (n={len(df_trio)})")

    # ---------- counts summary ----------
    stat_f = out_dir / "sample_counts.txt"
    with open(stat_f, "w") as w:
        w.write("=== STR project sample counts (v2) ===\n")
        w.write(f"Created: {ts()}\n")
        w.write(f"SampleInfo: {sample_info}\n\n")

        w.write("[EHdn profiling target]\n")
        w.write(df_ehdn_out["Group"].value_counts().to_string() + "\n\n")

        w.write("[Case-control target (controls=Healthy only)]\n")
        w.write(df_cc_out["Group"].value_counts().to_string() + "\n\n")

        if trio_links_file.exists():
            w.write("[Trio]\n")
            w.write(f"trio_links.tsv: {trio_links_file}\n")
            w.write(f"trio_samples.tsv: {trio_samples_file}\n")

    print(f"[{ts()}] [INFO] Wrote {stat_f}")

    end = datetime.now()
    print(f"[{ts()}] [DONE] Elapsed={(end-start).total_seconds():.1f}s")


if __name__ == "__main__":
    main()
