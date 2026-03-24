#!/usr/bin/env python3
# ============================================================================
# NOTE: This script contains hardcoded default paths specific to the original
# analysis environment (NIG supercomputer). To run in a different environment,
# update the paths below or set the corresponding environment variables
# (e.g. SAMPLE_INFO, PCA_EIGENVEC, CRAM_BASE_DIR1, CRAM_BASE_DIR2).
# ============================================================================
# 14_outlier_burden_rare_casecontrol_crossfit_v19.py
#
# v19: observed_clusters_total 追加 + ヘッダコメント修正
# - v18 High-PPV Edition をベースに以下を追加:
#   (1) per_sample出力に observed_clusters_total 列を追加
#       定義: motif filter + blacklist filter 後に、そのサンプルで≥1 IRRが検出された
#       ユニークcluster数。recurrent-hit filter前のdf_scから計算。
#       EHdnはreference-freeのため、サンプルごとに観測clusterが異なり、
#       これがそのサンプルの"observation opportunity"を反映する。
#   (2) rare_any 列を per_sample に直接追加
#   (3) ヘッダコメントの MIN_IRR_RAW デフォルト値を 10.0 に修正（v18で stale だった）
#   (4) observed_clusters_total の group 別要約をログ出力
# - v18 からの変更なし: motif filter, blacklist, clustering, Z-score, recurrent-hit filter
# - 実行時間を記録

from __future__ import annotations

import hashlib
import math
import os
import sys
import time
import subprocess
import tempfile
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

# -------------------------
# High-PPV Filtering Logic
# -------------------------

def is_artifact_motif(motif: str) -> bool:
    """
    【強化版】ノイズになりやすいモチーフを徹底的に除外する
    """
    m = str(motif).upper()
    
    # 1. Telomere / Centromere / Known Artifacts
    if "AACCCT" in m: return True  # Telomere
    if "AATGG" in m:  return True  # Centromere / Satellite
    if "ACATCC" in m: return True  # Specific Artifact
    
    # 2. Noisy Trinucleotides (PCR stutter prone)
    if m in ["AAC", "AAG", "ACT", "ATC"]: return True
    
    return False

def filter_with_blacklist(df: pd.DataFrame, blacklist_bed: Path) -> pd.DataFrame:
    """bedtools intersect -v でブラックリスト領域を除外"""
    if df.empty: return df
    req_cols = ["contig", "start", "end"]
    if not set(req_cols).issubset(df.columns): return df

    df["_tmp_idx"] = df.index
    tmp_bed_in = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".bed")
    try:
        sub = df[req_cols + ["_tmp_idx"]].copy()
        sub["start"] = sub["start"].astype(int)
        sub["end"] = sub["end"].astype(int)
        sub.to_csv(tmp_bed_in, sep="\t", header=False, index=False)
        tmp_bed_in.close()

        cmd = ["bedtools", "intersect", "-v", "-a", tmp_bed_in.name, "-b", str(blacklist_bed)]
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        keep_indices = set()
        for line in res.stdout.strip().split("\n"):
            if not line: continue
            parts = line.split("\t")
            if len(parts) >= 4:
                keep_indices.add(int(parts[3]))
        
        df_out = df[df["_tmp_idx"].isin(keep_indices)].copy()
        df_out.drop(columns=["_tmp_idx"], inplace=True)
        return df_out
    except Exception as e:
        print(f"[WARN] bedtools filter failed: {e}")
        return df
    finally:
        if os.path.exists(tmp_bed_in.name): os.remove(tmp_bed_in.name)

# -------------------------
# Utilities
# -------------------------
def ts() -> str: return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def detect_col(cols, candidates, required=False):
    low = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in low: return low[cand.lower()]
    if required: raise KeyError(f"Missing col: {candidates}")
    return None

def normalize_sex_to_male(x):
    s = str(x).strip().upper()
    if s in {"M", "MALE", "1"}: return 1
    if s in {"F", "FEMALE", "2"}: return 0
    try:
        if int(float(s)) == 1: return 1
        if int(float(s)) == 2: return 0
    except: pass
    return None

def canonical_repeat_unit(motif: str) -> str:
    m = str(motif).upper()
    if not m: return m
    rots = [m[i:] + m[:i] for i in range(len(m))]
    best = min(rots)
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    rc = "".join(comp.get(b, "N") for b in m[::-1])
    rots_rc = [rc[i:] + rc[:i] for i in range(len(rc))]
    return min(best, min(rots_rc))

def stable_case_fold(sample_id, k, seed):
    h = hashlib.md5(f"{seed}:{sample_id}".encode("utf-8")).hexdigest()
    return int(h[:8], 16) % k

# -------------------------
# IO & Clustering
# -------------------------
def read_casecontrol(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    c_id = detect_col(df.columns, ["SampleID", "sample"], True)
    c_grp = detect_col(df.columns, ["Group"], True)
    return df[[c_id, c_grp]].rename(columns={c_id: "SampleID", c_grp: "Group"})[lambda d: d["Group"].isin(["Healthy", "ASD", "SZ"])]

def read_sampleinfo(path):
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    c_id = detect_col(df.columns, ["SampleID"], True)
    ped = df.get("Pedigree_No", df[c_id]).astype(str).replace({"": "NA", "nan": "NA", ".": "NA"})
    ped[ped=="NA"] = df[c_id][ped=="NA"]
    out = pd.DataFrame({"SampleID": df[c_id], "Pedigree_No": ped})
    c_sex = detect_col(df.columns, ["Sex"], False)
    out["Sex_M"] = df[c_sex].map(normalize_sex_to_male) if c_sex else np.nan
    return out

def read_depths(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    c_id = detect_col(df.columns, ["sample", "SampleID"], True)
    c_dp = detect_col(df.columns, ["depth"], True)
    return dict(zip(df[c_id], pd.to_numeric(df[c_dp], errors="coerce")))

def read_eigenvec(path):
    try: df = pd.read_csv(path, sep=r"\s+", dtype=str, engine="python")
    except: df = pd.read_csv(path, sep="\t", dtype=str)
    if "PC1" not in df.columns:
        df = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
        df.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, df.shape[1]-1)]
    c_id = detect_col(df.columns, ["IID", "SampleID"], True)
    cols = [c_id] + [f"PC{i}" for i in range(1, 11)]
    return df[cols].set_axis(["SampleID"] + [f"PC{i}" for i in range(1, 11)], axis=1)

def choose_one(df, cols):
    tmp = df.assign(_ok=df[cols].notna().all(axis=1).astype(int), _prio=df["Group"].isin(["ASD", "SZ"]).map({True:0, False:1}))
    return tmp.sort_values(["Pedigree_No", "_ok", "_prio", "SampleID"], ascending=[True, False, True, True]).groupby("Pedigree_No", as_index=False).head(1).drop(columns=["_ok", "_prio"])

def build_clusters(locus_df, cluster_bp):
    df = locus_df.sort_values(["contig", "canon", "start", "end"]).reset_index(drop=True)
    rows, map_dict = [], {}
    c_idx, c_chr, c_can, c_s, c_e = 0, None, None, None, None
    for r in df.itertuples(index=False):
        new = (c_chr is None) or (r.contig != c_chr) or (r.canon != c_can) or ((r.start - c_s) > cluster_bp)
        if new:
            if c_chr: 
                rows.append((c_idx, c_chr, c_can, c_s, c_e))
                c_idx += 1
            c_chr, c_can, c_s, c_e = r.contig, r.canon, r.start, r.end
        else:
            c_s, c_e = min(c_s, r.start), max(c_e, r.end)
        map_dict[f"{r.contig}:{r.start}-{r.end}:{r.motif}"] = c_idx
    if c_chr: rows.append((c_idx, c_chr, c_can, c_s, c_e))
    return map_dict, pd.DataFrame(rows, columns=["cluster_idx", "contig", "canon_motif", "cluster_start", "cluster_end"])

# -------------------------
# Main
# -------------------------
def main():
    t0 = time.time()
    print(f"[{ts()}] [INFO] Start v19 (High-PPV + observed_clusters_total)")
    root = Path(__file__).resolve().parents[2]
    
    # Configs for High Precision (v18と同一パラメータ)
    cluster_bp = int(os.environ.get("CLUSTER_BP", "1000"))
    min_irr_raw = float(os.environ.get("MIN_IRR_RAW", "10.0"))  # v19: コメント修正 (実値は v18 から 10.0)
    max_hits = int(os.environ.get("MAX_OUTLIER_HITS", "10"))
    
    print(f"[{ts()}] [CONFIG] CLUSTER_BP = {cluster_bp}")
    print(f"[{ts()}] [CONFIG] MIN_IRR_RAW = {min_irr_raw}")
    print(f"[{ts()}] [CONFIG] MAX_OUTLIER_HITS = {max_hits}")

    merged_tsv = Path(os.environ.get("MERGED_NOVEL_TSV", str(root / "merged_results_novel" / "novel_loci_genic_norm.tsv")))
    blacklist = Path(os.environ.get("BLACKLIST_BED", str(root / "resources" / "blacklist" / "hg38_ehdn_blacklist.bed")))
    outdir = Path(os.environ.get("ANALYSIS_NOVEL_DIR", str(root / "analysis_results_novel")))
    outdir.mkdir(parents=True, exist_ok=True)

    # Load Data
    df_cc = read_casecontrol(Path(os.environ.get("CASECONTROL_SAMPLES", str(root / "sample_lists" / "casecontrol_samples.tsv"))))
    df_si = read_sampleinfo(Path(os.environ.get("SAMPLE_INFO", "/lustre12/home/kushima-pg/sampleInfo/GRIFIN_srWGS_SampleInfo_11242025.txt")))  # CONFIGURE
    dmap = read_depths(Path(os.environ.get("DEPTHS_ALL_TSV", str(root / "depth" / "depths_all.tsv"))))
    pca = read_eigenvec(Path(os.environ.get("PCA_EIGENVEC", "/lustre12/home/kushima-pg/PRS/population_stratfication_09012025/results_popstrat_20251006_v7/pca_jpn/pca.eigenvec")))  # CONFIGURE

    df = df_cc.merge(df_si, on="SampleID", how="left").merge(pca, on="SampleID", how="left")
    df["Depth"] = df["SampleID"].map(dmap)
    cols = ["Depth", "Sex_M"] + [f"PC{i}" for i in range(1, 11)]
    for c in cols: df[c] = pd.to_numeric(df[c], errors="coerce")
    df_ped = choose_one(df, cols)
    print(f"[{ts()}] [INFO] Samples after one-per-pedigree: {len(df_ped)}")

    # Folds
    rng = np.random.default_rng(20251228)
    h_ids = df_ped.loc[df_ped["Group"]=="Healthy", "SampleID"].tolist()
    rng.shuffle(h_ids)
    fold_map = {sid: i%5 for i, sid in enumerate(h_ids)}
    for sid in df_ped.loc[df_ped["Group"].isin(["ASD", "SZ"]), "SampleID"]:
        fold_map[sid] = stable_case_fold(sid, 5, 20251228)
    df_ped["fold"] = df_ped["SampleID"].map(fold_map)
    keep = set(df_ped["SampleID"])
    s2g, s2f = dict(zip(df_ped["SampleID"], df_ped["Group"])), dict(zip(df_ped["SampleID"], df_ped["fold"]))

    # Pass A: Filtering & Clustering (v18と同一)
    print(f"[{ts()}] [INFO] Pass A: Filtering & Clustering...")
    loci_set = set()
    h = pd.read_csv(merged_tsv, sep="\t", nrows=1)
    c_sid = detect_col(h.columns, ["SampleID"], True)
    c_c, c_s, c_e, c_m = detect_col(h.columns, ["contig"], True), detect_col(h.columns, ["start"], True), detect_col(h.columns, ["end"], True), detect_col(h.columns, ["motif"], True)
    c_r, c_n = detect_col(h.columns, ["num_anc_irrs"], True), detect_col(h.columns, ["num_anc_irrs_norm"], True)
    
    chunksize = int(os.environ.get("MERGED_CHUNKSIZE", "500000"))
    for chunk in pd.read_csv(merged_tsv, sep="\t", dtype=str, usecols=[c_sid,c_c,c_s,c_e,c_m,c_r,c_n], chunksize=chunksize):
        chunk = chunk[chunk[c_sid].isin(keep)].copy()
        if chunk.empty: continue
        chunk = chunk[~chunk[c_m].apply(is_artifact_motif)].copy()
        if blacklist.exists() and not chunk.empty:
            chunk = chunk.loc[filter_with_blacklist(chunk.rename(columns={c_c:"contig",c_s:"start",c_e:"end"}), blacklist).index].copy()
        if chunk.empty: continue
        for r in chunk[[c_c, c_s, c_e, c_m]].drop_duplicates().itertuples(index=False):
            loci_set.add(f"{r[0]}:{r[1]}-{r[2]}:{r[3]}")
    
    locus_df = pd.DataFrame([(k.split(":")[0], int(k.split(":")[1].split("-")[0]), int(k.split(":")[1].split("-")[1]), k.split(":")[2], canonical_repeat_unit(k.split(":")[2])) for k in loci_set], columns=["contig", "start", "end", "motif", "canon"])
    l2c, c_tbl = build_clusters(locus_df, cluster_bp)
    C = len(c_tbl)
    print(f"[{ts()}] [INFO] Unique Clusters (post motif+blacklist filter): {C:,}")

    # Pass B: Aggregating (v18と同一)
    print(f"[{ts()}] [INFO] Pass B: Aggregating...")
    rows = []
    for chunk in pd.read_csv(merged_tsv, sep="\t", dtype=str, usecols=[c_sid,c_c,c_s,c_e,c_m,c_r,c_n], chunksize=chunksize):
        chunk = chunk[chunk[c_sid].isin(keep)].copy()
        if chunk.empty: continue
        chunk = chunk[~chunk[c_m].apply(is_artifact_motif)].copy()
        if blacklist.exists() and not chunk.empty:
            chunk = chunk.loc[filter_with_blacklist(chunk.rename(columns={c_c:"contig",c_s:"start",c_e:"end"}), blacklist).index].copy()
        
        chunk["cluster_idx"] = (chunk[c_c]+":"+chunk[c_s]+"-"+chunk[c_e]+":"+chunk[c_m]).map(l2c)
        chunk = chunk.dropna(subset=["cluster_idx"]).assign(
            cluster_idx=lambda x: x["cluster_idx"].astype(int),
            raw=lambda x: pd.to_numeric(x[c_r], errors="coerce").fillna(0),
            norm=lambda x: pd.to_numeric(x[c_n], errors="coerce").fillna(0)
        )
        rows.append(chunk.groupby([c_sid, "cluster_idx"], as_index=False).agg({"raw":"sum", "norm":"sum"}))
    
    df_sc = pd.concat(rows).groupby([c_sid, "cluster_idx"], as_index=False).sum()
    df_sc = df_sc.rename(columns={c_sid: "SampleID", "raw": "raw_sum", "norm": "norm_sum"}).assign(
        Group=lambda x: x["SampleID"].map(s2g), fold=lambda x: x["SampleID"].map(s2f)
    )

    # ---- v19 NEW: observed_clusters_total を df_sc から計算 ----
    # recurrent-hit filter 前の df_sc から計算する。
    # recurrent-hit filter は outlier pruning であり、observation opportunity とは別。
    obs_clusters = df_sc.groupby("SampleID")["cluster_idx"].nunique().reset_index()
    obs_clusters.columns = ["SampleID", "observed_clusters_total"]
    print(f"[{ts()}] [INFO] observed_clusters_total: computed from df_sc (pre-recurrent-hit filter)")

    # Z-scores & Outlier detection (v18と同一)
    print(f"[{ts()}] [INFO] Calculating Z-scores...")
    sum_all, sumsq_all = np.zeros(C), np.zeros(C)
    sum_fold, sumsq_fold = np.zeros((5, C)), np.zeros((5, C))
    n_h = np.zeros(5, dtype=int)
    
    df_h = df_sc[df_sc["Group"]=="Healthy"]
    np.add.at(sum_all, df_h["cluster_idx"].values, df_h["norm_sum"].values)
    np.add.at(sumsq_all, df_h["cluster_idx"].values, df_h["norm_sum"].values**2)
    for f in range(5):
        n_h[f] = ((df_ped["Group"]=="Healthy")&(df_ped["fold"]==f)).sum()
        mask = (df_h["fold"]==f)
        np.add.at(sum_fold[f], df_h.loc[mask, "cluster_idx"].values, df_h.loc[mask, "norm_sum"].values)
        np.add.at(sumsq_fold[f], df_h.loc[mask, "cluster_idx"].values, df_h.loc[mask, "norm_sum"].values**2)

    details = []
    n_total = n_h.sum()
    zthr, rare_cut = float(os.environ.get("Z_THR", "5.0")), float(os.environ.get("RARE_CUT", "0.001"))

    sc_sid, sc_cid, sc_raw, sc_nrm, sc_fld, sc_grp = df_sc["SampleID"].values, df_sc["cluster_idx"].values, df_sc["raw_sum"].values, df_sc["norm_sum"].values, df_sc["fold"].values.astype(int), df_sc["Group"].values
    
    for f in range(5):
        nt = n_total - n_h[f]
        if nt < 2: continue
        mean = (sum_all - sum_fold[f])/nt
        var = ((sumsq_all - sumsq_fold[f])/nt) - mean**2
        sd = np.sqrt(np.maximum(var, 1e-9))
        
        mask_tr = (sc_grp=="Healthy")&(sc_fld!=f)
        z_tr = (sc_nrm[mask_tr] - mean[sc_cid[mask_tr]]) / sd[sc_cid[mask_tr]]
        is_out_tr = (z_tr > zthr) & (sc_raw[mask_tr] >= min_irr_raw)
        freq = np.bincount(sc_cid[mask_tr][is_out_tr], minlength=C) / nt
        
        mask_te = (sc_fld==f)
        if not mask_te.any(): continue
        z_te = (sc_nrm[mask_te] - mean[sc_cid[mask_te]]) / sd[sc_cid[mask_te]]
        is_out = (z_te > zthr) & (sc_raw[mask_te] >= min_irr_raw) & (freq[sc_cid[mask_te]] < rare_cut)
        
        if is_out.any():
            for i in np.where(is_out)[0]:
                details.append({
                    "SampleID": sc_sid[mask_te][i], "Group": s2g[sc_sid[mask_te][i]], "fold": f,
                    "cluster_idx": sc_cid[mask_te][i], "raw_irr": sc_raw[mask_te][i],
                    "norm_irr": sc_nrm[mask_te][i], "Z_score": z_te[i], "rare_freq": freq[sc_cid[mask_te][i]]
                })

    # Recurrent Filter (v18と同一: global post-hoc pruning)
    print(f"[{ts()}] [INFO] Recurrent Hit Filter (> {max_hits} samples)...")
    if details:
        cnt = Counter(d["cluster_idx"] for d in details)
        bad = {c for c, n in cnt.items() if n > max_hits}
        print(f"[{ts()}] [INFO] Removing {len(bad)} noisy clusters.")
        details = [d for d in details if d["cluster_idx"] not in bad]
    
    # Save per_sample
    s2i = {s: i for i, s in enumerate(df_ped["SampleID"])}
    counts = np.zeros(len(df_ped), dtype=int)
    for d in details: 
        if d["SampleID"] in s2i: counts[s2i[d["SampleID"]]] += 1
    
    df_out = df_ped.assign(rare_outlier_count=counts)
    df_out["rare_any"] = (df_out["rare_outlier_count"] >= 1).astype(int)
    
    # ---- v19 NEW: observed_clusters_total を per_sample に結合 ----
    df_out = df_out.merge(obs_clusters, on="SampleID", how="left")
    df_out["observed_clusters_total"] = df_out["observed_clusters_total"].fillna(0).astype(int)
    
    # ---- v19 NEW: group別 observed_clusters_total 要約 ----
    print(f"\n[{ts()}] [INFO] === observed_clusters_total summary (per Group) ===")
    for grp in ["Healthy", "ASD", "SZ"]:
        sub = df_out.loc[df_out["Group"] == grp, "observed_clusters_total"]
        if sub.empty:
            continue
        print(f"  {grp}: N={len(sub)}, "
              f"mean={sub.mean():.1f}, median={sub.median():.0f}, "
              f"min={sub.min()}, max={sub.max()}, "
              f"N_zero={int((sub == 0).sum())}")
    
    # exposure=0 のサンプルを警告（CRAMの処理失敗の可能性）
    zero_exp = df_out.loc[df_out["observed_clusters_total"] == 0, "SampleID"].tolist()
    if zero_exp:
        print(f"\n[{ts()}] [WARN] {len(zero_exp)} samples have observed_clusters_total=0 "
              f"(may indicate EHdn processing failure):")
        for sid in zero_exp[:20]:
            print(f"    {sid}")
        if len(zero_exp) > 20:
            print(f"    ... and {len(zero_exp) - 20} more")

    out_per_sample = outdir / "outlier_burden_rare_crossfit_v19.per_sample.tsv"
    df_out.to_csv(out_per_sample, sep="\t", index=False)
    print(f"\n[{ts()}] [INFO] Wrote {out_per_sample} (n={len(df_out)}, cols={list(df_out.columns)})")
    
    # Save outlier_details
    out_details = outdir / "outlier_burden_rare_crossfit_v19.outlier_details.tsv"
    if details:
        pd.DataFrame(details).merge(c_tbl, on="cluster_idx", how="left").to_csv(
            out_details, sep="\t", index=False, float_format="%.4f")
    else:
        pd.DataFrame().to_csv(out_details, sep="\t")
    print(f"[{ts()}] [INFO] Wrote {out_details} (n_outlier_events={len(details)})")

    elapsed = time.time() - t0
    print(f"\n[{ts()}] [DONE] v19 completed. Elapsed={elapsed:.1f}s")

if __name__ == "__main__":
    main()
