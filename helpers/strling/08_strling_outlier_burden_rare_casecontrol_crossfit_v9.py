#!/usr/bin/env python3
# ============================================================================
# NOTE: This script contains hardcoded default paths specific to the original
# analysis environment (NIG supercomputer). To run in a different environment,
# update the paths below or set the corresponding environment variables
# (e.g. SAMPLE_INFO, PCA_EIGENVEC, CRAM_BASE_DIR1, CRAM_BASE_DIR2).
# ============================================================================
# 08_strling_outlier_burden_rare_casecontrol_crossfit_v9.py
# - v9変更点: Poisson GLMにexposure>0ガードを追加（exposure=0サンプルの除外＋ログ出力）
#
# 処理内容:
#   - STRling outliers の STRs.tsv（in-bounds loci）を入力に、rare outlier burden を case-control で推定
#   - Healthy を K-fold に分割し、各サンプルは「自分のfoldを除いた train Healthy」で rare 定義（freq<rare_cut）される
#   - event 定義: outlier(z) > zthr かつ p_adj <= p_adj_thr かつ sum_str_counts >= min_sum_str_counts
#   - 近傍マージ (Clustering): 同一サンプル内で同一chr/motif/距離 < merge_dist の外れ値を1イベントにまとめる
#   - 出力: per_sample.tsv / group_summary.long.tsv / thresholds.tsv / model_summary.txt / burden_stats_results.tsv
#           / outlier_details.tsv (v8で追加)
#
# v8 変更点 (v7 → v8):
#   - merge_close_hits を改修: マージ後の代表hitの詳細情報（chrom, left, right, repeatunit, z,
#     p_adj, sum_str_counts, locus, rare_freq）を返す
#   - Pass 2 で各outlierの追加属性（p_adj, sum_str_counts, locus）も保持
#   - outlier_details.tsv を出力: 全サンプルの全rare outlier loci（マージ後）を1行1 locusで記録
#   - EHdn outlier_details.tsv と比較可能な形式
#   - burden計算・統計テストなど他の処理はv7と完全に同一
#
# v7 変更点 (v6 → v7):
#   - min_sum_str_counts のデフォルトを 0.0 → 1.0 に変更
#
# 使い方:
#   python 08_strling_outlier_burden_rare_casecontrol_crossfit_v9.py --merge_dist 1000

from __future__ import annotations

import argparse
import hashlib
import math
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from statsmodels.discrete.discrete_model import NegativeBinomial as SMNegBin
except Exception:
    sm = None
    smf = None
    SMNegBin = None

try:
    from scipy.stats import mannwhitneyu
except Exception:
    mannwhitneyu = None


def ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def detect_col(cols: Iterable[str], candidates: List[str]) -> Optional[str]:
    low = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    return None


def normalize_sex_to_male(x: object) -> Optional[int]:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    s = str(x).strip()
    if s == "":
        return None
    s_up = s.upper()
    if s_up in {"M", "MALE"}:
        return 1
    if s_up in {"F", "FEMALE"}:
        return 0
    if s in {"1"}:
        return 1
    if s in {"2"}:
        return 0
    try:
        v = int(float(s))
        if v == 1:
            return 1
        if v == 2:
            return 0
    except Exception:
        pass
    return None


def read_depths(depth_tsv: Path) -> Dict[str, float]:
    df = pd.read_csv(depth_tsv, sep="\t", dtype=str)
    c_sample = detect_col(df.columns, ["sample", "SampleID", "iid", "IID"])
    c_depth = detect_col(df.columns, ["depth", "Depth"])
    if c_sample is None or c_depth is None:
        raise ValueError(f"Depth file must contain sample/depth columns: {depth_tsv}")
    df = df[[c_sample, c_depth]].copy()
    df.columns = ["SampleID", "Depth"]
    df["SampleID"] = df["SampleID"].astype(str)
    df["Depth"] = pd.to_numeric(df["Depth"], errors="coerce")
    df = df.dropna(subset=["SampleID", "Depth"])
    return dict(zip(df["SampleID"], df["Depth"].astype(float)))


def read_casecontrol_samples(samples_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(samples_tsv, sep="\t", dtype=str)
    c_id = detect_col(df.columns, ["SampleID", "sample", "IID"])
    c_grp = detect_col(df.columns, ["Group", "Diagnosis", "group"])
    if c_id is None or c_grp is None:
        raise ValueError(f"casecontrol_samples.tsv must have SampleID and Group columns: {samples_tsv}")
    out = df[[c_id, c_grp]].copy()
    out.columns = ["SampleID", "Group"]
    out["SampleID"] = out["SampleID"].astype(str).str.strip()
    out["Group"] = out["Group"].astype(str).str.strip()
    out = out[out["Group"].isin(["Healthy", "ASD", "SZ"])].copy()
    return out


def read_sampleinfo(sampleinfo_tsv: Path, pedigree_col: str = "Pedigree_No") -> pd.DataFrame:
    df = pd.read_csv(sampleinfo_tsv, sep="\t", dtype=str, low_memory=False)
    c_id = detect_col(df.columns, ["SampleID", "sample", "id", "IID"])
    if c_id is None:
        raise ValueError(f"SampleInfo needs SampleID-like column: {sampleinfo_tsv}")
    sid = df[c_id].astype(str).str.strip()
    if pedigree_col in df.columns:
        ped_raw = df[pedigree_col].astype(str).str.strip()
    else:
        ped_raw = sid.copy()
    missing_tokens = {"", "NA", "N/A", "nan", "NaN", ".", "None", "NULL", "null", "Unknown", "unknown"}
    ped = ped_raw.copy()
    mask_missing = ped.isna() | ped.isin(missing_tokens)
    if mask_missing.any():
        ped.loc[mask_missing] = sid.loc[mask_missing]
    c_sex = detect_col(df.columns, ["Sex", "sex", "Gender", "gender"])
    out = pd.DataFrame({"SampleID": sid})
    out["Pedigree_No"] = ped
    if c_sex is not None:
        out["Sex_M"] = df[c_sex].map(normalize_sex_to_male)
    else:
        out["Sex_M"] = np.nan
    return out


def read_eigenvec(eigenvec_path: Path, n_pcs: int = 10) -> pd.DataFrame:
    try:
        df = pd.read_csv(eigenvec_path, sep=r"\s+", dtype=str, engine="python")
    except Exception:
        df = pd.read_csv(eigenvec_path, sep="\t", dtype=str)
    if detect_col(df.columns, ["PC1"]) is None:
        df2 = pd.read_csv(eigenvec_path, sep=r"\s+", header=None, engine="python")
        if df2.shape[1] < 2 + n_pcs:
            raise ValueError(f"eigenvec seems too few columns: {eigenvec_path}")
        cols = ["FID", "IID"] + [f"PC{i}" for i in range(1, df2.shape[1] - 1)]
        df2.columns = cols[: df2.shape[1]]
        df = df2
    c_iid = detect_col(df.columns, ["IID", "SampleID", "sample"])
    if c_iid is None:
        c_iid = detect_col(df.columns, ["iid"])
    if c_iid is None:
        raise ValueError(f"Could not detect IID column in eigenvec: {eigenvec_path}")
    pcs = []
    for i in range(1, n_pcs + 1):
        c = detect_col(df.columns, [f"PC{i}"])
        if c is None:
            raise ValueError(f"Missing PC{i} in eigenvec: {eigenvec_path}")
        pcs.append(c)
    out = df[[c_iid] + pcs].copy()
    out.columns = ["SampleID"] + [f"PC{i}" for i in range(1, n_pcs + 1)]
    out["SampleID"] = out["SampleID"].astype(str).str.strip()
    for i in range(1, n_pcs + 1):
        out[f"PC{i}"] = pd.to_numeric(out[f"PC{i}"], errors="coerce")
    return out


def choose_one_per_pedigree_prefer_complete(df: pd.DataFrame, need_cols: List[str]) -> pd.DataFrame:
    tmp = df.copy()
    for c in need_cols:
        if c not in tmp.columns:
            tmp[c] = np.nan
    tmp["_complete_ok"] = tmp[need_cols].notna().all(axis=1).astype(int)
    tmp["_case_prio"] = tmp["Group"].map(lambda g: 0 if g in ["ASD", "SZ"] else 1).astype(int)
    tmp = tmp.sort_values(
        ["Pedigree_No", "_complete_ok", "_case_prio", "SampleID"],
        ascending=[True, False, True, True]
    )
    out = tmp.groupby("Pedigree_No", as_index=False).head(1).drop(columns=["_complete_ok", "_case_prio"])
    return out


def assign_folds(healthy_ids: List[str], k: int, seed: int) -> Dict[str, int]:
    rng = np.random.default_rng(seed)
    ids = healthy_ids.copy()
    rng.shuffle(ids)
    return {sid: (i % k) for i, sid in enumerate(ids)}


def fold_for_case(sample_id: str, k: int, seed: int) -> int:
    h = hashlib.md5(f"{seed}:{sample_id}".encode("utf-8")).hexdigest()
    return int(h[:8], 16) % k


def load_inbounds_loci(bounds_file: Path) -> Tuple[List[str], Dict[str, int]]:
    df = pd.read_csv(bounds_file, sep="\t", dtype=str)
    c_chr = detect_col(df.columns, ["#chrom", "chrom", "chr", "contig"])
    c_left = detect_col(df.columns, ["left", "start"])
    c_right = detect_col(df.columns, ["right", "end", "stop"])
    c_rep = detect_col(df.columns, ["repeat", "repeatunit", "repeat_unit"])
    if c_chr is None or c_left is None or c_right is None or c_rep is None:
        raise ValueError(f"Bounds must contain chrom/left/right/repeat columns: {bounds_file}")
    loci_list = (df[c_chr].astype(str) + "-" + df[c_left].astype(str) + "-" +
                 df[c_right].astype(str) + "-" + df[c_rep].astype(str)).tolist()
    locus2idx = {l: i for i, l in enumerate(loci_list)}
    return loci_list, locus2idx


def merge_close_hits_with_details(hits_df: pd.DataFrame, dist_limit: int) -> Tuple[int, List[dict]]:
    """
    同一サンプル内で、chromが同じ、motifが同じ、かつ距離が dist_limit 以内のヒットをまとめる。
    グループ内で最大のZスコアを持つ行を代表とする。
    返り値: (マージ後のヒット数, 代表hitの詳細リスト)

    v8追加: 代表hitの詳細情報も返す
    """
    if hits_df.empty:
        return 0, []

    df_sorted = hits_df.sort_values(by=["chrom", "left"])
    merged_count = 0
    representative_hits = []

    for _, group in df_sorted.groupby(["chrom", "repeatunit"]):
        if len(group) == 1:
            merged_count += 1
            representative_hits.append(group.iloc[0].to_dict())
            continue

        group = group.sort_values("left")
        clusters = []
        current_cluster = [group.iloc[0]]
        current_cluster_end = int(group.iloc[0]["right"])

        for i in range(1, len(group)):
            row = group.iloc[i]
            left = int(row["left"])
            right = int(row["right"])
            dist = left - current_cluster_end

            if dist < dist_limit:
                current_cluster.append(row)
                current_cluster_end = max(current_cluster_end, right)
            else:
                clusters.append(current_cluster)
                current_cluster = [row]
                current_cluster_end = right

        clusters.append(current_cluster)

        for cluster in clusters:
            merged_count += 1
            best = max(cluster, key=lambda r: r["z"])
            representative_hits.append(best.to_dict())

    return merged_count, representative_hits


# ============================================================
# Per-comparison 統計テスト関数
# ============================================================

def run_per_comparison_tests(df_m: pd.DataFrame, case_grps: List[str],
                             exposure_col: str) -> pd.DataFrame:
    cov_formula = "Depth + Sex_M + " + " + ".join([f"PC{i}" for i in range(1, 11)])
    df_ctrl = df_m[df_m["Group"] == "Healthy"].copy()
    results = []

    for case_grp in case_grps:
        df_case = df_m[df_m["Group"] == case_grp].copy()
        df_sub = pd.concat([df_ctrl, df_case]).copy()
        df_sub["IsCase"] = (df_sub["Group"] == case_grp).astype(int)

        row = {
            "Comparison": f"{case_grp} vs Healthy",
            "N_Case": len(df_case),
            "N_Ctrl": len(df_ctrl),
            "rare_any_Case": df_case["rare_any"].mean(),
            "rare_any_Ctrl": df_ctrl["rare_any"].mean(),
            "mean_count_Case": df_case["rare_outlier_count"].mean(),
            "mean_count_Ctrl": df_ctrl["rare_outlier_count"].mean(),
        }

        if mannwhitneyu is not None:
            try:
                stat, p_mwu = mannwhitneyu(
                    df_case["rare_outlier_count"].values,
                    df_ctrl["rare_outlier_count"].values,
                    alternative="greater")
                row["MWU_P"] = p_mwu
            except Exception as e:
                print(f"[WARN] MWU failed for {case_grp}: {e}")
                row["MWU_P"] = np.nan
        else:
            row["MWU_P"] = np.nan

        formula_logit = f"rare_any ~ IsCase + {cov_formula}"
        try:
            model_logit = smf.logit(formula=formula_logit, data=df_sub).fit(
                disp=0, method="bfgs", maxiter=200)
            coef = model_logit.params["IsCase"]
            ci = model_logit.conf_int().loc["IsCase"]
            row["Logit_OR"] = np.exp(coef)
            row["Logit_OR_CI_low"] = np.exp(ci[0])
            row["Logit_OR_CI_high"] = np.exp(ci[1])
            row["Logit_P"] = model_logit.pvalues["IsCase"]
        except Exception as e:
            print(f"[WARN] Logistic regression failed for {case_grp}: {e}")
            row["Logit_OR"] = np.nan
            row["Logit_OR_CI_low"] = np.nan
            row["Logit_OR_CI_high"] = np.nan
            row["Logit_P"] = np.nan

        formula_pois = f"rare_outlier_count ~ IsCase + {cov_formula}"
        try:
            # exposure>0 ガード: log(0) 未定義を回避
            mask_exp = df_sub[exposure_col].astype(float) > 0
            n_excluded = (~mask_exp).sum()
            if n_excluded > 0:
                excluded_ids = df_sub.loc[~mask_exp, "SampleID"].tolist()
                print(f"[WARN] {case_grp}: excluded {n_excluded} samples with {exposure_col}=0: {excluded_ids}")
                df_sub = df_sub[mask_exp].copy()
            offset_vals = np.log(df_sub[exposure_col].astype(float).values)
            model_pois = smf.glm(
                formula=formula_pois, data=df_sub,
                family=sm.families.Poisson(),
                offset=offset_vals).fit()
            coef_p = model_pois.params["IsCase"]
            ci_p = model_pois.conf_int().loc["IsCase"]
            row["Poisson_RR"] = np.exp(coef_p)
            row["Poisson_RR_CI_low"] = np.exp(ci_p[0])
            row["Poisson_RR_CI_high"] = np.exp(ci_p[1])
            row["Poisson_P"] = model_pois.pvalues["IsCase"]
        except Exception as e:
            print(f"[WARN] Poisson GLM failed for {case_grp}: {e}")
            row["Poisson_RR"] = np.nan
            row["Poisson_RR_CI_low"] = np.nan
            row["Poisson_RR_CI_high"] = np.nan
            row["Poisson_P"] = np.nan

        results.append(row)

    return pd.DataFrame(results)


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    out_root = project_root / "strling_output_genomewide"

    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("--strs_tsv", default=str(out_root / "outliers_casecontrol_inbounds_v1" / "links" / "STRs.tsv"))
    ap.add_argument("--bounds_file", default=str(out_root / "str-results" / "joint-bounds.genic_1kbpad.len3_8.txt"))
    ap.add_argument("--casecontrol_samples", default=str(project_root / "sample_lists" / "casecontrol_samples.tsv"))
    ap.add_argument("--sampleinfo", default="/lustre12/home/kushima-pg/sampleInfo/GRIFIN_srWGS_SampleInfo_11242025.txt")  # CONFIGURE
    ap.add_argument("--pedigree_col", default="Pedigree_No")
    ap.add_argument("--pca_eigenvec", default="/lustre12/home/kushima-pg/PRS/population_stratfication_09012025/results_popstrat_20251006_v7/pca_jpn/pca.eigenvec")  # CONFIGURE
    ap.add_argument("--depths_tsv", default=str(project_root / "depth" / "depths_all.tsv"))
    ap.add_argument("--kfold", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20251228)
    ap.add_argument("--zthr", type=float, default=5.0)
    ap.add_argument("--p_adj_thr", type=float, default=0.05)
    ap.add_argument("--rare_cut", type=float, default=0.001)
    ap.add_argument("--min_sum_str_counts", type=float, default=1.0,
                        help="Minimum sum_str_counts to qualify as an outlier event. "
                             "sum_str_counts=0 entries lack STR read evidence and are "
                             "depth-driven artifacts (default: 1.0)")
    ap.add_argument("--merge_dist", type=int, default=1000, help="Distance to merge adjacent outliers (bp)")

    mx = ap.add_mutually_exclusive_group()
    mx.add_argument("--enforce_one_per_pedigree", action="store_true", default=True)
    mx.add_argument("--no_enforce_one_per_pedigree", action="store_true", default=False)

    ap.add_argument("--chunksize", type=int, default=500_000)
    ap.add_argument("--outdir", default=str(project_root / "analysis_results_strling"))
    ap.add_argument("--prefix", default="strling_outlier_burden_rare_crossfit_inbounds_v9")
    args = ap.parse_args()

    t0 = time.time()
    print(f"[{ts()}] [INFO] Start {Path(__file__).name}")
    print(f"[{ts()}] [INFO] zthr={args.zthr} rare_cut={args.rare_cut} merge_dist={args.merge_dist} "
          f"min_sum_str_counts={args.min_sum_str_counts}")

    strs_path = Path(args.strs_tsv)
    if not strs_path.exists():
        raise SystemExit(f"[ERROR] STRs.tsv not found: {strs_path}")

    inbounds_loci, locus2idx = load_inbounds_loci(Path(args.bounds_file))
    loci_set = set(inbounds_loci)
    L = len(inbounds_loci)
    print(f"[{ts()}] [INFO] IN_BOUNDS loci loaded: {L:,}")

    df_cc = read_casecontrol_samples(Path(args.casecontrol_samples))
    df_si = read_sampleinfo(Path(args.sampleinfo), pedigree_col=args.pedigree_col)
    df_pca = read_eigenvec(Path(args.pca_eigenvec), n_pcs=10)
    depth_map = read_depths(Path(args.depths_tsv))

    df = df_cc.merge(df_si[["SampleID", "Pedigree_No", "Sex_M"]], on="SampleID", how="left") \
              .merge(df_pca, on="SampleID", how="left")
    df["Depth"] = df["SampleID"].map(depth_map)

    need_cols = ["Depth", "Sex_M"] + [f"PC{i}" for i in range(1, 11)]
    for c in need_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    enforce = args.enforce_one_per_pedigree and (not args.no_enforce_one_per_pedigree)
    if enforce:
        df = choose_one_per_pedigree_prefer_complete(df, need_cols=need_cols)

    df["tested_loci_total"] = float(L)
    df_keep = df.dropna(subset=need_cols + ["tested_loci_total"]).copy()

    if len(df_keep) == 0:
        raise SystemExit("[ERROR] Complete-case is empty.")

    for grp in ["Healthy", "ASD", "SZ"]:
        n = (df_keep["Group"] == grp).sum()
        print(f"[{ts()}] [INFO] {grp}: N={n}")

    healthy_ids = df_keep.loc[df_keep["Group"] == "Healthy", "SampleID"].tolist()
    fold_map_h = assign_folds(healthy_ids, args.kfold, args.seed)

    fold_map: Dict[str, int] = {}
    for sid, grp in zip(df_keep["SampleID"], df_keep["Group"]):
        if grp == "Healthy":
            fold_map[sid] = fold_map_h[sid]
        else:
            fold_map[sid] = fold_for_case(sid, args.kfold, args.seed)

    df_keep["fold"] = df_keep["SampleID"].map(fold_map).astype(int)

    n_healthy_total = len(healthy_ids)
    n_healthy_fold = df_keep[df_keep["Group"] == "Healthy"].groupby("fold")["SampleID"].count().to_dict()
    for f in range(args.kfold):
        n_healthy_fold.setdefault(f, 0)

    keep_set = set(df_keep["SampleID"].tolist())
    sample2group = dict(zip(df_keep["SampleID"], df_keep["Group"]))
    sample2fold = dict(zip(df_keep["SampleID"], df_keep["fold"]))

    head = pd.read_csv(strs_path, sep="\t", nrows=5, dtype=str)
    c_sample = detect_col(head.columns, ["sample", "SampleID"])
    c_locus = detect_col(head.columns, ["locus"])
    c_outlier = detect_col(head.columns, ["outlier", "z", "zscore"])
    c_padj = detect_col(head.columns, ["p_adj", "padj", "p.adj"])
    c_sum = detect_col(head.columns, ["sum_str_counts", "sum_str_count"])
    c_chrom = detect_col(head.columns, ["#chrom", "chrom", "chr"])
    c_left = detect_col(head.columns, ["left", "start"])
    c_right = detect_col(head.columns, ["right", "end"])
    c_repeat = detect_col(head.columns, ["repeatunit", "repeat_unit", "motif"])

    if c_sample is None or c_locus is None or c_outlier is None or c_padj is None:
        raise SystemExit(f"[ERROR] Required columns (sample, locus, z, padj) not found.")

    usecols = [c_sample, c_locus, c_outlier, c_padj]
    if c_sum:
        usecols.append(c_sum)
    for c in [c_chrom, c_left, c_right, c_repeat]:
        if c and c not in usecols:
            usecols.append(c)

    print(f"[{ts()}] [INFO] Reading columns: {usecols}")

    def is_event(df_chunk: pd.DataFrame) -> pd.Series:
        z = pd.to_numeric(df_chunk[c_outlier], errors="coerce")
        padj = pd.to_numeric(df_chunk[c_padj], errors="coerce")
        ok = (z > args.zthr) & (padj <= args.p_adj_thr)
        if c_sum and args.min_sum_str_counts > 0:
            ss = pd.to_numeric(df_chunk[c_sum], errors="coerce").fillna(0.0)
            ok = ok & (ss >= args.min_sum_str_counts)
        return ok

    # --- Pass 1: Frequency Counting ---
    healthy_total = np.zeros(L, dtype=np.int32)
    healthy_fold = np.zeros((args.kfold, L), dtype=np.int32)

    print(f"[{ts()}] [INFO] Pass1: Counting Healthy carriers per locus...")
    for chunk in pd.read_csv(strs_path, sep="\t", dtype=str, usecols=usecols, chunksize=args.chunksize):
        chunk = chunk[chunk[c_sample].astype(str).isin(keep_set)]
        if chunk.empty:
            continue
        chunk = chunk[is_event(chunk)]
        if chunk.empty:
            continue
        chunk = chunk[chunk[c_locus].astype(str).isin(loci_set)]
        if chunk.empty:
            continue
        sids = chunk[c_sample].astype(str)
        ch = chunk[sids.map(lambda x: sample2group.get(x, "NA") == "Healthy")].copy()
        if ch.empty:
            continue
        ch["fold"] = ch[c_sample].astype(str).map(sample2fold).astype(int)
        ch["idx"] = ch[c_locus].astype(str).map(locus2idx).astype(int)
        vc = ch["idx"].value_counts()
        healthy_total[vc.index.to_numpy()] += vc.to_numpy(dtype=np.int32)
        g2 = ch.groupby(["fold", "idx"]).size()
        for (f, idx), cnt in g2.items():
            healthy_fold[int(f), int(idx)] += int(cnt)

    print(f"[{ts()}] [INFO] Pass1 done")

    # --- Pass 2: Collecting Rare Outliers ---
    # v8変更: p_adj, sum_str_counts, locusも保持
    print(f"[{ts()}] [INFO] Pass2: Collecting rare outliers per sample...")

    sample_rare_hits: Dict[str, List[dict]] = {}

    for chunk in pd.read_csv(strs_path, sep="\t", dtype=str, usecols=usecols, chunksize=args.chunksize):
        chunk = chunk[chunk[c_sample].astype(str).isin(keep_set)]
        if chunk.empty:
            continue
        chunk = chunk[is_event(chunk)]
        if chunk.empty:
            continue
        chunk = chunk[chunk[c_locus].astype(str).isin(loci_set)]
        if chunk.empty:
            continue

        sids = chunk[c_sample].astype(str)
        folds = sids.map(sample2fold).astype(int).to_numpy()
        idxs = chunk[c_locus].astype(str).map(locus2idx).astype(int).to_numpy()

        train_n = np.array([n_healthy_total - n_healthy_fold[int(f)] for f in folds], dtype=np.float64)
        train_c = healthy_total[idxs].astype(np.float64) - healthy_fold[folds, idxs].astype(np.float64)
        train_c[train_c < 0] = 0.0

        freq = np.divide(train_c, train_n, out=np.ones_like(train_c), where=train_n > 0)
        rare = freq < float(args.rare_cut)

        if rare.any():
            rare_chunk = chunk[rare].copy()
            rare_freq_vals = freq[rare]

            if c_chrom and c_chrom in rare_chunk.columns:
                rare_chunk["_chrom"] = rare_chunk[c_chrom]
            else:
                rare_chunk["_chrom"] = rare_chunk[c_locus].apply(
                    lambda x: x.split('-')[0] if '-' in str(x) else 'NA')
            if c_left and c_left in rare_chunk.columns:
                rare_chunk["_left"] = pd.to_numeric(rare_chunk[c_left], errors='coerce')
            else:
                rare_chunk["_left"] = rare_chunk[c_locus].apply(
                    lambda x: int(x.split('-')[1]) if '-' in str(x) else 0)
            if c_right and c_right in rare_chunk.columns:
                rare_chunk["_right"] = pd.to_numeric(rare_chunk[c_right], errors='coerce')
            else:
                rare_chunk["_right"] = rare_chunk[c_locus].apply(
                    lambda x: int(x.split('-')[2]) if '-' in str(x) else 0)
            if c_repeat and c_repeat in rare_chunk.columns:
                rare_chunk["_repeat"] = rare_chunk[c_repeat]
            else:
                rare_chunk["_repeat"] = "NA"

            rare_chunk["_z"] = pd.to_numeric(rare_chunk[c_outlier], errors='coerce')
            rare_chunk["_p_adj"] = pd.to_numeric(rare_chunk[c_padj], errors='coerce')
            rare_chunk["_sum_str_counts"] = pd.to_numeric(
                rare_chunk[c_sum], errors='coerce') if c_sum else 0.0
            rare_chunk["_locus"] = rare_chunk[c_locus].astype(str)
            rare_chunk["_rare_freq"] = rare_freq_vals

            for _, row in rare_chunk.iterrows():
                sid = str(row[c_sample])
                hit = {
                    "chrom": str(row["_chrom"]),
                    "left": int(row["_left"]),
                    "right": int(row["_right"]),
                    "repeatunit": str(row["_repeat"]),
                    "z": float(row["_z"]),
                    "p_adj": float(row["_p_adj"]),
                    "sum_str_counts": float(row["_sum_str_counts"]) if c_sum else 0.0,
                    "locus": str(row["_locus"]),
                    "rare_freq": float(row["_rare_freq"]),
                }
                if sid not in sample_rare_hits:
                    sample_rare_hits[sid] = []
                sample_rare_hits[sid].append(hit)

    print(f"[{ts()}] [INFO] Pass2 done. Start merging clusters (dist<{args.merge_dist}bp)...")

    # --- Merge Logic & Count ---
    # v8変更: merge_close_hits_with_details を使い、代表hitの詳細を保持
    final_rare_counts: Dict[str, int] = {}
    all_outlier_details: List[dict] = []

    for sid, hits_list in sample_rare_hits.items():
        if not hits_list:
            continue
        hits_df = pd.DataFrame(hits_list)
        cnt, rep_hits = merge_close_hits_with_details(hits_df, args.merge_dist)
        final_rare_counts[sid] = cnt

        grp = sample2group.get(sid, "NA")
        fld = sample2fold.get(sid, -1)
        for rh in rep_hits:
            all_outlier_details.append({
                "SampleID": sid,
                "Group": grp,
                "fold": fld,
                "chrom": rh["chrom"],
                "left": int(rh["left"]),
                "right": int(rh["right"]),
                "repeatunit": rh["repeatunit"],
                "locus": rh["locus"],
                "Z_score": rh["z"],
                "p_adj": rh["p_adj"],
                "sum_str_counts": rh["sum_str_counts"],
                "rare_freq": rh["rare_freq"],
            })

    # --- v8追加: outlier_details.tsv 出力 ---
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if all_outlier_details:
        df_details = pd.DataFrame(all_outlier_details)
        df_details = df_details.sort_values(["Group", "SampleID", "chrom", "left"]).reset_index(drop=True)
        details_tsv = outdir / f"{args.prefix}.outlier_details.tsv"
        df_details.to_csv(details_tsv, sep="\t", index=False)
        n_details = len(df_details)
        n_samples_with_outliers = df_details["SampleID"].nunique()
        print(f"[{ts()}] [INFO] Wrote {details_tsv} "
              f"({n_details} loci from {n_samples_with_outliers} samples)")
        for grp in ["Healthy", "ASD", "SZ"]:
            g_sub = df_details[df_details["Group"] == grp]
            print(f"  {grp}: {len(g_sub)} outlier loci from {g_sub['SampleID'].nunique()} samples")
    else:
        print(f"[{ts()}] [WARN] No rare outlier details to write")

    # --- Output Results ---
    df_out = df_keep.copy()
    df_out["rare_outlier_count"] = df_out["SampleID"].map(
        lambda x: int(final_rare_counts.get(x, 0))).astype(int)
    df_out["rare_any"] = (df_out["rare_outlier_count"] > 0).astype(int)
    df_out["rare_outlier_rate"] = df_out["rare_outlier_count"] / df_out["tested_loci_total"].astype(float)

    per_sample_tsv = outdir / f"{args.prefix}.per_sample.tsv"
    df_out.to_csv(per_sample_tsv, sep="\t", index=False)
    print(f"[{ts()}] [INFO] Wrote {per_sample_tsv} (n={len(df_out)})")

    g = df_out.groupby("Group", dropna=False)
    gs = pd.DataFrame({
        "N": g.size(),
        "mean_rare_outlier_count": g["rare_outlier_count"].mean(),
        "median_rare_outlier_count": g["rare_outlier_count"].median(),
        "mean_rare_outlier_rate": g["rare_outlier_rate"].mean(),
        "prop_rare_any": g["rare_any"].mean(),
    }).reset_index()
    group_tsv = outdir / f"{args.prefix}.group_summary.long.tsv"
    gs.to_csv(group_tsv, sep="\t", index=False, float_format="%.6g")
    print(f"[{ts()}] [INFO] Group summary:")
    for _, r in gs.iterrows():
        print(f"  {r['Group']}: N={r['N']}, prop_rare_any={r['prop_rare_any']:.3f}, "
              f"mean_count={r['mean_rare_outlier_count']:.4f}")

    thr = []
    h = df_out[df_out["Group"] == "Healthy"]["rare_outlier_count"].astype(float)
    if len(h) > 0:
        for q in [0.95, 0.99]:
            thr.append((q, float(h.quantile(q))))
    thr_df = pd.DataFrame(thr, columns=["quantile", "threshold"])
    thr_tsv = outdir / f"{args.prefix}.thresholds.tsv"
    thr_df.to_csv(thr_tsv, sep="\t", index=False, float_format="%.6g")

    # --- Statistical Tests ---
    print(f"\n[{ts()}] [INFO] === Statistical Testing ===")

    df_m = df_out.copy()
    exposure_col = "tested_loci_total"

    res_df = run_per_comparison_tests(df_m, case_grps=["ASD", "SZ"],
                                       exposure_col=exposure_col)

    stats_tsv = outdir / f"{args.prefix}.burden_stats_results.tsv"
    res_df.to_csv(stats_tsv, sep="\t", index=False)
    print(f"[{ts()}] [INFO] Wrote {stats_tsv}")

    print(f"\n=== Primary: Logistic Regression (rare_any -> OR) ===")
    for _, r in res_df.iterrows():
        or_val = r.get("Logit_OR", np.nan)
        ci_lo = r.get("Logit_OR_CI_low", np.nan)
        ci_hi = r.get("Logit_OR_CI_high", np.nan)
        p_val = r.get("Logit_P", np.nan)
        print(f"  {r['Comparison']}: OR={or_val:.2f} "
              f"(95%CI {ci_lo:.2f}-{ci_hi:.2f}), P={p_val:.2e}")

    print(f"\n=== Secondary: Poisson GLM (count -> RR, with offset) ===")
    for _, r in res_df.iterrows():
        rr_val = r.get("Poisson_RR", np.nan)
        ci_lo = r.get("Poisson_RR_CI_low", np.nan)
        ci_hi = r.get("Poisson_RR_CI_high", np.nan)
        p_val = r.get("Poisson_P", np.nan)
        print(f"  {r['Comparison']}: RR={rr_val:.2f} "
              f"(95%CI {ci_lo:.2f}-{ci_hi:.2f}), P={p_val:.2e}")

    print(f"\n=== Mann-Whitney U Test ===")
    for _, r in res_df.iterrows():
        print(f"  {r['Comparison']}: P={r['MWU_P']:.2e}")

    model_txt = outdir / f"{args.prefix}.model_summary.txt"
    with model_txt.open("w") as w:
        w.write(f"[{ts()}] Model fitting (STRling v9)\n")
        w.write(f"merge_dist: {args.merge_dist} bp\n")
        w.write(f"zthr: {args.zthr}, rare_cut: {args.rare_cut}, min_sum_str_counts: {args.min_sum_str_counts}\n")
        w.write(f"tested_loci_total: {L}\n\n")
        w.write(f"N used (complete-case): {len(df_m)}\n")
        w.write(f"Group counts:\n{df_m['Group'].value_counts().to_string()}\n\n")
        if all_outlier_details:
            w.write(f"Outlier details: {len(all_outlier_details)} loci from "
                    f"{len(set(d['SampleID'] for d in all_outlier_details))} samples\n\n")
        w.write("=" * 60 + "\n")
        w.write("Primary: Logistic Regression on rare_any -> OR\n")
        w.write("=" * 60 + "\n")
        for _, r in res_df.iterrows():
            w.write(f"\n  {r['Comparison']}:\n")
            w.write(f"    N_Case={r['N_Case']}, N_Ctrl={r['N_Ctrl']}\n")
            w.write(f"    rare_any: Case={r['rare_any_Case']:.3f}, Ctrl={r['rare_any_Ctrl']:.3f}\n")
            w.write(f"    OR={r.get('Logit_OR', np.nan):.4f} "
                    f"(95%CI {r.get('Logit_OR_CI_low', np.nan):.4f}-"
                    f"{r.get('Logit_OR_CI_high', np.nan):.4f})\n")
            w.write(f"    P={r.get('Logit_P', np.nan):.6e}\n")
        w.write("\n" + "=" * 60 + "\n")
        w.write("Secondary: Poisson GLM on rare_outlier_count (offset) -> RR\n")
        w.write("=" * 60 + "\n")
        for _, r in res_df.iterrows():
            w.write(f"\n  {r['Comparison']}:\n")
            w.write(f"    mean_count: Case={r['mean_count_Case']:.4f}, Ctrl={r['mean_count_Ctrl']:.4f}\n")
            w.write(f"    RR={r.get('Poisson_RR', np.nan):.4f} "
                    f"(95%CI {r.get('Poisson_RR_CI_low', np.nan):.4f}-"
                    f"{r.get('Poisson_RR_CI_high', np.nan):.4f})\n")
            w.write(f"    P={r.get('Poisson_P', np.nan):.6e}\n")
        w.write("\n" + "=" * 60 + "\n")
        w.write("Mann-Whitney U Test\n")
        w.write("=" * 60 + "\n")
        for _, r in res_df.iterrows():
            w.write(f"  {r['Comparison']}: P={r['MWU_P']:.6e}\n")

    print(f"[{ts()}] [INFO] Wrote {model_txt}")

    elapsed = time.time() - t0
    print(f"\n[{ts()}] [DONE] Total elapsed time: {elapsed:.1f}s ({elapsed/60:.1f}min)")


if __name__ == "__main__":
    main()
