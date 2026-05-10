#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ファイル名: 33_build_heffel_boundary_master_v9.py
# - 処理内容:
#   - BrainDev_raw.boundary.h5ad から raw support=2（高信頼 boundary）bin を抽出する
#   - BrainDev_impute.boundary.h5ad を重ねて、impute support 値と consensus support=2 を付与する
#   - L2_diff_domain_boundaries/*.bed.gz を読み込み、developmental boundary 変化を group_key × bin 単位で注釈する
#   - Heffel 2024 Supplementary Table 3 に基づく L3/L4 → L2 lineage 集約マッピングを使い、
#     h5ad obs（fine cluster 名）を diffbound 列（L2 pseudobulk 名）に一致させる
#   - v7 からの主な変更点:
#     (1) diff annotation の merge を 2-pass から 3-pass に拡張
#         Pass 1: exact match（group_key == group_key_from_diff）
#         Pass 2: normalized match（eCGE/eMGE → CGE/MGE）
#         Pass 3: L2 aggregated match（Supp3 由来の L3→L2 辞書で fine cluster を L2 に集約）
#     (2) obs_to_l2_group_key() ヘルパを追加
#     (3) master に group_key_L2_from_obs, matching_level 列を追加（pass1/pass2/pass3/none 診断）
#     (4) L2 集約マッピングは Supp Table 3 由来の訂正を反映:
#         - Exc-L4-5-FOXP2, Exc-L4-5-TOX, Exc-L4-6-LRRK1, Exc-L5-6-PDZRN4 → Exc-UL（v2 では DL と誤分類）
#         - Exc-Subiculum-FN1 → Exc-CA（v2 では lineage から除外）
#         - Exc-NP-TSHZ2 → Exc-DL
#     (5) 出力 suffix を v8 に統一
#   - v9 からの主な変更点:
#     (1) tad04212026/ パイプラインに移行
#     (2) パスを common/paths_v1.py で一元管理（Path.cwd() 依存を廃止）
#     (3) 出力 suffix を v9 に統一、出力先を OUT_01_BOUNDARY_MASTER に固定
#   - excitatory 用 flags（RG-anchor gain/loss）と inhibitory 用 flags（CGE/MGE × gain/loss）を付与する
#   - region（PFC/HPC）や lineage 情報も列として保持する
#   - 後続の SV overlap 解析に使う heffel_boundary_master_v9.tsv.gz と summary を出力する
#   - 実行時間を記録する

import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

# ----- common paths -----
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.paths_v1 import (
    HEFFEL_DOMAIN_BOUNDARIES,
    HEFFEL_L2_DIFF_DOMAIN_BOUNDARIES,
    OUT_01_BOUNDARY_MASTER,
    ensure_output_dirs,
)


def now_str() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    print(f"[{now_str()}] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Supp Table 3 由来の L3/L4 → L2 集約マッピング
# - Heffel 2024 Nature (DOI 10.1038/s41586-024-08030-7) Supplementary Table 3
# - 全 L3/L4 細胞型を L2 pseudobulk（diffbound 列名）に集約
# - v2 で誤っていた以下の分類を Supp3 L2 列に基づき訂正:
#     * Exc-L4-5-FOXP2, Exc-L4-5-TOX, Exc-L4-6-LRRK1, Exc-L5-6-PDZRN4 → Exc-UL
#     * Exc-Subiculum-FN1 → Exc-CA
#     * Exc-NP-TSHZ2 → Exc-DL (変更なし)
# - RG / Non-Neuron は diffbound file が存在しないため identity mapping
# ---------------------------------------------------------------------------
L3_TO_L2: Dict[str, str] = {
    # --- Astro (HPC/PFC 両方に diffbound file が存在) ---
    "Astro": "Astro",

    # --- Exc-CA (HPC) ---
    "Exc-CA":                "Exc-CA",  # L2-pseudobulk identity (2T stage)
    "Exc-CA1":               "Exc-CA",
    "Exc-CA1-ADCY1":         "Exc-CA",
    "Exc-CA1-RELN":          "Exc-CA",
    "Exc-CA3":               "Exc-CA",
    "Exc-CA-MossyCell":      "Exc-CA",
    "Exc-CA-NTNG2":          "Exc-CA",
    "Exc-CA-TSHZ2-PCP4":     "Exc-CA",
    "Exc-CA-TSHZ2-MEIS2":    "Exc-CA",
    "Exc-Subiculum-FN1":     "Exc-CA",  # Supp3 L2=Exc-CA (v2 訂正)

    # --- Exc-DG (HPC) ---
    "Exc-DG": "Exc-DG",

    # --- Exc-ENT (HPC) ---
    "Exc-ENT":               "Exc-ENT",
    "Exc-ENT-GALNT17":       "Exc-ENT",
    "Exc-ENT-GLIS1":         "Exc-ENT",
    "Exc-ENT-GRIK4":         "Exc-ENT",
    "Exc-ENT-GRIK4-TSHZ2":   "Exc-ENT",
    "Exc-ENT-TSHZ2":         "Exc-ENT",

    # --- Exc-DL (PFC) ---
    "Exc-DL-1":              "Exc-DL",
    "Exc-DL-2":              "Exc-DL",
    "Exc-DL-ASTN1":          "Exc-DL",
    "Exc-DL-GRIK4":          "Exc-DL",
    "Exc-DL-GRIK4-TSHZ2":    "Exc-DL",
    "Exc-NP-TSHZ2":          "Exc-DL",  # Near-Projecting も Supp3 L2=Exc-DL

    # --- Exc-UL (PFC) ---
    "Exc-UL-1":              "Exc-UL",  # 2T stage pseudobulk
    "Exc-L1-3-CUX2":         "Exc-UL",
    "Exc-L1-3-NRXN2":        "Exc-UL",
    "Exc-L4-RORB":           "Exc-UL",
    "Exc-L4-PLCH1":          "Exc-UL",
    "Exc-L4-5-FOXP2":        "Exc-UL",  # Supp3 L2=Exc-UL (v2 訂正: DL → UL)
    "Exc-L4-5-TOX":          "Exc-UL",  # Supp3 L2=Exc-UL (v2 訂正: DL → UL)
    "Exc-L4-6-LRRK1":        "Exc-UL",  # Supp3 L2=Exc-UL (v2 訂正: DL → UL)
    "Exc-L5-6-PDZRN4":       "Exc-UL",  # Supp3 L2=Exc-UL (v2 訂正: DL → UL)

    # --- Inh-CGE ---
    "Inh-CGE":               "Inh-CGE",
    "Inh-eCGE":              "Inh-CGE",  # early CGE は Pass 2 (normalize) でも対応
    "Inh-CGE-ALK":           "Inh-CGE",
    "Inh-CGE-CHRNA2":        "Inh-CGE",
    "Inh-CGE-CLMP":          "Inh-CGE",
    "Inh-CGE-EPHA4":         "Inh-CGE",
    "Inh-CGE-FRAS1":         "Inh-CGE",
    "Inh-CGE-MN1":           "Inh-CGE",
    "Inh-CGE-OXR1":          "Inh-CGE",
    "Inh-CGE-SOX13":         "Inh-CGE",
    "Inh-CGE-TOX":           "Inh-CGE",

    # --- Inh-MGE ---
    "Inh-MGE":               "Inh-MGE",
    "Inh-eMGE":              "Inh-MGE",  # early MGE は Pass 2 (normalize) でも対応
    "Inh-MGE-ALCAM":         "Inh-MGE",
    "Inh-MGE-ChC-UNC5B":     "Inh-MGE",
    "Inh-MGE-CNTNAP2":       "Inh-MGE",
    "Inh-MGE-CNTNAP4":       "Inh-MGE",
    "Inh-MGE-ERBB4":         "Inh-MGE",
    "Inh-MGE-MAN1A1":        "Inh-MGE",
    "Inh-MGE-RBFOX1":        "Inh-MGE",

    # --- RG (diffbound file なし; Pass 1 で RG anchor column として matches する obs のみ対応) ---
    "RG-1": "RG-1",
    "RG-2": "RG-2",

    # --- Non-Neuron (diffbound file なし; 集約しても match しない) ---
    "MGC-1": "MGC-1",
    "MGC-2": "MGC-2",
    "PC":    "PC",
    "EC":    "EC",
    "ODC":   "ODC",
    "OPC":   "OPC",
    "VLMC":  "VLMC",
}


def obs_to_l2_group_key(group_key: str) -> str:
    """
    h5ad obs_name（fine cluster）を L2 pseudobulk group_key に集約する。
    obs = "<region>_<stage>_<rest>" 形式（例: HPC_adult_Exc-CA1）
    未知の rest はそのまま返す（identity fallback）
    """
    gk = str(group_key)
    parts = gk.split("_", 2)
    if len(parts) != 3:
        return gk
    region, stage, rest = parts
    l2 = L3_TO_L2.get(rest, rest)
    return f"{region}_{stage}_{l2}"


# ---------------------------------------------------------------------------
# key uniqueness assertion
# ---------------------------------------------------------------------------
def assert_unique_keys(df: pd.DataFrame, key_cols: List[str], label: str, n_show: int = 10) -> None:
    """
    key_cols の組み合わせで重複があれば ValueError を投げる。
    エラーメッセージには重複件数と先頭 n_show 件の例を含む。
    """
    dup_mask = df.duplicated(subset=key_cols, keep=False)
    n_dup = int(dup_mask.sum())
    if n_dup > 0:
        dup_examples = (
            df.loc[dup_mask, key_cols]
            .drop_duplicates()
            .head(n_show)
            .to_string(index=False)
        )
        raise ValueError(
            f"[{label}] key uniqueness violation: "
            f"{n_dup} duplicated rows on {key_cols}.\n"
            f"First {min(n_show, n_dup)} duplicate key examples:\n{dup_examples}"
        )


def normalize_group_key(group_key: str) -> str:
    """
    group_key の表記ゆれを軽く正規化する。
    inhibitory early subtype の eCGE/eMGE を CGE/MGE に寄せる。
    ※ diffbound 側との対応付け専用。raw/impute の内部整合性には使わない。
    """
    x = str(group_key)
    x = x.replace("Inh-eCGE", "Inh-CGE")
    x = x.replace("Inh-eMGE", "Inh-MGE")
    return x


def stage_rank(colname: str) -> int:
    """
    developmental stage の順序を返す。
    小さいほど早期。
    """
    c = str(colname)
    if "2T" in c:
        return 1
    if "3T" in c:
        return 2
    if "infant" in c:
        return 3
    if "adult" in c:
        return 4
    return 99


def infer_trajectory_class_from_filename(basename: str) -> str:
    b = str(basename)
    if "Exc-" in b:
        return "neuronal_exc"
    if "Astro" in b:
        return "glial_astro"
    if "Inh-" in b:
        return "neuronal_inh"
    return "other"


def infer_lineage_key_from_filename(basename: str) -> str:
    b = str(basename)
    b = re.sub(r"_diffbound\.bed\.gz$", "", b)
    return b


def infer_region_from_lineage_key(lineage_key: str) -> str:
    lk = str(lineage_key)
    if lk.startswith("PFC_"):
        return "PFC"
    if lk.startswith("HPC_"):
        return "HPC"
    return "OTHER"


def extract_nonzero_sparse_long(
    adata, matrix_name: str, value_filter: Optional[float] = None
) -> pd.DataFrame:
    """
    sparse matrix の非ゼロ要素を long 形式に展開する。
    value_filter を指定した場合はその値のみ残す。
    """
    X = adata.X.tocoo() if sparse.issparse(adata.X) else sparse.coo_matrix(adata.X)
    df = pd.DataFrame({
        "row_idx": X.row,
        "col_idx": X.col,
        "value": X.data
    })
    if value_filter is not None:
        df = df.loc[df["value"] == value_filter].copy()
    obs_names = np.array(adata.obs_names)
    var_df = adata.var.reset_index(drop=True).copy()
    required_var_cols = ["chrom", "start", "end"]
    missing = [c for c in required_var_cols if c not in var_df.columns]
    if missing:
        raise ValueError(f"{matrix_name}: var に必要列がありません: {missing}")
    df["group_key"] = obs_names[df["row_idx"].to_numpy()]
    df["group_key_norm"] = df["group_key"].map(normalize_group_key)
    df["group_key_L2_from_obs"] = df["group_key"].map(obs_to_l2_group_key)
    df["chrom"] = var_df.loc[df["col_idx"], "chrom"].to_numpy()
    df["start0"] = var_df.loc[df["col_idx"], "start"].to_numpy()
    df["end"] = var_df.loc[df["col_idx"], "end"].to_numpy()
    df["bin_id"] = (
        df["chrom"].astype(str) + ":"
        + df["start0"].astype(str) + "-"
        + df["end"].astype(str)
    )
    df["matrix_name"] = matrix_name
    return df[
        ["matrix_name", "group_key", "group_key_norm", "group_key_L2_from_obs",
         "chrom", "start0", "end", "bin_id", "value"]
    ].reset_index(drop=True)


def read_diffbound_file(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="gzip")
    expected = {"chrom", "start", "end", "chi2_sc"}
    missing = expected - set(df.columns)
    if missing:
        raise ValueError(f"{path}: 必須列がありません: {missing}")
    return df


def summarize_diffbound_row_structure(df: pd.DataFrame, basename: str) -> Tuple[List[str], List[str], str]:
    """
    value columns を early / late に分ける。
    """
    meta_cols = {"chrom", "start", "end", "chi2_sc"}
    value_cols = [c for c in df.columns if c not in meta_cols]
    rg_cols = [c for c in value_cols if "RG" in c]
    if len(rg_cols) > 0:
        early_cols = sorted(rg_cols, key=stage_rank)
        late_cols = sorted([c for c in value_cols if c not in rg_cols], key=stage_rank)
        anchor_type = "RG_anchor"
    else:
        ordered = sorted(value_cols, key=stage_rank)
        if len(ordered) >= 2:
            early_cols = [ordered[0]]
            late_cols = ordered[1:]
        else:
            early_cols = ordered[:]
            late_cols = ordered[:]
        anchor_type = "temporal_noRG"
    if len(early_cols) == 0 or len(late_cols) == 0:
        ordered = sorted(value_cols, key=stage_rank)
        early_cols = ordered[:1]
        late_cols = ordered[:]
        anchor_type = "fallback_single"
    return early_cols, late_cols, anchor_type


def build_diffbound_long(path: Path) -> pd.DataFrame:
    """
    1つの diffbound file を long 形式に変換する。
    """
    basename = path.name
    lineage_key = infer_lineage_key_from_filename(basename)
    trajectory_class = infer_trajectory_class_from_filename(basename)
    region = infer_region_from_lineage_key(lineage_key)
    df = read_diffbound_file(path)
    early_cols, late_cols, anchor_type = summarize_diffbound_row_structure(df, basename)
    df["early_mean_prob"] = df[early_cols].mean(axis=1)
    df["late_mean_prob"] = df[late_cols].mean(axis=1)
    df["late_max_prob"] = df[late_cols].max(axis=1)
    df["dev_delta_prob"] = df["late_mean_prob"] - df["early_mean_prob"]

    def direction(x: float) -> str:
        if pd.isna(x):
            return "NA"
        if x > 0:
            return "gain"
        if x < 0:
            return "loss"
        return "flat"

    df["dev_direction"] = df["dev_delta_prob"].map(direction)
    meta_cols = [
        "chrom", "start", "end", "chi2_sc",
        "early_mean_prob", "late_mean_prob", "late_max_prob",
        "dev_delta_prob", "dev_direction"
    ]
    exclude_cols = {
        "early_mean_prob", "late_mean_prob", "late_max_prob",
        "dev_delta_prob", "dev_direction"
    }
    value_cols = [c for c in df.columns if c not in set(meta_cols)]
    group_cols = [c for c in value_cols if c not in exclude_cols and c != "chi2_sc"]
    long_list = []
    for group_col in group_cols:
        sub = df[meta_cols].copy()
        sub["group_key_from_diff"] = group_col
        sub["group_key_norm"] = sub["group_key_from_diff"].map(normalize_group_key)
        sub["group_prob_in_diff"] = df[group_col].to_numpy()
        sub["lineage_key"] = lineage_key
        sub["trajectory_class"] = trajectory_class
        sub["anchor_type"] = anchor_type
        sub["region"] = region
        sub["diff_file"] = basename
        sub["bin_id"] = (
            sub["chrom"].astype(str) + ":"
            + sub["start"].astype(str) + "-"
            + sub["end"].astype(str)
        )
        sub = sub.rename(columns={"start": "start0"})
        long_list.append(sub)
    out = pd.concat(long_list, axis=0, ignore_index=True)
    return out[
        [
            "group_key_from_diff", "group_key_norm", "chrom", "start0", "end", "bin_id",
            "group_prob_in_diff", "chi2_sc",
            "early_mean_prob", "late_mean_prob", "late_max_prob",
            "dev_delta_prob", "dev_direction",
            "trajectory_class", "anchor_type", "region", "lineage_key", "diff_file"
        ]
    ]


def deduplicate_diff_annotations(df: pd.DataFrame, key_cols: List[str]) -> pd.DataFrame:
    """
    key_cols で重複があれば、|dev_delta_prob| の大きいものを残す。
    """
    if df.empty:
        return df.copy()
    tmp = df.copy()
    tmp["abs_dev_delta_prob"] = tmp["dev_delta_prob"].abs()
    sort_cols = key_cols + ["abs_dev_delta_prob", "chi2_sc"]
    sort_asc = [True] * len(key_cols) + [False, False]
    tmp = tmp.sort_values(sort_cols, ascending=sort_asc)
    tmp = tmp.drop_duplicates(subset=key_cols, keep="first").copy()
    tmp = tmp.drop(columns=["abs_dev_delta_prob"])
    return tmp


# ---------------------------------------------------------------------------
# v8: three-pass merge for diff annotation
# ---------------------------------------------------------------------------
def merge_diff_annotation_three_pass(
    master: pd.DataFrame, diff_ann: pd.DataFrame
) -> pd.DataFrame:
    """
    diff annotation を master に merge する。
    Pass 1: exact match (master.group_key == diff_ann.group_key_from_diff AND bin_id)
    Pass 2: normalized match (master.group_key_norm == diff_ann.group_key_norm AND bin_id)
            ※ Pass 1 で matched 済みの master 行は除外
    Pass 3: L2 aggregated match (master.group_key_L2_from_obs == diff_ann.group_key_from_diff AND bin_id)
            ※ Pass 1/2 で matched 済みの master 行は除外
            ※ Supp Table 3 由来の L3→L2 lineage mapping を使った集約マッチ
    最後に結合して返し、matching_level 列を付与する。
    """
    diff_cols = [
        "group_key_from_diff", "group_prob_in_diff", "chi2_sc",
        "early_mean_prob", "late_mean_prob", "late_max_prob",
        "dev_delta_prob", "dev_direction",
        "trajectory_class", "anchor_type", "region", "lineage_key", "diff_file",
    ]

    # --- diff_ann を Pass 3 向けに group_key_from_diff × bin_id でも重複除去 ---
    # (Pass 1/3 は exact なので group_key_from_diff × bin_id でキー化)
    diff_ann_exact = deduplicate_diff_annotations(
        diff_ann, ["group_key_from_diff", "bin_id"]
    )

    # --- Pass 1: exact match ---
    pass1 = master[["group_key", "bin_id"]].merge(
        diff_ann_exact[["group_key_from_diff", "bin_id"] + [c for c in diff_cols if c != "group_key_from_diff"]],
        left_on=["group_key", "bin_id"],
        right_on=["group_key_from_diff", "bin_id"],
        how="inner"
    )
    pass1["matching_level"] = "exact"
    pass1_keys = set(zip(pass1["group_key"], pass1["bin_id"]))
    log(f"  diff merge Pass 1 (exact):        {len(pass1_keys)} master rows matched")

    # --- Pass 2: normalized match (unmatched rows only) ---
    master_keys = pd.Series(list(zip(master["group_key"], master["bin_id"])))
    master_unmatched_mask = ~master_keys.isin(pass1_keys).values
    master_unmatched = master.loc[
        master_unmatched_mask,
        ["group_key", "group_key_norm", "group_key_L2_from_obs", "bin_id"]
    ].copy()

    pass2 = master_unmatched.merge(
        diff_ann[["group_key_norm", "bin_id"] + diff_cols],
        on=["group_key_norm", "bin_id"],
        how="inner"
    )
    pass2["matching_level"] = "normalized"
    pass2_keys = set(zip(pass2["group_key"], pass2["bin_id"]))
    log(f"  diff merge Pass 2 (normalized):   {len(pass2_keys)} additional master rows matched")

    # --- Pass 3: L2 aggregated match (Pass 1/2 unmatched rows only) ---
    master_unmatched_2 = master_unmatched.loc[
        ~pd.Series(list(zip(master_unmatched["group_key"], master_unmatched["bin_id"]))).isin(pass2_keys).values,
        :
    ].copy()

    pass3 = master_unmatched_2.merge(
        diff_ann_exact[["group_key_from_diff", "bin_id"] + [c for c in diff_cols if c != "group_key_from_diff"]],
        left_on=["group_key_L2_from_obs", "bin_id"],
        right_on=["group_key_from_diff", "bin_id"],
        how="inner"
    )
    pass3["matching_level"] = "l2_aggregated"
    pass3_keys = set(zip(pass3["group_key"], pass3["bin_id"]))
    log(f"  diff merge Pass 3 (L2-aggregated): {len(pass3_keys)} additional master rows matched")

    # --- 結合 ---
    matched = pd.concat([pass1, pass2, pass3], axis=0, ignore_index=True)

    # many-to-many 検出: master 側の (group_key, bin_id) に重複がないことを確認
    assert_unique_keys(matched, ["group_key", "bin_id"], "diff_merge_matched")

    # master に left join
    join_cols = ["group_key", "bin_id"] + diff_cols + ["matching_level"]
    matched_for_join = matched[join_cols].copy()
    result = master.merge(
        matched_for_join,
        on=["group_key", "bin_id"],
        how="left",
        validate="1:1"
    )
    result["matching_level"] = result["matching_level"].fillna("none")
    n_matched = result["group_prob_in_diff"].notna().sum()
    n_total = result.shape[0]
    log(f"  diff merge total: {n_matched}/{n_total} master rows annotated "
        f"({n_matched/n_total*100:.1f}%)")

    # --- 診断サマリ: matching_level 別 / group_key 別 ---
    level_counts = result["matching_level"].value_counts()
    log(f"  matching_level distribution:")
    for level in ["exact", "normalized", "l2_aggregated", "none"]:
        cnt = int(level_counts.get(level, 0))
        log(f"    {level:<15s}: {cnt} rows")

    n_groups_total = result["group_key"].nunique()
    n_groups_any_match = result.loc[result["matching_level"] != "none", "group_key"].nunique()
    log(f"  unique group_keys with any diff annotation: {n_groups_any_match}/{n_groups_total}")

    return result


def main():
    t0 = time.time()
    log("Start 33_build_heffel_boundary_master_v9.py")

    ensure_output_dirs()

    raw_h5ad    = HEFFEL_DOMAIN_BOUNDARIES / "BrainDev_raw.boundary.h5ad"
    impute_h5ad = HEFFEL_DOMAIN_BOUNDARIES / "BrainDev_impute.boundary.h5ad"
    diff_dir    = HEFFEL_L2_DIFF_DOMAIN_BOUNDARIES
    for p in [raw_h5ad, impute_h5ad, diff_dir]:
        if not p.exists():
            raise FileNotFoundError(f"Not found: {p}")

    outdir = OUT_01_BOUNDARY_MASTER

    # 1) load h5ad
    log(f"Loading raw h5ad: {raw_h5ad}")
    adata_raw = ad.read_h5ad(raw_h5ad)
    log(f"Loading impute h5ad: {impute_h5ad}")
    adata_impute = ad.read_h5ad(impute_h5ad)

    # 2) raw support=2 を anchor にする
    log("Extracting raw support=2 bins")
    raw_support2 = extract_nonzero_sparse_long(adata_raw, matrix_name="raw", value_filter=2.0)
    raw_support2 = raw_support2.rename(columns={"value": "raw_value"})
    raw_support2["is_raw_support2"] = 1

    # uniqueness check on [group_key, bin_id]
    assert_unique_keys(raw_support2, ["group_key", "bin_id"], "raw_support2")
    log(f"raw_support2: {raw_support2.shape[0]} rows, unique keys [group_key, bin_id] verified")

    # normalize / L2 集約後キーの collision 診断
    n_norm_dup = raw_support2.duplicated(subset=["group_key_norm", "bin_id"], keep=False).sum()
    if n_norm_dup > 0:
        log(f"  INFO: {n_norm_dup} rows have duplicated [group_key_norm, bin_id] "
            "(expected when eCGE/eMGE and CGE/MGE coexist in raw h5ad)")
    n_l2_dup = raw_support2.duplicated(subset=["group_key_L2_from_obs", "bin_id"], keep=False).sum()
    if n_l2_dup > 0:
        log(f"  INFO: {n_l2_dup} rows have duplicated [group_key_L2_from_obs, bin_id] "
            "(expected when multiple fine clusters share the same L2 pseudobulk)")

    # group_key → L2 mapping の診断（ユニーク obs name 単位）
    unique_gk_to_l2 = (
        raw_support2[["group_key", "group_key_L2_from_obs"]].drop_duplicates()
        .sort_values(["group_key_L2_from_obs", "group_key"])
    )
    n_unique_obs = unique_gk_to_l2["group_key"].nunique()
    n_unique_l2 = unique_gk_to_l2["group_key_L2_from_obs"].nunique()
    n_identity = int((unique_gk_to_l2["group_key"] == unique_gk_to_l2["group_key_L2_from_obs"]).sum())
    log(f"  L3→L2 mapping: {n_unique_obs} unique obs → {n_unique_l2} unique L2 group_keys "
        f"(identity={n_identity}, aggregated={n_unique_obs - n_identity})")

    # 3) impute 非ゼロを join
    log("Extracting impute nonzero bins")
    impute_nonzero = extract_nonzero_sparse_long(adata_impute, matrix_name="impute", value_filter=None)
    impute_nonzero = impute_nonzero.rename(columns={"value": "impute_value"})
    impute_nonzero = impute_nonzero[
        ["group_key", "bin_id", "impute_value"]
    ].copy()

    assert_unique_keys(impute_nonzero, ["group_key", "bin_id"], "impute_nonzero")
    log(f"impute_nonzero: {impute_nonzero.shape[0]} rows, unique keys [group_key, bin_id] verified")

    master = raw_support2.merge(
        impute_nonzero,
        on=["group_key", "bin_id"],
        how="left",
        validate="1:1"
    )
    master["impute_value"] = master["impute_value"].fillna(0.0)
    master["is_impute_support2"] = (master["impute_value"] == 2.0).astype(int)
    master["is_consensus_support2"] = (
        (master["raw_value"] == 2.0) & (master["impute_value"] == 2.0)
    ).astype(int)
    log(f"master after raw×impute merge: {master.shape[0]} rows")

    # 4) diffbound developmental annotation
    log("Parsing diffbound files")
    diff_files = sorted(diff_dir.glob("*_diffbound.bed.gz"))
    if len(diff_files) == 0:
        raise FileNotFoundError(f"No diffbound files found in: {diff_dir}")

    diff_long_list = []
    for f in diff_files:
        log(f"Reading diffbound: {f.name}")
        diff_long = build_diffbound_long(f)
        diff_long_list.append(diff_long)
    diff_long_all = pd.concat(diff_long_list, axis=0, ignore_index=True)

    # group_key_norm × bin_id で dedup (Pass 2 用)
    diff_ann = deduplicate_diff_annotations(diff_long_all, ["group_key_norm", "bin_id"])

    diff_ann = diff_ann[
        [
            "group_key_from_diff", "group_key_norm", "bin_id",
            "group_prob_in_diff", "chi2_sc",
            "early_mean_prob", "late_mean_prob", "late_max_prob",
            "dev_delta_prob", "dev_direction",
            "trajectory_class", "anchor_type", "region", "lineage_key", "diff_file",
        ]
    ].copy()

    assert_unique_keys(diff_ann, ["group_key_norm", "bin_id"], "diff_ann")
    log(f"diff_ann: {diff_ann.shape[0]} rows, unique keys verified")

    # 5) v8: three-pass merge
    log("Merging developmental annotations (three-pass)")
    master = merge_diff_annotation_three_pass(master, diff_ann)
    master["overlaps_diffbound"] = master["group_prob_in_diff"].notna().astype(int)

    # final uniqueness check on master
    assert_unique_keys(master, ["group_key", "bin_id"], "master_final")
    log(f"master final: {master.shape[0]} rows, unique keys [group_key, bin_id] verified")

    # 6) generic flags
    master["is_rg_anchored"] = (master["anchor_type"] == "RG_anchor").astype(int)
    master["is_temporal_norg"] = (master["anchor_type"] == "temporal_noRG").astype(int)
    master["is_dev_gain"] = (master["dev_direction"] == "gain").astype(int)
    master["is_dev_loss"] = (master["dev_direction"] == "loss").astype(int)
    master["is_dev_exc"] = (master["trajectory_class"] == "neuronal_exc").astype(int)
    master["is_dev_astro"] = (master["trajectory_class"] == "glial_astro").astype(int)
    master["is_dev_inh"] = (master["trajectory_class"] == "neuronal_inh").astype(int)
    master["is_region_pfc"] = (master["region"] == "PFC").astype(int)
    master["is_region_hpc"] = (master["region"] == "HPC").astype(int)
    master["is_lineage_inh_cge"] = (
        master["lineage_key"].astype(str).str.contains("Inh-CGE", regex=False, na=False)
    ).astype(int)
    master["is_lineage_inh_mge"] = (
        master["lineage_key"].astype(str).str.contains("Inh-MGE", regex=False, na=False)
    ).astype(int)

    # 7) excitatory flags
    master["is_exc_rg_gain"] = (
        (master["overlaps_diffbound"] == 1)
        & (master["is_dev_exc"] == 1)
        & (master["is_rg_anchored"] == 1)
        & (master["is_dev_gain"] == 1)
    ).astype(int)
    master["is_exc_rg_loss"] = (
        (master["overlaps_diffbound"] == 1)
        & (master["is_dev_exc"] == 1)
        & (master["is_rg_anchored"] == 1)
        & (master["is_dev_loss"] == 1)
    ).astype(int)
    master["is_astro_rg_loss"] = (
        (master["overlaps_diffbound"] == 1)
        & (master["is_dev_astro"] == 1)
        & (master["is_rg_anchored"] == 1)
        & (master["is_dev_loss"] == 1)
    ).astype(int)

    # 8) inhibitory flags
    master["is_inh_cge_gain"] = (
        (master["overlaps_diffbound"] == 1)
        & (master["is_dev_inh"] == 1)
        & (master["is_lineage_inh_cge"] == 1)
        & (master["is_dev_gain"] == 1)
    ).astype(int)
    master["is_inh_cge_loss"] = (
        (master["overlaps_diffbound"] == 1)
        & (master["is_dev_inh"] == 1)
        & (master["is_lineage_inh_cge"] == 1)
        & (master["is_dev_loss"] == 1)
    ).astype(int)
    master["is_inh_mge_gain"] = (
        (master["overlaps_diffbound"] == 1)
        & (master["is_dev_inh"] == 1)
        & (master["is_lineage_inh_mge"] == 1)
        & (master["is_dev_gain"] == 1)
    ).astype(int)
    master["is_inh_mge_loss"] = (
        (master["overlaps_diffbound"] == 1)
        & (master["is_dev_inh"] == 1)
        & (master["is_lineage_inh_mge"] == 1)
        & (master["is_dev_loss"] == 1)
    ).astype(int)

    # region-specific inhibitory exploratory flags
    master["is_pfc_inh_cge_gain"] = (
        (master["is_inh_cge_gain"] == 1) & (master["is_region_pfc"] == 1)
    ).astype(int)
    master["is_hpc_inh_cge_gain"] = (
        (master["is_inh_cge_gain"] == 1) & (master["is_region_hpc"] == 1)
    ).astype(int)
    master["is_pfc_inh_mge_gain"] = (
        (master["is_inh_mge_gain"] == 1) & (master["is_region_pfc"] == 1)
    ).astype(int)
    master["is_hpc_inh_mge_gain"] = (
        (master["is_inh_mge_gain"] == 1) & (master["is_region_hpc"] == 1)
    ).astype(int)

    # 9) 必須列チェック
    ordered_cols = [
        "group_key", "group_key_norm", "group_key_L2_from_obs",
        "chrom", "start0", "end", "bin_id",
        "raw_value", "is_raw_support2", "impute_value", "is_impute_support2",
        "is_consensus_support2",
        "overlaps_diffbound", "matching_level",
        "group_key_from_diff", "group_prob_in_diff", "chi2_sc",
        "early_mean_prob", "late_mean_prob", "late_max_prob",
        "dev_delta_prob", "dev_direction",
        "trajectory_class", "anchor_type", "region", "lineage_key", "diff_file",
        "is_rg_anchored", "is_temporal_norg",
        "is_dev_gain", "is_dev_loss",
        "is_dev_exc", "is_dev_astro", "is_dev_inh",
        "is_region_pfc", "is_region_hpc",
        "is_lineage_inh_cge", "is_lineage_inh_mge",
        "is_exc_rg_gain", "is_exc_rg_loss", "is_astro_rg_loss",
        "is_inh_cge_gain", "is_inh_cge_loss",
        "is_inh_mge_gain", "is_inh_mge_loss",
        "is_pfc_inh_cge_gain", "is_hpc_inh_cge_gain",
        "is_pfc_inh_mge_gain", "is_hpc_inh_mge_gain",
    ]
    missing_cols = [c for c in ordered_cols if c not in master.columns]
    if missing_cols:
        raise KeyError(
            "master に必要列がありません: " + ", ".join(missing_cols)
            + "\n現在の列: " + ", ".join(master.columns.astype(str).tolist())
        )
    master = master[ordered_cols].copy()

    # 10) summary
    log("Building summary tables")
    summary_overall = pd.DataFrame([{
        "n_raw_support2_total": int(master.shape[0]),
        "n_consensus_support2": int(master["is_consensus_support2"].sum()),
        "n_overlaps_diffbound": int(master["overlaps_diffbound"].sum()),
        "n_exc_rg_gain": int(master["is_exc_rg_gain"].sum()),
        "n_astro_rg_loss": int(master["is_astro_rg_loss"].sum()),
        "n_inh_cge_gain": int(master["is_inh_cge_gain"].sum()),
        "n_inh_mge_gain": int(master["is_inh_mge_gain"].sum()),
        "n_exc_rg_loss": int(master["is_exc_rg_loss"].sum()),
        "n_inh_cge_loss": int(master["is_inh_cge_loss"].sum()),
        "n_inh_mge_loss": int(master["is_inh_mge_loss"].sum()),
        "n_match_exact": int((master["matching_level"] == "exact").sum()),
        "n_match_normalized": int((master["matching_level"] == "normalized").sum()),
        "n_match_l2_aggregated": int((master["matching_level"] == "l2_aggregated").sum()),
        "n_match_none": int((master["matching_level"] == "none").sum()),
    }])

    group_summary = (
        master.groupby(["group_key", "group_key_L2_from_obs"], as_index=False)
        .agg(
            n_raw_support2=("bin_id", "size"),
            n_consensus_support2=("is_consensus_support2", "sum"),
            n_overlaps_diffbound=("overlaps_diffbound", "sum"),
            n_exc_rg_gain=("is_exc_rg_gain", "sum"),
            n_astro_rg_loss=("is_astro_rg_loss", "sum"),
            n_inh_cge_gain=("is_inh_cge_gain", "sum"),
            n_inh_mge_gain=("is_inh_mge_gain", "sum"),
        )
        .sort_values(
            ["n_exc_rg_gain", "n_inh_cge_gain", "n_inh_mge_gain", "n_raw_support2"],
            ascending=[False, False, False, False]
        )
    )

    # matching_level × group_key の診断表（v8 新規）
    matching_summary = (
        master.groupby(["group_key", "group_key_L2_from_obs", "matching_level"], as_index=False)
        .size()
        .rename(columns={"size": "n_rows"})
        .sort_values(["group_key", "matching_level"])
    )

    lineage_summary = (
        master.groupby(["trajectory_class", "lineage_key", "anchor_type", "dev_direction"], as_index=False)
        .agg(
            n_bins=("bin_id", "size"),
            n_consensus_support2=("is_consensus_support2", "sum"),
            mean_group_prob=("group_prob_in_diff", "mean"),
            median_group_prob=("group_prob_in_diff", "median"),
            mean_dev_delta=("dev_delta_prob", "mean"),
            median_dev_delta=("dev_delta_prob", "median"),
        )
        .sort_values(["trajectory_class", "lineage_key", "anchor_type", "dev_direction"])
    )

    inhibitory_summary = (
        master.groupby(["region", "lineage_key"], as_index=False)
        .agg(
            n_inh_cge_gain=("is_inh_cge_gain", "sum"),
            n_inh_mge_gain=("is_inh_mge_gain", "sum"),
            n_inh_cge_loss=("is_inh_cge_loss", "sum"),
            n_inh_mge_loss=("is_inh_mge_loss", "sum"),
        )
        .sort_values(["region", "lineage_key"])
    )

    # 11) save
    master_path = outdir / "heffel_boundary_master_v9.tsv.gz"
    summary_path = outdir / "heffel_boundary_master_summary_v9.tsv"
    group_summary_path = outdir / "heffel_boundary_master_group_summary_v9.tsv"
    matching_summary_path = outdir / "heffel_boundary_master_matching_summary_v9.tsv"
    lineage_summary_path = outdir / "heffel_boundary_master_lineage_summary_v9.tsv"
    inhibitory_summary_path = outdir / "heffel_boundary_master_inhibitory_summary_v9.tsv"

    log(f"Writing master: {master_path}")
    master.to_csv(master_path, sep="\t", index=False, compression="gzip")
    log(f"Writing summary: {summary_path}")
    summary_overall.to_csv(summary_path, sep="\t", index=False)
    log(f"Writing group summary: {group_summary_path}")
    group_summary.to_csv(group_summary_path, sep="\t", index=False)
    log(f"Writing matching summary: {matching_summary_path}")
    matching_summary.to_csv(matching_summary_path, sep="\t", index=False)
    log(f"Writing lineage summary: {lineage_summary_path}")
    lineage_summary.to_csv(lineage_summary_path, sep="\t", index=False)
    log(f"Writing inhibitory summary: {inhibitory_summary_path}")
    inhibitory_summary.to_csv(inhibitory_summary_path, sep="\t", index=False)

    elapsed = time.time() - t0
    log(f"Done. elapsed_sec={elapsed:.2f}")


if __name__ == "__main__":
    main()
