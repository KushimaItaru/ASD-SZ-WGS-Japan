#!/usr/bin/env python3
"""
99_verify_pipeline_consistency_v1.py  (Script 99 Part A スケルトン)

- 処理内容（箇条書き）:
  * tad04212026 パイプライン Scripts 1-7 の出力について、論文記載の数値に
    一切依存せず、パイプライン内部の整合性だけを検証する（Pipeline-as-truth）。
  * 目的は「論文の数値と一致するか」ではなく「パイプラインが内部的に
    正しいか」の検出。論文の typo / 古い数値を根拠にしたくないため。
  * 9 種類のチェックを順次実行し、PASS / FAIL / SKIP / ERROR を判定。
  * 結果を TSV レポート + 詳細ログに書き出す。
  * 実装済みチェック:
      A-7: manifest (JSON) ↔ paths_v1.py 間の完全一致検証
      A-8: Step 6 出力の NaN / Inf 検出（OR / P値 / CI）
  * スタブ（TODO）:
      A-1: サンプル数保存則 (Scripts 3 → 4 → 5 → 6)
      A-2: case / control 合計整合
      A-3: bin 数分解則 (diff_all = static_only + differential)
      A-4: L2 バーデン同値性 (sum over L2 >= n_boundary_any)
      A-5: 共変量健全性 (gene covariate nonzero ratio)
      A-6: Overlap 閾値 0.10 フィルタ効果
      A-9: R glm convergence flag

- 実行例:
    python 99_verify_pipeline_consistency_v1.py
    python 99_verify_pipeline_consistency_v1.py --only A-7,A-8
    python 99_verify_pipeline_consistency_v1.py --outdir /path/to/out

- SBATCH (オプション、遺伝研 ncbn-cpu):
    #SBATCH -p ncbn-cpu
    #SBATCH --account=ncbn-cpu
    #SBATCH --cpus-per-task=4
    #SBATCH --mem=16G
    #SBATCH --time=1:00:00
    #SBATCH --output=verify_pipeline_consistency_v1_%j.log
"""

import argparse
import json
import sys
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# common.paths_v1 を import（PIPELINE_ROOT/common をパスに追加）
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
try:
    from common.paths_v1 import (
        PIPELINE_ROOT,
        OUT_99_VERIFY,
        ensure_output_dirs,
        # Inputs
        F_03_WGS_SV_EVENTS,
        F_04_EVENT_OVERLAP,
        F_04_SAMPLE_BURDEN_COV,
        F_04_SAMPLE_BURDEN_COV_SUM,
        F_05_SAMPLE_BURDEN_L2,
        F_05_SAMPLE_BURDEN_L2_SUM,
        F_06_B_PRIME_L2_RESULTS,
        F_06_COVARIATES,
        F_07_MATCHED_STATIC_MAIN,
        R_MANIFEST_JSON,
        # Part A outputs
        F_99_PIPELINE_CONSISTENCY_REPORT,
        F_99_PIPELINE_CONSISTENCY_LOG,
        # Constants used in manifest cross-check
        L2_CLASSES,
        BIN_SIZE_BP,
    )
except Exception as e:
    sys.stderr.write(
        "[ERROR] common.paths_v1 の import に失敗しました。\n"
        f"PIPELINE_ROOT/common/paths_v1.py を確認してください: {e}\n"
    )
    raise


# ===========================================================================
# フレームワーク: CheckResult dataclass + orchestrator
# ===========================================================================

# ステータス定数
STATUS_PASS = "PASS"
STATUS_FAIL = "FAIL"
STATUS_SKIP = "SKIP"
STATUS_ERROR = "ERROR"
STATUS_TODO = "TODO"


@dataclass
class CheckResult:
    check_id: str            # 例 "A-7"
    title: str               # 短い日本語タイトル
    status: str              # PASS / FAIL / SKIP / ERROR / TODO
    detail: str = ""         # 詳細メッセージ（失敗時の原因など）
    metrics: Dict[str, Any] = field(default_factory=dict)
    elapsed_sec: float = 0.0


def run_check(
    check_id: str,
    title: str,
    fn: Callable[[], Tuple[str, str, Dict[str, Any]]],
) -> CheckResult:
    """
    1 つのチェック関数を安全に実行し CheckResult を返す。
    fn は (status, detail, metrics) を返すこと。
    """
    t0 = time.time()
    try:
        status, detail, metrics = fn()
    except Exception:
        tb = traceback.format_exc()
        return CheckResult(
            check_id=check_id,
            title=title,
            status=STATUS_ERROR,
            detail=f"Unhandled exception:\n{tb}",
            elapsed_sec=time.time() - t0,
        )
    return CheckResult(
        check_id=check_id,
        title=title,
        status=status,
        detail=detail,
        metrics=metrics,
        elapsed_sec=time.time() - t0,
    )


def require_file(path: Path, check_id: str) -> Optional[CheckResult]:
    """入力ファイルが存在しない時は SKIP の CheckResult を返す。"""
    if not path.exists():
        return CheckResult(
            check_id=check_id,
            title="(required file missing)",
            status=STATUS_SKIP,
            detail=f"Required input missing: {path}",
        )
    return None


def log_line(fh, msg: str) -> None:
    """stdout と log ファイルの両方に 1 行書き出す。"""
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


# ===========================================================================
# 各チェックの実装 / スタブ
# ===========================================================================

def check_A1_sample_conservation() -> Tuple[str, str, Dict[str, Any]]:
    """
    Check A-1: サンプル数保存則
      Scripts 3 → 4 → 5 → 6 を通してサンプル集合が保存されているか、
      各ステップでの drop 数とログ上の理由が整合するかを検証。

    期待動作:
      - Step 3 (F_03_WGS_SV_EVENTS) 中のユニーク sample_id 集合 S3
      - Step 4 (F_04_SAMPLE_BURDEN_COV) 中の sample_id 集合 S4
      - Step 5 (F_05_SAMPLE_BURDEN_L2) 中の sample_id 集合 S5
      - Step 6 (F_06_COVARIATES) 中の sample_id 集合 S6
      について S6 ⊆ S5 ⊆ S4 ⊆ S3（逆向きではない）を検証。
      差集合があれば drop 理由が各 Step のログに記録されていることを期待。

    TODO:
      実装時は sample_id を含む列名を dynamic に検出（preference に従う）。
    """
    return (STATUS_TODO, "A-1 未実装。実装時は各ファイルから sample_id 列を"
            "動的に検出し、集合包含関係 (S6 ⊆ S5 ⊆ S4 ⊆ S3) を検証。", {})


def check_A2_case_control_sum() -> Tuple[str, str, Dict[str, Any]]:
    """
    Check A-2: case / control 合計整合
      Step 6 covariates で N_ASD + N_SZ + N_HC == N_total が成立するか。
      Diagnosis 列の値を数えて合計と比較。

    TODO:
      Step 6 の F_06_COVARIATES を読み、Diagnosis 列から case/control カウント。
      他の Diagnosis カテゴリ（例: Other）がある場合は別途 flag。
    """
    return (STATUS_TODO, "A-2 未実装。F_06_COVARIATES の Diagnosis 列カウント実装。",
            {})


def check_A3_bin_decomposition() -> Tuple[str, str, Dict[str, Any]]:
    """
    Check A-3: bin 数分解則
      Step 2 (bin_l2_annotation) において
        N_diff_all = N_static_only + N_differential
      が成立するか。

    TODO:
      F_02_BIN_L2_ANNOTATION を読み、overlaps_diffbound_any / diffbound カテゴリ
      列を dynamic に検出して bin 数を集計。
    """
    return (STATUS_TODO, "A-3 未実装。F_02_BIN_L2_ANNOTATION の bin 分解検証。",
            {})


def check_A4_L2_burden_identity() -> Tuple[str, str, Dict[str, Any]]:
    """
    Check A-4: L2 バーデン同値性
      各サンプルで sum over L2 classes of n_boundary_{class}_{SVTYPE}
      が全 L2 union (n_boundary_any_{SVTYPE} 相当) 以上か。
      bin が複数 L2 に属しうるため >= の関係（= ではない）。

    TODO:
      F_05_SAMPLE_BURDEN_L2 を読み、n_boundary_{class}_DEL / _DUP 列を列名から
      動的に抽出して合算、n_boundary_any_DEL / _DUP と比較。
    """
    return (STATUS_TODO, "A-4 未実装。F_05_SAMPLE_BURDEN_L2 の L2 sum vs any 検証。",
            {})


def check_A5_gene_covariate_sanity() -> Tuple[str, str, Dict[str, Any]]:
    """
    Check A-5: 共変量健全性
      F_04_SAMPLE_BURDEN_COV 中の total_gene_DEL / total_gene_DUP が
      全 case/control サンプルの 95% 以上で > 0 か。
      gene interval index の構築ミスや GTF 読込失敗を早期検出。

    TODO:
      F_04_SAMPLE_BURDEN_COV を読み、total_gene_DEL/DUP の nonzero 比率を算出。
      閾値は 0.95 をデフォルト、argparse で可変。
    """
    return (STATUS_TODO, "A-5 未実装。total_gene_DEL/DUP の nonzero ratio >= 0.95.",
            {})


def check_A6_overlap_threshold_filter() -> Tuple[str, str, Dict[str, Any]]:
    """
    Check A-6: Overlap 閾値 0.10 フィルタ効果
      Step 4 overlap table を --min-overlap-frac-boundary=0.10 で再 filter
      した結果と、最終 burden の event 数が一致するか。
      実際に 0.10 閾値が反映されているかの sanity。

    TODO:
      F_04_EVENT_OVERLAP を読み、overlap_frac_boundary 列を dynamic に検出、
      0.10 filter で残る event 数と F_04_SAMPLE_BURDEN_COV_SUM のカウントを比較。
    """
    return (STATUS_TODO, "A-6 未実装。overlap_frac_boundary >= 0.10 フィルタの再計算。",
            {})


# -----------------------------------------------------------
# A-7: manifest ↔ paths_v1.py 完全一致 (実装済み)
# -----------------------------------------------------------
def check_A7_manifest_consistency() -> Tuple[str, str, Dict[str, Any]]:
    """
    Check A-7: paths_v1.py で計算された値と manifest JSON の値が一致するか。
      - BIN_SIZE_BP (int)
      - L2_CLASSES (list, ソートして比較)
      - 主要ファイルパス (F_05_SAMPLE_BURDEN_L2, OUT_06_B_PRIME_L2 等)
      が完全一致することを確認。R 側が誤ったバージョンの manifest を読むリスクを排除。
    """
    if not R_MANIFEST_JSON.exists():
        return (STATUS_SKIP,
                f"manifest JSON が存在しない: {R_MANIFEST_JSON}\n"
                "まず `python3 -c 'from common.paths_v1 import export_r_manifest; "
                "export_r_manifest()'` で生成してください。",
                {})

    with open(R_MANIFEST_JSON, "r") as fh:
        manifest = json.load(fh)

    errors: List[str] = []
    metrics: Dict[str, Any] = {}

    # BIN_SIZE_BP check
    m_bin = manifest.get("BIN_SIZE_BP")
    metrics["py_BIN_SIZE_BP"] = BIN_SIZE_BP
    metrics["manifest_BIN_SIZE_BP"] = m_bin
    if int(m_bin) != int(BIN_SIZE_BP):
        errors.append(
            f"BIN_SIZE_BP mismatch: paths_v1={BIN_SIZE_BP}, manifest={m_bin}"
        )

    # L2_CLASSES check (順序非依存)
    m_l2 = manifest.get("L2_CLASSES", [])
    metrics["py_L2_CLASSES_n"] = len(L2_CLASSES)
    metrics["manifest_L2_CLASSES_n"] = len(m_l2)
    if sorted(L2_CLASSES) != sorted(m_l2):
        only_py = sorted(set(L2_CLASSES) - set(m_l2))
        only_mf = sorted(set(m_l2) - set(L2_CLASSES))
        errors.append(
            f"L2_CLASSES mismatch: "
            f"only_in_paths_v1={only_py}, only_in_manifest={only_mf}"
        )

    # パスの一致 (manifest は str を持つので等価性比較)
    path_keys = [
        ("F_05_SAMPLE_BURDEN_L2", F_05_SAMPLE_BURDEN_L2),
        ("OUT_06_B_PRIME_L2", F_06_B_PRIME_L2_RESULTS.parent),
        ("F_06_B_PRIME_L2_RESULTS", F_06_B_PRIME_L2_RESULTS),
    ]
    for key, py_path in path_keys:
        m_val = manifest.get(key)
        if m_val is None:
            errors.append(f"manifest key 欠落: {key}")
            continue
        if Path(m_val) != Path(py_path):
            errors.append(
                f"{key} mismatch: paths_v1={py_path}, manifest={m_val}"
            )

    if errors:
        return (STATUS_FAIL,
                "以下の manifest 不一致が検出されました:\n  - " +
                "\n  - ".join(errors),
                metrics)
    return (STATUS_PASS,
            f"manifest ↔ paths_v1 完全一致を確認 "
            f"(BIN_SIZE_BP={BIN_SIZE_BP}, L2_CLASSES n={len(L2_CLASSES)})",
            metrics)


# -----------------------------------------------------------
# A-8: Step 6 B' 出力の NaN / Inf 検出 (実装済み)
# -----------------------------------------------------------
def check_A8_nan_inf_in_results() -> Tuple[str, str, Dict[str, Any]]:
    """
    Check A-8: Step 6 F_06_B_PRIME_L2_RESULTS の OR / CI / P値に
      NaN / Inf が含まれていないか検出。
      ある場合は行を特定しログに記録（収束失敗・separation 等の検出）。
    """
    skip = require_file(F_06_B_PRIME_L2_RESULTS, "A-8")
    if skip:
        return (STATUS_SKIP, skip.detail, {})

    df = pd.read_csv(F_06_B_PRIME_L2_RESULTS, sep="\t", low_memory=False)

    # 動的列検出: OR 系 / P 系 / CI 系の列を柔軟に掴む
    def find_cols(patterns: List[str]) -> List[str]:
        return [c for c in df.columns
                if any(p.lower() in c.lower() for p in patterns)]

    or_cols = find_cols(["OR", "odds_ratio"])
    p_cols  = find_cols(["pvalue", "p_value", "p_val", "pval"])
    ci_cols = find_cols(["ci_lower", "ci_upper", "lower_ci", "upper_ci",
                         "CI_L", "CI_U", "lcl", "ucl"])

    # フィルタ: 数値 dtype のみを対象
    def numeric_only(cols: List[str]) -> List[str]:
        return [c for c in cols if pd.api.types.is_numeric_dtype(df[c])]

    or_cols = numeric_only(or_cols)
    p_cols  = numeric_only(p_cols)
    ci_cols = numeric_only(ci_cols)

    target_cols = or_cols + p_cols + ci_cols
    if not target_cols:
        return (STATUS_FAIL,
                "OR / P値 / CI の候補列が検出できませんでした。"
                f"列名: {list(df.columns)[:15]}...",
                {"n_rows": len(df)})

    nan_counts = {c: int(df[c].isna().sum()) for c in target_cols}
    inf_counts = {c: int(np.isinf(df[c]).sum()) for c in target_cols}

    bad_cols = [c for c in target_cols
                if nan_counts[c] > 0 or inf_counts[c] > 0]

    metrics = {
        "n_rows": len(df),
        "n_or_cols": len(or_cols),
        "n_p_cols": len(p_cols),
        "n_ci_cols": len(ci_cols),
        "nan_by_col": nan_counts,
        "inf_by_col": inf_counts,
    }

    if bad_cols:
        # FAIL 詳細: どの行に NaN/Inf があるか (最大 10 行)
        detail_lines = []
        for c in bad_cols:
            bad_mask = df[c].isna() | np.isinf(df[c])
            n_bad = int(bad_mask.sum())
            detail_lines.append(f"  {c}: NaN/Inf 行 = {n_bad}")
            bad_rows_head = df.index[bad_mask].tolist()[:10]
            if bad_rows_head:
                detail_lines.append(f"    行 index (先頭 10): {bad_rows_head}")
        return (STATUS_FAIL,
                f"OR / P値 / CI に NaN または Inf を検出\n" +
                "\n".join(detail_lines),
                metrics)

    return (STATUS_PASS,
            f"NaN / Inf なし (OR 列 {len(or_cols)}, P 列 {len(p_cols)}, "
            f"CI 列 {len(ci_cols)}, 全 {len(df)} 行)",
            metrics)


def check_A9_glm_convergence() -> Tuple[str, str, Dict[str, Any]]:
    """
    Check A-9: Step 6 R glm の .converged フラグが全行 TRUE か。
      R 側で glm(..., family=binomial) の結果に .converged == FALSE が
      混じっていないかを結果 TSV の該当列から確認。

    TODO:
      R が現在 converged 列を出力しているか未確認。出力していなければ
      R 側に列追加が必要。列名を動的検出してカウント。
    """
    return (STATUS_TODO, "A-9 未実装。R 側で .converged 列を出力しているか"
            "未確認。実装時に列を動的検出。",
            {})


# ===========================================================================
# メイン orchestrator
# ===========================================================================

# (check_id, title, callable) のリスト。順序はレポート出力順。
CHECK_REGISTRY: List[Tuple[str, str, Callable]] = [
    ("A-1", "サンプル数保存則 (Scripts 3 → 4 → 5 → 6)",
     check_A1_sample_conservation),
    ("A-2", "case / control 合計整合 (Step 6 共変量)",
     check_A2_case_control_sum),
    ("A-3", "bin 数分解則 (Step 2: diff_all = static + differential)",
     check_A3_bin_decomposition),
    ("A-4", "L2 バーデン同値性 (Step 5: sum over L2 >= any)",
     check_A4_L2_burden_identity),
    ("A-5", "共変量健全性 (Step 4: gene covariate nonzero ratio >= 0.95)",
     check_A5_gene_covariate_sanity),
    ("A-6", "Overlap 閾値フィルタ効果 (Step 4: 0.10 再計算一致)",
     check_A6_overlap_threshold_filter),
    ("A-7", "manifest ↔ paths_v1.py 完全一致",
     check_A7_manifest_consistency),
    ("A-8", "Step 6 B' 出力の NaN / Inf 検出",
     check_A8_nan_inf_in_results),
    ("A-9", "Step 6 R glm convergence flag",
     check_A9_glm_convergence),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Script 99 Part A: pipeline internal consistency checks"
    )
    p.add_argument("--outdir", type=Path, default=OUT_99_VERIFY,
                   help="Output directory (default: OUT_99_VERIFY)")
    p.add_argument("--report", type=Path, default=F_99_PIPELINE_CONSISTENCY_REPORT,
                   help="TSV report path")
    p.add_argument("--log", type=Path, default=F_99_PIPELINE_CONSISTENCY_LOG,
                   help="Log file path")
    p.add_argument("--only", type=str, default=None,
                   help="カンマ区切りの check_id リスト (例: 'A-7,A-8') で "
                        "実行するチェックを絞る。")
    p.add_argument("--skip", type=str, default=None,
                   help="カンマ区切りの check_id で除外するチェック。")
    return p.parse_args()


def select_checks(args: argparse.Namespace
                  ) -> List[Tuple[str, str, Callable]]:
    if args.only:
        keep = set(c.strip() for c in args.only.split(","))
        return [x for x in CHECK_REGISTRY if x[0] in keep]
    if args.skip:
        drop = set(c.strip() for c in args.skip.split(","))
        return [x for x in CHECK_REGISTRY if x[0] not in drop]
    return list(CHECK_REGISTRY)


def render_report_row(r: CheckResult) -> Dict[str, Any]:
    # metrics は JSON 文字列化（TSV 読みやすさのため）
    return {
        "check_id": r.check_id,
        "title": r.title,
        "status": r.status,
        "elapsed_sec": round(r.elapsed_sec, 3),
        "detail_head": r.detail.splitlines()[0] if r.detail else "",
        "metrics_json": json.dumps(r.metrics, ensure_ascii=False,
                                   default=str),
    }


def main() -> int:
    t0 = time.time()
    args = parse_args()

    ensure_output_dirs()
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.log.parent.mkdir(parents=True, exist_ok=True)

    checks = select_checks(args)

    with open(args.log, "w") as fh:
        log_line(fh, "=" * 70)
        log_line(fh, "99_verify_pipeline_consistency_v1.py (Part A)")
        log_line(fh, f"PIPELINE_ROOT: {PIPELINE_ROOT}")
        log_line(fh, f"N checks scheduled: {len(checks)}")
        log_line(fh, "=" * 70)

        results: List[CheckResult] = []
        for check_id, title, fn in checks:
            log_line(fh, f"[{check_id}] {title} ...")
            r = run_check(check_id, title, fn)
            results.append(r)
            log_line(fh, f"    -> {r.status} ({r.elapsed_sec:.2f}s)")
            if r.detail:
                # detail は多行のことがある
                for line in r.detail.splitlines():
                    log_line(fh, f"      | {line}")

        # -------------------------------------------------------------
        # レポート TSV 出力
        # -------------------------------------------------------------
        report_df = pd.DataFrame([render_report_row(r) for r in results])
        args.report.parent.mkdir(parents=True, exist_ok=True)
        report_df.to_csv(args.report, sep="\t", index=False)
        log_line(fh, "-" * 70)
        log_line(fh, f"Report TSV written: {args.report}")

        # -------------------------------------------------------------
        # サマリ: ステータス別件数
        # -------------------------------------------------------------
        from collections import Counter
        status_counts = Counter(r.status for r in results)
        log_line(fh, "Status summary:")
        for s in [STATUS_PASS, STATUS_FAIL, STATUS_ERROR,
                  STATUS_SKIP, STATUS_TODO]:
            log_line(fh, f"    {s:<6s}: {status_counts.get(s, 0)}")

        elapsed = time.time() - t0
        log_line(fh, f"Total elapsed: {elapsed:.2f} s")
        log_line(fh, "=" * 70)

    # exit code: FAIL or ERROR がある場合に非ゼロで返す
    # (SKIP / TODO はゼロのまま = CI で warning 扱い)
    n_bad = sum(1 for r in results
                if r.status in (STATUS_FAIL, STATUS_ERROR))
    return 0 if n_bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
