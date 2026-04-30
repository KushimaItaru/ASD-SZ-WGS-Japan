"""
tad04292026/common/naming_v1.py

Column-name normalization helpers used across the burden/specificity pipeline.

Design rationale
----------------
Step 5 (``38_compute_sample_burden_L2_and_specificity_v2.py``) emits
sample-burden columns using the CamelCase-dash Heffel labels, e.g.
``n_boundary_HPC_Exc-DG_DEL``. Downstream scripts (Step 7 matched-static
and Scripts 8-11 replication) were written against the legacy snake_case
token list (``hpc_exc_dg`` etc.) inherited from the unified pipeline
outputs. ``normalize_l2_burden_columns`` performs the rename in one place
so producers and consumers can keep their preferred conventions.

What it normalizes
------------------
Exact L2-class columns only:
    n_boundary_HPC_Exc-DG_DEL         -> n_boundary_hpc_exc_dg_DEL
    n_events_HPC_Exc-DG_DUP           -> n_events_hpc_exc_dg_DUP
    carrier_boundary_PFC_Inh-MGE_DEL  -> carrier_boundary_pfc_inh_mge_DEL

What it leaves untouched
------------------------
- Group/specificity columns (double underscore separator):
    n_boundary_group_primary__Diff_specific_n1_DEL   (unchanged)
- Non-L2 burden columns (log1p_total_del_bases, PC1, Sex, etc.)
- Any column whose body (between prefix and trailing SV type) does not
  exactly match a known L2 class label.

The SV type suffix (``_DEL`` / ``_DUP``) is always preserved verbatim —
this avoids the earlier bug in v6 where a naive ``.lower()`` converted
``_DEL`` into ``_del`` and broke ``get_observed_exposure_column()``.
"""

from __future__ import annotations
from typing import Dict, Iterable, Tuple

import pandas as pd

from common.paths_v1 import L2_CLASSES


_BURDEN_PREFIXES = ("n_boundary_", "n_events_", "carrier_boundary_")
_SV_SUFFIXES     = ("_DEL", "_DUP")


def l2_class_snake(cls: str) -> str:
    """Convert a single CamelCase-dash L2 class label to snake_case.

    ``"HPC_Exc-DG"`` -> ``"hpc_exc_dg"``
    ``"PFC_Inh-MGE"`` -> ``"pfc_inh_mge"``
    """
    return cls.replace("-", "_").lower()


def build_l2_rename_map(classes: Iterable[str] = L2_CLASSES) -> Dict[str, str]:
    """Return {CamelCaseDash: snake_case} mapping for all L2 classes."""
    return {c: l2_class_snake(c) for c in classes}


def normalize_l2_burden_columns(
    df: pd.DataFrame,
    classes: Iterable[str] = L2_CLASSES,
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """Rename only L2-class-specific burden columns; preserve SV suffix.

    Returns
    -------
    (df_renamed, rename_map)
        df_renamed: new DataFrame with renamed columns (copy if changes; else same object)
        rename_map: {old_name: new_name} for columns that were renamed.
    """
    cls_map = build_l2_rename_map(classes)
    cls_camel_set = set(cls_map.keys())

    rename_map: Dict[str, str] = {}
    for col in df.columns:
        if not isinstance(col, str):
            continue
        pref = next((p for p in _BURDEN_PREFIXES if col.startswith(p)), None)
        if pref is None:
            continue
        rest = col[len(pref):]
        suf = next((s for s in _SV_SUFFIXES if rest.endswith(s)), None)
        if suf is None:
            continue
        body = rest[:-len(suf)]
        # skip specificity-group columns (they use '__' after a 'group_xxx' prefix)
        if "__" in body:
            continue
        if body not in cls_camel_set:
            continue
        new_col = f"{pref}{cls_map[body]}{suf}"
        if new_col != col:
            rename_map[col] = new_col

    if rename_map:
        df = df.rename(columns=rename_map)
    return df, rename_map


__all__ = [
    "l2_class_snake",
    "build_l2_rename_map",
    "normalize_l2_burden_columns",
]
