# V18 → V19 Patch v3: Robust Sex mapping + Dual dump

## v2 → v3 changes (ChatGPT v4 review Priority 1-3 fix)
- Sex mapping を robust 化 (M/F, 1/0, 1/2, Male/Female, true/false 対応)
- v2 の dual dump (event-bin + sample covariates) は維持
- ancestry を string で保存する設計は維持 (downstream の `pd.get_dummies` で sanitize 化)

## Insertion location (v18 line ~1037)

Same as v1/v2 patches: **immediately before `del disrupted_bins_a, disrupted_bins_c`** in `run_phase4_bin_counts()`.

## Step 1: Copy v18 → v19

```bash
cd /lustre12/home/kushima-pg/tad04292026/09_mssng_sample_burden/
cp tad_replication_mssng_v18.py tad_replication_mssng_v19.py
cp 09a_mssng_replication_v18.sbatch 09a_mssng_replication_v19.sbatch
# Edit 09a_..._v19.sbatch to invoke v19.py
```

## Step 2: Modify `run_phase4_bin_counts` signature

Find:
```python
def run_phase4_bin_counts(sv, sample_status_dict, bin_index):
```

Replace with:
```python
def run_phase4_bin_counts(sv, sample_status_dict, bin_index, fid_map=None,
                            get_sex=None, get_platform_norm=None, get_ancestry=None):
```

## Step 3: Update the caller (around line 1501)

```python
# Before:
phase4_out = run_phase4_bin_counts(sv, sample_status_dict, bin_index)

# After:
phase4_out = run_phase4_bin_counts(
    sv, sample_status_dict, bin_index,
    fid_map=fid_map,
    get_sex=get_sex,
    get_platform_norm=get_platform_norm,
    get_ancestry=get_ancestry,
)
```

## Step 4: Insert dump block in PHASE 4 (v3: robust Sex mapping)

**INSERT THE FOLLOWING BLOCK IMMEDIATELY BEFORE `del disrupted_bins_a, disrupted_bins_c`**:

```python
    # ========================================================================
    # === V19 ADDITION (patch v3): Dual dump with robust Sex mapping
    # ========================================================================
    log("=" * 70)
    log("V19 ADDITION (patch v3): Dumping event-bin + sample covariate")
    log("=" * 70)
    import os as _os
    _OUT_DIR = "/lustre12/home/kushima-pg/tad04292026/14_constraint_enrichment_v1/output_v1"
    _os.makedirs(_OUT_DIR, exist_ok=True)

    # ---- v3 Sex mapping (robust to M/F, 1/0, 1/2, Male/Female, true/false) ----
    def _sex_to_numeric(x):
        """Map various Sex encodings to 1=Male, 0=Female, NaN=unknown.
        Verify against MSSNG subject_table.csv to confirm 1/2 coding direction."""
        s = str(x).strip().upper()
        if s in ("M", "MALE", "1", "TRUE"):
            return 1
        if s in ("F", "FEMALE", "0", "2", "FALSE"):
            return 0
        return float("nan")

    # ---- Dump A: event-bin records (sample, bin_id, sv_type, pattern, Diagnosis) ----
    _seen_pat_a = {svt: {} for svt in SV_TYPES}
    for _bclass in ALL_BCLASSES:
        for _svt in SV_TYPES:
            for _si, _bins in enumerate(disrupted_bins_a[_bclass][_svt]):
                _sample = all_samples[_si]
                if _sample not in _seen_pat_a[_svt]:
                    _seen_pat_a[_svt][_sample] = set()
                _seen_pat_a[_svt][_sample].update(_bins)

    def _diag_label(_status):
        if _status in ("Affected", "ASD", "Proband"):
            return "ASD"
        elif _status in ("Sibling", "Unaffected"):
            return "Sibling"
        return _status

    _evbin_records = []
    for _svt in SV_TYPES:
        for _sample, _bin_set in _seen_pat_a[_svt].items():
            _diag = _diag_label(sample_status_dict.get(_sample, ""))
            for _bid in _bin_set:
                _evbin_records.append({
                    "sample_id": _sample,
                    "bin_id": _bid,
                    "sv_type": _svt,
                    "pattern": "A",
                    "Diagnosis": _diag,
                })

    _evbin_path = _os.path.join(_OUT_DIR, "mssng_event_bins_dumped_v1.tsv.gz")
    pd.DataFrame(_evbin_records).to_csv(_evbin_path, sep="\t", index=False, compression="gzip")
    log(f"  Saved {len(_evbin_records)} event-bin records: {_evbin_path}")

    # ---- Dump B: sample-level covariate file (REQUIRED by 14_fit v5) ----
    _cov_records = []
    for _si, _sample in enumerate(all_samples):
        _diag = _diag_label(sample_status_dict.get(_sample, ""))
        _rec = {
            "sample_id": _sample,
            "Diagnosis": _diag,
            "log1p_total_del_bases": float(np.log1p(sample_total_base["DEL"][_si])),
            "log1p_total_gene_DEL": float(np.log1p(sample_total_gene["DEL"][_si])),
        }
        # FAMILYID (NaN → _solo_<sample>)
        if fid_map is not None:
            _rec["FAMILYID"] = fid_map.get(_sample, f"_solo_{_sample}")
        else:
            _rec["FAMILYID"] = f"_solo_{_sample}"
        # Sex (v3: robust mapping)
        if get_sex is not None:
            _rec["Sex"] = _sex_to_numeric(get_sex(_sample))
        else:
            _rec["Sex"] = float("nan")
        # Platform (string; downstream sanitizes)
        if get_platform_norm is not None:
            _rec["Platform"] = get_platform_norm(_sample) or ""
        else:
            _rec["Platform"] = ""
        # Ancestry (string; downstream sanitizes + dummies)
        if get_ancestry is not None:
            _rec["ancestry"] = get_ancestry(_sample) or ""
        else:
            _rec["ancestry"] = ""
        _cov_records.append(_rec)

    _cov_path = _os.path.join(_OUT_DIR, "mssng_sample_covariates_dumped_v1.tsv.gz")
    pd.DataFrame(_cov_records).to_csv(_cov_path, sep="\t", index=False, compression="gzip")
    log(f"  Saved {len(_cov_records)} sample covariate records: {_cov_path}")

    # ---- v3 Sanity: Sex NaN count ----
    _df_chk = pd.DataFrame(_cov_records)
    _n_sex_na = _df_chk["Sex"].isna().sum()
    log(f"  Sex NaN count: {_n_sex_na} / {len(_df_chk)}")
    if _n_sex_na > 0:
        log(f"  WARNING: {_n_sex_na} samples have Sex=NaN. "
            f"Verify subject_table.csv encoding before downstream pipeline.")

    del _seen_pat_a, _evbin_records, _cov_records, _df_chk
    # === END V19 ADDITION (patch v3) ===
```

## Step 5: Pre-flight check on subject_table.csv Sex encoding

**重要**: Patch を流す前に MSSNG `subject_table.csv` の Sex 列の unique values を確認:

```bash
cd /lustre12/home/kushima-pg/resource/mssng/
awk -F',' 'NR==1{for(i=1;i<=NF;i++) if($i=="SEX"||$i=="Sex"||$i=="sex") c=i; next} {print $c}' subject_table.csv | sort -u
```

期待される結果と patch v3 の対応:

| Encoding | Mapping (v3) | Verification needed |
|---|---|---|
| M / F | M=1, F=0 | 直接対応 |
| Male / Female | Male=1, Female=0 | 直接対応 |
| 1 / 0 | 1=1, 0=0 | 直接対応 |
| 1 / 2 | **1=Male, 2=Female** ← 確認必須 | MSSNG が PLINK 1/2 coding を採用しているか確認 |
| その他 | NaN | warning + 14_compute_per_sample_burden_v5 で hard fail |

**1/2 encoding の場合**: ABRF/MSSNG メタデータが PLINK 標準（1=Male, 2=Female）を採用しているなら patch v3 のまま OK。逆方向の場合は `_sex_to_numeric` の `("F", "FEMALE", "0", "2", "FALSE")` を見直してください。

## Step 6: Re-submit v19

```bash
cd /lustre12/home/kushima-pg/tad04292026/09_mssng_sample_burden/
sbatch 09a_mssng_replication_v19.sbatch
```

## Step 7: Verify both dumps

```bash
ls -lh /lustre12/home/kushima-pg/tad04292026/14_constraint_enrichment_v1/output_v1/mssng_*

zcat .../mssng_event_bins_dumped_v1.tsv.gz | head -3
zcat .../mssng_sample_covariates_dumped_v1.tsv.gz | head -3

# v3 critical check: Sex NaN count
zcat .../mssng_sample_covariates_dumped_v1.tsv.gz | \
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Sex") c=i; next}
                $c==""||$c=="nan"{n++} END{print "Sex NaN:", n+0}'
```

## Verification checklist

After v19 v3 runs:
- [ ] v18 results (Pattern A + Pattern C primary replication) reproduced
- [ ] `mssng_event_bins_dumped_v1.tsv.gz` exists with expected columns
- [ ] `mssng_sample_covariates_dumped_v1.tsv.gz` exists with all required B' covariates
- [ ] Sex NaN count == 0 (or fully accounted by missing subject_table entries)
- [ ] Number of unique samples in covariate dump == number of MSSNG analytic samples
- [ ] Diagnosis distribution matches v18 (ASD probands + Siblings)
- [ ] FAMILYID is non-empty for non-`_solo_` rows
- [ ] log1p_total_del_bases / log1p_total_gene_DEL are non-zero for samples with rare DEL
