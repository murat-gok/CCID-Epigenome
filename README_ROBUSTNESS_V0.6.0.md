# CCID-Epigenome v0.6.0 — O3 robustness and time-series validation

This stage is run **after** the frozen O1/O2 signature and the first held-out O3 validation. It does not modify the frozen 500-CpG candidate list.

## Analyses implemented

1. **Direct successful-versus-failed comparison**
   - O3 successful day-10 sample minus O3 failed day-10 sample.
   - Uses the frozen day/batch/stage inferred from the O1/O2 matched-pair file.

2. **Multi-baseline sensitivity**
   - Successful minus exact same-day baseline.
   - Successful minus the mean of other-day O3 baselines.
   - Successful minus the mean of all O3 baselines in the frozen batch/stage.
   - Successful minus failed reprogramming in the same stratum.
   - Reports raw and genome-wide-median-centred metrics.

3. **Direction- and effect-matched empirical null**
   - Random CpG sets preserve both discovery direction and absolute-effect bin.
   - Runs on raw and global-median-centred effects.
   - This is the preferred null for the strongly one-directional candidate list.

4. **Leave-one-donor-out robustness**
   - O2+O3 → O1, O1+O3 → O2, O1+O2 → O3.
   - Each rotation rebuilds its candidate list using only the two non-held-out donors.
   - Rotation overlap and three-way consensus are reported.

5. **Time-series consistency**
   - Uses all available O1/O2/O3 `reprogrammed`, `reprogrammed_intermediate`, and `failed` samples in the frozen batch.
   - Every sample is compared with an exact donor/day/batch/stage baseline.
   - Reports per-timepoint metrics and per-CpG success/failed consistency.

## Install the patch

```bash
cd ~/ccid_epigenome
unzip -o /mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.6.0.zip -d ~/ccid_epigenome
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
python -c "import ccid; print(ccid.__version__)"
```

Expected:

```text
36 passed
0.6.0
```

## Run the complete robustness stage

```bash
mkdir -p results/real/o3_robustness

python scripts/04b_validate_signature_robustness.py \
  --beta data/processed/fibroblast_beta.csv.gz \
  --metadata data/processed/fibroblast_methylation_metadata.csv \
  --signature results/real/safe_osk_target_signature_O1O2.csv.gz \
  --discovery-pairs results/real/safe_osk_target_signature_O1O2_matched_pairs.csv \
  --all-loci results/real/safe_osk_target_signature_O1O2_all_loci.csv.gz \
  --output-dir results/real/o3_robustness \
  --validation-donor O3 \
  --donors O1,O2,O3 \
  --lodo-top-n 500 \
  --lodo-min-abs-delta 0.03 \
  --null-iterations 1000 \
  --null-bins 10 \
  --seed 20260721 \
  --chunksize 10000
```

## Main outputs

```text
results/real/o3_robustness/
├── robustness_summary.json
├── robustness_analysis_manifest.json
├── o3_successful_vs_failed.csv.gz
├── o3_baseline_sensitivity_metrics.csv
├── o3_baseline_sensitivity_cpg.csv.gz
├── o3_direction_balanced_null.csv.gz
├── lodo_summary.csv
├── lodo_candidates.csv.gz
├── lodo_overlap.csv
├── time_series_metrics.csv
├── time_series_cpg_consistency.csv.gz
├── time_series_sample_manifest.csv
└── figures/
    ├── baseline_sensitivity.png
    ├── direction_balanced_null.png
    ├── lodo_sign_concordance.png
    ├── successful_vs_failed.png
    └── time_series_sign_concordance.png
```

## First inspection commands

```bash
cat results/real/o3_robustness/robustness_summary.json
column -s, -t < results/real/o3_robustness/o3_baseline_sensitivity_metrics.csv
column -s, -t < results/real/o3_robustness/lodo_summary.csv
column -s, -t < results/real/o3_robustness/lodo_overlap.csv
column -s, -t < results/real/o3_robustness/time_series_metrics.csv | less -S
```

## Interpretation safeguards

- O3 remains one independent held-out donor; population-level confirmation is not claimed.
- LODO rotations are sensitivity analyses, not three independent confirmatory studies.
- Raw and globally centred results must be reported together when global methylation shift is plausible.
- The frozen O1/O2 500-CpG list and all thresholds remain unchanged after seeing O3.
