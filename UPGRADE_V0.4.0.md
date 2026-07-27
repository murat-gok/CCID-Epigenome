# CCID-Epigenome v0.4.0 upgrade

## Purpose

This release implements the next real-data stage:

1. leakage-free safe-window selection on discovery donors;
2. exact negative-control matching by donor, day, experimental batch and stage;
3. donor-level methylation target signature construction;
4. chunked processing of EPIC-scale beta matrices;
5. pairing and QC reports for auditability.

## Upgrade

```bash
cd ~/ccid_epigenome
unzip -o /mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.4.0.zip -d ~/ccid_epigenome
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
```

Expected: `23 passed`.

## Step 1: define the discovery safe window

```bash
mkdir -p results/real

python scripts/02_define_safe_window.py \
  --expression data/processed/fibroblast_expression.csv.gz \
  --metadata data/processed/fibroblast_expression_metadata.csv \
  --discovery-donors O1,O2 \
  --identity-floor -0.50 \
  --pluripotency-ceiling 0.50 \
  --min-rejuvenation 0.0 \
  --min-donors-per-day 2 \
  --allowed-conditions reprogrammed,reprogrammed_intermediate \
  --output results/real/safe_window_scores_O1O2.csv \
  --summary-output results/real/safe_window_summary_O1O2.json
```

Inspect:

```bash
cat results/real/safe_window_summary_O1O2.json
```

## Step 2: build the exact-matched methylation target signature

```bash
python scripts/03_build_target_signature.py \
  --beta data/processed/fibroblast_beta.csv.gz \
  --metadata data/processed/fibroblast_methylation_metadata.csv \
  --safe-window results/real/safe_window_scores_O1O2.csv \
  --discovery-donors O1,O2 \
  --match-columns donor_id,day,batch,stage \
  --min-abs-delta 0.03 \
  --fdr 0.10 \
  --min-sign-consistency 1.0 \
  --min-matched-donors 2 \
  --chunksize 5000 \
  --output results/real/safe_osk_target_signature_O1O2.csv.gz
```

Generated files:

- `safe_osk_target_signature_O1O2.csv.gz`: strict signature, or ranked fallback if underpowered;
- `safe_osk_target_signature_O1O2_all_loci.csv.gz`: all CpG statistics;
- `safe_osk_target_signature_O1O2_matched_pairs.csv`: exact safe/control matches;
- `safe_osk_target_signature_O1O2_donor_effects.csv.gz`: O1/O2 donor-level effects;
- `safe_osk_target_signature_O1O2_qc.json`: matching and selection QC.

## Required checks

```bash
cat results/real/safe_osk_target_signature_O1O2_qc.json
column -s, -t < results/real/safe_osk_target_signature_O1O2_matched_pairs.csv | less -S
```

The QC report should show both O1 and O2 under `matched_donors`. Unmatched samples are not silently assigned to controls from another day or batch.
