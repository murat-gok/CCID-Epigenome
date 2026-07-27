# CCID-Epigenome v0.4.1 upgrade

## What this release adds

1. Safe-window scores are calibrated using discovery donors only; held-out O3 cannot influence O1/O2 score scaling.
2. Safe samples are matched exactly to negative controls by donor, day, batch and stage.
3. Pair effects are averaged within donor before CpG ranking.
4. With only O1 and O2, inferential FDR selection is explicitly marked underpowered; the output is an exploratory ranked target list for O3 validation.
5. EPIC-scale beta matrices are processed in chunks.

## Install

```bash
cd ~/ccid_epigenome
unzip -o /mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.4.1.zip -d ~/ccid_epigenome
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
```

Expected: `25 passed`.

## Step 1 — Discovery-only safe window

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
  --utility-tolerance 0.10 \
  --allowed-conditions reprogrammed,reprogrammed_intermediate \
  --output results/real/safe_window_scores_O1O2.csv \
  --summary-output results/real/safe_window_summary_O1O2.json
```

Inspect:

```bash
cat results/real/safe_window_summary_O1O2.json
```

## Step 2 — Exact-matched methylation target signature

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
  --min-inferential-donors 3 \
  --fallback-top-n 500 \
  --chunksize 5000 \
  --output results/real/safe_osk_target_signature_O1O2.csv.gz
```

Generated files:

- `safe_osk_target_signature_O1O2.csv.gz` — ranked target CpGs; with two discovery donors this is exploratory.
- `safe_osk_target_signature_O1O2_all_loci.csv.gz` — statistics for all CpGs.
- `safe_osk_target_signature_O1O2_matched_pairs.csv` — exact safe/control matches.
- `safe_osk_target_signature_O1O2_donor_effects.csv.gz` — independent O1/O2 effects.
- `safe_osk_target_signature_O1O2_qc.json` — matching and inferential status.

## Required checks

```bash
cat results/real/safe_osk_target_signature_O1O2_qc.json
column -s, -t < results/real/safe_osk_target_signature_O1O2_matched_pairs.csv | less -S
```

Expected QC logic:

- `matched_donors` contains O1 and O2.
- `inferential_status` is `exploratory_underpowered` because only two discovery donors are used.
- `used_underpowered_fallback` is normally `true`.
- O3 is not used for safe-window selection or target-CpG ranking; it remains held out for validation.
