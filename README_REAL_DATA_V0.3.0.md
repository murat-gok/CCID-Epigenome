# Real-data preprocessing (v0.3.0)

This release adds streaming preprocessing for the public MPTR and epigenome-editing datasets.

## Critical completeness check

Two files were missing from the original download list and are required before preprocessing:

```bash
python scripts/01_download_public_data.py \
  --datasets gse165177_rna_part2,gse269758_tcell_bystander_r3 \
  --output data/raw
```

The required files are:

- `GSE165177_Log2_RPM_Transient_reprogramming_part2_170621.txt.gz`
- `GSE269758_T_cell_R3_Hyperlibrary_beta_detp.csv.gz`

The first completes the GSE165177 expression trajectory. The second adds the third T-cell targeted/control replicate.

## Install the updated package

```bash
cd ~/ccid_epigenome
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
```

Expected result: `18 passed` (plus the known precision-loss warning in one synthetic test).

## Run all preprocessing

```bash
python scripts/02_preprocess_all.py \
  --raw-dir data/raw \
  --processed-dir data/processed \
  --chunksize 5000 \
  --detection-p-threshold 0.01 \
  --max-locus-missing 0.20
```

The large GSE165179 matrix is processed in chunks rather than loaded into memory at once.

## Expected processed files

```text
data/processed/
├── fibroblast_expression.csv.gz
├── fibroblast_expression_metadata.csv
├── fibroblast_beta.csv.gz
├── fibroblast_methylation_metadata.csv
├── tcell_beta.csv.gz
├── tcell_metadata.csv
├── msc_beta.csv.gz
├── msc_metadata.csv
└── qc/
    ├── fibroblast_expression_qc.json
    ├── fibroblast_methylation_sample_qc.csv
    ├── fibroblast_methylation_summary.json
    ├── bystander_preprocessing_summary.json
    └── fibroblast_sample_overlap.json
```

## What the preprocessing does

- Merges both GSE165177 expression parts.
- Extracts donor, day, experimental batch and biological condition from sample names.
- Maps negative controls and fibroblast samples to `baseline`.
- Maps successful MPTR samples to `reprogrammed` or `reprogrammed_intermediate`.
- Maps failed reprogramming samples to `failed`, preventing accidental inclusion in the safe-window search.
- Reads beta/detection-p column pairs from GSE165179 in chunks.
- Masks beta values with detection p-value above 0.01 or values outside [0, 1].
- Removes loci exceeding the prespecified missingness threshold.
- Merges T-cell replicates 1–2 with replicate 3.
- Normalizes EPIC-v2 probe suffixes such as `cg21870274_BC21` to base CpG IDs and reports duplicate collapsing.
- Writes sample-level and dataset-level QC reports.
- Checks expression/methylation sample overlap and requires at least two paired donors.

## Run the next analysis stages

Inspect the overlap report first:

```bash
cat data/processed/qc/fibroblast_sample_overlap.json
```

Then define discovery donors. Do not use held-out donors when selecting thresholds or days.

```bash
python scripts/02_define_safe_window.py \
  --expression data/processed/fibroblast_expression.csv.gz \
  --metadata data/processed/fibroblast_expression_metadata.csv \
  --discovery-donors O1,O2 \
  --identity-floor -0.50 \
  --pluripotency-ceiling 0.50 \
  --output results/real/safe_window_scores.csv
```

The exact discovery/validation split must be chosen after inspecting donor coverage in the overlap report.
