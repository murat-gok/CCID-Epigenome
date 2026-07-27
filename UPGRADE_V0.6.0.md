# Upgrade to CCID-Epigenome v0.6.0

## 1. Download and apply the patch

```bash
cd ~/ccid_epigenome
unzip -o "/mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.6.0.zip" -d ~/ccid_epigenome
```

## 2. Reinstall and verify

```bash
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

## 3. Run the next scientific test package

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

The beta matrix is read in chunks. The direction-balanced null is run twice (raw and global-centred), so the null output contains 2,000 rows with the default settings.
