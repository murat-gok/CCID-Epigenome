# Upgrade CCID-Epigenome v0.2.0 → v0.3.0

## 1. Apply the patch

Extract the patch ZIP into the existing project directory, allowing files to be overwritten.

```bash
cd ~/ccid_epigenome
unzip -o /mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.3.0.zip -d .
```

## 2. Reinstall and test

```bash
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
```

Expected: `18 passed`.

## 3. Download the two newly required files

```bash
python scripts/01_download_public_data.py \
  --datasets gse165177_rna_part2,gse269758_tcell_bystander_r3 \
  --output data/raw
```

## 4. Run preprocessing

```bash
python scripts/02_preprocess_all.py \
  --raw-dir data/raw \
  --processed-dir data/processed \
  --chunksize 5000
```

## 5. Inspect QC

```bash
find data/processed -maxdepth 2 -type f -printf "%p\t%s bytes\n" | sort
cat data/processed/qc/fibroblast_sample_overlap.json
```
