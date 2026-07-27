# Upgrade CCID-Epigenome v0.3.0 → v0.3.1

This patch fixes a pandas 3.x compatibility error in real-data preprocessing:

```text
ValueError: assignment destination is read-only
```

## 1. Apply the patch

```bash
cd ~/ccid_epigenome
unzip -o /mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.3.1.zip -d .
```

## 2. Reinstall and test

```bash
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
```

Expected result:

```text
................... [100%]
19 passed
```

## 3. Resume preprocessing

You do not need to redownload data. Remove partial processed outputs, then rerun:

```bash
rm -rf data/processed/*
python scripts/02_preprocess_all.py \
  --raw-dir data/raw \
  --processed-dir data/processed \
  --chunksize 5000 \
  --detection-p-threshold 0.01 \
  --max-locus-missing 0.20
```
