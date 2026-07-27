# Upgrade to v0.7.1

## Apply the patch

```bash
cd ~/ccid_epigenome
unzip -o \
  "/mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.7.1.zip" \
  -d ~/ccid_epigenome
```

## Reinstall and test

```bash
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
python -c "import ccid; print(ccid.__version__)"
```

Expected:

```text
45 passed
0.7.1
```

SciPy precision-loss warnings in synthetic tests are non-fatal.

## Run the analysis

Use the command in `README_BYSTANDER_VALIDATION_V0.7.1.md`. The new run writes to:

```text
results/real/bystander_validation_v071
```

The existing v0.7.0 directory is not overwritten.
