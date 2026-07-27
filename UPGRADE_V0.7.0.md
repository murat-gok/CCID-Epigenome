# Upgrade to v0.7.0

## Apply the patch

```bash
cd ~/ccid_epigenome
unzip -o \
  "/mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.7.0.zip" \
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
40 passed
0.7.0
```

SciPy precision-loss warnings can occur in deterministic synthetic tests with nearly identical pair effects; they do not indicate a failed test.

## New files

```text
src/ccid/bystander_validation.py
scripts/04c_validate_bystander_editability.py
scripts/04c_export_epicv2_annotation.R
tests/test_bystander_validation.py
README_BYSTANDER_VALIDATION_V0.7.0.md
UPGRADE_V0.7.0.md
MANIFEST_V0.7.0.txt
```

## Prepare annotation and run

See `README_BYSTANDER_VALIDATION_V0.7.0.md` for the exact commands. Probe-level validation works without annotation; regional clustering requires the optional EPICv2 hg38 annotation export.
