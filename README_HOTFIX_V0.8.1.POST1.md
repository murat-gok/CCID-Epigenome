# CCID-Epigenome v0.8.1.post1 hotfix

This post-release hotfix fixes pandas 3.x missing-value handling in linked-gene and region-gene aggregation.

## Symptom

`python -m pytest -q` failed in three `test_regulatory_safety.py` tests with:

```text
TypeError: sequence item 0: expected str instance, float found
```

## Cause

Under pandas 3.x, missing values may remain scalar `NaN` objects in string-oriented columns. The v0.8.1 region aggregation attempted to join the resulting mixed string/float sequence.

## Fix

- Parse linked genes one value at a time.
- Exclude `NaN`, `pd.NA`, `None`, `null`, and related missing markers.
- Deduplicate and sort only validated gene symbols.
- Add a regression test for pandas 3.x missing scalar handling.
- Bump the package version to `0.8.1.post1`.

Scientific scores, weights, reference accessions, and interpretation guardrails are unchanged.

## Install over v0.8.1

```bash
cd ~/ccid_epigenome
unzip -o "/mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.8.1.post1.zip" -d ~/ccid_epigenome
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
python -c "import ccid; print(ccid.__version__)"
```

Expected:

```text
54 passed
0.8.1.post1
```

The regulatory reference files already downloaded into `data/reference/v081` do not need to be downloaded again.
