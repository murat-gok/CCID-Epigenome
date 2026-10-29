# Restoration record — 2026-10-29

## What was wrong

The working tree used to produce the manuscript's post5 results lived in WSL at
`~/ccid_epigenome` and was lost. Two folders survived on the Windows side:

- `C:\Users\user\CCID-Epigenome`       — git repo, source, no results
- `C:\Users\user\CCID-Epigenome_build` — a dead staging copy, strict subset, post4

## The finding that matters for the manuscript

`C:\Users\user\CCID-Epigenome` is tagged **v0.9.0.post5** in git and its HEAD
commit is "Add CCID-Epigenome v0.9.0.post5 reproducibility package". Despite
that:

- `src/ccid/__init__.py` read `__version__ = "0.9.0.post4"`
- `pyproject.toml` read `version = "0.9.0.post4"`
- `grep -c Sample_description src/ccid/orthogonal_validation.py` returned **0**,
  and preserving `!Sample_description` is the entire purpose of post5
- `apply_post5.py` was committed but had never been executed
- **the post5 acceptance tests failed**: 2 failed, 70 passed

So the public tag `github.com/murat-gok/CCID-Epigenome/tree/v0.9.0.post5`
contains post4 code plus an unapplied patch script. The manuscript's data
availability statement points at that tag. A reviewer who clones it cannot
reproduce the post5 GSE179847 results, and the post5 tests in that tag fail
against its own source. **This must be corrected before submission.**

## What was done

1. Ran the official `apply_post5.py` against the tree. No `--force` was needed
   because the tree is genuine post4, which is what the patch expects.
   Backups: `*.post4.bak`. Result: version strings now post5,
   `Sample_description` handling present, post5 acceptance tests pass.

2. Two post4-era tests then still failed:
   `test_gse179847_processed_csv_streaming_with_direct_titles` and
   `test_natural_aging_accepts_gse179847_processed_csv`.

   Cause: in `read_gse179847_processed_matrix`, when no `" beta"`-suffixed
   columns are found the post5 code unconditionally delegates to
   `read_geo_series_matrix`, which requires a **tab-delimited** file containing
   `!series_matrix_table_begin`. A conventional **comma-delimited** processed CSV
   whose sample titles sit directly in the header therefore raised
   `ValueError: No series matrix table was found`.

   The handoff notes referred to a "CSV compatibility hotfix" that
   "distinguishes a genuine tab-delimited GEO Series Matrix metadata header from
   a conventional processed CSV whose sample titles appear directly in the
   header". That fix existed only in the lost WSL tree and is in no archive, so
   it was rewritten.

   The fix inspects the header before delegating: a genuine series matrix
   collapses to a single comma-field beginning with `!`, and only that case is
   delegated. A conventional processed CSV is handled inline, treating every
   column after the first as a sample column.
   Backup: `src/ccid/orthogonal_validation.py.pre_csvfix.bak`.

3. Result: **72 passed**, which matches the originally reported run exactly.

## Verification

| Check | Before | After |
|---|---|---|
| Full suite | 2 failed, 70 passed | **72 passed** |
| `tests/test_orthogonal_validation_post5.py` | 2 failed | 2 passed |
| `__version__` | 0.9.0.post4 | 0.9.0.post5 |
| `Sample_description` in orthogonal_validation.py | 0 | 4 |

Environment used for verification: Linux, Python 3.10, pandas 2.3.3, numpy 2.2.6.
The SciPy precision-loss warnings in `test_signature.py` and the bystander tests
are the same ones reported originally and are not failures.

## One decision you need to make

The tree now contains a fix that is **not** part of post5, so the version string
`0.9.0.post5` is no longer strictly accurate — which is the same class of problem
as the mislabelled tag above. The version was deliberately left unchanged because
it is what the manuscript cites.

Recommended: bump to `0.9.1`, tag it, push it, and mint a Zenodo DOI, then update
the manuscript's data availability statement to the new tag and DOI. That is a
one-line change in `src/ccid/__init__.py` and `pyproject.toml`.

## What is still missing

Neither folder contains any results or raw data:

- `outputs/`, `candidates/`, `manifests/`, `config/` are all empty
- no `*_all_loci.csv.gz`, no `*_donor_effects.csv.gz`
- no raw GEO matrices
- `git log --all --diff-filter=A` confirms results were **never** committed, so
  there is nothing to recover from git history

The FAZ 0 diagnostic therefore still requires rerunning the discovery chain.
See `START_HERE.md` section 3.
