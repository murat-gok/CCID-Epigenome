# Upgrade to v0.8.1

## Apply the patch

```bash
cd ~/ccid_epigenome
unzip -o "/mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.8.1.zip" -d ~/ccid_epigenome
```

## Reinstall and verify

```bash
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
python -c "import ccid; print(ccid.__version__)"
```

Expected version: `0.8.1`.

## Preserve earlier results

v0.8.1 reads the frozen v0.8.0 output directory and writes to a new directory. It does not overwrite v0.8.0 results.

```text
input:  results/real/genomic_causality_v080
output: results/real/regulatory_safety_v081
```

## Download regulatory references

```bash
python scripts/00_fetch_v081_regulatory_references.py \
  --output-dir data/reference/v081
```

## Run

```bash
python scripts/04f_complete_regulatory_safety.py \
  --v080-dir results/real/genomic_causality_v080 \
  --gencode-gtf data/reference/v081/gencode.v50.annotation.gtf.gz \
  --literature-evidence data/reference/literature_evidence_seed_v081.csv \
  --config-json configs/regulatory_safety_v081.json \
  --shortlist-size 100 \
  --output-dir results/real/regulatory_safety_v081
```

## Breaking semantic correction

In v0.8.0, absent biological or safety evidence could appear as a zero. v0.8.1 replaces this with explicit `not_reviewed`/`not_supplied` statuses and missing scores. Downstream scripts must not fill these missing values with zero.
