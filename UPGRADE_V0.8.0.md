# Upgrade to v0.8.0

v0.8.0 is applied on top of a working v0.7.1 project. It does not overwrite the v0.7.1 results.

## Apply patch

```bash
cd ~/ccid_epigenome
unzip -o \
  "/mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.8.0.zip" \
  -d ~/ccid_epigenome
```

## Install and verify

```bash
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
python -c "import ccid; print(ccid.__version__)"
```

Expected release:

```text
0.8.0
```

## Run

```bash
python scripts/04e_score_genomic_causality.py \
  --v071-dir results/real/bystander_validation_v071 \
  --probe-annotation data/reference/epicv2_probe_annotation_hg38.csv.gz \
  --gene-evidence-seed data/reference/biological_evidence_seed_v080.csv \
  --successful-vs-failed results/real/o3_robustness/o3_successful_vs_failed.csv.gz \
  --time-series-consistency results/real/o3_robustness/time_series_cpg_consistency.csv.gz \
  --weights-json configs/causality_scoring_v080.json \
  --output-dir results/real/genomic_causality_v080
```

See `README_GENOMIC_CAUSALITY_V0.8.0.md` for evidence-table and GMT extensions.
