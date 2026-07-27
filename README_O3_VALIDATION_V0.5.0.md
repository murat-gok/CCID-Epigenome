# O3 held-out validation — v0.5.0

This stage validates the frozen O1/O2 safe-OSK CpG ranking in donor O3 without reselecting candidates or retuning thresholds from O3.

## Primary design

The script reads the frozen O1/O2 matched-pair table and requires exactly one discovery stratum. For the current analysis this is:

- day: 10
- batch: exp2
- stage: endpoint
- O3 validation sample: `O3_transiently_reprogrammed_10days_exp2`
- matched baseline: `O3_negative_control_10days_exp2`

Wrong-day baselines, cross-donor baselines and failed-reprogramming samples are used only as diagnostic null controls.

## Install/update

```bash
cd ~/ccid_epigenome
unzip -o /mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.5.0.zip -d ~/ccid_epigenome
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
python -c "import ccid; print(ccid.__version__)"
```

Expected: 32 tests passed and version `0.5.0`.

## Run O3 validation

```bash
mkdir -p results/real/o3_validation

python scripts/04_validate_signature_o3.py \
  --beta data/processed/fibroblast_beta.csv.gz \
  --metadata data/processed/fibroblast_methylation_metadata.csv \
  --signature results/real/safe_osk_target_signature_O1O2.csv.gz \
  --discovery-pairs results/real/safe_osk_target_signature_O1O2_matched_pairs.csv \
  --all-loci results/real/safe_osk_target_signature_O1O2_all_loci.csv.gz \
  --output-dir results/real/o3_validation \
  --validation-donor O3 \
  --tier1-min-abs-delta 0.03 \
  --tier1-min-retention 0.25 \
  --tier1-max-retention 2.0 \
  --tier2-min-abs-delta 0.01 \
  --tier2-min-retention 0.10 \
  --permutations 1000 \
  --random-sets 1000 \
  --random-effect-bins 10 \
  --seed 20260721 \
  --chunksize 10000
```

The `--all-loci` argument can be omitted; the script will try to infer the sibling `_all_loci.csv.gz` file.

## Main outputs

- `o3_cpg_validation.csv.gz`: all 500 frozen candidates with O3 effect and tier
- `o3_tier1_candidates.csv.gz`: same direction, |O3 Δβ| ≥ 0.03 and effect retention 0.25–2.0
- `o3_tier2_candidates.csv.gz`: same direction with weaker or overshooting effect, but |O3 Δβ| ≥ 0.01 and retention ≥ 0.10
- `o3_tier3_candidates.csv.gz`: reversed, missing or weak candidates
- `o3_validation_summary.json`: primary metrics and interpretation safeguards
- `o3_null_tests.csv`: permutation, effect-matched random CpG sets, wrong-day, cross-donor and failed-reprogramming controls
- `o3_topk_metrics.csv`: top-50/100/250/500 sign agreement and correlations
- `o3_sample_qc.csv`: genome-wide beta-distribution QC
- `o3_primary_pair.csv`: frozen O3 primary comparison
- `o3_auxiliary_controls.csv`: diagnostic wrong-day, cross-donor and failed-reprogramming controls
- `o3_frozen_analysis_manifest.json`: hashes and frozen thresholds
- `figures/`: correlation, distribution, top-k, sample-QC and null figures

## Interpretation rule

O3 is one independent held-out donor. Tier labels and locus-level correlations are descriptive prioritisation results, not population-level confirmation. Correlated CpGs invalidate naïve independence assumptions; empirical permutation and effect-matched random-set null tests are preferred over exact CpG-level p-values.
