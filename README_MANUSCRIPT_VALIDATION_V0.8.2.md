# CCID-Epigenome v0.8.2: Manuscript-Readiness Sensitivity and Negative-Control Validation

v0.8.2 adds a final computational robustness layer before manuscript drafting is frozen. It tests whether the candidate ranking is stable under plausible changes in evidence weights and measurement uncertainty, whether any single evidence component dominates the result, and whether the observed aggregation of evidence exceeds a marginal-preserving permutation null.

## What the suite tests

1. **Weight sensitivity:** 1,000 log-normal perturbations of the evidence weights, renormalized on every iteration.
2. **Component uncertainty:** 1,000 Monte Carlo perturbations of component scores with bounded noise.
3. **Leave-one-component-out ablation:** every evidence component is removed once and ranks are recomputed.
4. **Ranking-specificity null:** each component is independently permuted across CpGs while preserving its marginal distribution and missingness pattern.
5. **Evidence redundancy:** pairwise Spearman correlations, high-correlation pairs and effective rank of the component correlation matrix.
6. **Candidate-set enrichment:** Fisher exact tests for stable-core, consensus, O3, frozen-signature and multi-CpG-region membership among robust top-k candidates.

Missing evidence remains unavailable. It is never imputed as zero.

## Inputs

The default input is the v0.8.1 candidate table:

```text
results/real/regulatory_safety_v081/candidate_cpg_regulatory_safety_v081.csv.gz
```

When `results/real/regulatory_safety_v081_reviewed/` exists, the script automatically prefers the reviewed candidate table. A specific file can always be supplied with `--candidate-file`.

## Run

```bash
python scripts/04g_validate_manuscript_readiness.py \
  --candidate-file results/real/regulatory_safety_v081/candidate_cpg_regulatory_safety_v081.csv.gz \
  --config-json configs/manuscript_validation_v082.json \
  --output-dir results/real/manuscript_validation_v082
```

For a faster diagnostic run:

```bash
python scripts/04g_validate_manuscript_readiness.py \
  --candidate-file results/real/regulatory_safety_v081/candidate_cpg_regulatory_safety_v081.csv.gz \
  --iterations 200 \
  --output-dir results/real/manuscript_validation_v082_quick
```

Use at least 1,000 iterations for manuscript values. Re-run with 5,000 iterations for final empirical p-values if runtime permits.

## Main outputs

```text
candidate_rank_stability_v082.csv.gz
weight_perturbation_iterations_v082.csv.gz
component_uncertainty_iterations_v082.csv.gz
ablation_summary_v082.csv
ablation_candidate_rank_shifts_v082.csv.gz
ranking_specificity_null_summary_v082.csv
ranking_specificity_null_iterations_v082.csv.gz
component_correlation_v082.csv
component_redundancy_summary_v082.csv
candidate_set_enrichment_v082.csv
manuscript_ready_robust_shortlist_v082.csv
manuscript_validation_v082_summary.json
manuscript_validation_v082_manifest.json
```

Figures:

```text
figures/top_robust_candidates.png
figures/weight_rank_stability.png
figures/ablation_rank_impact.png
figures/ranking_specificity_null.png
figures/component_correlation.png
```

## Interpretation guide

A candidate is labelled `stable_primary_shortlist` when it appears in the primary top-k in at least 70% of both weight-perturbation and component-uncertainty iterations. `stable_extended_shortlist` indicates stability only within the extended top-k. Other candidates are `sensitivity_dependent`.

Recommended manuscript evidence:

- median Spearman correlation between perturbed and baseline rankings;
- median top-k Jaccard overlap;
- empirical p-value for the top-k permutation-null statistic;
- the three most influential ablations;
- component effective rank and any pairs with |rho| >= 0.80;
- robust shortlist inclusion frequencies and rank IQRs.

These tests establish ranking robustness, not biological causality, clinical safety or sequence-level editability. Literature/safety completion and wet-lab validation remain separate requirements.
