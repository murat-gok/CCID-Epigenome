# v0.7.0 — T-cell/MSC Bystander and Regional Responsiveness Validation

This stage tests the frozen safe-reprogramming candidate hierarchy in two public primary-cell dCas9-DNMT3A/3L datasets:

- T cells: three targeted-versus-scramble replicate pairs.
- Mesenchymal stromal cells (MSCs): three targeted-versus-scramble donor pairs.

The public experiments multiplexed guides against five age-hypermethylated regions (`PDE4C`, `FHL2`, `ELOVL2`, `KLF14`, `TEAD1`). The package therefore evaluates whether the frozen candidate CpGs participate in reproducible genome-wide bystander responses. It does **not** relabel a DNMT3A response as proof of the desired TET-like hypomethylating intervention.

## Scientific guardrails

1. Candidate sets are frozen before this analysis.
2. Targeted and scramble samples are paired by replicate/donor.
3. The independent unit is the replicate/donor pair, not each CpG or sample column.
4. Raw locus-level t-test/FDR values are descriptive with only three pairs.
5. Matched empirical nulls are preferred for candidate-set enrichment.
6. Regional scores mean **regional bystander responsiveness**, not causal editability.
7. A direct-target-gene flag identifies clusters annotated to the five published target genes.
8. A future dCas9-TET1 dataset or wet-lab experiment is still required for directional validation of hypomethylating protocols.

## Prerequisite files

The v0.3 preprocessing stage should already have produced:

```text
data/processed/tcell_beta.csv.gz
data/processed/tcell_metadata.csv
data/processed/msc_beta.csv.gz
data/processed/msc_metadata.csv
```

Check them:

```bash
ls -lh \
  data/processed/tcell_beta.csv.gz \
  data/processed/tcell_metadata.csv \
  data/processed/msc_beta.csv.gz \
  data/processed/msc_metadata.csv
```

If missing, run:

```bash
python scripts/02c_preprocess_bystander_methylation.py
```

Candidate hierarchy:

```text
results/real/o3_robustness/lodo_consensus_3of3_ids.csv
results/real/o3_robustness/lodo_consensus_at_least_2of3.csv
results/real/o3_validation/o3_tier1_candidates.csv.gz
results/real/safe_osk_target_signature_O1O2.csv.gz
results/real/safe_osk_target_signature_O1O2_all_loci.csv.gz
```

## Recommended probe annotation

Regional clustering requires genomic coordinates. Export the EPICv2 hg38 annotation through Bioconductor:

```bash
mkdir -p data/reference
Rscript scripts/04c_export_epicv2_annotation.R \
  data/reference/epicv2_probe_annotation_hg38.csv.gz
```

The T-cell preprocessing retains loci shared across its EPIC/EPICv2 files, so the EPICv2 annotation provides coordinates for the cross-platform loci that are also represented on EPICv2. Candidate probes not present in the annotation remain in probe-level analyses but are omitted from regional clustering.

The main analysis can run without `--probe-annotation`; regional tables will then be empty and the summary will state that annotation was absent.

## Main command

```bash
mkdir -p results/real/bystander_validation

python scripts/04c_validate_bystander_editability.py \
  --tcell-beta data/processed/tcell_beta.csv.gz \
  --tcell-metadata data/processed/tcell_metadata.csv \
  --msc-beta data/processed/msc_beta.csv.gz \
  --msc-metadata data/processed/msc_metadata.csv \
  --stable-core results/real/o3_robustness/lodo_consensus_3of3_ids.csv \
  --consensus-2of3 results/real/o3_robustness/lodo_consensus_at_least_2of3.csv \
  --tier1 results/real/o3_validation/o3_tier1_candidates.csv.gz \
  --frozen-signature results/real/safe_osk_target_signature_O1O2.csv.gz \
  --discovery-all-loci results/real/safe_osk_target_signature_O1O2_all_loci.csv.gz \
  --probe-annotation data/reference/epicv2_probe_annotation_hg38.csv.gz \
  --output-dir results/real/bystander_validation \
  --min-abs-delta 0.02 \
  --fdr 0.10 \
  --null-iterations 1000 \
  --null-bins 10 \
  --cluster-gap-bp 2000 \
  --seed 20260721
```

## Candidate sets

| Set | Interpretation |
|---|---|
| `stable_core_3of3` | 50 CpGs selected in all three LODO rotations; primary set |
| `consensus_at_least_2of3` | 101 CpGs selected in at least two rotations |
| `o3_tier1` | O3-supported Tier 1 sensitivity set |
| `frozen_500` | Complete frozen O1/O2 discovery signature |

## Statistical outputs

For each cell type, pair-level effects are:

```text
editor Δβ = mean(targeted samples in pair) − mean(scramble controls in pair)
```

The package reports:

- mean and median editor Δβ;
- pair-level SD/SE and sign consistency;
- descriptive one-sample t-test and BH-FDR;
- significant-effect fraction at the chosen Δβ/FDR thresholds;
- coverage of each frozen candidate set;
- cross-cell-type sign concordance, Pearson and Spearman replication;
- relation between DNMT3A direction and the frozen safe-reprogramming target direction;
- matched random-set empirical tests.

Matched null sets use baseline-methylation bins and, when `--discovery-all-loci` is provided, absolute discovery-effect bins. This reduces bias from different initial methylation levels and target-effect magnitudes.

## Regional scoring

Annotated candidate CpGs on the same chromosome are grouped when consecutive probes are no more than `--cluster-gap-bp` apart. A regional responsiveness score combines:

- mean absolute T-cell/MSC editor response;
- replicate/donor sign consistency;
- cross-cell-type sign concordance;
- measured-probe fraction.

The file retains the explicit label:

```text
regional_bystander_responsiveness_not_direct_editability
```

Regions overlapping `PDE4C`, `FHL2`, `ELOVL2`, `KLF14`, or `TEAD1` are flagged because they may represent direct multiplex targets rather than distal bystanders.

## Output files

```text
results/real/bystander_validation/
├── candidate_probe_coverage.csv
├── candidate_set_summary.csv
├── candidate_bystander_effects.csv.gz
├── direction_compatibility.csv
├── cross_celltype_candidate_effects.csv.gz
├── cross_celltype_replication.csv
├── empirical_null_tests.csv.gz
├── regional_candidate_clusters.csv.gz
├── regional_editability_scores.csv
├── sample_pair_manifest.csv
├── sample_qc.csv
├── tcell_cpg_effects.csv.gz
├── msc_cpg_effects.csv.gz
├── tcell_pair_effects.csv.gz
├── msc_pair_effects.csv.gz
├── bystander_validation_summary.json
├── bystander_analysis_manifest.json
└── figures/
    ├── candidate_coverage.png
    ├── tcell_msc_effect_correlation.png
    ├── candidate_set_comparison.png
    ├── donor_consistency.png
    ├── regional_editability.png
    └── matched_null_comparison.png
```

## First result checks

```bash
cat results/real/bystander_validation/bystander_validation_summary.json

column -s, -t \
  < results/real/bystander_validation/candidate_probe_coverage.csv

column -s, -t \
  < results/real/bystander_validation/cross_celltype_replication.csv

column -s, -t \
  < results/real/bystander_validation/candidate_set_summary.csv
```

Compact empirical-null summary:

```bash
python - <<'PY'
import pandas as pd

path = "results/real/bystander_validation/empirical_null_tests.csv.gz"
df = pd.read_csv(path)
summary = (
    df.groupby(["candidate_set", "cell_type", "metric"], as_index=False)
      .agg(
          observed=("observed", "first"),
          null_mean=("null_value", "mean"),
          null_q95=("null_value", lambda x: x.quantile(0.95)),
          empirical_p=("empirical_p", "first"),
      )
)
print(summary.to_string(index=False))
PY
```

## Success criteria before moving forward

The phase supports advancement when the stable core or 2/3 consensus shows several of the following:

- high probe coverage in both cell types;
- larger median absolute editor effects than matched random sets;
- reproducible pair/donor direction;
- cross-cell-type sign concordance above null expectation;
- at least a subset of regions with high regional responsiveness;
- clear separation of direct-target-gene regions from distal candidate regions.

A negative result is still informative: it would indicate that the safe-reprogramming signature is not broadly responsive to this particular hypermethylating perturbation context and should not be used as editability evidence.
