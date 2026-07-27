# v0.7.1 — Positive-Control, Platform-Aware and Unique-Region Validation

This stage closes the main methodological gaps left by v0.7.0. It keeps the frozen O1/O2/O3 candidate hierarchy unchanged and adds six prespecified analyses:

1. **Published-target positive controls** for `PDE4C`, `FHL2`, `ELOVL2`, `KLF14`, and `TEAD1`.
2. **Platform-aware T-cell analysis** in which R1/R2 and R3 remain on their native probe spaces.
3. **Available-pair meta-analysis** so a CpG can contribute with one, two, or three T-cell pairs.
4. **Unique genomic regions** clustered once across the candidate union, with candidate-set membership flags.
5. **Coverage-bias diagnostics** comparing measured and missing candidates.
6. **Multiplicity-corrected and cross-cell-type empirical nulls**.

## Scientific guardrails

- The dCas9-DNMT3A/3L datasets support perturbation responsiveness, not direct proof of the desired TET-like hypomethylating intervention.
- The five positive-control genes are identified through EPICv2 gene annotation.
- Distance-decay plots use the **median coordinate of gene-annotated probes as an annotation-derived anchor**. This is not the unpublished sgRNA cut site.
- T-cell fixed/random-effect quantities use empirical pair-variance weighting because each pair has one targeted-minus-control delta per CpG. They are sensitivity summaries, not confirmatory clinical meta-analysis.
- Unique-region scores are prioritisation evidence, not causal proof.
- The primary regional table contains only regions with at least two candidate CpGs. Singletons remain available separately.

## Required files

```text
data/raw/GSE269758_T_cell_hyper_beta_detp.csv.gz
data/raw/GSE269758_T_cell_R3_Hyperlibrary_beta_detp.csv.gz
data/processed/msc_beta.csv.gz
data/processed/msc_metadata.csv
data/reference/epicv2_probe_annotation_hg38.csv.gz
results/real/o3_robustness/lodo_consensus_3of3_ids.csv
results/real/o3_robustness/lodo_consensus_at_least_2of3.csv
results/real/o3_validation/o3_tier1_candidates.csv.gz
results/real/safe_osk_target_signature_O1O2.csv.gz
results/real/safe_osk_target_signature_O1O2_all_loci.csv.gz
```

Check them:

```bash
ls -lh \
  data/raw/GSE269758_T_cell_hyper_beta_detp.csv.gz \
  data/raw/GSE269758_T_cell_R3_Hyperlibrary_beta_detp.csv.gz \
  data/processed/msc_beta.csv.gz \
  data/processed/msc_metadata.csv \
  data/reference/epicv2_probe_annotation_hg38.csv.gz
```

## Main command

```bash
mkdir -p results/real/bystander_validation_v071

python scripts/04d_validate_bystander_v071.py \
  --tcell-raw-inputs \
data/raw/GSE269758_T_cell_hyper_beta_detp.csv.gz,data/raw/GSE269758_T_cell_R3_Hyperlibrary_beta_detp.csv.gz \
  --msc-beta data/processed/msc_beta.csv.gz \
  --msc-metadata data/processed/msc_metadata.csv \
  --stable-core results/real/o3_robustness/lodo_consensus_3of3_ids.csv \
  --consensus-2of3 results/real/o3_robustness/lodo_consensus_at_least_2of3.csv \
  --tier1 results/real/o3_validation/o3_tier1_candidates.csv.gz \
  --frozen-signature results/real/safe_osk_target_signature_O1O2.csv.gz \
  --discovery-all-loci results/real/safe_osk_target_signature_O1O2_all_loci.csv.gz \
  --probe-annotation data/reference/epicv2_probe_annotation_hg38.csv.gz \
  --output-dir results/real/bystander_validation_v071 \
  --min-abs-delta 0.02 \
  --fdr 0.10 \
  --null-iterations 1000 \
  --null-bins 10 \
  --cluster-gap-bp 2000 \
  --seed 20260721
```

## Platform-aware T-cell model

Each raw T-cell source is read separately. Pair-level effects are computed on the native platform:

```text
pair Δβ = targeted pair − scramble pair
```

CpGs are then combined across all available R1, R2 and R3 pairs. A CpG present only on the R3 platform is retained with `n_pairs=1`; it is not removed by an R1/R2/R3 intersection filter.

Outputs contain:

- available-case mean and median;
- empirical-variance weighted fixed effect;
- DerSimonian–Laird random-effect sensitivity estimate;
- `tau2`, Cochran Q, available-pair count and sign consistency;
- platform and pair manifests.

## Positive controls

The package locates EPICv2 probes annotated to the five published multiplex target genes and creates:

- gene-associated probe summaries;
- ±5 kb annotation-anchor neighbourhoods;
- 0–250, 251–500, 501–1000, 1001–2000 and 2001–5000 bp distance bins;
- T-cell and MSC effect and consistency summaries.

These analyses verify whether the public perturbation datasets exhibit expected local responsiveness. They do not reconstruct the exact guide-centred decay curve because guide cut coordinates are not provided in the current inputs.

## Unique regions

Candidate CpGs from all four sets are merged before clustering. Each coordinate appears once and receives Boolean membership columns such as:

```text
in_stable_core_3of3
in_consensus_at_least_2of3
in_o3_tier1
in_frozen_500
```

Main regions require at least two candidate CpGs. Singleton loci are retained in the all-region table but are not counted as multi-CpG regions.

## Multiple testing

The matched-null summary reports:

```text
empirical_p
bh_q_within_celltype
bh_q_global
```

The cross-cell-type null matches random CpG sets on common measurement, baseline methylation and discovery-effect magnitude. It tests:

- sign concordance;
- absolute Spearman correlation;
- mean absolute T-cell response;
- mean absolute MSC response.

## Output files

```text
results/real/bystander_validation_v071/
├── bystander_validation_v071_summary.json
├── bystander_validation_v071_manifest.json
├── candidate_probe_coverage_platform_aware.csv
├── candidate_set_summary_platform_aware.csv
├── candidate_bystander_effects_platform_aware.csv.gz
├── cross_celltype_candidate_effects_platform_aware.csv.gz
├── cross_celltype_replication_platform_aware.csv
├── tcell_platform_meta_effects.csv.gz
├── tcell_platform_pair_effects.csv.gz
├── tcell_platform_pair_manifest.csv
├── tcell_platform_qc.csv
├── msc_effects_v071.csv.gz
├── positive_control_probe_effects.csv.gz
├── positive_control_gene_summary.csv
├── positive_control_distance_decay.csv
├── unique_region_membership.csv.gz
├── unique_region_scores_all.csv
├── unique_region_scores_multicpg.csv
├── coverage_bias_per_candidate.csv.gz
├── coverage_bias_summary.csv
├── matched_null_tests_v071.csv.gz
├── matched_null_summary_bh.csv
├── cross_celltype_null_tests.csv.gz
├── cross_celltype_null_summary.csv
└── figures/
    ├── positive_control_decay_tcell.png
    ├── positive_control_decay_msc.png
    ├── matched_null_bh_summary.png
    ├── unique_multicpg_regions.png
    ├── coverage_bias.png
    └── cross_celltype_empirical_null.png
```

## First checks

```bash
cat results/real/bystander_validation_v071/bystander_validation_v071_summary.json

column -s, -t \
  < results/real/bystander_validation_v071/positive_control_gene_summary.csv

column -s, -t \
  < results/real/bystander_validation_v071/matched_null_summary_bh.csv

column -s, -t \
  < results/real/bystander_validation_v071/cross_celltype_null_summary.csv

column -s, -t \
  < results/real/bystander_validation_v071/coverage_bias_summary.csv
```

Top multi-CpG regions:

```bash
python - <<'PY'
import pandas as pd

path = "results/real/bystander_validation_v071/unique_region_scores_multicpg.csv"
df = pd.read_csv(path)
print("Unique multi-CpG regions:", len(df))
print(df.head(40).to_string(index=False))
PY
```

## Decision rule

The bystander/editability phase can be frozen after this run when:

- positive-control genes show measurable local editor response;
- the MSC stable-core enrichment remains significant after BH correction or is clearly labelled nominal/supportive;
- platform-aware T-cell results establish whether the earlier negative result persists;
- coverage-bias tests do not reveal a severe rank/effect imbalance, or the imbalance is explicitly adjusted for;
- cross-cell-type concordance is interpreted against its empirical matched null;
- the manuscript reports unique multi-CpG regions, not duplicated candidate-set-specific singleton records.
