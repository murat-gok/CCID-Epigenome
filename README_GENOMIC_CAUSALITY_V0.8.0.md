# v0.8.0 — Genomic Annotation, Biological Relevance and Causality Evidence Scoring

This stage freezes the v0.7.1 bystander analysis and converts the candidate CpGs and unique regions into a transparent experimental-prioritisation table.

## Scientific scope

The package integrates nine evidence components:

1. discovery effect/rank strength;
2. donor stability from LODO selection frequency;
3. held-out O3 support;
4. successful-versus-failed specificity;
5. successful/failed time-series consistency;
6. T-cell/MSC editor-responsiveness evidence;
7. unique-region support;
8. EPICv2 genomic context;
9. biological-relevance evidence.

It also reports a separate proximity-based safety-risk score. The combined outputs are **experimental priorities, not proof that editing a CpG causes rejuvenation or is safe**.

## Required inputs

```text
data/reference/epicv2_probe_annotation_hg38.csv.gz
data/reference/biological_evidence_seed_v080.csv
results/real/bystander_validation_v071/unique_region_membership.csv.gz
results/real/bystander_validation_v071/unique_region_scores_all.csv
results/real/o3_robustness/o3_successful_vs_failed.csv.gz
results/real/o3_robustness/time_series_cpg_consistency.csv.gz
```

The two O3 robustness files are optional at runtime. If they are absent, their components remain missing, the available weights are renormalised, and `evidence_completeness` decreases. Missing evidence is not silently replaced by zero.

## Main command

```bash
rm -rf results/real/genomic_causality_v080
mkdir -p results/real/genomic_causality_v080

python scripts/04e_score_genomic_causality.py \
  --v071-dir results/real/bystander_validation_v071 \
  --probe-annotation data/reference/epicv2_probe_annotation_hg38.csv.gz \
  --gene-evidence-seed data/reference/biological_evidence_seed_v080.csv \
  --successful-vs-failed results/real/o3_robustness/o3_successful_vs_failed.csv.gz \
  --time-series-consistency results/real/o3_robustness/time_series_cpg_consistency.csv.gz \
  --weights-json configs/causality_scoring_v080.json \
  --shortlist-size 100 \
  --output-dir results/real/genomic_causality_v080
```

## Optional literature-backed evidence

The bundled seed panel is intentionally labelled `canonical_seed` and is only a reproducible computational triage aid. Manuscript-level biological claims should use a traceable external table:

```bash
python scripts/04e_score_genomic_causality.py \
  ... \
  --external-gene-evidence data/reference/literature_gene_evidence.csv
```

Required columns:

```text
gene,category,relevance_weight,risk_weight
```

Recommended provenance columns:

```text
evidence_level,source_id,source_note
```

Weights must be between 0 and 1. Multiple files can be comma separated. The strongest evidence is used for scoring while all source rows remain available for audit.

Generic GMT files can be added for enrichment:

```bash
--gmt data/reference/go_terms.gmt,data/reference/reactome_terms.gmt
```

No proprietary GMT collection is bundled.

## Prespecified score

Default weights are stored in `configs/causality_scoring_v080.json`:

```text
discovery strength                 0.10
donor stability                   0.18
O3 support                        0.14
successful-vs-failed specificity  0.12
temporal consistency              0.12
editor responsiveness             0.10
regional support                  0.08
genomic context                   0.06
biological relevance              0.10
```

The score is renormalised over available components. The package always reports `evidence_completeness` so candidates with missing evidence can be distinguished from candidates with negative evidence.

The safety-adjusted priority is:

```text
causal evidence × proximity-risk penalty × completeness modifier
```

Risk flags are conservative proximity flags for manual review; they are not clinical safety classifications.

## Output files

```text
results/real/genomic_causality_v080/
├── candidate_cpg_annotation_causality_scores.csv.gz
├── candidate_region_annotation_causality_scores.csv
├── gene_evidence_scores.csv
├── biological_category_enrichment.csv
├── safety_flags.csv
├── evidence_component_summary.csv
├── top_candidate_shortlist.csv
├── genomic_causality_v080_summary.json
├── genomic_causality_v080_manifest.json
└── figures/
    ├── cpg_score_distribution.png
    ├── top_cpg_evidence_heatmap.png
    ├── evidence_vs_safety_risk.png
    ├── top_candidate_regions.png
    └── biological_category_enrichment.png
```

The enrichment figure is produced only when at least one category is testable.

## Main CpG columns

```text
causal_evidence_score
evidence_completeness
safety_risk_score
safety_adjusted_priority_score
causal_evidence_grade
experimental_priority_tier
directional_editor_hypothesis
```

`directional_editor_hypothesis` is based only on the sign of the frozen safe-reprogramming target:

- negative target Δβ → TET-like hypomethylating editor hypothesis;
- positive target Δβ → DNMT3A-like hypermethylating editor hypothesis.

It is not sequence-level editor selection.

## Priority tiers

- `P1`: highest experimental priority after safety penalty;
- `P2`: supported priority;
- `P3`: exploratory;
- `P4`: insufficient current evidence.

Region tiers use `R1`–`R4`. A high region score does not establish direct editability or causality.

## First checks

```bash
cat results/real/genomic_causality_v080/genomic_causality_v080_summary.json

column -s, -t \
  < results/real/genomic_causality_v080/evidence_component_summary.csv

python - <<'PY'
import pandas as pd

base = "results/real/genomic_causality_v080"
short = pd.read_csv(f"{base}/top_candidate_shortlist.csv")
regions = pd.read_csv(f"{base}/candidate_region_annotation_causality_scores.csv")

print("Top CpGs:")
print(short.head(30).to_string(index=False))
print("\nTop multi-CpG regions:")
print(regions[regions["n_candidate_cpgs"] >= 2].head(20).to_string(index=False))
PY
```

## Manuscript-safe interpretation

Safe wording:

> Candidates were prioritised by prespecified triangulation of donor stability, held-out validation, temporal specificity, cell-context editor responsiveness, genomic context and biological evidence, with separate safety-proximity penalties.

Unsafe wording:

> The highest-scoring CpGs are proven causal rejuvenation targets.
