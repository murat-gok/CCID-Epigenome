# CCID-Epigenome

> **v0.9.0.post4 hotfix:** GSE179847 extended validation now downloads the deposited 2-GB processed beta matrix, streams only the frozen CpGs, uses the small Series Matrix only for metadata, resumes interrupted downloads, and fails clearly on metadata-only input. See `README_HOTFIX_V0.9.0.POST4.md`.

> **v0.9.0.post3 hotfix:** ENCODE rE2G files are now streamed line by line, preventing WSL out-of-memory termination and correctly using the `TargetGene` and `Score` columns. See `README_HOTFIX_V0.9.0.POST3.md`.

> **v0.9.0.post2 hotfix:** GPL6883 files returned by NCBI GEO are now parsed as SOFT documents by extracting the delimited platform table. See `README_HOTFIX_V0.9.0.POST2.md`.

> **v0.9.0.post1 hotfix:** the retired GPL6883 `annot/` URL was replaced with NCBI GEO's supported Accession Display platform-table endpoint. See `README_HOTFIX_V0.9.0.POST1.md`.

**Current release: v0.9.0.post4**

**Causality-Constrained Inverse Design of Minimal Epigenome Edits that Recapitulate Safe Partial Reprogramming**

## v0.9.0 orthogonal validation

The v0.9.0 workflow evaluates the frozen 16-CpG panel in independent transient-reprogramming, fibroblast-aging and matched methylation–expression datasets. It also provides optional fibroblast rE2G overlap and SpCas9 sequence-feasibility analyses. Original ranking weights are not refit. See `README_ORTHOGONAL_VALIDATION_V0.9.0.md`.

Core run:

```bash
bash scripts/run_orthogonal_validation_v090.sh
```

Optional 479-sample multidonor lifespan run:

```bash
RUN_EXTENDED=1 bash scripts/run_orthogonal_validation_v090.sh
```


This repository is a runnable research scaffold for the revised Fikir B. It separates biological target construction, causal/feasibility filtering, editor-effect modeling, constrained multi-objective optimization and paper-level validation.

## What is implemented

- Leakage-aware safe-window scoring from expression data.
- Donor-paired methylation target-signature construction.
- Multi-evidence causal scoring that does not treat missing evidence as absence.
- sgRNA/editor physical-feasibility checks.
- Editor/bystander effect surrogate with group-aware split and conformal uncertainty.
- Four-objective constrained optimization: target mismatch, safety, edit count and uncertainty.
- NSGA-II, MOEA/D, binary differential evolution, random, greedy and exact-small baselines.
- Seed-independent hypervolume, IGD and IGD+ with frozen per-instance normalization and a common reference front.
- Separate no-intervention comparator; zero-edit solutions are excluded from candidate protocol fronts.
- Prespecified minimum edit count and minimum relative target-improvement feasibility constraints.
- Objective-evaluation counts and runtime for optimizer benchmarking.
- Planted-set recovery, seed stability, permutation, multiplex sensitivity and ablation tests.
- Friedman, Wilcoxon-Holm and Vargha-Delaney A12 utilities for full paper benchmarks.

The included demonstration is **synthetic**. It tests software and experimental logic; it is not evidence of biological rejuvenation.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate       # Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
export PYTHONPATH="$PWD/src"    # Windows PowerShell: $env:PYTHONPATH="$PWD/src"
```

## 1. Verify the package

```bash
pytest -q
python -m compileall -q src scripts tests
```

## 2. Run the synthetic end-to-end validation

```bash
python scripts/run_synthetic_demo.py --output results --seed 7
```

Outputs include benchmark/ablation/stability CSV files and Pareto/effect figures under `results/`.

## 3. Download public data

List the configured datasets:

```bash
python scripts/01_download_public_data.py --list
```

Download the smaller starting set:

```bash
python scripts/01_download_public_data.py \
  --datasets gse165177_rna,gse269758_tcell_bystander,gse285251_msc_bystander \
  --output data/raw
```

The large MPTR methylation matrix can be downloaded separately:

```bash
python scripts/01_download_public_data.py \
  --datasets gse165179_mptr_methylation \
  --output data/raw
```

## 4. Define the safe partial-reprogramming window

Prepare an expression matrix (genes × samples) and metadata with `sample_id, donor_id, day, condition`. Choose discovery donors before running:

```bash
python scripts/02_define_safe_window.py \
  --expression data/processed/fibroblast_expression.csv.gz \
  --metadata data/processed/fibroblast_expression_metadata.csv \
  --discovery-donors D1,D2,D3 \
  --identity-floor -0.50 \
  --pluripotency-ceiling 0.50 \
  --output results/real/safe_window_scores.csv
```

Freeze thresholds and selected days before held-out-donor validation.

## 5. Build the donor-paired methylation target signature

```bash
python scripts/03_build_target_signature.py \
  --beta data/processed/fibroblast_beta.csv.gz \
  --metadata data/processed/fibroblast_methylation_metadata.csv \
  --safe-window results/real/safe_window_scores.csv \
  --baseline-condition baseline \
  --min-abs-delta 0.03 \
  --fdr 0.10 \
  --output results/real/safe_OSK_target_signature.csv
```

The manuscript primary endpoint should be weighted mismatch to this frozen signature. Clock changes are secondary outcomes.


## 6A. Robust held-out-donor and time-series validation

After `04_validate_signature_o3.py`, run the prespecified robustness package:

```bash
python scripts/04b_validate_signature_robustness.py \
  --beta data/processed/fibroblast_beta.csv.gz \
  --metadata data/processed/fibroblast_methylation_metadata.csv \
  --signature results/real/safe_osk_target_signature_O1O2.csv.gz \
  --discovery-pairs results/real/safe_osk_target_signature_O1O2_matched_pairs.csv \
  --all-loci results/real/safe_osk_target_signature_O1O2_all_loci.csv.gz \
  --output-dir results/real/o3_robustness \
  --null-iterations 1000 \
  --seed 20260721
```

This stage adds successful-versus-failed, multi-baseline, global-centred, direction-balanced null, leave-one-donor-out and time-series analyses. See `README_ROBUSTNESS_V0.6.0.md`.


## 6C. Positive-control, platform-aware and unique-region validation

After v0.7.0, run the closing bystander validation stage:

```bash
python scripts/04d_validate_bystander_v071.py \
  --probe-annotation data/reference/epicv2_probe_annotation_hg38.csv.gz \
  --output-dir results/real/bystander_validation_v071 \
  --null-iterations 1000 \
  --seed 20260721
```

This stage preserves the native T-cell R1/R2 and R3 probe spaces, adds published-target positive controls, unique multi-CpG regions, coverage-bias diagnostics, BH-adjusted matched-null results and cross-cell-type empirical nulls. See `README_BYSTANDER_VALIDATION_V0.7.1.md`.

## 6. Quantify reproducible bystander effects

Standardize T-cell or MSC methylation matrices and metadata with `sample_id, treatment, replicate`:

```bash
python scripts/04_analyze_bystander_data.py \
  --beta data/processed/tcell_beta.csv.gz \
  --metadata data/processed/tcell_metadata.csv \
  --age-loci data/processed/age_cpgs.txt \
  --permutations 1000 \
  --output-dir results/real/tcell_bystander
```

Repeat for MSC and any editor-specific dataset. Keep cell type and editor as explicit model features.


## 6B. T-cell/MSC bystander and regional responsiveness validation

After freezing the O3/LODO candidate hierarchy, validate the 50/101/422/500 CpG sets in the paired dCas9-DNMT3A T-cell and MSC datasets. First export EPICv2 hg38 probe annotation (recommended for regional clustering):

```bash
Rscript scripts/04c_export_epicv2_annotation.R \
  data/reference/epicv2_probe_annotation_hg38.csv.gz
```

Then run the paired cross-cell-type analysis:

```bash
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
  --null-iterations 1000 \
  --seed 20260721
```

The analysis uses replicate/donor-matched targeted-minus-scramble effects, cross-cell-type replication, baseline/target-effect matched null sets, target-direction compatibility, and optional regional clustering. Because these public datasets use a hypermethylating dCas9-DNMT3A/3L editor, the primary claim is **bystander responsiveness and regional prioritisation**, not direct proof that dCas9-TET1 can recreate the hypomethylating safe-reprogramming target. See `README_BYSTANDER_VALIDATION_V0.7.0.md`.


## 6C. Genomic annotation and causality-evidence prioritisation

After the v0.7.1 platform-aware bystander validation, run:

```bash
python scripts/04e_score_genomic_causality.py \
  --v071-dir results/real/bystander_validation_v071 \
  --probe-annotation data/reference/epicv2_probe_annotation_hg38.csv.gz \
  --gene-evidence-seed data/reference/biological_evidence_seed_v080.csv \
  --output-dir results/real/genomic_causality_v080
```

The score is a transparent experimental-prioritisation device, not causal proof. See `README_GENOMIC_CAUSALITY_V0.8.0.md`.

## 7. Build the causal and feasible candidate library

The guide table must contain guide/editor identity, PAM/direction/chromatin/amplicon feasibility and risk scores. The evidence table should contain MR, colocalization, mediation, perturbation, damage-clock and safe-reprogramming support where available.

```bash
python scripts/04_build_candidate_library.py \
  --guides data/processed/guide_candidates.csv \
  --evidence data/processed/region_evidence.csv \
  --merge-key region_id \
  --output results/real/candidate_library.csv
```

Tier-1 and exploratory candidates should be reported separately.

## 8. Train and calibrate the effect surrogate

Prepare observation features, regional effects, grouping labels and candidate features as NumPy arrays. Groups must encode study/donor batches, not random rows.

```bash
python scripts/05_train_effect_surrogate.py \
  --features data/processed/surrogate_X.npy \
  --effects data/processed/surrogate_Y.npy \
  --groups data/processed/surrogate_groups.npy \
  --candidate-features data/processed/candidate_X.npy \
  --test-groups 6,7 \
  --calibration-groups 5 \
  --output-dir results/real/surrogate
```

Report random-CV versus grouped-CV error and 80/90/95% empirical coverage. A model that is accurate only under random CV is not acceptable.

## 9. Optimize protocols

```bash
python scripts/06_optimize_protocols.py \
  --candidates results/real/candidate_library.csv \
  --effect-mean results/real/surrogate/candidate_effect_mean.npy \
  --effect-std results/real/surrogate/candidate_effect_std.npy \
  --target-signature results/real/safe_OSK_target_signature.csv \
  --min-edits 1 \
  --max-edits 4 \
  --min-relative-improvement 0.10 \
  --population 160 \
  --generations 250 \
  --seed 0 \
  --output-dir results/real/optimization_seed_0
```

Run 30 seeds with an identical objective-evaluation budget for NSGA-II, MOEA/D, binary DE and Random. Greedy and Exact-small are deterministic/structural baselines and their own evaluation counts are reported separately. The scaffold supplies all six implementations. Add a true MILP solver for the linearized subset if a licensed/open solver is available; Exact-small already provides a solver-independent correctness reference. The zero-edit state is written separately as `no_intervention_baseline.csv` and is never treated as a candidate protocol.

## 10. Statistical comparison

Create a long-format CSV with `instance, algorithm, hypervolume` (or another metric):

```bash
python scripts/07_statistical_tests.py \
  --results results/real/all_optimizer_runs.csv \
  --metric hypervolume \
  --higher-is-better \
  --output-dir results/real/statistics
```

Instances should be donor/study bootstraps, target definitions or candidate subsets—not repeated seeds alone.

## Paper-ready minimum test matrix

Use `configs/full_paper.json` and `REVIEW_AND_TEST_PLAN.md`. The minimum includes grouped and region holdouts, calibration, negative controls, 1,000 permutations, six ablations, guide-count sensitivity, 30 optimizer seeds, exact or MILP validation on a tractable subset and external T-cell/MSC transfer tests.

## Current sandbox result

The synthetic run recovered two of three planted edits. Grouped CV error was higher than random CV error, demonstrating why grouped validation is necessary. Nominal 90% uncertainty coverage fell below 90% under synthetic domain shift, correctly flagging calibration as a central risk. See `RUN_LOG.md` and `results/demo_summary.json`.

## 11. Equal-budget multi-optimizer benchmark

For each bootstrap/study-target instance, run the same candidate and target files with a unique `--instance` label:

```bash
python scripts/08_benchmark_optimizers.py \
  --candidates results/real/candidate_library.csv \
  --effect-mean results/real/surrogate/candidate_effect_mean.npy \
  --effect-std results/real/surrogate/candidate_effect_std.npy \
  --target-signature results/real/safe_OSK_target_signature.csv \
  --instance donor_bootstrap_001 \
  --seeds 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29 \
  --min-edits 1 \
  --max-edits 4 \
  --min-relative-improvement 0.10 \
  --population 160 \
  --generations 250 \
  --output results/real/benchmarks/donor_bootstrap_001.csv
```

The benchmark freezes one normalization range, one pooled nondominated reference front, one hypervolume reference point and one Monte Carlo sample cloud per instance. `Exact-small` is computed once and reused unchanged across seed rows. Random automatically receives the evolutionary objective-evaluation budget unless `--random-samples` is explicitly supplied. Each run produces the main CSV plus `_no_intervention.csv` and `_metric_metadata.json` sidecars.

Concatenate all instance CSVs and apply `07_statistical_tests.py`. Do not treat repeated random seeds as independent biological instances in the Friedman test. See `UPGRADE_V0.2.0.md` for replacement and verification commands.


## Real-data safe-window and target-signature stage

See `README_NEXT_STEP_V0.4.1.md` and `UPGRADE_V0.4.1.md`.

## v0.5.0: held-out O3 validation

Run `scripts/04_validate_signature_o3.py` after freezing the O1/O2 signature. See `README_O3_VALIDATION_V0.5.0.md` for the exact command and interpretation rules.

## O3 held-out validation (v0.5.0)

After building the frozen O1/O2 signature, run `scripts/04_validate_signature_o3.py` to test it in held-out donor O3. The script infers the frozen day/batch/stage from the matched-pair file, generates Tier 1–3 candidate lists, empirical null controls, sample QC and publication-ready figures. See `README_O3_VALIDATION_V0.5.0.md`.

## v0.8.1 regulatory, literature and safety completion

v0.8.1 corrects the distinction between unknown and negative evidence, adds GENCODE nearest-TSS and ENCODE cCRE context, supports traceable candidate-specific literature and authorised Cancer Gene Census inputs, and penalises unresolved safety uncertainty.

See [`README_REGULATORY_SAFETY_V0.8.1.md`](README_REGULATORY_SAFETY_V0.8.1.md).


## v0.8.2 manuscript-readiness sensitivity validation

v0.8.2 adds weight perturbation, component-uncertainty Monte Carlo, leave-one-component-out ablation, marginal-preserving permutation nulls, evidence-redundancy diagnostics and candidate-set enrichment. See [`README_MANUSCRIPT_VALIDATION_V0.8.2.md`](README_MANUSCRIPT_VALIDATION_V0.8.2.md).
