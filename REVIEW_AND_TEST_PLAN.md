# Revised Fikir B: final scientific additions and test plan

## Required additions before manuscript execution

1. **Cell-type alignment.** Use human dermal fibroblasts as the primary discovery system because the transient OSK/OSKM methylation and expression trajectories used to define the safe target are fibroblast based. Use T-cell and MSC perturbation datasets as external transfer tests, not as interchangeable training samples.
2. **Leakage-free safe-window definition.** Select identity/pluripotency/rejuvenation thresholds and safe days on discovery donors only. Freeze the rule before evaluating held-out donors or studies.
3. **Four primary Pareto objectives.** Optimize target-signature mismatch, safety risk, edit burden and predictive uncertainty. Treat identity retention, PAM/direction feasibility, chromatin accessibility, causal-evidence minimum and multiplex delivery limits as constraints.
4. **Editor terminology.** Call dCas9-TET1, dCas9-DNMT3A and CRISPRoff sequence-preserving epigenome editors. Do not call all of them transient or intrinsically reversible; persistence is locus/editor dependent and CRISPRoff requires a counter-editor for reversal.
5. **Evidence tiers rather than a binary causal filter.** Combine MR, meQTL/colocalization, mediation and perturbation support. Missing evidence is unknown, not negative. Report Tier-1 confirmatory and Tier-2/3 exploratory solutions separately.
6. **Robust decision rule.** Select protocols using calibrated uncertainty intervals or a worst-case/CVaR criterion; report interval coverage and width under donor/study/region shift.
7. **Physical feasibility layer.** The final unit is region + sgRNA + editor + direction + expected regional effect. Require PAM, chromatin, amplicon, off-target and multiplex compatibility before optimization.
8. **Claim boundary.** The computational endpoint is prioritization of protocols that recapitulate a safe partial-reprogramming signature. Biological rejuvenation remains a wet-lab hypothesis.

## Primary analysis

- Define a safe reprogramming window from expression data while retaining fibroblast identity and avoiding pluripotency activation.
- Derive a donor-paired methylation target vector from the safe window.
- Estimate editor-specific target and bystander effects with grouped train/calibration/test splits.
- Build an evidence-weighted, physically feasible sgRNA-editor candidate library.
- Solve the four-objective constrained inverse-design problem and return a Pareto front.
- Select a predeclared target-priority protocol plus safety-priority alternatives.

## Required tests

### Data and target validity

- Donor-held-out safe-window reproducibility.
- Study-held-out target-signature concordance.
- Sign consistency, effect-size confidence intervals and probe-quality filtering.
- Cell-composition/batch sensitivity where applicable.

### Surrogate validity

- Random CV versus grouped donor/study CV.
- Leave-one-editor-out and leave-one-cell-type-out validation.
- Genomic-region or chromosome-block holdout.
- RMSE/MAE and 80/90/95% interval coverage plus interval width.
- Calibration before and after conformal correction.

### Biological and negative controls

- Shuffled target signature.
- Direction-preserving matched random CpGs.
- Age-associated but non-safe-OSK CpGs.
- Non-targeting/scramble guides and inactive-editor controls where data allow.
- Alternative clock panels used only as secondary endpoints.

### Optimization benchmarking

- Exact enumeration or MILP on a prespecified tractable subset.
- NSGA-II, MOEA/D, binary differential evolution, greedy and random baselines; NSGA-III is optional if the final design expands beyond four primary objectives.
- Equal objective-evaluation budget for NSGA-II, MOEA/D, binary DE and Random, with Greedy/Exact reported using their actual structural counts; 30 independent stochastic seeds.
- Freeze one normalization range, one hypervolume reference point/sample cloud and one pooled nondominated IGD/IGD+ reference front per biological problem instance.
- Compute Exact-small once per instance and reuse the identical front across seed rows.
- Keep the zero-edit/no-intervention state outside candidate protocol fronts; report it in a separate comparator table.
- Hypervolume, IGD+, epsilon indicator, feasibility rate, runtime, objective-evaluation count and selected-set stability.
- Friedman omnibus followed by paired Wilcoxon-Holm and Vargha-Delaney A12 across resampled problem instances.

### Ablation and sensitivity

- Remove causal evidence, bystander model, uncertainty, identity constraint, functional annotation and nonlinear interactions one at a time.
- Vary maximum guides from 1 to 8.
- Vary evidence, safety and safe-window thresholds.
- Repeat with alternative safe-window definitions and target-construction rules.
- Bootstrap donors/studies and report guide-selection frequency/Jaccard stability.
- At least 1,000 target permutations for the final manuscript.

### External transfer

- Fibroblast discovery and validation across independent transient/full reprogramming datasets.
- T-cell and MSC epigenome-editing datasets as domain-shift tests.
- Do not pool cell types before quantifying transfer error and recalibrating uncertainty.
