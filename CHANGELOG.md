# CCID-Epigenome — changelog

- Corrects GSE179847 extended-data acquisition: the 63-KB Series Matrix is retained only for sample metadata, while the deposited ~2-GB processed beta CSV is downloaded for methylation values.
- Adds a streaming GSE179847 CSV reader that retains only the frozen CpG panel.
- Maps GSM columns to deposited sample titles when necessary and accepts direct unique-variable-name columns.
- Adds resumable `.part` downloads for large supplementary files.
- Replaces the previous `KeyError` with a clear error when a metadata-only matrix is supplied.
- Adds three regression tests for processed-matrix streaming, multidonor aging analysis, and metadata-only rejection.

# CCID-Epigenome v0.9.0.post3

## Memory-safe ENCODE rE2G overlap hotfix

- Streams large gzip-compressed ENCODE element-gene-link BED files instead of loading them into pandas.
- Restricts overlap testing to the frozen 16-CpG panel, preventing WSL out-of-memory `Killed` failures.
- Parses the header-aware `TargetGene` and combined `Score` columns; earlier generic BED assumptions incorrectly treated the element name/class columns as gene/score fields.
- Uses correct zero-based, half-open BED overlap semantics against one-based CpG coordinates.
- Adds regression coverage for ENCODE-style rE2G headers and gzip streaming.
- No candidate ranks, thresholds, or prior validation results are changed.

# CCID-Epigenome v0.9.0.post2

## GEO SOFT platform-table parser hotfix

- Parses NCBI GEO ``view=data&form=text`` platform files only inside the
  ``!platform_table_begin`` / ``!platform_table_end`` section.
- Retains support for ordinary plain or gzip-compressed TSV annotation files.
- Prevents pandas variable-field ``ParserError`` failures caused by SOFT metadata.
- Adds regression coverage using a realistic GPL6883-style SOFT document.

## 0.9.0.post1 — GEO platform download hotfix

- Replaced the retired GEO `GPL6883/annot/GPL6883.annot.gz` endpoint, which returns HTTP 403, with NCBI GEO's supported Accession Display `view=data&form=text` endpoint.
- The HumanRef-8 v3 platform table is now stored as `GPL6883_data.txt`.
- Updated the eQTM script and end-to-end runner to use the corrected platform-table path.
- No scientific scoring, frozen candidates, or analysis thresholds were changed.

# Changelog

## 0.9.0 — Orthogonal in-silico validation

- Freezes the 16-CpG v0.8.2 primary shortlist before external testing.
- Adds independent paired transient-reprogramming direction validation in GSE142439.
- Adds replicative, chronological and optional multidonor fibroblast-aging direction analyses.
- Adds matched primary-fibroblast CpG–expression association testing using GSE52025/GSE53243.
- Adds optional ENCODE-rE2G overlap analysis and SpCas9-NGG sequence-feasibility assessment.
- Integrates orthogonal evidence without refitting the original nine-component ranking.
- Adds six tests; complete suite now contains 64 passing tests.

## 0.8.2 — manuscript-readiness sensitivity and negative-control validation

- Adds Monte Carlo evidence-weight perturbation and candidate rank inclusion frequencies.
- Adds bounded component-uncertainty simulations with rank IQR and score uncertainty.
- Adds leave-one-component-out ablations and candidate-specific worst rank drops.
- Adds marginal-preserving independent-component permutation nulls for top-k score aggregation.
- Adds component Spearman correlations, effective-rank redundancy diagnostics and high-correlation flags.
- Adds robust-set and multi-CpG enrichment tests with BH correction.
- Produces a manuscript-ready robust shortlist, audit manifest, ten tables and five figures.
- Adds four regression/integration tests; full suite: 58 tests.

## 0.8.1.post1 — pandas 3.x linked-gene hotfix

- Fixes a `TypeError` when missing linked-gene values remain numeric `NaN` scalars under pandas 3.x.
- Aggregates region genes value-by-value and excludes `NaN`, `pd.NA`, `None` and null-like strings.
- Adds a regression test for pandas 3.x missing scalar handling.
- Scientific scoring logic and reference inputs are unchanged.

## 0.8.0 — genomic annotation, biological relevance and causality-evidence scoring

- Adds transparent CpG-, region- and gene-level evidence scoring over frozen v0.7.1 candidates.
- Integrates discovery strength, LODO stability, held-out O3 specificity, time-series consistency, T-cell/MSC responsiveness, unique-region support and EPICv2 context.
- Separates biological relevance, evidence completeness and conservative safety-proximity penalties.
- Adds a configurable heuristic seed panel and support for traceable external evidence CSVs and generic GMT enrichment.
- Adds explicit TET-like versus DNMT3A-like directional hypotheses without claiming sequence-level feasibility.
- Adds P1–P4 CpG and R1–R4 region priorities, audit manifests, seven tables and up to five figures.
- Adds three regression/integration tests; full suite: 48 tests.

## 0.7.1 — positive controls, platform-aware T-cell meta-analysis and unique regions

- Preserves T-cell R1/R2 and R3 native probe spaces instead of enforcing a three-pair intersection.
- Adds available-pair empirical-variance fixed/random-effect sensitivity summaries.
- Adds annotation-derived positive controls for PDE4C, FHL2, ELOVL2, KLF14 and TEAD1.
- Adds annotation-anchor distance-decay summaries with an explicit guide-coordinate guardrail.
- Clusters the candidate union once, attaches candidate-set membership flags, and separates singleton loci from multi-CpG regions.
- Adds measured-versus-missing candidate coverage-bias tests.
- Adds BH q-values within cell type and globally for matched-null results.
- Adds common-measured, baseline/effect-matched cross-cell-type empirical nulls.
- Adds six figures, machine-readable manifests and five regression/integration tests.
- Full suite: 45 tests.

## 0.7.0 — paired T-cell/MSC bystander responsiveness validation

- Adds paired targeted-minus-scramble analysis for three T-cell and three MSC replicate/donor pairs.
- Evaluates frozen 50/101/422/500 CpG candidate tiers separately.
- Adds probe coverage, pair-level effects, descriptive FDR, sign consistency and cross-cell-type replication.
- Adds explicit target-direction versus DNMT3A-direction compatibility categories.
- Adds baseline-methylation and discovery-effect matched empirical random-set tests.
- Adds optional EPICv2 hg38 annotation export and regional candidate clustering.
- Flags the five published multiplex target genes to distinguish potential direct targets from distal bystanders.
- Labels regional scores as bystander responsiveness rather than causal/direct editability.
- Adds six paper-ready figures, machine-readable summaries and a frozen analysis manifest.
- Adds four regression/integration tests; the full suite now contains 40 tests.

## 0.6.0 — robustness, LODO and time-series validation

- Adds direct O3 successful-versus-failed comparison in the frozen stratum.
- Adds exact-, other-day-, all-baseline-, and failed-control sensitivity analyses.
- Reports raw and genome-wide-median-centred metrics.
- Adds direction-and-effect-matched empirical nulls for raw and centred effects.
- Adds leakage-safe leave-one-donor-out rotations that rebuild candidates using only non-held-out donors.
- Adds rotation overlap and three-way consensus reporting.
- Adds exact-baseline time-series consistency for successful, intermediate and failed samples.
- Adds five paper-ready figures and a frozen robustness manifest.
- Adds four regression/integration tests; the full suite now contains 36 tests.


## 0.4.1 — leakage-free target-signature stage

- Fits safe-window score standardization on O1/O2 discovery samples only.
- Adds a regression test proving that extreme changes in held-out O3 do not change discovery scores or selected samples.
- Matches safe samples to controls exactly by donor, day, batch and stage.
- Adds donor-level effects and chunked target-signature construction.
- Requires at least three independent donors for confirmatory t-test/FDR selection.
- With O1/O2 only, writes a clearly labelled exploratory ranked fallback for later O3/external validation.
- Adds matched-pair, donor-effect and QC outputs.
- 25 tests pass.


## 0.3.1 — pandas 3.x preprocessing compatibility fix

- Fixes `ValueError: assignment destination is read-only` in chunked beta-matrix preprocessing.
- Requests an explicit writable NumPy copy before masking invalid beta values.
- Adds a regression test that simulates read-only pandas/NumPy views.
- Keeps detection-p arrays read-only because they are never modified.

## 0.2.0 — benchmark validity fix

- Uses one seed-independent normalization range per benchmark instance.
- Uses one pooled nondominated IGD/IGD+ reference front for all algorithms and seeds.
- Uses one fixed hypervolume reference point and one reusable Monte Carlo sample cloud.
- Computes Exact-small once and reuses the identical front for every seed row.
- Adds IGD+ to benchmark output.
- Excludes zero-edit solutions from candidate protocol fronts.
- Writes no intervention as a separate comparator CSV.
- Adds `min_edits` and prespecified `min_relative_improvement` feasibility constraints.
- Adds metric metadata JSON containing normalization bounds and reference policy.

## 0.4.0

- Added exact safe/control matching on donor, day, batch and stage.
- Added donor-level independent-unit target signature inference.
- Added chunked EPIC beta-matrix processing to reduce memory use.
- Safe-window conditions are now matched exactly and require cross-donor day support.
- Added safe-window summary, matched-pair table, donor effects and QC outputs.
- Added four regression/integration tests; 23 tests pass.

## 0.5.0

- Added leakage-free held-out O3 validation (`scripts/04_validate_signature_o3.py`).
- Added exact O3 donor/day/batch/stage matching and diagnostic wrong-day, cross-donor, and failed-reprogramming controls.
- Added Pearson/Spearman, sign-concordance, effect-retention, top-k, permutation, and optional effect-matched random-CpG null tests.
- Added descriptive Tier 1–3 candidate prioritisation without re-fitting the frozen O1/O2 signature.
- Added four standalone figures and machine-readable QC/summary tables.

## v0.5.0

- Added leakage-free held-out O3 validation module and CLI.
- Frozen day/batch/stage are inferred from the O1/O2 matched-pair file.
- Added Tier 1/2/3 candidate prioritisation with maximum Tier-1 effect-retention guard.
- Added top-k sign concordance, Pearson/Spearman, effect retention and slope metrics.
- Added CpG-label permutation and effect-matched random-CpG-set null tests.
- Added wrong-day, cross-donor baseline and failed-reprogramming diagnostic controls.
- Added sample-level methylation QC, figures and a frozen-analysis manifest.
- Added seven O3-specific regression tests; full suite now contains 32 tests.

## 0.8.1 — literature-grounded regulatory and safety evidence completion

- Separates unreviewed or unsupplied evidence from reviewed negative evidence.
- Adds GENCODE nearest-gene/TSS annotation with signed and absolute distance.
- Adds ENCODE/SCREEN-compatible cCRE overlap for agnostic, T-cell and fibroblast contexts.
- Adds optional enhancer–gene link integration with confidence and provenance.
- Adds traceable literature relevance and safety evidence schemas using source IDs, PMIDs and DOIs.
- Adds an optional user-supplied COSMIC Cancer Gene Census adapter without redistributing licensed COSMIC data.
- Adds regulatory completeness, literature-review completeness, safety-review completeness and safety-uncertainty scores.
- Renames the v0.8.0 score in this stage to `integrated_experimental_evidence_score` to avoid causal overstatement.
- Adds risk- and uncertainty-adjusted CpG and region priorities, evidence-readiness grades and a mandatory manual-review queue.
- Adds five paper-ready figures, machine-readable summaries and a frozen v0.8.1 manifest.
- Adds five regression/integration tests; the full suite now contains 53 tests.

## 0.9.0.post5
- Maps GSE179847 processed beta columns through GEO !Sample_description.
- Excludes detection-p columns from the methylation matrix.
- Estimates 21% O2 slopes per donor x study-part series before donor aggregation.
- Adds separate 3% O2 sensitivity fields when available.

## 0.9.1 — 2026-10-29

First release whose source actually contains the post5 changes.

### Fixed

- **Applied the official v0.9.0.post5 hotfix, which had been committed but never
  executed.** The v0.9.0.post5 tag carried post4 source plus an unapplied
  `apply_post5.py`: both version strings read `0.9.0.post4`,
  `src/ccid/orthogonal_validation.py` contained no `!Sample_description`
  handling, and the post5 acceptance tests failed against it
  (2 failed, 70 passed). Anyone cloning that tag could not reproduce the
  GSE179847 results.
- **`read_gse179847_processed_matrix` CSV compatibility.** When no
  `" beta"`-suffixed columns are present the reader no longer delegates
  unconditionally to the tab-delimited `read_geo_series_matrix`, which requires
  `!series_matrix_table_begin` and therefore raised
  `ValueError: No series matrix table was found` on a conventional
  comma-delimited processed CSV. The header is now inspected first: a genuine
  GEO Series Matrix (a single comma-field beginning with `!`) is still
  delegated, while a processed CSV whose sample titles sit directly in the
  header is handled inline, treating every column after the first as a sample
  column. Restores `test_gse179847_processed_csv_streaming_with_direct_titles`
  and `test_natural_aging_accepts_gse179847_processed_csv`.

### Added

- `scripts/faz0_diagnostics.py` — discovery-ranking diagnostic gate: SE-floor
  anatomy of the frozen panel, ranking sensitivity across six variance
  treatments, genome-wide delta distribution and the true null direction
  probability, baseline-variance percentiles, and probe-mask / solo-WCGW overlap.
- `START_HERE.md`, `RESTORATION.md`, `setup.ps1`, `setup.sh`, `run_tests.ps1`,
  `run_faz0.ps1`.
- `data/derived/` post5 summaries and the candidate-level aging audit, and
  `data/metadata/dataset_manifest.csv`.

### Not distributed

The manuscript DOCX files and figure images are present in the working project
folder but are excluded from version control by `.gitignore` (`manuscript/`,
`*.docx`), so they are deliberately absent from this release and from its
archived snapshot. The repository distributes code, configuration, reference
definitions and derived result summaries only.

### Verification

Full suite: **72 passed**, matching the originally reported run.
Verified on Python 3.10 with pandas 2.3.3 and numpy 2.2.6.

### Known limitations

No `results/` tree or raw GEO matrices are included; results were never
committed to version control. Regenerating the discovery outputs is documented
in `START_HERE.md` section 3.
