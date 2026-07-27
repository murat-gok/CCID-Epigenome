# CCID-Epigenome v0.9.0: orthogonal in-silico validation

## Purpose

This release addresses the absence of wet-lab validation by evaluating the previously frozen 16-CpG panel in data sources that were not used to fit the original nine-component ranking. It does **not** claim that any CpG causally rejuvenates a cell.

The analysis answers six narrower questions:

1. Does the frozen CpG direction replicate in an independent transient-reprogramming experiment?
2. Is the proposed reprogramming direction opposite to natural or replicative fibroblast aging?
3. Is candidate methylation associated with expression of an annotated nearby gene in primary fibroblasts?
4. Does a candidate overlap a fibroblast enhancer–gene prediction when an rE2G file is available?
5. Is the locus sequence-accessible to at least one basic SpCas9-NGG guide near the CpG?
6. Which candidates remain supported when these orthogonal evidence layers are integrated without refitting the original model?

## Frozen hypothesis panel

`data/reference/frozen_primary_shortlist_v090.csv` contains the 16 v0.8.2 stable-primary CpGs. The file is fixed before external validation. Original ranks and evidence weights are not re-estimated from the new datasets.

## Public data

### Core set

- `GSE142439`: four paired human skin-fibroblast donors, normal versus four-day OSKMLN treatment; independent direction validation.
- `GSE131280`: EPIC methylation across the replicative lifespan of primary fibroblasts in normal glucose; technical replicates at the same day are averaged.
- `GSE22595`: HumanMethylation27 profiles from young and older dermal-fibroblast donors; passage-3 samples only. Anatomical-site matched old-minus-young differences are calculated where possible.
- `GSE52025` and `GSE53243`: matched methylation and expression data from primary fibroblast lines; used for candidate CpG–gene expression association tests.
- `GPL6883`: expression-probe to gene-symbol annotation.

### Extended set

- `GSE179847`: 479 EPIC profiles across fibroblast lifespan. The analysis uses untreated healthy-control samples and estimates within-donor trends. This file is large and optional.

### Optional regulatory links

- ENCODE-rE2G annotation accessions `ENCSR360BQF` and `ENCSR886NXW` represent skin-fibroblast enhancer–gene predictions. Availability of downloadable association-graph files should be checked through the ENCODE API.

## Core workflow

```bash
cd ~/ccid_epigenome
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q

bash scripts/run_orthogonal_validation_v090.sh
```

Expected software test result: `64 passed`.

## Extended multidonor lifespan workflow

```bash
RUN_EXTENDED=1 bash scripts/run_orthogonal_validation_v090.sh
```

The extended matrix is large. The GEO Series Matrix parser streams the file and retains only the 16 frozen CpGs, but the compressed file still must be downloaded and decompressed sequentially.

## Individual stages

### 1. Download core data

```bash
python scripts/00_download_orthogonal_validation_data.py \
  --datasets core \
  --output data/raw/orthogonal_validation
```

A SHA-256 download manifest is written to the same directory.

### 2. Independent transient-reprogramming validation

```bash
python scripts/04h_validate_external_reprogramming.py
```

Primary outputs:

- `external_reprogramming_validation_v090.csv.gz`
- `external_reprogramming_summary_v090.json`
- effect forest and discovery-versus-external figures.

This is a directional generalization test across a different reprogramming cocktail and protocol. It is not an exact experimental replication of MPTR.

### 3. Natural-aging direction validation

```bash
python scripts/04i_validate_natural_aging_direction.py
```

Primary outputs:

- `natural_aging_direction_validation_v090.csv.gz`
- `natural_aging_summary_v090.json`
- aging-direction matrix and quadrant figures.

A candidate is directionally supported when the MPTR-associated methylation change opposes the aging trend. GSE131280 is a longitudinal single-line time course; GSE22595 is cross-sectional and has limited HM27 coverage. These limitations are retained in interpretation.

### 4. Fibroblast eQTM analysis

```bash
python scripts/04j_validate_fibroblast_eqtm.py
```

Primary outputs:

- `fibroblast_eqtm_validation_v090.csv`
- `fibroblast_eqtm_summary_v090.json`
- eQTM effect figure.

Spearman association and Benjamini–Hochberg FDR are reported. Association does not prove mediation or causality. Novel Ensembl transcripts absent from the older expression array cannot be evaluated in this stage.

### 5. Regulatory links and guide feasibility

```bash
python scripts/04k_assess_regulatory_links_and_editability.py
```

This stage retrieves hg38 sequence from the UCSC API and enumerates SpCas9-NGG protospacers within the configured distance from each CpG. It reports only basic sequence feasibility. It does not replace genome-wide off-target scoring, chromatin accessibility assays, or wet-lab guide validation.

Optional ENCODE files can be fetched and supplied as follows:

```bash
python scripts/00_fetch_encode_re2g.py

python scripts/04k_assess_regulatory_links_and_editability.py \
  --re2g-files data/raw/orthogonal_validation/encode_re2g/<file1> \
               data/raw/orthogonal_validation/encode_re2g/<file2> \
  --biosample-names fibroblast_skin_scalp fibroblast_skin_back
```

Inspect the ENCODE manifest before using the files because formats may change between portal releases.

### 6. Evidence integration

```bash
python scripts/04l_integrate_orthogonal_validation.py
```

The integrated file is:

`results/real/orthogonal_validation_v090/integrated/orthogonally_validated_candidate_panel_v090.csv`

Tiers:

- `V1`: independent reprogramming direction plus at least four total orthogonal supports;
- `V2`: independent direction plus at least two total supports;
- `V3`: at least two mechanistic/aging/editability supports but incomplete external direction support;
- `V4`: insufficient orthogonal support;
- `S`: pre-specified safety sentinel.

These tiers are validation summaries, not a newly fitted prediction model.

## Interpretation rules

- Do not use “experimentally validated target.”
- Use “orthogonally supported computational candidate.”
- Do not interpret non-significance as absence of biological relevance when a platform lacks the CpG or transcript.
- Do not equate OSKMLN direction concordance with proof of OSK/OSKM mechanism.
- Do not claim rejuvenation from aging-direction reversal alone.
- Preserve all safety sentinels regardless of external direction.
