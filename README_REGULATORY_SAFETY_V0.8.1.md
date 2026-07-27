# CCID-Epigenome v0.8.1

## Literature-Grounded Regulatory Annotation and Safety Evidence Completion

v0.8.1 extends the frozen v0.8.0 candidate scores without changing the discovery, O3, LODO, bystander or v0.8.0 ranking inputs.

The release fixes the principal v0.8.0 evidence-semantics problem:

- **not reviewed / not supplied** is represented as `NaN` and an explicit status;
- **reviewed negative / no known signal** may be represented as a true score of `0`;
- absence of an overlap or literature row is never interpreted as biological irrelevance or safety.

The output is an experimental-prioritisation layer. It is not causal proof, a clinical safety assessment or an sgRNA design.

## Added evidence layers

1. Nearest GENCODE gene and transcription start site (TSS), with signed and absolute distance.
2. ENCODE/SREEN-compatible candidate cis-regulatory-element overlap.
3. Optional enhancer–gene links with confidence and source provenance.
4. Direct EPIC gene, nearest-TSS gene and enhancer-linked gene relationships with separate relationship weights.
5. Traceable literature evidence with PMID/DOI/source identifiers.
6. User-supplied COSMIC Cancer Gene Census adapter; COSMIC data are not redistributed.
7. Separate relevance, risk, regulatory completeness and safety-uncertainty scores.
8. A mandatory manual-review queue for high experimental candidates lacking reviewed literature or safety evidence.
9. Candidate and region rankings that penalise known risk and unresolved safety uncertainty.

## 1. Upgrade and test

```bash
cd ~/ccid_epigenome
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
python -c "import ccid; print(ccid.__version__)"
```

Expected:

```text
53 passed
0.8.1
```

The exact pytest dot display may vary. The existing SciPy precision-loss warnings come from nearly identical synthetic test values and do not indicate a failed test.

## 2. Download official regulatory references

```bash
cd ~/ccid_epigenome
source .venv/bin/activate

python scripts/00_fetch_v081_regulatory_references.py \
  --output-dir data/reference/v081
```

The fetcher retrieves:

- GENCODE release 50 gene annotation for GRCh38.p14;
- the ENCODE agnostic GRCh38 cCRE registry;
- activated T-cell cCREs;
- WI38 fibroblast cCREs.

WI38 is a fibroblast regulatory-context proxy. It must not be described as an MSC-specific regulatory map.

The fetcher records URLs, ENCODE accessions, file sizes and SHA-256 values in:

```text
data/reference/v081/reference_manifest_v081.json
```

Downloaded reference files are intentionally not bundled in the release ZIP.

## 3. Run v0.8.1

The main command automatically discovers `encode_ccre_*.bed.gz` files beside the GENCODE GTF.

```bash
cd ~/ccid_epigenome
source .venv/bin/activate

rm -rf results/real/regulatory_safety_v081
mkdir -p results/real/regulatory_safety_v081

python scripts/04f_complete_regulatory_safety.py \
  --v080-dir results/real/genomic_causality_v080 \
  --gencode-gtf data/reference/v081/gencode.v50.annotation.gtf.gz \
  --literature-evidence data/reference/literature_evidence_seed_v081.csv \
  --config-json configs/regulatory_safety_v081.json \
  --shortlist-size 100 \
  --output-dir results/real/regulatory_safety_v081
```

## 4. Optional enhancer–gene links

Supply a traceable local CSV/TSV using this schema:

```text
chromosome,start,end,target_gene,confidence,source_id
```

Run:

```bash
python scripts/04f_complete_regulatory_safety.py \
  ... \
  --enhancer-gene-links data/reference/enhancer_gene_links.csv
```

An input template is provided at:

```text
data/reference/enhancer_gene_links_template_v081.csv
```

If no enhancer-link file is supplied, the enhancer component remains unknown. It is not set to zero.

## 5. Candidate-specific literature completion

The bundled literature table is deliberately small and traceable. It is not an exhaustive literature review.

After the first run, inspect:

```text
results/real/regulatory_safety_v081/manual_literature_safety_review_queue.csv
```

This file includes candidate genes, recommended review actions and PubMed search links. Record reviewed evidence in a new CSV using:

```text
data/reference/literature_evidence_template_v081.csv
```

Required fields:

```text
gene
evidence_role
evidence_domain
evidence_status
evidence_strength
source_id
reviewed
```

Recommended provenance fields:

```text
evidence_level
source_type
source_year
source_url_or_doi
source_note
```

Allowed roles:

- `relevance`
- `risk`
- `protective`

Recommended statuses:

- relevance: `positive`, `negative`, `mixed`
- risk: `positive`, `risk`, `mixed`, `no_known_signal`, `negative`

Only `reviewed=true` rows with a populated `source_id` contribute to scoring.

Rerun with both files:

```bash
python scripts/04f_complete_regulatory_safety.py \
  ... \
  --literature-evidence \
data/reference/literature_evidence_seed_v081.csv,data/reference/candidate_literature_evidence_v081.csv
```

## 6. Optional Cancer Gene Census safety evidence

COSMIC files are not bundled because access and reuse are governed by the COSMIC licence. After obtaining an authorised Cancer Gene Census table, provide it locally:

```bash
python scripts/04f_complete_regulatory_safety.py \
  ... \
  --cancer-gene-census data/reference/Cosmic_CancerGeneCensus_Tsv_v102_GRCh38.tsv
```

The adapter detects common `Gene Symbol`, `Role in Cancer` and `Tier` columns and converts them into traceable proximity-risk evidence.

A Cancer Gene Census match is a **manual safety flag**. It does not prove that editing the CpG changes the gene or causes cancer.

## Outputs

```text
results/real/regulatory_safety_v081/
├── candidate_cpg_regulatory_safety_v081.csv.gz
├── candidate_region_regulatory_safety_v081.csv
├── candidate_gene_links_v081.csv.gz
├── literature_evidence_audit.csv
├── safety_evidence_audit.csv
├── evidence_status_summary.csv
├── manual_literature_safety_review_queue.csv
├── top_candidate_shortlist_v081.csv
├── regulatory_safety_v081_summary.json
├── regulatory_safety_v081_manifest.json
└── figures/
    ├── nearest_tss_distance.png
    ├── regulatory_annotation_status.png
    ├── experimental_evidence_vs_safety_uncertainty.png
    ├── top_candidate_cpgs_v081.png
    └── top_candidate_regions_v081.png
```

## Main score fields

### Experimental evidence

```text
integrated_experimental_evidence_score
```

This is the v0.8.0 experimental score renamed to avoid overstating causality.

### Regulatory evidence

```text
nearest_tss_gene
nearest_tss_distance_bp
regulatory_classes
regulatory_biosamples
regulatory_evidence_score
regulatory_evidence_completeness
regulatory_annotation_status
```

### Literature and safety

```text
literature_relevance_score
literature_evidence_status
safety_risk_score_v081
safety_evidence_status
safety_uncertainty_score
literature_source_ids
safety_source_ids
```

### Final priority

```text
literature_grounded_priority_score
risk_uncertainty_adjusted_priority_score
overall_evidence_completeness_v081
evidence_readiness_grade
experimental_priority_tier_v081
manual_review_required
```

## Evidence-readiness grades

- `E1_design_ready`: high score, sufficient regulatory coverage, reviewed safety evidence and no strong risk signal.
- `E2_high_priority_review_required`: strong candidate but literature or safety review remains incomplete.
- `E3_exploratory`: moderate support.
- `E4_insufficient`: low priority or insufficient evidence.

A lack of E1 candidates in the first run is expected when candidate-specific safety evidence has not yet been reviewed.

## Interpretation guardrails

Do not write:

- “No safety risk was found” when `safety_evidence_status=not_reviewed`.
- “The CpG regulates the nearest gene” from distance alone.
- “ENCODE overlap proves enhancer function in MSCs.”
- “The score proves causal rejuvenation.”

Appropriate wording:

> Candidate CpGs were prioritised by integrating frozen experimental evidence with nearest-TSS and cCRE context, traceable literature evidence where available, known-risk proximity flags and explicit penalties for unresolved safety uncertainty.
