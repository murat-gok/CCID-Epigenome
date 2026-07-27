# v0.8.1 reference evidence

This directory contains only small, auditable configuration and evidence files. Large regulatory references are downloaded locally and are not redistributed in the package.

## Evidence semantics

- `not_reviewed` is represented as a missing score (`NaN`), never as zero.
- `0` is permitted only when a traceable manual review explicitly records negative evidence or `no_known_signal`.
- Literature rows contribute only when `reviewed=true` and `source_id` is populated.
- A gene or regulatory-element overlap is proximity evidence, not proof that editing the CpG changes that gene.

## Recommended references

Run `scripts/00_fetch_v081_regulatory_references.py` to retrieve:

- GENCODE release 50 GTF for GRCh38.p14.
- ENCODE cCRE annotations for the agnostic registry, activated T cells, and WI38 fibroblasts.

WI38 is used only as a fibroblast regulatory-context proxy. It must not be described as an MSC-specific map.

## Candidate-specific literature review

After the first v0.8.1 run, complete `manual_literature_safety_review_queue.csv` using the schema in `literature_evidence_template_v081.csv`. Add reviewed rows to a new CSV and rerun the main command with both evidence files separated by a comma.
