# Gene evidence schema for v0.8.0

The bundled `biological_evidence_seed_v080.csv` is a heuristic triage panel. It is not an exhaustive literature review and must not be cited as one.

External evidence tables may contain:

```text
gene                 uppercase or standard gene symbol
category             biological or risk category
relevance_weight     0–1 biological relevance weight
risk_weight          0–1 conservative safety-proximity weight
evidence_level       e.g. experimental, review, canonical_seed
source_id            DOI, PMID, database release, or internal evidence ID
source_note          short auditable note
```

A gene may occur in several categories and sources. Use positive relevance weights for biological-priority evidence and positive risk weights only for conservative manual-review flags.
