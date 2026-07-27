# Data layout

Public data are not redistributed in this package.

- `raw/`: downloaded GEO/CRISPRepi files.
- `processed/`: standardized expression, methylation, metadata, guide, evidence and surrogate matrices.

Expected matrix convention: rows are genes/CpGs/features and columns are sample IDs unless a script explicitly requests NumPy arrays.

Core metadata columns:

- Reprogramming: `sample_id, donor_id, day, condition`
- Bystander: `sample_id, treatment, replicate`
- Target signature: `cpg_id, target_delta_beta, weight`
- Guide library: see `scripts/04_build_candidate_library.py` and `src/ccid/candidates.py`.

Large public datasets should be downloaded with `scripts/01_download_public_data.py`. Verify GEO filenames before a production run because repository-side filenames may change.
