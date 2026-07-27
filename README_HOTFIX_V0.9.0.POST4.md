# CCID-Epigenome v0.9.0.post4 hotfix

## Problem corrected

The previous extended-data definition downloaded `GSE179847_series_matrix.txt.gz` and described it as the two-gigabyte methylation matrix. That file is a small GEO Series Matrix used mainly for sample metadata and does not contain the deposited 866k-by-479 beta matrix. It therefore produced no frozen CpG rows and the downstream code attempted to access a missing multidonor slope column.

## Correct files

- Metadata: `GSE179847_series_matrix.txt.gz`
- Normalized beta values: `GSE179847_Cell_lifespan_DNAm_processed_matrix.csv.gz`

The processed beta file is read as a stream. Only the prespecified frozen CpGs are retained in memory. The metadata Series Matrix is used to translate GSM accessions to sample titles when needed.

## Install

```bash
cd ~/ccid_epigenome
unzip -o /mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.9.0.post4.zip -d ~/ccid_epigenome
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
python -c "import ccid; print(ccid.__version__)"
```

Expected version: `0.9.0.post4`.

## Remove the obsolete mistaken input

The small Series Matrix is still useful as metadata and should remain under its original name. Do not rename it as the processed beta file. Remove any stale partial processed download only when intentionally restarting:

```bash
rm -f data/raw/orthogonal_validation/GSE179847_Cell_lifespan_DNAm_processed_matrix.csv.gz.part
```

## Download extended inputs

```bash
python scripts/00_download_orthogonal_validation_data.py \
  --datasets extended \
  --output data/raw/orthogonal_validation
```

Interrupted downloads retain a `.part` file and resume on the next invocation when the server supports HTTP Range requests.

Verify:

```bash
ls -lh data/raw/orthogonal_validation/GSE179847*
gzip -t data/raw/orthogonal_validation/GSE179847_Cell_lifespan_DNAm_processed_matrix.csv.gz \
  && echo "GSE179847 processed matrix is intact"
```

## Run extended validation

```bash
RUN_EXTENDED=1 bash scripts/run_orthogonal_validation_v090.sh
```

The analysis includes untreated healthy-control trajectories only and computes donor-specific Theil-Sen slopes before panel-level directional assessment.
