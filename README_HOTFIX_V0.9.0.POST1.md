# CCID-Epigenome v0.9.0.post1 hotfix

## Problem

The legacy GEO platform annotation URL
`https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6883/annot/GPL6883.annot.gz`
returns HTTP 403.

## Fix

The downloader now retrieves the full GPL6883 platform data table from NCBI GEO's supported Accession Display endpoint:

`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6883&targ=self&view=data&form=text`

and writes it to:

`data/raw/orthogonal_validation/GPL6883_data.txt`

## Resume

```bash
cd ~/ccid_epigenome
source .venv/bin/activate
python scripts/00_download_orthogonal_validation_data.py --datasets core --output data/raw/orthogonal_validation
bash scripts/run_orthogonal_validation_v090.sh
```

Existing GEO files are detected and skipped, so only the missing platform table is downloaded.
