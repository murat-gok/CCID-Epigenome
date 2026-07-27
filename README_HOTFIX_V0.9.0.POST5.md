# CCID-Epigenome v0.9.0.post5 hotfix

This hotfix corrects GSE179847 sample linkage and longitudinal analysis.

## Fixed

1. Preserves GEO `!Sample_description`, which contains the Sentrix array/position identifier used by the processed matrix.
2. Retains only the 479 beta-value columns and excludes the 479 detection-p columns.
3. Maps each beta column to its GEO title and accession through `!Sample_description`.
4. Uses healthy `Normal_Control` samples at 21% O2 for the primary analysis.
5. Estimates slopes within each donor × study-part series before aggregating within and across donors.
6. Reports 3% O2 results separately as sensitivity fields when those series are available.

The patch is applied in place to a v0.9.0.post4 project so unrelated local changes remain intact. Backup files are written with the suffix `.post4.bak`.
