# Next real-data step: safe window and matched OSK target signature

Version 0.4.0 prevents day and batch confounding when constructing the methylation target.

## Scientific design

- O1 and O2 are discovery donors.
- O3 is held out and is not used to choose the safe window or target CpGs.
- Eligible safe-window conditions are matched exactly: `reprogrammed` and `reprogrammed_intermediate`.
- A selected day must be represented by both discovery donors.
- Each safe sample is matched to negative controls with the same donor, day, batch and stage.
- Pair effects are averaged within donor before inference, so technical samples are not treated as independent donors.
- The EPIC matrix is read in chunks to reduce memory requirements.

See `UPGRADE_V0.4.0.md` for commands.
