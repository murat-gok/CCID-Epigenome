# Real-data stage: leakage-free safe window and matched OSK target signature

Version 0.4.1 implements the next analysis stage.

- O1 and O2 are discovery donors.
- O3 is held out from score calibration, window selection and CpG ranking.
- Eligible conditions are exactly `reprogrammed` and `reprogrammed_intermediate`.
- A selected day must be represented by both discovery donors.
- Each selected sample is matched to negative controls with the same donor, day, batch and stage.
- Technical/temporal pairs are averaged within donor before ranking CpGs.
- Two donors are insufficient for confirmatory genome-wide FDR inference; the target list is explicitly labelled exploratory and must be validated on O3 and external perturbation data.
- The large beta matrix is processed in chunks.

Run the commands in `UPGRADE_V0.4.1.md`.
