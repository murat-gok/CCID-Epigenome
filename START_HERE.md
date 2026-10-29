# CCID-Epigenome — project state and what to do next

**Version:** 0.9.0.post5 + CSV-compatibility fix — see `RESTORATION.md`
**Restored:** 2026-10-29 from the recovered post4 tree
**Manuscript:** *Evidence-convergent prioritization of partial-reprogramming CpGs reveals limited aging concordance* (single author, M. Gök)
**Current task:** FAZ 0 diagnostic gate, before any manuscript revision.

---

## 1. What this folder is

The recovered CCID-Epigenome git repository, brought up to a verified working
state. It contains:

- the genuine source tree, now at post5 **plus** a rewritten CSV-compatibility
  fix — **72 tests pass**, matching the original run exactly (see §5)
- all reference data, plus the post5 derived summaries and candidate audit
- both manuscript DOCX versions and all seven figures, in `manuscript/`
  (present in this working folder only: `.gitignore` excludes `manuscript/`
  and `*.docx`, so they are not in the git repository or any released archive)
- the FAZ 0 diagnostic script, `scripts/faz0_diagnostics.py`
- `RESTORATION.md` — what was broken, what was fixed, and the one open decision

**It does not contain raw GEO data or any results.** `outputs/`, `candidates/`
and `manifests/` are empty and results were never committed to git, so there is
nothing to recover. Both must be regenerated — §3 and §4 say how.

Read `RESTORATION.md` before touching the code. It documents a reproducibility
problem with the published v0.9.0.post5 git tag that has to be fixed before
submission.

---

## 2. Set up the environment

Windows PowerShell, from this folder:

```powershell
.\setup.ps1
```

That creates `.venv`, installs `requirements.txt`, and runs the test suite.
**Expect "72 passed".** If you prefer WSL, this same folder is reachable at
`/mnt/c/Users/user/CCID-Epigenome` and `setup.sh` does the equivalent.

Dependencies are pure-Python wheels (numpy, pandas, scipy, scikit-learn,
statsmodels, matplotlib, pytest), so Windows-native works without WSL. The only
R file in the tree, `scripts/04c_export_epicv2_annotation.R`, is an optional
annotation exporter and is not on the FAZ 0 path.

---

## 3. Regenerate the FAZ 0 inputs

FAZ 0 needs two files that the discovery stage writes as side outputs:

```
results/real/safe_osk_target_signature_O1O2_all_loci.csv.gz
results/real/safe_osk_target_signature_O1O2_donor_effects.csv.gz
```

Neither survived. Producing them means rerunning the discovery chain. Only
steps 1 and 3 are slow; nothing here depends on the orthogonal-validation
datasets.

```powershell
# 1. download the two MPTR datasets (large; the methylation matrix is the big one)
python scripts\01_download_public_data.py --datasets gse165177_rna --output data\raw
python scripts\01_download_public_data.py --datasets gse165179_mptr_methylation --output data\raw

# 2. preprocess expression
python scripts\02a_preprocess_fibroblast_expression.py

# 3. preprocess methylation  (slow; writes data\processed\fibroblast_beta.csv.gz)
python scripts\02b_preprocess_fibroblast_methylation.py

# 4. operational day-10 window, calibrated on the two discovery donors
python scripts\02_define_safe_window.py `
  --expression data\processed\fibroblast_expression.csv.gz `
  --metadata   data\processed\fibroblast_expression_metadata.csv `
  --discovery-donors O1,O2 `
  --identity-floor -0.50 `
  --pluripotency-ceiling 0.50 `
  --output results\real\safe_window_scores.csv

# 5. donor-paired methylation signature  -> produces the FAZ 0 inputs
python scripts\03_build_target_signature.py `
  --beta        data\processed\fibroblast_beta.csv.gz `
  --metadata    data\processed\fibroblast_methylation_metadata.csv `
  --safe-window results\real\safe_window_scores.csv `
  --baseline-condition baseline `
  --discovery-donors O1,O2 `
  --min-abs-delta 0.03 --fdr 0.10 --fallback-top-n 500 `
  --output results\real\safe_osk_target_signature_O1O2.csv.gz
```

Check the exact `--discovery-donors` labels and `--baseline-condition` value
against `data\processed\fibroblast_methylation_metadata.csv` before running
step 5; the donor identifiers in GSE165179 must match what the script is given.

Then run the gate:

```powershell
.\run_faz0.ps1
```

or explicitly:

```powershell
python scripts\faz0_diagnostics.py `
  --all-loci      results\real\safe_osk_target_signature_O1O2_all_loci.csv.gz `
  --donor-effects results\real\safe_osk_target_signature_O1O2_donor_effects.csv.gz `
  --signature     results\real\safe_osk_target_signature_O1O2.csv.gz `
  --shortlist     data\reference\frozen_primary_shortlist_v090.csv `
  --outdir        results\faz0
```

`--signature` is not optional in practice: it triggers a sanity check that the
recomputed top-500 reproduces the stored frozen list. If that Jaccard is not
1.0, the `--all-loci` file does not come from the run that froze the panel.

---

## 4. Full reproduction (only if needed later)

The orthogonal-validation stages additionally need
`results/real/regulatory_safety_v081/candidate_cpg_regulatory_safety_v081.csv.gz`,
which is produced by stages `04` through `04f`. `scripts/run_orthogonal_validation_v090.sh`
drives the v0.9.0 orthogonal stages once that file exists, and expects the raw
matrices under `data/raw/orthogonal_validation/`:

```powershell
python scripts\00_download_orthogonal_validation_data.py --datasets core     --output data\raw\orthogonal_validation
python scripts\00_download_orthogonal_validation_data.py --datasets extended --output data\raw\orthogonal_validation
```

Accessions, roles and official record URLs are in
`data/metadata/dataset_manifest.csv`. GSE179847 needs both the processed matrix
and the series matrix; the exact URLs are in `data/README_DATA_HANDOFF.md`.

This is not on the critical path. FAZ 0 decides the revision route first.

---

## 5. Verification status of this rebuild

| Check | Result |
|---|---|
| `tests/test_orthogonal_validation_post5.py` | **2 passed** |
| Full suite | **72 passed** — matches the original run exactly |
| SciPy precision-loss warnings | `test_signature.py` and the bystander tests, as originally reported |
| `__version__` / `pyproject.toml` | `0.9.0.post5` |
| `Sample_description` handling present | yes (4 occurrences; was 0 before restoration) |

Verified on Linux / Python 3.10 / pandas 2.3.3 / numpy 2.2.6. Rerun on your own
machine with `.\run_tests.ps1` and expect 72.

---

## 6. Open items carried over from the handoff

These were flagged before the rebuild and are still open:

1. **The manuscript quotes PRE-POST5 integration tier counts** (12 V2 / 2 sentinels
   / 2 V3). `data/derived/orthogonal_validation_summary_PRE_POST5.json` marks them
   historical. Rerun `scripts/04l_integrate_orthogonal_validation.py` against the
   post5 natural-aging output before quoting final tiers.
2. **`data/derived/natural_aging_candidate_audit_v090.csv` is a reconstruction**
   from pasted console output, not an independent re-analysis. Verify against
   the original pipeline output before submission.
3. **`biological_shortlist_rank` in `data/reference/frozen_primary_shortlist_v090.csv`
   has duplicate values** (5, 6 and 8 each appear twice) and gaps (11, 14 missing),
   while manuscript Table 4 shows a clean 1–10. The two disagree.
4. **Manuscript Table 4 shows only 10 of the 16 candidates**; the six secondary
   candidates are omitted without the table saying so.
5. **Lifespan slopes are per passage, not per population doubling.** The 21% vs
   3% oxygen comparison is not valid until this is corrected, because low oxygen
   extends replicative lifespan and therefore changes doublings per passage.
6. **Do not regenerate publication DOCX/PDF through Google Docs.** Earlier
   conversions broke formatting. Generate DOCX natively (`python-docx`) and
   render the PDF from a controlled layout, then inspect every page.

---

## 7. Claims discipline

The manuscript must not assert that any CpG causes rejuvenation, is a validated
therapeutic target, is experimentally validated, or is safe to edit. The
defensible framing is "orthogonally supported computational candidate" and
"context-aware experimental hypothesis". Keep the two safety sentinels
(`cg13203228` at SSH1/AC087893.2 and `cg08575408` at GLIS3-AS1) visibly separated
from proposed first-line experimental targets.

Keep these five distinct and never interchangeable: (a) MPTR direction
concordance, (b) aging-direction opposition, (c) score stability, (d) sequence
feasibility, (e) causal or functional validation.

---

## 8. Map of this folder

```
START_HERE.md            this file — authoritative current state
RESTORATION.md           what was broken, what was fixed, the open decision
README.md                original v0.9.0 pipeline documentation (full reference)
CHANGELOG.md             version history, post5 entry appended
setup.ps1 / setup.sh     environment bootstrap
run_tests.ps1            test suite
run_faz0.ps1             FAZ 0 diagnostic gate
src/ccid/                library
scripts/                 pipeline stages; faz0_diagnostics.py added
configs/                 run profiles
data/reference/          frozen 16-CpG shortlist and evidence seeds
data/derived/            post5 summaries and the candidate audit
data/metadata/           dataset_manifest.csv
data/raw/                EMPTY — download target
data/processed/          EMPTY — preprocessing target
results/                 EMPTY — regenerate per §3
manuscript/              both DOCX versions
manuscript/figures/      Figures 1-7 as PNG
docs/handoff/            original handoff documents, verbatim
.git/                    history; tag v0.9.0.post5 (see RESTORATION.md)
*.post4.bak              pre-post5 backups
*.pre_csvfix.bak         pre-CSV-fix backup
tests/                   72 tests including the 2 post5 acceptance tests
```
