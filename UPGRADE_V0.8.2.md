# Upgrade to v0.8.2

## Install

Apply the patch over v0.8.1.post1, then reinstall editable mode:

```bash
cd ~/ccid_epigenome
unzip -o /path/to/ccid_epigenome_patch_v0.8.2.zip -d ~/ccid_epigenome
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
python -c "import ccid; print(ccid.__version__)"
```

Expected:

```text
58 passed
0.8.2
```

## Run manuscript-readiness validation

```bash
rm -rf results/real/manuscript_validation_v082
mkdir -p results/real/manuscript_validation_v082

python scripts/04g_validate_manuscript_readiness.py \
  --candidate-file results/real/regulatory_safety_v081/candidate_cpg_regulatory_safety_v081.csv.gz \
  --config-json configs/manuscript_validation_v082.json \
  --output-dir results/real/manuscript_validation_v082
```

If a reviewed v0.8.1 table has been produced, use:

```bash
--candidate-file results/real/regulatory_safety_v081_reviewed/candidate_cpg_regulatory_safety_v081.csv.gz
```

## Verify

```bash
cat results/real/manuscript_validation_v082/manuscript_validation_v082_summary.json

python - <<'PY'
import pandas as pd
p = "results/real/manuscript_validation_v082/manuscript_ready_robust_shortlist_v082.csv"
df = pd.read_csv(p)
cols = [
    "manuscript_shortlist_rank", "cpg_id", "gene", "nearest_tss_gene",
    "v082_baseline_score", "robustness_composite_score",
    "weight_top_k_frequency", "uncertainty_top_k_frequency",
    "weight_rank_iqr", "uncertainty_rank_iqr",
    "worst_ablation_rank_drop", "manuscript_stability_status",
]
print(df[[c for c in cols if c in df.columns]].head(30).to_string(index=False))
PY
```
