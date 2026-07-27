# Upgrade to CCID-Epigenome v0.2.0

This release fixes the benchmark-validity issues identified after the initial smoke test.

## What changed

1. Hypervolume and IGD no longer use seed-specific normalization.
2. One frozen normalization range is calculated per benchmark instance.
3. One pooled nondominated reference front is used for IGD and IGD+.
4. One fixed hypervolume reference point and Monte Carlo sample cloud are reused for every method and seed.
5. Exact-small is computed once and reused unchanged for every seed row.
6. Zero-edit solutions are excluded from protocol Pareto fronts.
7. No intervention is written to a separate comparator file.
8. A prespecified minimum relative target improvement is enforced.
9. Objective-evaluation counts and runtime are reported.
10. Random search automatically receives the same objective-evaluation budget as the evolutionary methods unless explicitly overridden.

## Recommended upgrade on WSL

Back up the current project:

```bash
cd ~
cp -a ccid_epigenome ccid_epigenome_backup_v0.1.0
```

Extract the patch ZIP directly over the current project:

```bash
unzip -o /mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.2.0.zip \
  -d ~/ccid_epigenome
```

Activate the existing environment and reinstall the editable package:

```bash
cd ~/ccid_epigenome
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
```

Expected test result:

```text
13 passed
```

## Re-run the corrected smoke benchmark

Do not pass `--random-samples`; the script will automatically match Random search to the evolutionary objective-evaluation budget.

```bash
python scripts/08_benchmark_optimizers.py \
  --candidates results/benchmark_smoke_inputs/candidates.csv \
  --effect-mean results/benchmark_smoke_inputs/effect_mean.npy \
  --effect-std results/benchmark_smoke_inputs/effect_std.npy \
  --target-signature results/benchmark_smoke_inputs/target.csv \
  --instance local_smoke_test \
  --seeds 0,1,2 \
  --population 40 \
  --generations 30 \
  --output results/local_optimizer_test_v0.2.0.csv
```

The command creates three files:

```text
results/local_optimizer_test_v0.2.0.csv
results/local_optimizer_test_v0.2.0_no_intervention.csv
results/local_optimizer_test_v0.2.0_metric_metadata.json
```

Check the corrected invariants:

```bash
python - <<'PY'
import pandas as pd
p = "results/local_optimizer_test_v0.2.0.csv"
df = pd.read_csv(p)
exact = df[df["algorithm"] == "Exact-small"]
print("Exact HV:", exact["hypervolume"].tolist())
print("Exact IGD:", exact["igd"].tolist())
print("Zero-edit protocols:", int((df["edit_count"] < 1).sum()))
print("Equal-budget stochastic rows:",
      df[df["algorithm"].isin(["NSGA-II", "MOEA-D", "Binary-DE", "Random"])][
          "budget_equal_to_evolutionary"
      ].all())
PY
```

Expected properties:

- Exact-small hypervolume is identical for every seed.
- Exact-small IGD is identical for every seed.
- No protocol row has `edit_count = 0`.
- NSGA-II, MOEA-D, Binary-DE and Random have equal objective-evaluation budgets.
- No intervention appears only in the separate baseline CSV.

## Important scientific note

The default minimum relative improvement is 10%:

```text
--min-relative-improvement 0.10
```

This threshold must be prespecified before the main experiment and subjected to sensitivity analysis. It should not be changed after inspecting the final results merely to increase the number of feasible protocols.
