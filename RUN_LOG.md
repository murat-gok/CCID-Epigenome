# Sandbox run log — v0.2.0

## Scope

Only the synthetic demonstration, unit tests and a small multi-optimizer software benchmark were run. Public GEO matrices were not downloaded in the execution environment, so no biological conclusions are drawn.

## Synthetic command

```bash
PYTHONPATH=src python scripts/run_synthetic_demo.py --output results --seed 7
```

## Synthetic summary

- Planted candidates: 7, 23, 24
- Selected candidates: 7, 24
- Surrogate test RMSE: 0.03558
- Nominal 90% interval coverage: 0.82199
- Random-CV RMSE: 0.02391
- Group-CV RMSE: 0.03090
- Random/group CV ratio: 0.77386
- Mean selected-set Jaccard across quick seeds: 1.0
- Observed relative improvement: 0.37747
- Null mean relative improvement: 0.05281
- Permutation p-value: 0.05 (19 permutations; coarse quick-profile resolution)

Interpretation: the code can recover a sparse planted signal, while the deliberately grouped split exposes optimistic random-CV performance and undercoverage under domain shift. The final manuscript must use the full profile, independent public data and wet-lab-neutral claims.

## Unit tests

```text
13 passed
```

The remaining SciPy precision-loss warning occurs only in a deliberately near-identical synthetic signature test and is not a failed test.

## Corrected multi-optimizer smoke benchmark

The corrected benchmark used:

```text
population = 40
generations = 30
evolutionary objective-evaluation budget = 1240
Random objective-evaluation budget = 1240
seeds = 0,1,2
```

Validation checks:

- Exact-small hypervolume is identical across all three seed rows: 0.367474459.
- Exact-small IGD and IGD+ are 0 across all seed rows.
- No candidate protocol has zero edits.
- No intervention is stored separately.
- NSGA-II, MOEA-D, Binary-DE and Random use the same objective-evaluation budget.
- Greedy and Exact-small report their own structural evaluation counts rather than being falsely labeled equal-budget.

Outputs:

```text
results/tables/optimizer_benchmark_smoke.csv
results/tables/optimizer_benchmark_smoke_no_intervention.csv
results/tables/optimizer_benchmark_smoke_metric_metadata.json
```

These values verify execution and metric invariants only. They must not be interpreted as biological evidence or final optimizer rankings.
