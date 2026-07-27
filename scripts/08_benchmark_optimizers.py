#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import sys
import time
import warnings
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import numpy as np
import pandas as pd

from ccid.metrics import (
    apply_normalization,
    compute_normalization_bounds,
    igd,
    igd_plus,
    make_hypervolume_samples,
    monte_carlo_hypervolume,
    nondominated_points,
)
from ccid.objectives import no_intervention_summary, relative_target_improvement
from ccid.optimizers import (
    binary_differential_evolution,
    exact_enumeration,
    greedy_search,
    moead,
    nsga2,
    random_search,
    select_target_priority,
)
from ccid.schema import OptimizationContext, SolutionRecord, validate_candidate_table


def build_context(args: argparse.Namespace) -> OptimizationContext:
    candidates = pd.read_csv(args.candidates)
    validate_candidate_table(candidates)
    mean = np.load(args.effect_mean)
    std = np.load(args.effect_std)
    target = pd.read_csv(args.target_signature)
    if not {"target_delta_beta", "weight"}.issubset(target.columns):
        raise ValueError("Target signature requires target_delta_beta and weight columns")
    if mean.shape != std.shape or mean.shape != (len(candidates), len(target)):
        raise ValueError("Candidate/effect/target dimensions do not match")
    return OptimizationContext(
        predicted_effects=mean,
        predicted_std=std,
        target_effect=target["target_delta_beta"].to_numpy(float),
        locus_weights=target["weight"].to_numpy(float),
        candidate_table=candidates,
        max_edits=args.max_edits,
        min_edits=args.min_edits,
        min_relative_improvement=args.min_relative_improvement,
        uncertainty_quantile=args.uncertainty_quantile,
    )


def _protocol_front(front: list[SolutionRecord], context: OptimizationContext) -> list[SolutionRecord]:
    return [
        record
        for record in front
        if record.feasible and record.objectives[2] >= context.min_edits
    ]


def _objective_matrix(front: list[SolutionRecord], context: OptimizationContext) -> np.ndarray:
    protocol = _protocol_front(front, context)
    if not protocol:
        raise RuntimeError(
            "An optimizer returned no feasible protocol satisfying min_edits and "
            "min_relative_improvement. Relax only prespecified thresholds or increase its budget."
        )
    return np.vstack([record.objectives for record in protocol])


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Equal-budget optimizer benchmark with seed-independent metrics"
    )
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--effect-mean", required=True)
    parser.add_argument("--effect-std", required=True)
    parser.add_argument("--target-signature", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--instance", default="instance_0")
    parser.add_argument("--seeds", default=",".join(str(x) for x in range(30)))
    parser.add_argument("--max-edits", type=int, default=4)
    parser.add_argument("--min-edits", type=int, default=1)
    parser.add_argument("--min-relative-improvement", type=float, default=0.10)
    parser.add_argument("--population", type=int, default=160)
    parser.add_argument("--generations", type=int, default=250)
    parser.add_argument(
        "--random-samples",
        type=int,
        default=0,
        help="0 matches the evolutionary objective-evaluation budget",
    )
    parser.add_argument("--exact-pool-size", type=int, default=16)
    parser.add_argument("--uncertainty-quantile", type=float, default=1.645)
    parser.add_argument("--metric-seed", type=int, default=20260721)
    parser.add_argument("--hv-samples", type=int, default=100000)
    parser.add_argument("--reference-margin", type=float, default=1.10)
    args = parser.parse_args()

    if not 0.0 <= args.min_relative_improvement < 1.0:
        raise ValueError("--min-relative-improvement must be in [0, 1)")
    if args.reference_margin <= 1.0:
        raise ValueError("--reference-margin must be > 1.0")

    context = build_context(args)
    seeds = [int(value) for value in args.seeds.split(",") if value]
    if not seeds:
        raise ValueError("At least one seed is required")

    candidates = context.candidate_table
    priority = (
        candidates["causal_score"]
        + candidates["editability"]
        - candidates["sequence_offtarget"]
        - candidates["oncogenic_risk"]
        - candidates["identity_risk"]
    ).to_numpy()
    exact_pool_size = min(args.exact_pool_size, 18, len(candidates))
    exact_pool = np.argsort(priority)[-exact_pool_size:]

    evolutionary_budget = args.population * (args.generations + 1)
    random_budget = args.random_samples if args.random_samples > 0 else evolutionary_budget
    if random_budget != evolutionary_budget:
        warnings.warn(
            "Random search uses a different objective-evaluation budget. "
            "For paper comparisons omit --random-samples (or set it to "
            f"{evolutionary_budget}).",
            RuntimeWarning,
        )

    # Exact-small is deterministic and intentionally computed only once.
    exact_start = time.perf_counter()
    exact_front = exact_enumeration(context, candidate_indices=exact_pool, seed=0)
    exact_runtime = time.perf_counter() - exact_start
    exact_evaluations = sum(
        math.comb(exact_pool_size, k)
        for k in range(context.min_edits, context.max_edits + 1)
        if k <= exact_pool_size
    )

    all_runs: dict[int, dict[str, list[SolutionRecord]]] = {}
    run_meta: dict[tuple[int, str], dict[str, float | int | bool]] = {}
    for seed in seeds:
        all_runs[seed] = {}

        timed_calls = {
            "NSGA-II": lambda: nsga2(
                context,
                population_size=args.population,
                generations=args.generations,
                seed=seed,
            ),
            "MOEA-D": lambda: moead(
                context,
                population_size=args.population,
                generations=args.generations,
                seed=seed,
            ),
            "Binary-DE": lambda: binary_differential_evolution(
                context,
                population_size=args.population,
                generations=args.generations,
                seed=seed,
            ),
            "Random": lambda: random_search(
                context, n_samples=random_budget, seed=seed
            ),
            "Greedy": lambda: greedy_search(context, seed=seed),
        }
        for algorithm, call in timed_calls.items():
            started = time.perf_counter()
            all_runs[seed][algorithm] = call()
            elapsed = time.perf_counter() - started
            if algorithm in {"NSGA-II", "MOEA-D", "Binary-DE"}:
                evaluations = evolutionary_budget
                count_exact = True
            elif algorithm == "Random":
                evaluations = random_budget
                count_exact = True
            else:
                # One empty-mask evaluation plus all single-addition candidates
                # considered at each greedy step; marked estimated because a rare
                # early stop can reduce the actual count.
                evaluations = 1 + sum(
                    max(len(candidates) - step, 0)
                    for step in range(context.max_edits)
                )
                count_exact = False
            run_meta[(seed, algorithm)] = {
                "runtime_seconds": elapsed,
                "objective_evaluations": evaluations,
                "evaluation_count_exact": count_exact,
                "budget_equal_to_evolutionary": evaluations == evolutionary_budget,
            }

        all_runs[seed]["Exact-small"] = exact_front
        run_meta[(seed, "Exact-small")] = {
            "runtime_seconds": exact_runtime,
            "objective_evaluations": exact_evaluations,
            "evaluation_count_exact": True,
            "budget_equal_to_evolutionary": exact_evaluations == evolutionary_budget,
        }

    # Compute ONE global set of normalization bounds across all algorithms and
    # seeds. These bounds are then frozen for every row in this instance.
    raw_matrices: list[np.ndarray] = []
    for fronts in all_runs.values():
        raw_matrices.extend(_objective_matrix(front, context) for front in fronts.values())
    lower, upper = compute_normalization_bounds(raw_matrices)

    normalized_by_run: dict[tuple[int, str], np.ndarray] = {}
    for seed, fronts in all_runs.items():
        for algorithm, front in fronts.items():
            normalized_by_run[(seed, algorithm)] = apply_normalization(
                _objective_matrix(front, context), lower, upper
            )

    # ONE pooled nondominated reference front for IGD/IGD+ across all runs.
    pooled = np.vstack(list(normalized_by_run.values()))
    reference_front = nondominated_points(pooled)

    # ONE common reference point and ONE common Monte Carlo cloud for all rows.
    n_objectives = pooled.shape[1]
    hv_reference = np.ones(n_objectives, dtype=float) * args.reference_margin
    hv_samples = make_hypervolume_samples(
        n_objectives,
        hv_reference,
        n_samples=args.hv_samples,
        seed=args.metric_seed,
    )

    rows: list[dict[str, object]] = []
    for seed, fronts in all_runs.items():
        for algorithm, front in fronts.items():
            protocol_front = _protocol_front(front, context)
            points = normalized_by_run[(seed, algorithm)]
            selected = select_target_priority(protocol_front)
            rows.append(
                {
                    "instance": args.instance,
                    "seed": seed,
                    "algorithm": algorithm,
                    "front_size": len(protocol_front),
                    "hypervolume": monte_carlo_hypervolume(
                        points,
                        reference=hv_reference,
                        samples=hv_samples,
                    ),
                    "igd": igd(points, reference_front),
                    "igd_plus": igd_plus(points, reference_front),
                    "feasible": selected.feasible,
                    "target_mismatch": selected.objectives[0],
                    "relative_target_improvement": relative_target_improvement(
                        selected.objectives[0], context
                    ),
                    "safety_risk": selected.objectives[1],
                    "edit_count": selected.objectives[2],
                    "uncertainty": selected.objectives[3],
                    **run_meta[(seed, algorithm)],
                }
            )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output, index=False)

    no_intervention_path = output.with_name(f"{output.stem}_no_intervention.csv")
    pd.DataFrame([{"instance": args.instance, **no_intervention_summary(context)}]).to_csv(
        no_intervention_path, index=False
    )

    metadata_path = output.with_name(f"{output.stem}_metric_metadata.json")
    metadata = {
        "instance": args.instance,
        "normalization_lower": lower.tolist(),
        "normalization_upper": upper.tolist(),
        "hypervolume_reference": hv_reference.tolist(),
        "metric_seed": args.metric_seed,
        "hypervolume_samples": args.hv_samples,
        "evolutionary_objective_evaluation_budget": evolutionary_budget,
        "random_objective_evaluation_budget": random_budget,
        "random_budget_equal": random_budget == evolutionary_budget,
        "reference_front_size": int(len(reference_front)),
        "reference_front_policy": "global pooled nondominated union across all algorithms and seeds",
        "exact_small_policy": "computed once and reused unchanged for every seed",
        "min_edits": context.min_edits,
        "min_relative_improvement": context.min_relative_improvement,
        "no_intervention_file": no_intervention_path.name,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print(f"Wrote {len(rows)} benchmark rows to {output}")
    print(f"Wrote no-intervention comparator to {no_intervention_path}")
    print(f"Wrote fixed metric metadata to {metadata_path}")


if __name__ == "__main__":
    main()
