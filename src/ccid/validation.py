from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd

from .metrics import mean_pairwise_jaccard, set_recovery
from .objectives import scalarized_score
from .optimizers import nsga2, select_compromise, select_target_priority
from .schema import OptimizationContext, SolutionRecord


def run_seed_stability(
    context: OptimizationContext,
    seeds: list[int],
    population_size: int = 24,
    generations: int = 10,
) -> tuple[list[SolutionRecord], float]:
    compromises: list[SolutionRecord] = []
    for seed in seeds:
        front = nsga2(
            context,
            population_size=population_size,
            generations=generations,
            seed=seed,
        )
        compromises.append(select_target_priority(front))
    return compromises, mean_pairwise_jaccard(compromises)


def permutation_target_test(
    context: OptimizationContext,
    observed: SolutionRecord,
    n_permutations: int = 40,
    seed: int = 91,
) -> dict[str, float]:
    """Test whether optimization improves over the no-edit baseline more than shuffled targets."""
    from .objectives import evaluate_mask

    rng = np.random.default_rng(seed)
    empty = np.zeros(len(context.candidate_table), dtype=bool)
    baseline_obj, _, _ = evaluate_mask(empty, context)
    observed_improvement = float(
        (baseline_obj[0] - observed.objectives[0]) / max(baseline_obj[0], 1e-12)
    )
    null_improvements: list[float] = []

    for permutation in range(n_permutations):
        shuffled_target = rng.permutation(context.target_effect)
        null_context = replace(context, target_effect=shuffled_target)
        null_baseline, _, _ = evaluate_mask(empty, null_context)
        front = nsga2(
            null_context,
            population_size=22,
            generations=8,
            seed=seed + permutation + 1,
        )
        chosen = select_target_priority(front)
        improvement = (null_baseline[0] - chosen.objectives[0]) / max(null_baseline[0], 1e-12)
        null_improvements.append(float(improvement))

    null_array = np.asarray(null_improvements)
    p_value = (1 + np.sum(null_array >= observed_improvement)) / (len(null_array) + 1)
    return {
        "observed_relative_improvement": observed_improvement,
        "null_mean_relative_improvement": float(null_array.mean()),
        "null_sd_relative_improvement": float(null_array.std(ddof=1)) if len(null_array) > 1 else 0.0,
        "permutation_p_value": float(p_value),
    }


def ablation_study(
    context: OptimizationContext,
    planted_indices: tuple[int, ...],
    seed: int = 13,
) -> pd.DataFrame:
    variants: dict[str, OptimizationContext] = {
        "full": context,
        "no_causal_filter": replace(context, causal_threshold=0.0),
        "no_uncertainty_penalty": replace(
            context,
            predicted_std=np.zeros_like(context.predicted_std),
            uncertainty_quantile=0.0,
        ),
        "no_bystander_risk": replace(
            context,
            candidate_table=context.candidate_table.assign(bystander_risk=0.0),
        ),
        "relaxed_identity": replace(context, identity_risk_limit=10.0),
    }

    rows: list[dict[str, float | str | bool]] = []
    for offset, (name, variant) in enumerate(variants.items()):
        front = nsga2(
            variant,
            population_size=24,
            generations=10,
            seed=seed + offset,
        )
        compromise = select_target_priority(front)
        recovery = set_recovery(compromise, planted_indices)
        rows.append(
            {
                "variant": name,
                "feasible": compromise.feasible,
                "target_mismatch": compromise.objectives[0],
                "safety_risk": compromise.objectives[1],
                "edit_count": compromise.objectives[2],
                "uncertainty": compromise.objectives[3],
                "recovery_precision": recovery["precision"],
                "recovery_recall": recovery["recall"],
                "recovery_jaccard": recovery["jaccard"],
            }
        )
    return pd.DataFrame(rows)


def multiplex_sensitivity(
    context: OptimizationContext,
    max_edits_values: list[int],
    seed: int = 37,
) -> pd.DataFrame:
    rows: list[dict[str, float | int]] = []
    for offset, max_edits in enumerate(max_edits_values):
        variant = replace(context, max_edits=max_edits)
        front = nsga2(
            variant,
            population_size=24,
            generations=10,
            seed=seed + offset,
        )
        compromise = select_target_priority(front)
        rows.append(
            {
                "max_edits": max_edits,
                "target_mismatch": compromise.objectives[0],
                "safety_risk": compromise.objectives[1],
                "selected_edits": int(compromise.objectives[2]),
                "uncertainty": compromise.objectives[3],
            }
        )
    return pd.DataFrame(rows)
