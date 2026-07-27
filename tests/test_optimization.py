from __future__ import annotations

import numpy as np

from ccid.metrics import set_recovery
from ccid.optimizers import exact_enumeration, nsga2, select_target_priority
from ccid.schema import OptimizationContext
from ccid.synthetic import make_synthetic_problem


def test_optimizer_returns_feasible_protocol() -> None:
    problem = make_synthetic_problem(seed=7)
    # Use a controlled surrogate approximation in this unit test so the test
    # isolates optimization correctness rather than model generalization.
    context = OptimizationContext(
        predicted_effects=problem.true_effects,
        predicted_std=np.full_like(problem.true_effects, 0.004),
        target_effect=problem.target_effect,
        locus_weights=problem.locus_weights,
        candidate_table=problem.candidate_table,
        max_edits=4,
        uncertainty_quantile=1.0,
    )
    priority_pool = np.unique(
        np.concatenate([
            np.array(problem.planted_indices),
            np.argsort(problem.candidate_table["causal_score"].to_numpy())[-9:],
        ])
    )
    if len(priority_pool) > 16:
        priority_pool = priority_pool[-16:]
    front = exact_enumeration(context, candidate_indices=priority_pool, seed=1)
    selected = select_target_priority(front)
    recovery = set_recovery(selected, problem.planted_indices)
    assert selected.feasible
    assert selected.objectives[2] <= 4
    assert recovery["recall"] >= 2 / 3

    nsga_front = nsga2(context, population_size=28, generations=12, seed=2)
    assert any(record.feasible for record in nsga_front)
