from __future__ import annotations

import numpy as np

from ccid.metrics import (
    apply_normalization,
    igd,
    igd_plus,
    make_hypervolume_samples,
    monte_carlo_hypervolume,
    nondominated_points,
)
from ccid.objectives import evaluate_mask, no_intervention_summary
from ccid.optimizers import exact_enumeration, select_target_priority
from ccid.schema import OptimizationContext
from ccid.synthetic import make_synthetic_problem


def _context() -> OptimizationContext:
    problem = make_synthetic_problem(seed=19, n_candidates=14, n_loci=18)
    return OptimizationContext(
        predicted_effects=problem.true_effects,
        predicted_std=np.full_like(problem.true_effects, 0.003),
        target_effect=problem.target_effect,
        locus_weights=problem.locus_weights,
        candidate_table=problem.candidate_table,
        max_edits=4,
        min_edits=1,
        min_relative_improvement=0.05,
        uncertainty_quantile=1.0,
    )


def test_zero_edit_is_separate_nonprotocol_baseline() -> None:
    context = _context()
    empty = np.zeros(len(context.candidate_table), dtype=bool)
    objectives, feasible, violations = evaluate_mask(empty, context)
    baseline = no_intervention_summary(context)

    assert not feasible
    assert violations > 0
    assert objectives[2] == 0
    assert baseline["edit_count"] == 0
    assert baseline["is_protocol_candidate"] is False


def test_exact_front_has_no_zero_edit_protocol() -> None:
    context = _context()
    front = exact_enumeration(context, candidate_indices=np.arange(12), seed=0)
    selected = select_target_priority(front)
    assert selected.feasible
    assert selected.objectives[2] >= 1
    assert all(record.objectives[2] >= 1 for record in front if record.feasible)


def test_fixed_hypervolume_and_reference_metrics_are_reproducible() -> None:
    a = np.array([[0.10, 0.80, 0.40, 0.30], [0.30, 0.30, 0.50, 0.20]])
    b = np.array([[0.20, 0.60, 0.20, 0.40], [0.40, 0.20, 0.30, 0.10]])
    lower = np.zeros(4)
    upper = np.ones(4)
    a_norm = apply_normalization(a, lower, upper)
    b_norm = apply_normalization(b, lower, upper)
    reference_front = nondominated_points(np.vstack([a_norm, b_norm]))
    reference = np.ones(4) * 1.10
    samples = make_hypervolume_samples(4, reference, n_samples=20000, seed=123)

    hv_first = monte_carlo_hypervolume(a_norm, reference=reference, samples=samples)
    hv_second = monte_carlo_hypervolume(a_norm, reference=reference, samples=samples)
    assert hv_first == hv_second
    assert igd(reference_front, reference_front) == 0.0
    assert igd_plus(reference_front, reference_front) == 0.0
