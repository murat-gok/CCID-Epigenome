from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .metrics import (
    igd,
    igd_plus,
    make_hypervolume_samples,
    monte_carlo_hypervolume,
    nondominated_points,
    normalize_objectives,
    set_recovery,
)
from .objectives import no_intervention_summary, predict_combined_effect, relative_target_improvement
from .optimizers import (
    exact_enumeration,
    greedy_search,
    nsga2,
    random_search,
    select_compromise,
    select_target_priority,
)
from .reporting import records_to_frame, save_effect_plot, save_pareto_plot
from .schema import OptimizationContext
from .surrogate import (
    BootstrapEffectSurrogate,
    compare_random_vs_group_cv,
    grouped_train_cal_test_split,
)
from .synthetic import make_surrogate_training_data, make_synthetic_problem
from .validation import ablation_study, multiplex_sensitivity, permutation_target_test, run_seed_stability


def run_demo(output_dir: str | Path, seed: int = 7) -> dict[str, object]:
    output_dir = Path(output_dir)
    tables_dir = output_dir / "tables"
    figures_dir = output_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    problem = make_synthetic_problem(seed=seed)
    x, y, groups, candidate_indices = make_surrogate_training_data(problem, seed=seed + 4)
    train_idx, cal_idx, test_idx = grouped_train_cal_test_split(groups)

    surrogate = BootstrapEffectSurrogate(random_state=seed + 10)
    surrogate.fit(x[train_idx], y[train_idx])
    surrogate.calibrate(x[cal_idx], y[cal_idx], level=0.90)
    test_evaluation = surrogate.evaluate(x[test_idx], y[test_idx])
    candidate_mean, candidate_std = surrogate.predict_raw(problem.candidate_features)

    leakage = compare_random_vs_group_cv(x, y, groups, random_state=seed + 20)

    context = OptimizationContext(
        predicted_effects=candidate_mean,
        predicted_std=candidate_std,
        target_effect=problem.target_effect,
        locus_weights=problem.locus_weights,
        candidate_table=problem.candidate_table,
        max_edits=4,
        min_edits=1,
        min_relative_improvement=0.10,
        uncertainty_quantile=surrogate.conformal_q_,
    )

    # Exact baseline on a prespecified small candidate pool.
    priority = (
        problem.candidate_table["causal_score"]
        + problem.candidate_table["editability"]
        - problem.candidate_table["sequence_offtarget"]
        - problem.candidate_table["oncogenic_risk"]
        - problem.candidate_table["identity_risk"]
    ).to_numpy()
    exact_pool = np.argsort(priority)[-12:]

    fronts = {
        "NSGA-II": nsga2(context, population_size=32, generations=16, seed=seed),
        "Random": random_search(context, n_samples=300, seed=seed + 1),
        "Greedy": greedy_search(context, seed=seed + 2),
        "Exact-small": exact_enumeration(context, candidate_indices=exact_pool, seed=seed + 3),
    }

    protocol_fronts = {
        method: [r for r in front if r.feasible and r.objectives[2] >= context.min_edits]
        for method, front in fronts.items()
    }
    matrices = [np.vstack([r.objectives for r in front]) for front in protocol_fronts.values()]
    normalized, _, _ = normalize_objectives(matrices)
    union = np.vstack(normalized)
    reference_front = nondominated_points(union)
    hv_reference = np.ones(union.shape[1], dtype=float) * 1.10
    hv_samples = make_hypervolume_samples(
        union.shape[1], hv_reference, n_samples=60000, seed=20260721
    )

    benchmark_rows: list[dict[str, object]] = []
    for (method, front), norm_points in zip(protocol_fronts.items(), normalized):
        compromise = select_target_priority(front)
        recovery = set_recovery(compromise, problem.planted_indices)
        benchmark_rows.append(
            {
                "method": method,
                "front_size": len(front),
                "hypervolume_mc": monte_carlo_hypervolume(
                    norm_points, reference=hv_reference, samples=hv_samples
                ),
                "igd_to_union": igd(norm_points, reference_front),
                "igd_plus_to_union": igd_plus(norm_points, reference_front),
                "compromise_feasible": compromise.feasible,
                "target_mismatch": compromise.objectives[0],
                "relative_target_improvement": relative_target_improvement(
                    compromise.objectives[0], context
                ),
                "safety_risk": compromise.objectives[1],
                "edit_count": int(compromise.objectives[2]),
                "uncertainty": compromise.objectives[3],
                "recovery_precision": recovery["precision"],
                "recovery_recall": recovery["recall"],
                "recovery_jaccard": recovery["jaccard"],
            }
        )
    benchmark = pd.DataFrame(benchmark_rows).sort_values("hypervolume_mc", ascending=False)

    nsga_compromise = select_target_priority(fronts["NSGA-II"])
    combined_mean, combined_std = predict_combined_effect(nsga_compromise.mask, context)
    lower = combined_mean - surrogate.conformal_q_ * combined_std
    upper = combined_mean + surrogate.conformal_q_ * combined_std

    selected_protocol = problem.candidate_table.loc[
        list(nsga_compromise.selected_indices)
    ].copy()
    selected_protocol.insert(0, "candidate_index", list(nsga_compromise.selected_indices))
    selected_protocol["selected_by"] = "NSGA-II compromise"

    ablations = ablation_study(context, problem.planted_indices, seed=seed + 30)
    sensitivity = multiplex_sensitivity(context, [1, 2, 3, 4, 5], seed=seed + 40)
    stability_records, stability_jaccard = run_seed_stability(
        context,
        seeds=[0, 1, 2],
    )
    permutation = permutation_target_test(
        context,
        nsga_compromise,
        n_permutations=19,
        seed=seed + 50,
    )

    problem.candidate_table.to_csv(tables_dir / "candidate_library.csv", index=False)
    benchmark.to_csv(tables_dir / "benchmark_summary.csv", index=False)
    pd.DataFrame([no_intervention_summary(context)]).to_csv(
        tables_dir / "no_intervention_baseline.csv", index=False
    )
    selected_protocol.to_csv(tables_dir / "selected_protocol.csv", index=False)
    ablations.to_csv(tables_dir / "ablation_results.csv", index=False)
    sensitivity.to_csv(tables_dir / "multiplex_sensitivity.csv", index=False)
    records_to_frame(fronts["NSGA-II"]).to_csv(tables_dir / "nsga2_pareto_front.csv", index=False)
    records_to_frame(stability_records).to_csv(tables_dir / "seed_stability_solutions.csv", index=False)

    save_pareto_plot(fronts, figures_dir / "pareto_fronts.png")
    save_effect_plot(
        problem.target_effect,
        combined_mean,
        lower,
        upper,
        figures_dir / "selected_protocol_effect.png",
    )

    summary = {
        "synthetic_only": True,
        "planted_indices": problem.planted_indices,
        "selected_indices": nsga_compromise.selected_indices,
        "surrogate_test_rmse": test_evaluation.rmse,
        "surrogate_90pct_coverage": test_evaluation.empirical_coverage,
        "conformal_q": test_evaluation.conformal_q,
        "random_cv_rmse": leakage["random_cv_rmse"],
        "group_cv_rmse": leakage["group_cv_rmse"],
        "cv_optimism_ratio": leakage["optimism_ratio"],
        "seed_stability_mean_jaccard": stability_jaccard,
        **permutation,
    }
    (output_dir / "demo_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary
