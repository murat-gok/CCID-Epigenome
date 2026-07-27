from ccid.optimizers import binary_differential_evolution, moead
from ccid.schema import OptimizationContext
from ccid.synthetic import make_synthetic_problem
from ccid.surrogate import BootstrapEffectSurrogate, grouped_train_cal_test_split
from ccid.synthetic import make_surrogate_training_data


def _context():
    problem = make_synthetic_problem(seed=11, n_candidates=18, n_loci=16)
    x, y, groups, _ = make_surrogate_training_data(problem, seed=12)
    train, cal, _ = grouped_train_cal_test_split(groups)
    model = BootstrapEffectSurrogate(n_models=3, trees_per_model=8, random_state=13)
    model.fit(x[train], y[train])
    model.calibrate(x[cal], y[cal], level=0.90)
    mean, std = model.predict_raw(problem.candidate_features)
    return OptimizationContext(
        predicted_effects=mean,
        predicted_std=std,
        target_effect=problem.target_effect,
        locus_weights=problem.locus_weights,
        candidate_table=problem.candidate_table,
        max_edits=3,
        uncertainty_quantile=model.conformal_q_,
    )


def test_binary_de_runs_and_returns_front():
    front = binary_differential_evolution(_context(), population_size=12, generations=3, seed=1)
    assert front
    assert all(record.objectives.shape == (4,) for record in front)


def test_moead_runs_and_returns_front():
    front = moead(_context(), population_size=12, generations=3, neighborhood_size=4, seed=2)
    assert front
    assert all(record.objectives.shape == (4,) for record in front)
