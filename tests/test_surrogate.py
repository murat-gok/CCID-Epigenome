from __future__ import annotations

from ccid.surrogate import (
    BootstrapEffectSurrogate,
    compare_random_vs_group_cv,
    grouped_train_cal_test_split,
)
from ccid.synthetic import make_surrogate_training_data, make_synthetic_problem


def test_group_split_and_uncertainty() -> None:
    problem = make_synthetic_problem(seed=7)
    x, y, groups, _ = make_surrogate_training_data(problem, seed=11)
    train, cal, test = grouped_train_cal_test_split(groups)
    model = BootstrapEffectSurrogate(n_models=4, trees_per_model=12, random_state=4)
    model.fit(x[train], y[train])
    model.calibrate(x[cal], y[cal], level=0.90)
    evaluation = model.evaluate(x[test], y[test])
    assert evaluation.rmse < 0.08
    assert evaluation.empirical_coverage > 0.70

    leakage = compare_random_vs_group_cv(x, y, groups, random_state=3)
    assert leakage["random_cv_rmse"] < leakage["group_cv_rmse"]
