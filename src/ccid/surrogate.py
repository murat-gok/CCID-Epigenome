from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GroupKFold, KFold

from .schema import SurrogatePredictions


@dataclass
class BootstrapEffectSurrogate:
    """Bootstrap ensemble with split-conformal uncertainty calibration."""

    n_models: int = 6
    trees_per_model: int = 20
    min_samples_leaf: int = 2
    random_state: int = 19

    def __post_init__(self) -> None:
        self.models: list[ExtraTreesRegressor] = []
        self.conformal_q_: float = 1.645
        self.aleatoric_floor_: float = 0.005

    def fit(self, x: np.ndarray, y: np.ndarray) -> "BootstrapEffectSurrogate":
        rng = np.random.default_rng(self.random_state)
        self.models = []
        for model_id in range(self.n_models):
            indices = rng.integers(0, len(x), size=len(x))
            model = ExtraTreesRegressor(
                n_estimators=self.trees_per_model,
                min_samples_leaf=self.min_samples_leaf,
                max_features=0.85,
                n_jobs=-1,
                random_state=self.random_state + model_id,
            )
            model.fit(x[indices], y[indices])
            self.models.append(model)
        return self

    def predict_raw(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if not self.models:
            raise RuntimeError("The surrogate must be fitted before prediction")
        predictions = np.stack([model.predict(x) for model in self.models], axis=0)
        mean = predictions.mean(axis=0)
        std = predictions.std(axis=0, ddof=1)
        std = np.sqrt(np.square(std) + self.aleatoric_floor_ ** 2)
        std = np.maximum(std, 1e-4)
        return mean, std

    def calibrate(self, x_cal: np.ndarray, y_cal: np.ndarray, level: float = 0.90) -> float:
        mean, std = self.predict_raw(x_cal)
        residual = np.abs(y_cal - mean)
        self.aleatoric_floor_ = float(max(np.median(residual), 0.004))
        mean, std = self.predict_raw(x_cal)
        scores = residual / std
        self.conformal_q_ = float(np.quantile(scores, level, method="higher"))
        return self.conformal_q_

    def evaluate(self, x_test: np.ndarray, y_test: np.ndarray) -> SurrogatePredictions:
        mean, std = self.predict_raw(x_test)
        lower = mean - self.conformal_q_ * std
        upper = mean + self.conformal_q_ * std
        coverage = float(np.mean((y_test >= lower) & (y_test <= upper)))
        rmse = float(np.sqrt(mean_squared_error(y_test, mean)))
        return SurrogatePredictions(
            mean=mean,
            std=std,
            conformal_q=self.conformal_q_,
            empirical_coverage=coverage,
            rmse=rmse,
        )


def grouped_train_cal_test_split(
    groups: np.ndarray,
    test_groups: tuple[int, ...] = (6, 7),
    calibration_groups: tuple[int, ...] = (5,),
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    test_mask = np.isin(groups, test_groups)
    cal_mask = np.isin(groups, calibration_groups)
    train_mask = ~(test_mask | cal_mask)
    if not train_mask.any() or not cal_mask.any() or not test_mask.any():
        raise ValueError("Every split must contain observations")
    return np.flatnonzero(train_mask), np.flatnonzero(cal_mask), np.flatnonzero(test_mask)


def compare_random_vs_group_cv(
    x: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    random_state: int = 23,
) -> dict[str, float]:
    """Quantify optimistic bias caused by random study mixing."""

    def evaluate_splits(splits: list[tuple[np.ndarray, np.ndarray]]) -> float:
        fold_rmse: list[float] = []
        for fold, (train_idx, test_idx) in enumerate(splits):
            model = ExtraTreesRegressor(
                n_estimators=35,
                min_samples_leaf=2,
                max_features=0.85,
                n_jobs=-1,
                random_state=random_state + fold,
            )
            model.fit(x[train_idx], y[train_idx])
            pred = model.predict(x[test_idx])
            fold_rmse.append(float(np.sqrt(mean_squared_error(y[test_idx], pred))))
        return float(np.mean(fold_rmse))

    kfold = KFold(n_splits=5, shuffle=True, random_state=random_state)
    random_splits = list(kfold.split(x))
    group_splits = list(GroupKFold(n_splits=min(5, len(np.unique(groups)))).split(x, y, groups))
    random_rmse = evaluate_splits(random_splits)
    group_rmse = evaluate_splits(group_splits)
    return {
        "random_cv_rmse": random_rmse,
        "group_cv_rmse": group_rmse,
        "optimism_ratio": random_rmse / group_rmse,
    }
