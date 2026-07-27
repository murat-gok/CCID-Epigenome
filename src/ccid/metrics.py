from __future__ import annotations

import numpy as np

from .schema import SolutionRecord


def compute_normalization_bounds(
    matrices: list[np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    if not matrices:
        raise ValueError("At least one objective matrix is required")
    stacked = np.vstack(matrices)
    if stacked.ndim != 2 or stacked.shape[0] == 0:
        raise ValueError("Objective matrices must contain at least one row")
    lower = stacked.min(axis=0)
    upper = stacked.max(axis=0)
    return lower, upper


def apply_normalization(
    matrix: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=float)
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    span = np.maximum(upper - lower, 1e-12)
    return np.clip((matrix - lower) / span, 0.0, 1.0)


def normalize_objectives(
    matrices: list[np.ndarray],
    lower: np.ndarray | None = None,
    upper: np.ndarray | None = None,
) -> tuple[list[np.ndarray], np.ndarray, np.ndarray]:
    """Normalize all fronts with one shared, fixed set of bounds."""
    if lower is None or upper is None:
        lower, upper = compute_normalization_bounds(matrices)
    normalized = [apply_normalization(m, lower, upper) for m in matrices]
    return normalized, np.asarray(lower, dtype=float), np.asarray(upper, dtype=float)


def nondominated_points(points: np.ndarray) -> np.ndarray:
    """Return the unique minimization nondominated subset of a point cloud."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2:
        raise ValueError("points must be a 2-D matrix")
    if len(points) == 0:
        return points.copy()
    unique = np.unique(points, axis=0)
    keep = np.ones(len(unique), dtype=bool)
    for i, point in enumerate(unique):
        dominated = np.any(
            np.all(unique <= point, axis=1)
            & np.any(unique < point, axis=1)
        )
        keep[i] = not dominated
    return unique[keep]


def make_hypervolume_samples(
    n_objectives: int,
    reference: np.ndarray,
    n_samples: int = 60000,
    seed: int = 0,
) -> np.ndarray:
    reference = np.asarray(reference, dtype=float)
    if reference.shape != (n_objectives,):
        raise ValueError("reference must have one value per objective")
    if np.any(reference <= 0):
        raise ValueError("reference values must be positive")
    rng = np.random.default_rng(seed)
    return rng.random((n_samples, n_objectives)) * reference


def monte_carlo_hypervolume(
    points: np.ndarray,
    reference: np.ndarray | None = None,
    n_samples: int = 60000,
    seed: int = 0,
    samples: np.ndarray | None = None,
) -> float:
    """Estimate minimization hypervolume using a fixed reusable sample cloud.

    Passing the same ``samples`` for every algorithm/seed makes the estimator
    exactly reproducible and removes Monte Carlo noise from paired comparisons.
    """
    points = np.asarray(points, dtype=float)
    if len(points) == 0:
        return 0.0
    if reference is None:
        reference = np.ones(points.shape[1], dtype=float) * 1.10
    reference = np.asarray(reference, dtype=float)
    if samples is None:
        samples = make_hypervolume_samples(
            points.shape[1], reference, n_samples=n_samples, seed=seed
        )
    else:
        samples = np.asarray(samples, dtype=float)
        if samples.ndim != 2 or samples.shape[1] != points.shape[1]:
            raise ValueError("samples dimensionality does not match points")
    clipped = points[np.all(points <= reference, axis=1)]
    if len(clipped) == 0:
        return 0.0
    dominated = np.zeros(len(samples), dtype=bool)
    for point in clipped:
        dominated |= np.all(point <= samples, axis=1)
    return float(np.mean(dominated) * np.prod(reference))


def igd(approximation: np.ndarray, reference_front: np.ndarray) -> float:
    if len(approximation) == 0 or len(reference_front) == 0:
        return float("inf")
    distances = np.sqrt(
        np.sum(
            np.square(reference_front[:, None, :] - approximation[None, :, :]),
            axis=2,
        )
    )
    return float(np.mean(np.min(distances, axis=1)))


def igd_plus(approximation: np.ndarray, reference_front: np.ndarray) -> float:
    """IGD+ for minimization; lower is better and weak Pareto compliance holds."""
    approximation = np.asarray(approximation, dtype=float)
    reference_front = np.asarray(reference_front, dtype=float)
    if len(approximation) == 0 or len(reference_front) == 0:
        return float("inf")
    # For each reference point r and approximation point a, penalize only
    # objective components where a is worse than r.
    positive = np.maximum(
        approximation[None, :, :] - reference_front[:, None, :],
        0.0,
    )
    distances = np.sqrt(np.sum(np.square(positive), axis=2))
    return float(np.mean(np.min(distances, axis=1)))


def set_recovery(record: SolutionRecord, planted_indices: tuple[int, ...]) -> dict[str, float]:
    selected = set(record.selected_indices)
    planted = set(planted_indices)
    intersection = len(selected & planted)
    precision = intersection / max(len(selected), 1)
    recall = intersection / max(len(planted), 1)
    jaccard = intersection / max(len(selected | planted), 1)
    return {
        "precision": precision,
        "recall": recall,
        "jaccard": jaccard,
    }


def mean_pairwise_jaccard(records: list[SolutionRecord]) -> float:
    if len(records) < 2:
        return 1.0
    values: list[float] = []
    for i in range(len(records)):
        a = set(records[i].selected_indices)
        for j in range(i + 1, len(records)):
            b = set(records[j].selected_indices)
            values.append(len(a & b) / max(len(a | b), 1))
    return float(np.mean(values))
