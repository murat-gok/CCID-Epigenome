from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class SyntheticProblem:
    """Container for a complete inverse-design benchmark instance."""

    candidate_table: pd.DataFrame
    candidate_features: np.ndarray
    true_effects: np.ndarray
    target_effect: np.ndarray
    locus_weights: np.ndarray
    planted_indices: tuple[int, ...]
    locus_groups: np.ndarray
    editor_names: tuple[str, ...]


@dataclass(frozen=True)
class SurrogatePredictions:
    mean: np.ndarray
    std: np.ndarray
    conformal_q: float
    empirical_coverage: float
    rmse: float


@dataclass(frozen=True)
class OptimizationContext:
    predicted_effects: np.ndarray
    predicted_std: np.ndarray
    target_effect: np.ndarray
    locus_weights: np.ndarray
    candidate_table: pd.DataFrame
    max_edits: int = 4
    min_edits: int = 1
    min_relative_improvement: float = 0.10
    causal_threshold: float = 0.35
    identity_risk_limit: float = 1.10
    safety_risk_limit: float = 2.40
    uncertainty_quantile: float = 1.645


@dataclass(frozen=True)
class SolutionRecord:
    mask: np.ndarray
    objectives: np.ndarray
    feasible: bool
    violations: float
    method: str
    seed: int

    @property
    def selected_indices(self) -> tuple[int, ...]:
        return tuple(np.flatnonzero(self.mask).tolist())


def validate_candidate_table(table: pd.DataFrame) -> None:
    required = {
        "candidate_id",
        "region_id",
        "editor",
        "causal_score",
        "editability",
        "sequence_offtarget",
        "bystander_risk",
        "oncogenic_risk",
        "identity_risk",
        "pam_feasible",
        "direction_feasible",
    }
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"Candidate table is missing columns: {sorted(missing)}")

    if table["candidate_id"].duplicated().any():
        raise ValueError("candidate_id values must be unique")

    bounded = [
        "causal_score",
        "editability",
        "sequence_offtarget",
        "bystander_risk",
        "oncogenic_risk",
        "identity_risk",
    ]
    for col in bounded:
        values = table[col].to_numpy(dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{col} contains non-finite values")
        if np.any((values < 0) | (values > 1)):
            raise ValueError(f"{col} must lie in [0, 1]")


def masks_to_frame(
    masks: Sequence[np.ndarray],
    candidate_table: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for solution_id, mask in enumerate(masks):
        for idx in np.flatnonzero(mask):
            row = candidate_table.iloc[int(idx)]
            rows.append(
                {
                    "solution_id": solution_id,
                    "candidate_index": int(idx),
                    "candidate_id": row["candidate_id"],
                    "region_id": row["region_id"],
                    "editor": row["editor"],
                }
            )
    return pd.DataFrame(rows)
