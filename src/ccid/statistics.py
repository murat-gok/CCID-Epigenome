from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def vargha_delaney_a12(x: np.ndarray, y: np.ndarray) -> float:
    """Probability that a random x value is better (smaller) than random y."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    wins = sum(np.sum(value < y) for value in x)
    ties = sum(np.sum(value == y) for value in x)
    return float((wins + 0.5 * ties) / (len(x) * len(y)))


def compare_algorithms(
    results: pd.DataFrame,
    metric: str,
    lower_is_better: bool = True,
    block_column: str = "instance",
    algorithm_column: str = "algorithm",
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Friedman omnibus + paired Wilcoxon-Holm + A12 effect sizes."""

    pivot = results.pivot_table(
        index=block_column,
        columns=algorithm_column,
        values=metric,
        aggfunc="mean",
    ).dropna()
    if pivot.shape[0] < 3 or pivot.shape[1] < 3:
        raise ValueError("At least 3 complete blocks and 3 algorithms are required")

    friedman = stats.friedmanchisquare(*[pivot[col].to_numpy() for col in pivot.columns])
    rows: list[dict[str, float | str]] = []
    raw_p: list[float] = []
    pairs = list(combinations(pivot.columns, 2))
    for a, b in pairs:
        x = pivot[a].to_numpy()
        y = pivot[b].to_numpy()
        alternative = "less" if lower_is_better else "greater"
        try:
            test = stats.wilcoxon(x, y, alternative=alternative, zero_method="pratt")
            p_value = float(test.pvalue)
            statistic = float(test.statistic)
        except ValueError:
            p_value = 1.0
            statistic = 0.0
        raw_p.append(p_value)
        rows.append(
            {
                "algorithm_a": a,
                "algorithm_b": b,
                "wilcoxon_statistic": statistic,
                "p_raw": p_value,
                "a12_a_better": vargha_delaney_a12(x, y) if lower_is_better else vargha_delaney_a12(-x, -y),
                "median_difference_a_minus_b": float(np.median(x - y)),
            }
        )
    adjusted = multipletests(raw_p, method="holm")[1]
    for row, p_adj in zip(rows, adjusted):
        row["p_holm"] = float(p_adj)
    return pd.DataFrame(rows), {
        "friedman_statistic": float(friedman.statistic),
        "friedman_p_value": float(friedman.pvalue),
        "n_blocks": float(pivot.shape[0]),
        "n_algorithms": float(pivot.shape[1]),
    }
