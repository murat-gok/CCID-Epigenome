import numpy as np
import pandas as pd

from ccid.statistics import compare_algorithms, vargha_delaney_a12


def test_a12_direction():
    assert vargha_delaney_a12(np.array([1, 1]), np.array([2, 2])) == 1.0
    assert vargha_delaney_a12(np.array([2, 2]), np.array([1, 1])) == 0.0


def test_compare_algorithms_returns_pairwise_and_omnibus():
    rows = []
    for instance in range(5):
        rows.extend(
            [
                {"instance": instance, "algorithm": "A", "metric": 1.0 + instance * 0.01},
                {"instance": instance, "algorithm": "B", "metric": 1.5 + instance * 0.01},
                {"instance": instance, "algorithm": "C", "metric": 2.0 + instance * 0.01},
            ]
        )
    pairwise, omnibus = compare_algorithms(pd.DataFrame(rows), metric="metric")
    assert len(pairwise) == 3
    assert omnibus["n_blocks"] == 5
    assert omnibus["n_algorithms"] == 3
    assert set(pairwise.columns) >= {"p_raw", "p_holm", "a12_a_better"}
