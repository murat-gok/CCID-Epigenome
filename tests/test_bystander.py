import numpy as np
import pandas as pd

from ccid.bystander import analyze_bystander_effects


def test_bystander_detects_direction_and_replicates():
    beta = pd.DataFrame(
        {
            "c1": [0.10, 0.20, 0.30, 0.40],
            "t1": [0.30, 0.20, 0.32, 0.40],
            "c2": [0.11, 0.21, 0.29, 0.41],
            "t2": [0.31, 0.21, 0.31, 0.41],
            "c3": [0.09, 0.19, 0.30, 0.39],
            "t3": [0.29, 0.19, 0.33, 0.39],
        },
        index=["cg1", "cg2", "cg3", "cg4"],
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["c1", "t1", "c2", "t2", "c3", "t3"],
            "treatment": ["control", "targeted"] * 3,
            "replicate": [1, 1, 2, 2, 3, 3],
        }
    )
    result = analyze_bystander_effects(
        beta,
        metadata,
        age_loci={"cg1"},
        min_abs_delta=0.05,
        fdr=0.25,
        n_permutations=19,
    )
    cg1 = result.locus_results.set_index("locus_id").loc["cg1"]
    assert cg1["delta_beta"] > 0.15
    assert bool(cg1["age_associated"])
    assert len(result.replicate_correlations) == 3


def test_bystander_requires_paired_replication():
    beta = pd.DataFrame({"a": [0.1], "b": [0.2]}, index=["cg1"])
    metadata = pd.DataFrame(
        {"sample_id": ["a", "b"], "treatment": ["control", "targeted"], "replicate": [1, 1]}
    )
    try:
        analyze_bystander_effects(beta, metadata)
    except ValueError as exc:
        assert "At least two" in str(exc)
    else:
        raise AssertionError("Expected replication check to fail")
