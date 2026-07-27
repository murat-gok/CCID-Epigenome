from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ccid.o3_validation import (
    infer_validation_stratum,
    select_validation_pair,
    validate_o3_signature,
)


def _metadata() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_id": [
                "O1_safe_d10",
                "O1_ctrl_d10",
                "O2_safe_d10",
                "O2_ctrl_d10",
                "O3_safe_d10",
                "O3_ctrl_d10",
                "O3_failed_d10",
                "O3_ctrl_d13",
            ],
            "donor_id": ["O1", "O1", "O2", "O2", "O3", "O3", "O3", "O3"],
            "condition": [
                "reprogrammed",
                "baseline",
                "reprogrammed",
                "baseline",
                "reprogrammed",
                "baseline",
                "failed",
                "baseline",
            ],
            "day": [10, 10, 10, 10, 10, 10, 10, 13],
            "batch": ["exp2"] * 8,
            "stage": ["endpoint"] * 8,
        }
    )


def _discovery_pairs() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "pair_id": ["pair_001", "pair_002"],
            "safe_sample_id": ["O1_safe_d10", "O2_safe_d10"],
            "baseline_sample_ids": ["O1_ctrl_d10", "O2_ctrl_d10"],
            "n_baseline_samples": [1, 1],
            "donor_id": ["O1", "O2"],
            "condition": ["reprogrammed", "reprogrammed"],
            "day": [10, 10],
            "batch": ["exp2", "exp2"],
            "stage": ["endpoint", "endpoint"],
        }
    )


def test_infer_validation_stratum_from_frozen_pairs() -> None:
    stratum = infer_validation_stratum(_discovery_pairs())
    assert int(stratum["day"]) == 10
    assert stratum["batch"] == "exp2"
    assert stratum["stage"] == "endpoint"


def test_infer_validation_stratum_rejects_multiple_strata() -> None:
    pairs = _discovery_pairs().copy()
    pairs.loc[1, "day"] = 13
    try:
        infer_validation_stratum(pairs)
    except ValueError as exc:
        assert "exactly one frozen discovery stratum" in str(exc)
    else:
        raise AssertionError("Expected multiple-stratum failure")


def test_select_validation_pair_uses_only_o3_primary_pair() -> None:
    pair = select_validation_pair(
        _metadata(), validation_donor="O3", day=10, batch="exp2", stage="endpoint"
    )
    assert pair.validation_sample_ids == ("O3_safe_d10",)
    assert pair.baseline_sample_ids == ("O3_ctrl_d10",)
    assert pair.failed_sample_ids == ("O3_failed_d10",)
    assert pair.wrong_day_baseline_ids == ("O3_ctrl_d13",)
    assert set(pair.cross_donor_baseline_ids) == {"O1_ctrl_d10", "O2_ctrl_d10"}


def test_o3_validation_tiers_and_metrics(tmp_path: Path) -> None:
    beta = pd.DataFrame(
        {
            "O1_ctrl_d10": [0.70, 0.50, 0.50, 0.50, 0.20],
            "O2_ctrl_d10": [0.70, 0.50, 0.50, 0.50, 0.20],
            "O3_safe_d10": [0.40, 0.65, 0.51, 0.30, np.nan],
            "O3_ctrl_d10": [0.70, 0.50, 0.50, 0.50, 0.20],
            "O3_failed_d10": [0.68, 0.51, 0.50, 0.48, 0.20],
            "O3_ctrl_d13": [0.60, 0.55, 0.50, 0.45, 0.20],
        },
        index=["cg1", "cg2", "cg3", "cg4", "cg5"],
    )
    path = tmp_path / "beta.csv.gz"
    beta.to_csv(path, index=True, index_label="cpg_id")
    signature = pd.DataFrame(
        {
            "cpg_id": ["cg1", "cg2", "cg3", "cg4", "cg5"],
            "target_delta_beta": [-0.30, 0.10, 0.20, 0.20, -0.10],
            "rank_score": [5, 4, 3, 2, 1],
        }
    )
    result = validate_o3_signature(
        path,
        _metadata(),
        signature,
        validation_donor="O3",
        day=10,
        batch="exp2",
        stage="endpoint",
        tier1_min_abs_delta=0.03,
        tier1_min_retention=0.25,
        tier1_max_retention=2.0,
        tier2_min_abs_delta=0.005,
        tier2_min_retention=0.01,
        top_k=(2, 5),
        n_permutations=10,
        n_random_sets=0,
        random_seed=7,
        chunksize=2,
    )
    table = result.cpg_validation.set_index("cpg_id")
    assert table.loc["cg1", "tier"] == "Tier1"
    assert table.loc["cg2", "tier"] == "Tier1"
    assert table.loc["cg3", "tier"] == "Tier2"
    assert table.loc["cg4", "tier"] == "Tier3"
    assert table.loc["cg5", "tier"] == "Tier3"
    assert result.summary["validation_donor"] == "O3"
    assert result.summary["analysis_status"] == "held_out_single_donor_descriptive_validation"
    assert result.summary["tier_counts"]["Tier1"] == 2
    assert result.summary["tier_counts"]["Tier2"] == 1
    assert "cpg_label_permutation" in set(result.null_tests["null_type"])


def test_tier1_rejects_implausible_effect_overshoot(tmp_path: Path) -> None:
    beta = pd.DataFrame(
        {
            "O1_ctrl_d10": [0.5],
            "O2_ctrl_d10": [0.5],
            "O3_safe_d10": [0.9],
            "O3_ctrl_d10": [0.5],
            "O3_failed_d10": [0.5],
            "O3_ctrl_d13": [0.5],
        },
        index=["cg1"],
    )
    path = tmp_path / "beta.csv.gz"
    beta.to_csv(path, index=True, index_label="cpg_id")
    signature = pd.DataFrame({"cpg_id": ["cg1"], "target_delta_beta": [0.05]})
    result = validate_o3_signature(
        path,
        _metadata(),
        signature,
        tier1_min_abs_delta=0.03,
        tier1_min_retention=0.25,
        tier1_max_retention=2.0,
        tier2_min_abs_delta=0.01,
        tier2_min_retention=0.10,
        n_permutations=0,
        n_random_sets=0,
        chunksize=2,
    )
    assert result.cpg_validation.loc[0, "effect_retention"] > 2.0
    assert result.cpg_validation.loc[0, "tier"] == "Tier2"


def test_o3_validation_null_is_reproducible(tmp_path: Path) -> None:
    beta = pd.DataFrame(
        {
            "O1_ctrl_d10": [0.5, 0.5, 0.5, 0.5],
            "O2_ctrl_d10": [0.5, 0.5, 0.5, 0.5],
            "O3_safe_d10": [0.4, 0.6, 0.3, 0.8],
            "O3_ctrl_d10": [0.7, 0.5, 0.5, 0.5],
            "O3_failed_d10": [0.69, 0.51, 0.49, 0.51],
            "O3_ctrl_d13": [0.6, 0.55, 0.45, 0.55],
        },
        index=["cg1", "cg2", "cg3", "cg4"],
    )
    path = tmp_path / "beta.csv.gz"
    beta.to_csv(path, index=True, index_label="cpg_id")
    signature = pd.DataFrame(
        {"cpg_id": beta.index, "target_delta_beta": [-0.3, 0.1, -0.2, 0.3]}
    )
    kwargs = dict(
        beta_path=path,
        metadata=_metadata(),
        signature=signature,
        n_permutations=20,
        n_random_sets=0,
        random_seed=123,
        chunksize=2,
    )
    first = validate_o3_signature(**kwargs)
    second = validate_o3_signature(**kwargs)
    pd.testing.assert_frame_equal(first.null_tests, second.null_tests)


def test_o3_pair_fails_without_exact_baseline() -> None:
    metadata = _metadata()
    metadata = metadata[metadata["sample_id"] != "O3_ctrl_d10"]
    try:
        select_validation_pair(metadata, validation_donor="O3", day=10, batch="exp2", stage="endpoint")
    except ValueError as exc:
        assert "No held-out baseline sample" in str(exc)
    else:
        raise AssertionError("Expected exact-match failure")
