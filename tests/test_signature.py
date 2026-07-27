from __future__ import annotations

import numpy as np
import pandas as pd

from ccid.signature import build_methylation_target_signature, define_safe_window


def test_safe_window_and_signature() -> None:
    genes = ["COL1A1", "DCN", "POU5F1", "NANOG", "CDKN2A", "SERPINE1"]
    samples = ["d1_base", "d1_day10", "d2_base", "d2_day10"]
    expression = pd.DataFrame(
        [
            [8, 7.8, 8.2, 7.9],
            [7, 6.8, 7.2, 6.9],
            [0, 0.2, 0, 0.1],
            [0, 0.1, 0, 0.2],
            [5, 3.0, 5.2, 3.1],
            [4, 2.5, 4.1, 2.4],
        ],
        index=genes,
        columns=samples,
    )
    metadata = pd.DataFrame(
        {
            "sample_id": samples,
            "donor_id": ["d1", "d1", "d2", "d2"],
            "day": [0, 10, 0, 10],
            "condition": ["baseline", "reprogrammed", "baseline", "reprogrammed"],
        }
    )
    safe = define_safe_window(
        expression,
        metadata,
        identity_genes=["COL1A1", "DCN"],
        pluripotency_genes=["POU5F1", "NANOG"],
        rejuvenation_genes=["CDKN2A", "SERPINE1"],
        identity_floor=-2.0,
        pluripotency_ceiling=2.0,
    )
    assert set(safe.selected_samples) == {"d1_day10", "d2_day10"}

    beta = pd.DataFrame(
        {
            "d1_base": [0.7, 0.3, 0.5],
            "d1_day10": [0.5, 0.4, 0.5],
            "d2_base": [0.75, 0.25, 0.5],
            "d2_day10": [0.55, 0.36, 0.5],
        },
        index=["cg1", "cg2", "cg3"],
    )
    signature = build_methylation_target_signature(
        beta,
        metadata,
        safe_samples=list(safe.selected_samples),
        min_abs_delta=0.02,
        fdr=1.0,
    )
    assert "cg1" in set(signature["cpg_id"])
    assert np.isfinite(signature["weight"]).all()


def test_safe_window_requires_cross_donor_day_support() -> None:
    genes = ["COL1A1", "POU5F1", "CDKN2A"]
    samples = ["o1_d10", "o1_d13", "o2_d10", "o2_d13"]
    expression = pd.DataFrame(
        [
            [5.0, 5.0, 5.0, 5.0],
            [0.0, 0.0, 0.0, 0.0],
            [3.0, 1.0, 3.0, 8.0],
        ],
        index=genes,
        columns=samples,
    )
    metadata = pd.DataFrame(
        {
            "sample_id": samples,
            "donor_id": ["O1", "O1", "O2", "O2"],
            "day": [10, 13, 10, 13],
            "condition": ["reprogrammed"] * 4,
        }
    )
    result = define_safe_window(
        expression,
        metadata,
        identity_genes=["COL1A1"],
        pluripotency_genes=["POU5F1"],
        rejuvenation_genes=["CDKN2A"],
        identity_floor=-10,
        pluripotency_ceiling=10,
        min_rejuvenation=-10,
        discovery_donors=["O1", "O2"],
        min_donors_per_day=2,
    )
    assert set(result.sample_scores.loc[result.sample_scores["selected_safe_window"], "donor_id"]) == {
        "O1",
        "O2",
    }


def test_exact_matched_signature_uses_same_day_batch_and_stage() -> None:
    from ccid.signature import build_matched_methylation_target_signature

    metadata = pd.DataFrame(
        {
            "sample_id": [
                "O1_safe_d10",
                "O1_ctrl_d10",
                "O1_ctrl_d13",
                "O2_safe_d10",
                "O2_ctrl_d10",
                "O2_ctrl_d13",
            ],
            "donor_id": ["O1", "O1", "O1", "O2", "O2", "O2"],
            "day": [10, 10, 13, 10, 10, 13],
            "condition": ["reprogrammed", "baseline", "baseline", "reprogrammed", "baseline", "baseline"],
            "batch": ["exp2", "exp2", "exp2", "exp2", "exp2", "exp2"],
            "stage": ["endpoint"] * 6,
        }
    )
    beta = pd.DataFrame(
        {
            "O1_safe_d10": [0.40, 0.60],
            "O1_ctrl_d10": [0.70, 0.50],
            "O1_ctrl_d13": [0.10, 0.95],
            "O2_safe_d10": [0.35, 0.65],
            "O2_ctrl_d10": [0.65, 0.55],
            "O2_ctrl_d13": [0.05, 0.90],
        },
        index=["cg1", "cg2"],
    )
    result = build_matched_methylation_target_signature(
        beta,
        metadata,
        safe_samples=["O1_safe_d10", "O2_safe_d10"],
        discovery_donors=["O1", "O2"],
        min_abs_delta=0.01,
        fdr=1.0,
        min_sign_consistency=1.0,
    )
    cg1 = result.all_loci.set_index("cpg_id").loc["cg1"]
    assert np.isclose(cg1["target_delta_beta"], -0.30)
    assert result.qc["n_matched_safe_samples"] == 2
    assert set(result.pair_table["baseline_sample_ids"]) == {"O1_ctrl_d10", "O2_ctrl_d10"}


def test_exact_matched_signature_reports_unmatched_samples() -> None:
    from ccid.signature import build_matched_methylation_target_signature

    metadata = pd.DataFrame(
        {
            "sample_id": ["O1_safe", "O1_ctrl", "O2_safe", "O2_ctrl", "O2_unmatched"],
            "donor_id": ["O1", "O1", "O2", "O2", "O2"],
            "day": [10, 10, 10, 10, 13],
            "condition": ["reprogrammed", "baseline", "reprogrammed", "baseline", "reprogrammed"],
            "batch": ["exp2"] * 5,
            "stage": ["endpoint"] * 5,
        }
    )
    beta = pd.DataFrame(
        {
            "O1_safe": [0.4],
            "O1_ctrl": [0.7],
            "O2_safe": [0.4],
            "O2_ctrl": [0.7],
            "O2_unmatched": [0.2],
        },
        index=["cg1"],
    )
    result = build_matched_methylation_target_signature(
        beta,
        metadata,
        safe_samples=["O1_safe", "O2_safe", "O2_unmatched"],
        min_abs_delta=0.01,
        fdr=1.0,
        min_sign_consistency=1.0,
    )
    assert result.qc["n_unmatched_safe_samples"] == 1
    assert result.qc["unmatched_safe_samples"] == ["O2_unmatched"]


def test_chunked_file_signature_matches_in_memory(tmp_path) -> None:
    from ccid.signature import (
        build_matched_methylation_target_signature,
        build_matched_methylation_target_signature_file,
    )

    metadata = pd.DataFrame(
        {
            "sample_id": ["O1_safe", "O1_ctrl", "O2_safe", "O2_ctrl"],
            "donor_id": ["O1", "O1", "O2", "O2"],
            "day": [10, 10, 10, 10],
            "condition": ["reprogrammed", "baseline", "reprogrammed", "baseline"],
            "batch": ["exp2"] * 4,
            "stage": ["endpoint"] * 4,
        }
    )
    beta = pd.DataFrame(
        {
            "O1_safe": [0.4, 0.6, 0.2],
            "O1_ctrl": [0.7, 0.5, 0.2],
            "O2_safe": [0.35, 0.65, 0.2],
            "O2_ctrl": [0.65, 0.55, 0.2],
        },
        index=["cg1", "cg2", "cg3"],
    )
    path = tmp_path / "beta.csv.gz"
    beta.to_csv(path, index=True, index_label="cpg_id")
    kwargs = dict(
        metadata=metadata,
        safe_samples=["O1_safe", "O2_safe"],
        min_abs_delta=0.01,
        fdr=1.0,
        min_sign_consistency=1.0,
    )
    memory = build_matched_methylation_target_signature(beta=beta, **kwargs)
    chunked = build_matched_methylation_target_signature_file(
        beta_path=path, chunksize=2, **kwargs
    )
    left = memory.all_loci.set_index("cpg_id")["target_delta_beta"].sort_index()
    right = chunked.all_loci.set_index("cpg_id")["target_delta_beta"].sort_index()
    assert np.allclose(left, right, equal_nan=True)
    assert chunked.qc["chunksize"] == 2


def test_safe_window_calibration_excludes_held_out_donor() -> None:
    genes = ["COL1A1", "POU5F1", "CDKN2A"]
    samples = ["O1_d10", "O2_d10", "O3_d10"]
    metadata = pd.DataFrame(
        {
            "sample_id": samples,
            "donor_id": ["O1", "O2", "O3"],
            "day": [10, 10, 10],
            "condition": ["reprogrammed"] * 3,
        }
    )
    expression_a = pd.DataFrame(
        [[5.0, 5.2, 4.9], [0.1, 0.2, 0.3], [3.0, 2.8, 3.1]],
        index=genes,
        columns=samples,
    )
    expression_b = expression_a.copy()
    expression_b.loc[:, "O3_d10"] = [1000.0, -1000.0, 5000.0]

    kwargs = dict(
        metadata=metadata,
        identity_genes=["COL1A1"],
        pluripotency_genes=["POU5F1"],
        rejuvenation_genes=["CDKN2A"],
        identity_floor=-10,
        pluripotency_ceiling=10,
        min_rejuvenation=-10,
        discovery_donors=["O1", "O2"],
        min_donors_per_day=2,
    )
    result_a = define_safe_window(expression=expression_a, **kwargs)
    result_b = define_safe_window(expression=expression_b, **kwargs)
    cols = ["identity_score", "pluripotency_score", "rejuvenation_score"]
    scores_a = result_a.sample_scores.set_index("sample_id").loc[["O1_d10", "O2_d10"], cols]
    scores_b = result_b.sample_scores.set_index("sample_id").loc[["O1_d10", "O2_d10"], cols]
    assert np.allclose(scores_a, scores_b)
    assert set(result_a.selected_samples) == set(result_b.selected_samples)


def test_two_donor_signature_is_marked_exploratory() -> None:
    from ccid.signature import build_matched_methylation_target_signature

    metadata = pd.DataFrame(
        {
            "sample_id": ["O1_safe", "O1_ctrl", "O2_safe", "O2_ctrl"],
            "donor_id": ["O1", "O1", "O2", "O2"],
            "day": [10, 10, 10, 10],
            "condition": ["reprogrammed", "baseline", "reprogrammed", "baseline"],
            "batch": ["exp2"] * 4,
            "stage": ["endpoint"] * 4,
        }
    )
    beta = pd.DataFrame(
        {
            "O1_safe": [0.4, 0.6],
            "O1_ctrl": [0.7, 0.5],
            "O2_safe": [0.4, 0.6],
            "O2_ctrl": [0.7, 0.5],
        },
        index=["cg1", "cg2"],
    )
    result = build_matched_methylation_target_signature(
        beta,
        metadata,
        safe_samples=["O1_safe", "O2_safe"],
        min_abs_delta=0.01,
        fdr=1.0,
        min_sign_consistency=1.0,
        min_inferential_donors=3,
    )
    assert result.qc["inferential_status"] == "exploratory_underpowered"
    assert result.qc["n_strictly_selected_loci"] == 0
    assert result.qc["used_underpowered_fallback"] is True
    assert (result.all_loci["q_value"] == 1.0).all()
