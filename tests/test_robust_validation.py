from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ccid.robust_validation import (
    direction_balanced_effect_matched_null,
    leave_one_donor_out_analysis,
    run_robust_validation,
    time_series_consistency_analysis,
)


def _metadata() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for donor in ("O1", "O2", "O3"):
        for day in (10, 13):
            rows.extend(
                [
                    {
                        "sample_id": f"{donor}_success_d{day}",
                        "donor_id": donor,
                        "day": day,
                        "condition": "reprogrammed",
                        "batch": "exp2",
                        "stage": "endpoint",
                    },
                    {
                        "sample_id": f"{donor}_baseline_d{day}",
                        "donor_id": donor,
                        "day": day,
                        "condition": "baseline",
                        "batch": "exp2",
                        "stage": "endpoint",
                    },
                    {
                        "sample_id": f"{donor}_failed_d{day}",
                        "donor_id": donor,
                        "day": day,
                        "condition": "failed",
                        "batch": "exp2",
                        "stage": "endpoint",
                    },
                ]
            )
        rows.extend(
            [
                {
                    "sample_id": f"{donor}_intermediate_d10",
                    "donor_id": donor,
                    "day": 10,
                    "condition": "reprogrammed_intermediate",
                    "batch": "exp2",
                    "stage": "intermediate",
                },
                {
                    "sample_id": f"{donor}_baseline_intermediate_d10",
                    "donor_id": donor,
                    "day": 10,
                    "condition": "baseline",
                    "batch": "exp2",
                    "stage": "intermediate",
                },
            ]
        )
    return pd.DataFrame(rows)


def _pairs() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "pair_id": ["pair_001", "pair_002"],
            "safe_sample_id": ["O1_success_d10", "O2_success_d10"],
            "baseline_sample_ids": ["O1_baseline_d10", "O2_baseline_d10"],
            "n_baseline_samples": [1, 1],
            "donor_id": ["O1", "O2"],
            "condition": ["reprogrammed", "reprogrammed"],
            "day": [10, 10],
            "batch": ["exp2", "exp2"],
            "stage": ["endpoint", "endpoint"],
        }
    )


def _beta_and_signature(tmp_path: Path) -> tuple[Path, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(42)
    cpgs = [f"cg{i:04d}" for i in range(120)]
    base = rng.uniform(0.25, 0.75, size=len(cpgs))
    true_delta = np.concatenate(
        [np.linspace(-0.30, -0.08, 30), np.linspace(0.08, 0.20, 10), np.zeros(80)]
    )
    frame: dict[str, np.ndarray] = {}
    for donor_index, donor in enumerate(("O1", "O2", "O3")):
        donor_noise = rng.normal(0, 0.005, size=len(cpgs))
        for day, scale in ((10, 1.0), (13, 1.15)):
            baseline = np.clip(base + donor_noise + rng.normal(0, 0.003, len(cpgs)), 0, 1)
            success = np.clip(baseline + true_delta * scale + rng.normal(0, 0.006, len(cpgs)), 0, 1)
            failed = np.clip(baseline + rng.normal(0, 0.006, len(cpgs)), 0, 1)
            frame[f"{donor}_baseline_d{day}"] = baseline
            frame[f"{donor}_success_d{day}"] = success
            frame[f"{donor}_failed_d{day}"] = failed
        baseline_intermediate = np.clip(base + donor_noise, 0, 1)
        frame[f"{donor}_baseline_intermediate_d10"] = baseline_intermediate
        frame[f"{donor}_intermediate_d10"] = np.clip(
            baseline_intermediate + true_delta * 0.55 + rng.normal(0, 0.006, len(cpgs)), 0, 1
        )
    beta = pd.DataFrame(frame, index=cpgs)
    path = tmp_path / "beta.csv.gz"
    beta.to_csv(path, index=True, index_label="cpg_id")

    donor_deltas = []
    for donor in ("O1", "O2"):
        donor_deltas.append(beta[f"{donor}_success_d10"] - beta[f"{donor}_baseline_d10"])
    discovery = pd.concat(donor_deltas, axis=1).mean(axis=1)
    signature = pd.DataFrame(
        {
            "cpg_id": cpgs[:40],
            "target_delta_beta": discovery.iloc[:40].to_numpy(),
            "rank_score": np.arange(40, 0, -1),
            "weight": np.ones(40),
        }
    )
    all_loci = pd.DataFrame(
        {"cpg_id": cpgs, "target_delta_beta": discovery.to_numpy()}
    )
    return path, signature, all_loci


def test_direction_balanced_null_is_reproducible(tmp_path: Path) -> None:
    beta_path, signature, all_loci = _beta_and_signature(tmp_path)
    beta = pd.read_csv(beta_path, index_col=0)
    delta = beta["O3_success_d10"] - beta["O3_baseline_d10"]
    first = direction_balanced_effect_matched_null(
        signature, all_loci, delta, n_iterations=20, n_bins=4, seed=11
    )
    second = direction_balanced_effect_matched_null(
        signature, all_loci, delta, n_iterations=20, n_bins=4, seed=11
    )
    pd.testing.assert_frame_equal(first, second)
    assert set(first["analysis_scale"]) == {"raw"}
    assert (first["fallback_count"] >= 0).all()


def test_lodo_uses_only_non_holdout_donors(tmp_path: Path) -> None:
    beta_path, _, _ = _beta_and_signature(tmp_path)
    beta = pd.read_csv(beta_path, index_col=0)
    donor_deltas = pd.DataFrame(
        {
            donor: beta[f"{donor}_success_d10"] - beta[f"{donor}_baseline_d10"]
            for donor in ("O1", "O2", "O3")
        }
    )
    summary, candidates, overlap = leave_one_donor_out_analysis(
        donor_deltas, top_n=20, min_abs_delta=0.03
    )
    assert set(summary["holdout_donor"]) == {"O1", "O2", "O3"}
    assert len(candidates) == 60
    assert not candidates.apply(
        lambda row: row["holdout_donor"] in row["discovery_donors"].split(";"), axis=1
    ).any()
    assert "consensus" in set(overlap["holdout_a"])


def test_time_series_exact_baseline_matching(tmp_path: Path) -> None:
    beta_path, signature, _ = _beta_and_signature(tmp_path)
    metrics, consistency, manifest = time_series_consistency_analysis(
        beta_path,
        _metadata(),
        signature,
        donors=("O1", "O2", "O3"),
        batch="exp2",
        chunksize=25,
    )
    assert len(metrics) == 15
    assert len(manifest) == 15
    failed = metrics[metrics["condition"] == "failed"]
    successful = metrics[metrics["condition"] == "reprogrammed"]
    assert successful["sign_concordance"].mean() > failed["sign_concordance"].mean()
    assert len(consistency) == 40


def test_full_robust_validation_outputs(tmp_path: Path) -> None:
    beta_path, signature, all_loci = _beta_and_signature(tmp_path)
    result = run_robust_validation(
        beta_path=beta_path,
        metadata=_metadata(),
        signature=signature,
        discovery_pair_table=_pairs(),
        all_loci=all_loci,
        validation_donor="O3",
        donors=("O1", "O2", "O3"),
        lodo_top_n=20,
        lodo_min_abs_delta=0.03,
        null_iterations=20,
        null_bins=4,
        seed=123,
        chunksize=25,
    )
    assert not result.successful_vs_failed.empty
    assert set(result.baseline_sensitivity_metrics["comparison"]) == {
        "successful_vs_exact_baseline",
        "successful_vs_other_day_baseline_mean",
        "successful_vs_all_same_donor_baseline_mean",
        "successful_vs_failed_same_stratum",
    }
    assert set(result.direction_balanced_null["analysis_scale"]) == {
        "raw",
        "global_centered",
    }
    direct = result.baseline_sensitivity_metrics[
        result.baseline_sensitivity_metrics["comparison"]
        == "successful_vs_failed_same_stratum"
    ].iloc[0]
    assert direct["raw_sign_concordance"] > 0.8
    assert len(result.lodo_summary) == 3
    assert result.summary["candidate_list_frozen_before_validation"] is True
