from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ccid.manuscript_validation import run_manuscript_validation, write_manuscript_validation_outputs


COMPONENTS = [
    "discovery_strength_score",
    "donor_stability_score",
    "o3_support_score",
    "successful_failed_specificity_score",
    "temporal_consistency_score",
    "editor_responsiveness_score",
    "regional_support_score",
    "genomic_context_score",
    "regulatory_evidence_score",
]


def _candidates(n: int = 60) -> pd.DataFrame:
    rng = np.random.default_rng(11)
    latent = np.linspace(0.95, 0.15, n)
    data: dict[str, object] = {
        "cpg_id": [f"cg{i:08d}" for i in range(n)],
        "gene": ["TOPGENE" if i < 4 else "" for i in range(n)],
        "nearest_tss_gene": ["TOPGENE" if i < 4 else f"G{i}" for i in range(n)],
        "unique_region_id": [f"chr1:{1000 + i * 10}-{1000 + i * 10}" for i in range(n)],
        "in_stable_core_3of3": [i < 12 for i in range(n)],
        "in_consensus_at_least_2of3": [i < 25 for i in range(n)],
        "in_o3_tier1": [i < 45 for i in range(n)],
        "in_frozen_500": True,
        "n_candidate_cpgs": [2 if i < 10 else 1 for i in range(n)],
    }
    for j, col in enumerate(COMPONENTS):
        noise = rng.normal(0, 0.025 + j * 0.001, n)
        data[col] = np.clip(latent + noise, 0, 1)
    frame = pd.DataFrame(data)
    frame.loc[5:8, "donor_stability_score"] = np.nan
    return frame


def _config(iterations: int = 80) -> dict[str, object]:
    return {
        "top_k": 10,
        "secondary_top_k": 20,
        "weight_iterations": iterations,
        "uncertainty_iterations": iterations,
        "null_iterations": iterations,
        "shortlist_size": 20,
        "seed": 1234,
    }


def test_missing_components_are_unavailable_not_zero() -> None:
    frame = _candidates()
    result = run_manuscript_validation(frame, _config())
    row = result.candidate_stability.loc[result.candidate_stability["cpg_id"].eq("cg00000005")].iloc[0]
    assert row["v082_component_availability"] < 1.0
    assert np.isfinite(row["v082_baseline_score"])
    assert row["v082_baseline_score"] > 0


def test_reproducible_perturbations_and_stable_signal() -> None:
    frame = _candidates()
    a = run_manuscript_validation(frame, _config())
    b = run_manuscript_validation(frame, _config())
    pd.testing.assert_frame_equal(a.weight_iteration_summary, b.weight_iteration_summary)
    pd.testing.assert_frame_equal(a.uncertainty_iteration_summary, b.uncertainty_iteration_summary)
    first = a.candidate_stability.sort_values("robustness_rank").iloc[0]
    assert first["cpg_id"] in {"cg00000000", "cg00000001", "cg00000002", "cg00000003"}
    assert first["weight_top_k_frequency"] >= 0.9
    assert first["uncertainty_top_k_frequency"] >= 0.9


def test_ablation_null_redundancy_and_enrichment_outputs() -> None:
    result = run_manuscript_validation(_candidates(), _config(100))
    assert set(result.ablation_summary["removed_component"]) == set(COMPONENTS)
    assert len(result.ablation_candidate_shifts) == 60 * len(COMPONENTS)
    assert set(result.null_summary["statistic"]) == {"mean_top_k_score", "maximum_score"}
    assert result.null_summary["empirical_p"].between(0, 1).all()
    assert "effective_rank" in result.redundancy_summary.columns
    assert not result.candidate_set_enrichment.empty
    stable = result.candidate_stability["manuscript_stability_status"].eq("stable_primary_shortlist")
    assert stable.sum() > 0


def test_write_outputs(tmp_path: Path) -> None:
    result = run_manuscript_validation(_candidates(40), _config(30))
    write_manuscript_validation_outputs(result, tmp_path)
    expected = [
        "candidate_rank_stability_v082.csv.gz",
        "ablation_summary_v082.csv",
        "ranking_specificity_null_summary_v082.csv",
        "component_correlation_v082.csv",
        "candidate_set_enrichment_v082.csv",
        "manuscript_ready_robust_shortlist_v082.csv",
        "manuscript_validation_v082_summary.json",
        "manuscript_validation_v082_manifest.json",
        "figures/top_robust_candidates.png",
        "figures/weight_rank_stability.png",
        "figures/ablation_rank_impact.png",
        "figures/ranking_specificity_null.png",
        "figures/component_correlation.png",
    ]
    for rel in expected:
        assert (tmp_path / rel).exists(), rel
