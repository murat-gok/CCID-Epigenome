from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ccid.genomic_causality import (
    DEFAULT_WEIGHTS,
    category_enrichment,
    load_gene_evidence,
    run_genomic_causality_scoring,
    score_candidate_cpgs,
    write_genomic_causality_outputs,
)


def _membership() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cpg_id": [f"cg{i:08d}" for i in range(1, 9)],
            "chromosome": ["1"] * 4 + ["2"] * 4,
            "position": [1000, 1200, 5000, 9000, 1000, 1150, 8000, 12000],
            "gene": ["SIRT1", "SIRT1", "CDKN2A", "GENE1", "KLF4", "KLF4", "GENE2", ""],
            "gene_group": ["TSS200", "Body", "TSS1500", "Body", "TSS200", "Body", "Body", ""],
            "island_relation": ["Island", "N_Shore", "Island", "OpenSea", "Island", "S_Shore", "OpenSea", "OpenSea"],
            "in_stable_core_3of3": [True, True, False, False, True, True, False, False],
            "in_consensus_at_least_2of3": [True, True, True, False, True, True, True, False],
            "in_o3_tier1": [True, True, True, True, True, True, True, False],
            "in_frozen_500": [True] * 8,
            "target_delta_beta": [-0.6, -0.5, -0.4, -0.3, -0.5, -0.45, -0.2, -0.1],
            "rank_score": [8, 7, 6, 5, 4, 3, 2, 1],
            "n_lodo_rotations": [3, 3, 2, 1, 3, 3, 2, 1],
            "unique_region_id": ["chr1:1000-1200", "chr1:1000-1200", "chr1:5000-5000", "chr1:9000-9000", "chr2:1000-1150", "chr2:1000-1150", "chr2:8000-8000", "chr2:12000-12000"],
            "tcell_delta": [-0.04, -0.03, -0.01, 0.01, -0.05, -0.04, 0.005, np.nan],
            "tcell_consistency": [1.0, 1.0, 0.67, 0.67, 1.0, 1.0, 0.67, np.nan],
            "msc_delta": [-0.06, -0.05, -0.02, 0.02, -0.04, -0.03, 0.01, np.nan],
            "msc_consistency": [1.0, 1.0, 0.67, 0.67, 1.0, 1.0, 0.67, np.nan],
        }
    )


def _regions() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "unique_region_id": ["chr1:1000-1200", "chr1:5000-5000", "chr1:9000-9000", "chr2:1000-1150", "chr2:8000-8000", "chr2:12000-12000"],
            "n_candidate_cpgs": [2, 1, 1, 2, 1, 1],
            "region_class": ["multi_cpg_region", "singleton_locus", "singleton_locus", "multi_cpg_region", "singleton_locus", "singleton_locus"],
            "regional_responsiveness_score": [0.08, 0.02, 0.01, 0.07, 0.01, np.nan],
        }
    )


def _annotation() -> pd.DataFrame:
    frame = _membership()[["cpg_id", "chromosome", "position", "gene", "gene_group", "island_relation"]].copy()
    frame["genome_build"] = "hg38"
    return frame


def _evidence(tmp_path: Path) -> pd.DataFrame:
    path = tmp_path / "evidence.csv"
    pd.DataFrame(
        {
            "gene": ["SIRT1", "CDKN2A", "CDKN2A", "KLF4", "KLF4"],
            "category": ["aging_longevity", "senescence_sasp", "tumor_suppressor_risk", "pluripotency_reprogramming", "pluripotency_reprogramming_risk"],
            "relevance_weight": [0.9, 0.9, 0.0, 0.9, 0.0],
            "risk_weight": [0.0, 0.0, 0.95, 0.0, 0.7],
            "source_id": ["test"] * 5,
        }
    ).to_csv(path, index=False)
    return load_gene_evidence(path)


def _optional_evidence() -> tuple[pd.DataFrame, pd.DataFrame]:
    svf = pd.DataFrame(
        {
            "cpg_id": [f"cg{i:08d}" for i in range(1, 9)],
            "discovery_delta_beta": [-0.6, -0.5, -0.4, -0.3, -0.5, -0.45, -0.2, -0.1],
            "validation_delta_beta": [-0.5, -0.4, -0.3, -0.1, -0.4, -0.4, 0.02, 0.0],
            "same_direction": [True, True, True, True, True, True, False, False],
        }
    )
    ts = pd.DataFrame(
        {
            "cpg_id": [f"cg{i:08d}" for i in range(1, 9)],
            "success_sign_consistency": [1, 1, 0.9, 0.8, 1, 1, 0.6, 0.5],
            "failed_sign_consistency": [0.2, 0.2, 0.4, 0.5, 0.2, 0.2, 0.6, 0.5],
            "o3_success_sign_consistency": [1, 1, 1, 0.8, 1, 1, 0.5, 0.5],
            "earliest_o3_supported_day": [10, 10, 10, 13, 10, 10, 15, np.nan],
        }
    )
    return svf, ts


def test_candidate_scoring_rewards_stable_supported_loci(tmp_path: Path) -> None:
    svf, ts = _optional_evidence()
    scored, gene_map = score_candidate_cpgs(
        _membership(),
        _regions(),
        _evidence(tmp_path),
        successful_vs_failed=svf,
        time_series=ts,
    )
    assert len(scored) == 8
    assert not gene_map.empty
    top = scored.iloc[0]
    assert top["cpg_id"] in {"cg00000001", "cg00000002"}
    assert top["causal_evidence_score"] > 60
    assert scored["evidence_completeness"].between(0, 1).all()
    risky = scored.set_index("cpg_id").loc["cg00000005"]
    assert risky["safety_risk_score"] > 0
    assert risky["safety_adjusted_priority_score"] < risky["causal_evidence_score"]


def test_enrichment_and_full_outputs(tmp_path: Path) -> None:
    svf, ts = _optional_evidence()
    result = run_genomic_causality_scoring(
        membership=_membership(),
        region_scores=_regions(),
        annotation=_annotation(),
        gene_evidence=_evidence(tmp_path),
        successful_vs_failed=svf,
        time_series=ts,
        weights=DEFAULT_WEIGHTS,
        shortlist_size=5,
    )
    assert result.summary["n_candidate_cpgs"] == 8
    assert result.summary["n_unique_regions"] == 6
    assert len(result.shortlist) == 5
    assert "region_causal_priority_score" in result.region_scores.columns
    assert "bh_q_global" in result.category_enrichment.columns
    assert not result.safety_flags.empty

    out = tmp_path / "out"
    write_genomic_causality_outputs(result, out)
    expected = [
        "candidate_cpg_annotation_causality_scores.csv.gz",
        "candidate_region_annotation_causality_scores.csv",
        "gene_evidence_scores.csv",
        "biological_category_enrichment.csv",
        "safety_flags.csv",
        "evidence_component_summary.csv",
        "top_candidate_shortlist.csv",
        "genomic_causality_v080_summary.json",
        "genomic_causality_v080_manifest.json",
        "figures/cpg_score_distribution.png",
        "figures/top_cpg_evidence_heatmap.png",
        "figures/top_candidate_regions.png",
    ]
    for name in expected:
        assert (out / name).exists(), name


def test_missing_optional_evidence_is_reported_not_imputed(tmp_path: Path) -> None:
    scored, _ = score_candidate_cpgs(_membership(), _regions(), _evidence(tmp_path))
    assert scored["successful_failed_specificity_score"].isna().all()
    assert scored["temporal_consistency_score"].isna().all()
    assert (scored["evidence_completeness"] < 1).all()
