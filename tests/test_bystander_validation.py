from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ccid.bystander_validation import (
    build_paired_celltype_effects,
    load_candidate_sets,
    run_bystander_validation,
    write_bystander_validation_outputs,
)


def _matrix_and_meta(cell: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    loci = [f"cg{i}" for i in range(1, 21)]
    columns: dict[str, list[float]] = {}
    meta_rows = []
    for idx, pair in enumerate(["R1", "R2", "R3"], start=1):
        control = f"{cell}_control_{pair}"
        targeted = f"{cell}_targeted_{pair}"
        baseline = np.linspace(0.20, 0.70, len(loci)) + idx * 0.001
        effect = np.array([0.10, 0.08, 0.06, -0.04, -0.03] + [0.002] * 15)
        columns[control] = baseline.tolist()
        columns[targeted] = (baseline + effect).tolist()
        meta_rows.extend(
            [
                {"sample_id": control, "treatment": "control", "replicate": pair},
                {"sample_id": targeted, "treatment": "targeted", "replicate": pair},
            ]
        )
    return pd.DataFrame(columns, index=loci), pd.DataFrame(meta_rows)


def test_build_paired_celltype_effects() -> None:
    beta, meta = _matrix_and_meta("t")
    result = build_paired_celltype_effects(beta, meta, cell_type="t_cell")
    assert len(result.pair_manifest) == 3
    first = result.locus_effects.set_index("cpg_id").loc["cg1"]
    assert first["editor_delta_beta"] == pytest.approx(0.10)
    assert first["sign_consistency"] == pytest.approx(1.0)
    assert bool(first["significant_editor_effect"])


def test_load_candidate_sets_accepts_nested_files(tmp_path: Path) -> None:
    pd.DataFrame({"cpg_id": ["cg1", "cg2"]}).to_csv(tmp_path / "stable.csv", index=False)
    pd.DataFrame({"cpg_id": ["cg1", "cg2", "cg3"]}).to_csv(tmp_path / "consensus.csv", index=False)
    pd.DataFrame({"cpg_id": ["cg1", "cg2", "cg3", "cg4"]}).to_csv(tmp_path / "tier1.csv", index=False)
    pd.DataFrame(
        {"cpg_id": [f"cg{i}" for i in range(1, 6)], "target_delta_beta": [-0.1] * 5}
    ).to_csv(tmp_path / "frozen.csv", index=False)
    sets, target = load_candidate_sets(
        stable_core_path=tmp_path / "stable.csv",
        consensus_2of3_path=tmp_path / "consensus.csv",
        tier1_path=tmp_path / "tier1.csv",
        frozen_signature_path=tmp_path / "frozen.csv",
    )
    assert len(sets["stable_core_3of3"]) == 2
    assert target.shape[0] == 5


def test_full_bystander_validation_and_outputs(tmp_path: Path) -> None:
    t_beta, t_meta = _matrix_and_meta("t")
    m_beta, m_meta = _matrix_and_meta("m")
    # Slight cell-type shift while retaining direction.
    m_beta.loc[:, [c for c in m_beta if "targeted" in c]] += 0.01

    candidate_sets = {
        "stable_core_3of3": {"cg1", "cg2"},
        "consensus_at_least_2of3": {"cg1", "cg2", "cg3"},
        "o3_tier1": {"cg1", "cg2", "cg3", "cg4", "cg5"},
        "frozen_500": {f"cg{i}" for i in range(1, 11)},
    }
    target_map = pd.DataFrame(
        {"cpg_id": [f"cg{i}" for i in range(1, 21)], "target_delta_beta": [-0.2] * 20}
    )
    annotation = pd.DataFrame(
        {
            "cpg_id": [f"cg{i}" for i in range(1, 21)],
            "chromosome": ["1"] * 20,
            "position": np.arange(1000, 21000, 1000),
            "gene": ["GENE1"] * 20,
            "island_relation": ["OpenSea"] * 20,
        }
    )
    result = run_bystander_validation(
        tcell_beta=t_beta,
        tcell_metadata=t_meta,
        msc_beta=m_beta,
        msc_metadata=m_meta,
        candidate_sets=candidate_sets,
        target_map=target_map,
        discovery_all_loci=target_map,
        probe_annotation=annotation,
        null_iterations=10,
        null_bins=3,
        seed=7,
    )
    assert result.summary["n_tcell_pairs"] == 3
    assert result.summary["n_msc_pairs"] == 3
    assert len(result.regional_scores) > 0
    stable = result.cross_celltype_summary.set_index("candidate_set").loc["stable_core_3of3"]
    assert stable["cross_cell_sign_concordance"] == pytest.approx(1.0)
    assert not result.null_tests.empty

    output = tmp_path / "out"
    write_bystander_validation_outputs(result, output)
    expected = [
        "candidate_probe_coverage.csv",
        "candidate_bystander_effects.csv.gz",
        "cross_celltype_replication.csv",
        "empirical_null_tests.csv.gz",
        "regional_editability_scores.csv",
        "bystander_validation_summary.json",
        "figures/candidate_coverage.png",
        "figures/tcell_msc_effect_correlation.png",
    ]
    for relative in expected:
        assert (output / relative).exists(), relative


def test_candidate_null_is_deterministic() -> None:
    t_beta, t_meta = _matrix_and_meta("t")
    m_beta, m_meta = _matrix_and_meta("m")
    candidate_sets = {
        "stable_core_3of3": {"cg1", "cg2"},
        "consensus_at_least_2of3": {"cg1", "cg2", "cg3"},
        "o3_tier1": {"cg1", "cg2", "cg3", "cg4", "cg5"},
        "frozen_500": {f"cg{i}" for i in range(1, 11)},
    }
    target = pd.DataFrame({"cpg_id": [f"cg{i}" for i in range(1, 21)], "target_delta_beta": [-0.2] * 20})
    kwargs = dict(
        tcell_beta=t_beta,
        tcell_metadata=t_meta,
        msc_beta=m_beta,
        msc_metadata=m_meta,
        candidate_sets=candidate_sets,
        target_map=target,
        discovery_all_loci=target,
        null_iterations=5,
        null_bins=2,
        seed=17,
    )
    a = run_bystander_validation(**kwargs)
    b = run_bystander_validation(**kwargs)
    pd.testing.assert_frame_equal(a.null_tests, b.null_tests)
