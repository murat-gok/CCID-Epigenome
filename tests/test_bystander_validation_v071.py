from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ccid.bystander_validation_v071 import (
    build_candidate_metadata,
    build_platform_aware_tcell_effects,
    build_positive_control_analysis,
    build_unique_regions,
    run_bystander_validation_v071,
    summarise_matched_nulls,
    write_bystander_validation_v071_outputs,
)


def _write_raw(path: Path, pairs: list[str], loci: list[str], extra_locus: str | None = None) -> None:
    ids = loci + ([extra_locus] if extra_locus else [])
    data: dict[str, object] = {"probeID": ids}
    for idx, pair in enumerate(pairs, start=1):
        baseline = np.linspace(0.2, 0.6, len(ids)) + idx * 0.001
        effect = np.array([0.08, 0.06, -0.04, 0.02, -0.01] + [0.003] * max(0, len(ids) - 5))
        data[f"{pair}_Scramble.beta"] = baseline
        data[f"{pair}_Scramble.detp"] = np.zeros(len(ids))
        data[f"{pair}_Hyper.beta"] = baseline + effect
        data[f"{pair}_Hyper.detp"] = np.zeros(len(ids))
    pd.DataFrame(data).to_csv(path, index=False, compression="gzip")


def _msc() -> tuple[pd.DataFrame, pd.DataFrame]:
    loci = [f"cg{i:08d}" for i in range(1, 16)]
    cols = {}
    meta = []
    for idx, donor in enumerate(["MSC170", "MSC179", "MSC182"]):
        c = f"{donor}_Scramble"
        t = f"{donor}_Hyper"
        base = np.linspace(0.25, 0.65, len(loci)) + idx * 0.001
        effect = np.array([0.09, 0.07, -0.03, 0.03, -0.02] + [0.004] * 10)
        cols[c] = base
        cols[t] = base + effect
        meta += [
            {"sample_id": c, "treatment": "control", "replicate": donor},
            {"sample_id": t, "treatment": "targeted", "replicate": donor},
        ]
    return pd.DataFrame(cols, index=loci), pd.DataFrame(meta)


def _annotation() -> pd.DataFrame:
    loci = [f"cg{i:08d}" for i in range(1, 17)]
    return pd.DataFrame(
        {
            "cpg_id": loci,
            "chromosome": ["1"] * 8 + ["2"] * 8,
            "position": [1000, 1200, 1600, 5000, 5100, 10000, 13000, 17000] * 2,
            "gene": ["ELOVL2", "ELOVL2", "ELOVL2", "GENE1", "GENE1", "FHL2", "", ""] * 2,
            "island_relation": ["OpenSea"] * 16,
        }
    )


def test_platform_aware_union_preserves_r3_only_probe(tmp_path: Path) -> None:
    loci = [f"cg{i:08d}" for i in range(1, 11)]
    p12 = tmp_path / "GSE269758_R12.csv.gz"
    p3 = tmp_path / "GSE269758_R3.csv.gz"
    _write_raw(p12, ["R1", "R2"], loci)
    _write_raw(p3, ["R3"], loci, extra_locus="cg99999999")
    result = build_platform_aware_tcell_effects([p12, p3])
    row = result.locus_effects.set_index("cpg_id").loc["cg99999999"]
    assert row["n_pairs"] == 1
    assert len(result.pair_manifest) == 3


def test_positive_controls_and_unique_regions() -> None:
    loci = [f"cg{i:08d}" for i in range(1, 16)]
    effects = pd.DataFrame(
        {
            "cpg_id": loci,
            "editor_delta_beta": np.linspace(0.01, 0.10, len(loci)),
            "sign_consistency": [1.0] * len(loci),
            "n_pairs": [3] * len(loci),
        }
    )
    probes, genes, decay = build_positive_control_analysis(_annotation(), effects, effects)
    assert {"ELOVL2", "FHL2"}.issubset(set(genes["positive_control_gene"]))
    assert not probes.empty
    assert not decay.empty

    sets = {
        "stable_core_3of3": {loci[0], loci[1]},
        "consensus_at_least_2of3": {loci[0], loci[1], loci[2]},
        "o3_tier1": set(loci[:6]),
        "frozen_500": set(loci[:10]),
    }
    meta = pd.DataFrame({"cpg_id": loci, "target_delta_beta": [-0.2] * len(loci)})
    membership, all_regions, multi = build_unique_regions(
        sets, meta, _annotation(), effects, effects, cluster_gap_bp=500
    )
    assert membership["cpg_id"].nunique() == 10
    assert all_regions["unique_region_id"].is_unique
    assert (multi["n_candidate_cpgs"] >= 2).all()


def test_matched_null_summary_has_bh_columns() -> None:
    rows = []
    for cell in ["t_cell", "mesenchymal_stromal_cell"]:
        for metric, p in [("mean_abs_delta", 0.01), ("median_abs_delta", 0.20)]:
            for i in range(5):
                rows.append(
                    {
                        "candidate_set": "stable_core_3of3",
                        "cell_type": cell,
                        "metric": metric,
                        "observed": 0.2,
                        "null_value": 0.1 + i * 0.01,
                        "empirical_p": p,
                        "n_candidates_measured": 10,
                    }
                )
    summary = summarise_matched_nulls(pd.DataFrame(rows))
    assert {"bh_q_global", "bh_q_within_celltype"}.issubset(summary.columns)


def test_full_v071_pipeline_and_outputs(tmp_path: Path) -> None:
    loci = [f"cg{i:08d}" for i in range(1, 16)]
    p12 = tmp_path / "GSE269758_T_cell_hyper.csv.gz"
    p3 = tmp_path / "GSE269758_T_cell_R3_Hyperlibrary.csv.gz"
    _write_raw(p12, ["R1", "R2"], loci)
    _write_raw(p3, ["R3"], loci, extra_locus="cg00000016")
    msc_beta, msc_meta = _msc()
    sets = {
        "stable_core_3of3": set(loci[:2]),
        "consensus_at_least_2of3": set(loci[:4]),
        "o3_tier1": set(loci[:7]),
        "frozen_500": set(loci[:10]),
    }
    target = pd.DataFrame({"cpg_id": loci, "target_delta_beta": [-0.2] * len(loci)})
    candidate_meta = target.copy()
    candidate_meta["rank_score"] = np.arange(len(loci), 0, -1)
    candidate_meta["n_lodo_rotations"] = 3
    result = run_bystander_validation_v071(
        tcell_raw_paths=[p12, p3],
        msc_beta=msc_beta,
        msc_metadata=msc_meta,
        candidate_sets=sets,
        target_map=target,
        candidate_metadata=candidate_meta,
        discovery_all_loci=target,
        probe_annotation=_annotation(),
        null_iterations=5,
        null_bins=2,
        cluster_gap_bp=500,
        seed=9,
    )
    assert result.summary["n_tcell_pairs"] == 3
    assert result.summary["n_unique_regions_all"] > 0
    assert "bh_q_global" in result.matched_null_summary.columns
    output = tmp_path / "out"
    write_bystander_validation_v071_outputs(result, output)
    expected = [
        "bystander_validation_v071_summary.json",
        "tcell_platform_meta_effects.csv.gz",
        "positive_control_gene_summary.csv",
        "unique_region_scores_multicpg.csv",
        "matched_null_summary_bh.csv",
        "cross_celltype_null_summary.csv",
        "figures/matched_null_bh_summary.png",
    ]
    for name in expected:
        assert (output / name).exists(), name


def test_build_candidate_metadata(tmp_path: Path) -> None:
    pd.DataFrame({"cpg_id": ["cg1"]}).to_csv(tmp_path / "stable.csv", index=False)
    pd.DataFrame({"cpg_id": ["cg1", "cg2"], "n_lodo_rotations": [3, 2]}).to_csv(tmp_path / "cons.csv", index=False)
    pd.DataFrame({"cpg_id": ["cg1", "cg2", "cg3"]}).to_csv(tmp_path / "tier.csv", index=False)
    pd.DataFrame(
        {"cpg_id": ["cg1", "cg2", "cg3"], "target_delta_beta": [-0.1, -0.2, -0.3], "rank_score": [3, 2, 1]}
    ).to_csv(tmp_path / "frozen.csv", index=False)
    sets, target, meta = build_candidate_metadata(
        stable_core_path=tmp_path / "stable.csv",
        consensus_2of3_path=tmp_path / "cons.csv",
        tier1_path=tmp_path / "tier.csv",
        frozen_signature_path=tmp_path / "frozen.csv",
    )
    assert sets["stable_core_3of3"] == {"cg1"}
    assert target.shape[0] == 3
    assert "in_stable_core_3of3" in meta.columns
