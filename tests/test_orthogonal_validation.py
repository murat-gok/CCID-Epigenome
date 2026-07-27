from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from ccid.orthogonal_validation import (
    enumerate_spcas9_guides,
    integrate_orthogonal_evidence,
    read_geo_series_matrix,
    validate_external_reprogramming,
    validate_fibroblast_eqtm,
    validate_natural_aging,
)


def write_geo_matrix(path: Path, titles: list[str], rows: dict[str, list[float]]) -> None:
    sample_ids = [f"GSM{i+1}" for i in range(len(titles))]
    lines = [
        '!Series_title\t"synthetic"',
        "!Sample_title\t" + "\t".join(f'"{x}"' for x in titles),
        "!Sample_geo_accession\t" + "\t".join(f'"{x}"' for x in sample_ids),
        "!series_matrix_table_begin",
        '"ID_REF"\t' + "\t".join(f'"{x}"' for x in sample_ids),
    ]
    for feature, values in rows.items():
        lines.append(f'"{feature}"\t' + "\t".join(str(x) for x in values))
    lines.append("!series_matrix_table_end")
    text = "\n".join(lines) + "\n"
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)


def write_inputs(tmp_path: Path) -> tuple[Path, Path]:
    candidates = pd.DataFrame(
        {
            "cpg_id": ["cg00000001", "cg00000002"],
            "target_delta_beta": [0.10, -0.08],
            "candidate_rank": [1, 2],
            "chromosome": ["1", "2"],
            "position": [101, 201],
        }
    )
    candidate_path = tmp_path / "candidates.csv"
    candidates.to_csv(candidate_path, index=False)
    frozen = pd.DataFrame(
        {
            "robust_rank": [1, 2],
            "biological_shortlist_rank": [1, 2],
            "cpg_id": ["cg00000001", "cg00000002"],
            "direct_gene": ["GENEA", "GENEB"],
            "nearest_tss_gene": ["GENEA", "GENEB"],
            "candidate_identity": ["A", "B"],
            "panel_role": ["Experimental target", "Experimental target"],
        }
    )
    frozen_path = tmp_path / "frozen.csv"
    frozen.to_csv(frozen_path, index=False)
    return candidate_path, frozen_path


def test_geo_parser_filters_rows(tmp_path: Path) -> None:
    matrix = tmp_path / "test.txt.gz"
    write_geo_matrix(matrix, ["A", "B"], {"cg00000001": [0.1, 0.2], "cg00000002": [0.3, 0.4]})
    result = read_geo_series_matrix(matrix, target_ids=["cg00000002"])
    assert result.values.index.tolist() == ["cg00000002"]
    assert result.metadata["sample_title"].tolist() == ["A", "B"]


def test_external_reprogramming_direction_validation(tmp_path: Path) -> None:
    candidate_path, frozen_path = write_inputs(tmp_path)
    titles = []
    for person in range(1, 5):
        titles += [f"Person {person} SkinFibroblasts Normal", f"Person {person} SkinFibroblasts Treated"]
    matrix = tmp_path / "external.txt.gz"
    cpg1 = []
    cpg2 = []
    for i in range(4):
        cpg1 += [0.40 + i * 0.01, 0.50 + i * 0.01]
        cpg2 += [0.60 - i * 0.01, 0.52 - i * 0.01]
    write_geo_matrix(matrix, titles, {"cg00000001": cpg1, "cg00000002": cpg2})
    result, summary = validate_external_reprogramming(
        candidate_path, frozen_path, matrix, bootstrap_iterations=100, seed=1
    )
    assert result["external_mean_direction_concordant"].all()
    assert summary["direction_concordant_candidates"] == 2
    assert summary["n_paired_donors"] == 4


def test_natural_aging_three_datasets(tmp_path: Path) -> None:
    candidate_path, frozen_path = write_inputs(tmp_path)
    rep = tmp_path / "rep.txt.gz"
    rep_titles = ["HDF_NG_0days", "HDF_NG_10days", "HDF_NG_20days", "HDF_HG_20days"]
    write_geo_matrix(
        rep,
        rep_titles,
        {
            "cg00000001": [0.60, 0.55, 0.50, 0.9],
            "cg00000002": [0.30, 0.34, 0.38, 0.1],
        },
    )
    age = tmp_path / "age.txt.gz"
    age_titles = ["HDF-arm_20yrs_P3", "HDF-arm_40yrs_P3", "HDF-arm_70yrs_P3", "HDF-arm_70yrs_P21"]
    write_geo_matrix(
        age,
        age_titles,
        {
            "cg00000001": [0.65, 0.58, 0.48, 0.2],
            "cg00000002": [0.25, 0.31, 0.40, 0.8],
        },
    )
    multi = tmp_path / "multi.txt.gz"
    multi_titles = []
    vals1 = []
    vals2 = []
    for donor in ["HC1", "HC2", "HC3"]:
        for passage in [1, 5, 10, 15]:
            multi_titles.append(f"Methylation_{donor}_sp2_Normal_Control_ox21_p{passage}")
            vals1.append(0.70 - 0.01 * passage)
            vals2.append(0.20 + 0.01 * passage)
    write_geo_matrix(multi, multi_titles, {"cg00000001": vals1, "cg00000002": vals2})
    result, summary = validate_natural_aging(
        candidate_path,
        frozen_path,
        replicative_matrix_file=rep,
        chronological_matrix_file=age,
        multidonor_matrix_file=multi,
    )
    assert result["replicative_reversal_concordant"].all()
    assert result["chronological_reversal_concordant"].all()
    assert result["multidonor_reversal_concordant"].all()
    assert summary["multidonor_reversal_concordant_candidates"] == 2


def test_fibroblast_eqtm(tmp_path: Path) -> None:
    candidate_path, frozen_path = write_inputs(tmp_path)
    samples = [f"GM{i:05d}" for i in range(1, 13)]
    methyl = tmp_path / "methyl.txt.gz"
    expression = tmp_path / "expression.txt.gz"
    x = np.linspace(0.1, 0.9, len(samples))
    write_geo_matrix(methyl, samples, {"cg00000001": x.tolist(), "cg00000002": x[::-1].tolist()})
    write_geo_matrix(
        expression,
        [f"{s} [gene expression]" for s in samples],
        {"probeA": (2 - x).tolist(), "probeB": (1 + x).tolist()},
    )
    annotation = tmp_path / "GPL6883.annot.gz"
    with gzip.open(annotation, "wt") as handle:
        handle.write("ID\tGene Symbol\nprobeA\tGENEA\nprobeB\tGENEB\n")
    result, summary = validate_fibroblast_eqtm(
        candidate_path,
        frozen_path,
        methyl,
        expression,
        annotation,
        frozen_path,
    )
    assert summary["n_matched_fibroblast_samples"] == 12
    assert summary["n_tested_pairs"] == 2
    assert set(result["gene"]) == {"GENEA", "GENEB"}
    assert (result["eqtm_spearman_rho"].abs() > 0.9).all()


def test_spcas9_guide_enumeration() -> None:
    sequence = "A" * 30 + "GCGCGCGCGCGCGCGCGCGC" + "AGG" + "T" * 30
    guides = enumerate_spcas9_guides(
        sequence, sequence_start_1based=1, cpg_position_1based=42, max_distance=100
    )
    assert not guides.empty
    assert "basic_sequence_pass" in guides.columns


def test_integrated_evidence_tiers(tmp_path: Path) -> None:
    _, frozen_path = write_inputs(tmp_path)
    external = pd.DataFrame(
        {
            "cpg_id": ["cg00000001", "cg00000002"],
            "external_mean_direction_concordant": [True, False],
        }
    )
    aging = pd.DataFrame(
        {
            "cpg_id": ["cg00000001", "cg00000002"],
            "replicative_reversal_concordant": [True, True],
            "chronological_reversal_concordant": [True, False],
            "multidonor_reversal_concordant": [True, False],
        }
    )
    guides = pd.DataFrame(
        {"cpg_id": ["cg00000001", "cg00000002"], "guide_feasible": [True, True]}
    )
    external_path = tmp_path / "external.csv"
    aging_path = tmp_path / "aging.csv"
    guide_path = tmp_path / "guides.csv"
    external.to_csv(external_path, index=False)
    aging.to_csv(aging_path, index=False)
    guides.to_csv(guide_path, index=False)
    result, summary = integrate_orthogonal_evidence(
        frozen_path, external_path, aging_path, guide_summary_file=guide_path
    )
    tier = result.set_index("cpg_id")["orthogonal_validation_tier"]
    assert tier["cg00000001"].startswith("V1")
    assert summary["n_candidates"] == 2


def test_read_platform_annotation_plain_geo_data_table(tmp_path):
    from ccid.orthogonal_validation import read_platform_annotation

    path = tmp_path / "GPL6883_data.txt"
    path.write_text(
        "ID\tSymbol\tGene Title\n"
        "ILMN_1\tTMCC1\ttransmembrane and coiled-coil domain family 1\n"
        "ILMN_2\tSSH1; AMOTL2\tmultiple annotation\n",
        encoding="utf-8",
    )
    result = read_platform_annotation(path)
    observed = set(map(tuple, result[["probe_id", "gene_symbol"]].to_records(index=False)))
    assert ("ILMN_1", "TMCC1") in observed
    assert ("ILMN_2", "SSH1") in observed
    assert ("ILMN_2", "AMOTL2") in observed


def test_read_platform_annotation_geo_soft_document(tmp_path):
    from ccid.orthogonal_validation import read_platform_annotation

    path = tmp_path / "GPL6883_data.txt"
    path.write_text(
        "^PLATFORM = GPL6883\n"
        "#ID = Unique identifier for the probe\n"
        "#Symbol = Gene symbol\n"
        "!Platform_title = Illumina HumanRef-8 v3.0\n"
        "!platform_table_begin\n"
        "ID\tSpecies\tSymbol\tGene Title\n"
        "ILMN_1\tHomo sapiens\tTMCC1\ttransmembrane and coiled-coil domain family 1\n"
        "ILMN_2\tHomo sapiens\tSSH1; AMOTL2\tmultiple annotation\n"
        "!platform_table_end\n",
        encoding="utf-8",
    )
    result = read_platform_annotation(path)
    observed = set(map(tuple, result[["probe_id", "gene_symbol"]].to_records(index=False)))
    assert ("ILMN_1", "TMCC1") in observed
    assert ("ILMN_2", "SSH1") in observed
    assert ("ILMN_2", "AMOTL2") in observed


def test_streaming_encode_re2g_overlap_uses_target_gene_and_score(tmp_path: Path) -> None:
    from ccid.orthogonal_validation import overlap_candidates_with_re2g

    candidate_path, frozen_path = write_inputs(tmp_path)
    re2g = tmp_path / "encode_links.bed.gz"
    header = [
        "chr", "start", "end", "name", "class", "activity_base",
        "normalized_dhs", "activity_base_squared", "TargetGene",
        "TargetGeneTSS", "TargetGeneExpression", "TargetGenePromoterActivityQuantile",
        "TargetGeneIsExpressed", "TargetGeneEnsembl_ID", "normalized_dhs_b",
        "distance", "isSelfPromoter", "powerlaw_contact", "powerlaw_contact_reference",
        "hic_contact", "hic_contact_pl_scaled", "hic_pseudocount",
        "hic_contact_pl_scaled_adj", "ABC.Score.Numerator", "ABC.Score",
        "powerlaw.Score.Numerator", "powerlaw.Score", "CellType",
        "hic_contact_squared", "Score",
    ]
    row1 = [
        "chr1", "99", "150", "enhancer|chr1:99-150", "enhancer", "0.1", "0.1", "0.01",
        "GENEA", "1000", "1.2", "0.8", "True", "ENSG1", "0.2", "900", "False",
        "0.1", "0.1", "0", "0", "0.1", "0.1", "0.01", "0.02", "0.03", "0.04",
        "fibroblast of skin", "0", "0.99",
    ]
    row2 = [
        "chr3", "10", "20", "enhancer|chr3:10-20", "enhancer", "0.1", "0.1", "0.01",
        "OTHER", "1000", "1.2", "0.8", "True", "ENSG2", "0.2", "900", "False",
        "0.1", "0.1", "0", "0", "0.1", "0.1", "0.01", "0.02", "0.03", "0.04",
        "fibroblast of skin", "0", "0.50",
    ]
    with gzip.open(re2g, "wt") as handle:
        handle.write("#" + "\t".join(header) + "\n")
        handle.write("\t".join(row1) + "\n")
        handle.write("\t".join(row2) + "\n")

    result = overlap_candidates_with_re2g(
        candidate_path, frozen_path, [re2g], biosample_names=["skin"]
    )
    linked = result[result["re2g_linked"].fillna(False)]
    assert linked["cpg_id"].tolist() == ["cg00000001"]
    assert linked.iloc[0]["re2g_target_gene"] == "GENEA"
    assert float(linked.iloc[0]["re2g_score"]) == 0.99
    assert result.loc[result["cpg_id"] == "cg00000002", "re2g_linked"].eq(False).all()


def test_gse179847_processed_csv_streaming_with_direct_titles(tmp_path: Path) -> None:
    from ccid.orthogonal_validation import read_gse179847_processed_matrix

    matrix = tmp_path / "GSE179847_Cell_lifespan_DNAm_processed_matrix.csv.gz"
    titles = [
        "HC1_sp2_Normal_Control_ox21_p1",
        "HC1_sp2_Normal_Control_ox21_p5",
        "HC1_sp2_Normal_Control_ox21_p10",
        "HC2_sp2_Normal_Control_ox21_p1",
        "HC2_sp2_Normal_Control_ox21_p5",
        "HC2_sp2_Normal_Control_ox21_p10",
    ]
    frame = pd.DataFrame(
        [
            ["cg00000001", 0.70, 0.65, 0.60, 0.72, 0.67, 0.62],
            ["cg00000002", 0.20, 0.25, 0.30, 0.22, 0.27, 0.32],
            ["cg99999999", 0.40, 0.40, 0.40, 0.40, 0.40, 0.40],
        ],
        columns=["ID_REF", *titles],
    )
    frame.to_csv(matrix, index=False, compression="gzip")
    result = read_gse179847_processed_matrix(
        matrix, target_ids=["cg00000001", "cg00000002"]
    )
    assert result.values.index.tolist() == ["cg00000001", "cg00000002"]
    assert result.metadata["sample_title"].tolist() == titles


def test_natural_aging_accepts_gse179847_processed_csv(tmp_path: Path) -> None:
    candidate_path, frozen_path = write_inputs(tmp_path)
    matrix = tmp_path / "GSE179847_Cell_lifespan_DNAm_processed_matrix.csv.gz"
    titles: list[str] = []
    vals1: list[float] = []
    vals2: list[float] = []
    for donor in ["HC1", "HC2", "HC3"]:
        for passage in [1, 5, 10, 15]:
            titles.append(f"{donor}_sp2_Normal_Control_ox21_p{passage}")
            vals1.append(0.70 - 0.01 * passage)
            vals2.append(0.20 + 0.01 * passage)
    pd.DataFrame(
        [["cg00000001", *vals1], ["cg00000002", *vals2]],
        columns=["ID_REF", *titles],
    ).to_csv(matrix, index=False, compression="gzip")
    result, summary = validate_natural_aging(
        candidate_path,
        frozen_path,
        multidonor_matrix_file=matrix,
    )
    assert result["multidonor_reversal_concordant"].all()
    assert summary["multidonor_available_candidates"] == 2


def test_natural_aging_rejects_metadata_only_multidonor_matrix(tmp_path: Path) -> None:
    candidate_path, frozen_path = write_inputs(tmp_path)
    matrix = tmp_path / "GSE179847_series_matrix.txt.gz"
    write_geo_matrix(matrix, ["Methylation_HC1_sp2_Normal_Control_ox21_p1"], {})
    import pytest
    with pytest.raises(ValueError, match="processed_matrix"):
        validate_natural_aging(
            candidate_path,
            frozen_path,
            multidonor_matrix_file=matrix,
        )
