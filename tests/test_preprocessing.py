from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ccid.dataio import read_matrix
from ccid.preprocessing import (
    build_overlap_report,
    parse_bystander_sample_name,
    parse_reprogramming_sample_name,
    preprocess_bystander_matrices,
    preprocess_fibroblast_expression,
    preprocess_paired_beta_matrix,
)


def _write_gzip_text(path: Path, text: str) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def test_parse_reprogramming_names() -> None:
    baseline = parse_reprogramming_sample_name(
        "O1_negative_control_15days_exp1", source_series="GSE165179"
    )
    assert baseline["donor_id"] == "O1"
    assert baseline["day"] == 15
    assert baseline["condition"] == "baseline"

    successful = parse_reprogramming_sample_name(
        "O3_transiently_reprogrammed_15days_exp1", source_series="GSE165179"
    )
    assert successful["condition"] == "reprogrammed"

    failed = parse_reprogramming_sample_name(
        "O2_failed_to_transiently_reprogram_13days_exp1", source_series="GSE165179"
    )
    assert failed["condition"] == "failed"

    intermediate = parse_reprogramming_sample_name(
        "O3_transient_reprogramming_intermediate_10days_exp2", source_series="GSE165179"
    )
    assert intermediate["condition"] == "reprogrammed_intermediate"


def test_parse_bystander_names() -> None:
    target = parse_bystander_sample_name(
        "DNMT3A_HyperCpGs_D21_R2", cell_type="t_cell", source_series="GSE269758"
    )
    assert target["treatment"] == "targeted"
    assert target["replicate"] == "R2"
    assert target["day"] == 21

    control = parse_bystander_sample_name(
        "MSC170_scrDNMT", cell_type="mesenchymal_stromal_cell", source_series="GSE285251"
    )
    assert control["treatment"] == "control"
    assert control["replicate"] == "MSC170"


def test_expression_preprocessing_merges_parts(tmp_path: Path) -> None:
    annotation = "Probe\tChromosome\tStart\tEnd\tProbe Strand\tFeature\tID\tDescription\tFeature Strand\tType\tFeature Orientation\tDistance"
    part1 = tmp_path / "part1.txt.gz"
    part2 = tmp_path / "part2.txt.gz"
    _write_gzip_text(
        part1,
        annotation
        + "\tO1_negative_control_15days_exp1\tO1_transiently_reprogrammed_15days_exp1\n"
        + "COL1A1\t1\t1\t2\t+\tCOL1A1\tENSG1\tdesc\t+\tgene\tmatch\t0\t1.0\t2.0\n"
        + "POU5F1\t1\t3\t4\t+\tPOU5F1\tENSG2\tdesc\t+\tgene\tmatch\t0\t0.0\t1.0\n",
    )
    _write_gzip_text(
        part2,
        annotation
        + "\tO2_negative_control_13days_exp2\tO2_transient_reprogramming_intermediate_13days_exp2\n"
        + "COL1A1\t1\t1\t2\t+\tCOL1A1\tENSG1\tdesc\t+\tgene\tmatch\t0\t1.5\t2.5\n"
        + "POU5F1\t1\t3\t4\t+\tPOU5F1\tENSG2\tdesc\t+\tgene\tmatch\t0\t0.1\t0.8\n",
    )

    expression_path = tmp_path / "expression.csv.gz"
    metadata_path = tmp_path / "metadata.csv"
    qc_path = tmp_path / "qc.json"
    expression, metadata, qc = preprocess_fibroblast_expression(
        [part1, part2], expression_path, metadata_path, qc_path
    )
    assert expression.shape == (2, 4)
    assert set(metadata["condition"]) == {
        "baseline",
        "reprogrammed",
        "reprogrammed_intermediate",
    }
    assert qc["n_samples"] == 4
    reloaded = read_matrix(expression_path)
    assert reloaded.shape == (2, 4)


def test_paired_beta_streaming_preprocessing(tmp_path: Path) -> None:
    raw = tmp_path / "mptr.csv.gz"
    _write_gzip_text(
        raw,
        "ID_REF,O1_negative_control_15days_exp1,Detection Pval,O1_transiently_reprogrammed_15days_exp1,Detection Pval,O2_negative_control_15days_exp1,Detection Pval,O2_transiently_reprogrammed_15days_exp1,Detection Pval\n"
        "cg1,0.2,0.0,0.4,0.0,0.3,0.0,0.5,0.0\n"
        "cg2,0.2,0.2,0.4,0.0,0.3,0.0,0.5,0.0\n"
        "cg3,1.2,0.0,0.4,0.0,0.3,0.0,0.5,0.0\n",
    )
    result = preprocess_paired_beta_matrix(
        raw,
        tmp_path / "beta.csv.gz",
        tmp_path / "meta.csv",
        tmp_path / "sample_qc.csv",
        tmp_path / "summary.json",
        source_series="GSE165179",
        max_locus_missing=0.20,
        chunksize=2,
    )
    # cg2 and cg3 each have one invalid value out of four and fail the 20% threshold.
    assert result.n_input_loci == 3
    assert result.n_output_loci == 1
    beta = read_matrix(result.matrix_path)
    assert beta.index.tolist() == ["cg1"]
    meta = pd.read_csv(result.metadata_path)
    assert set(meta["condition"]) == {"baseline", "reprogrammed"}



def test_paired_beta_preprocessing_handles_read_only_numpy_views(
    tmp_path: Path, monkeypatch
) -> None:
    """Regression test for pandas 3.x read-only arrays returned by to_numpy()."""

    original_to_numpy = pd.DataFrame.to_numpy

    def guarded_to_numpy(self, *args, **kwargs):
        values = original_to_numpy(self, *args, **kwargs)
        if not kwargs.get("copy", False):
            values.setflags(write=False)
        return values

    monkeypatch.setattr(pd.DataFrame, "to_numpy", guarded_to_numpy)

    raw = tmp_path / "mptr_readonly.csv.gz"
    _write_gzip_text(
        raw,
        "ID_REF,O1_negative_control_15days_exp1,Detection Pval,O1_transiently_reprogrammed_15days_exp1,Detection Pval\n"
        "cg1,0.2,0.0,0.4,0.0\n"
        "cg2,1.2,0.0,0.4,0.0\n",
    )
    result = preprocess_paired_beta_matrix(
        raw,
        tmp_path / "beta.csv.gz",
        tmp_path / "meta.csv",
        tmp_path / "sample_qc.csv",
        tmp_path / "summary.json",
        source_series="GSE165179",
        max_locus_missing=0.50,
        chunksize=1,
    )

    beta = read_matrix(result.matrix_path)
    assert beta.index.tolist() == ["cg1", "cg2"]
    assert np.isnan(beta.loc["cg2", "O1_negative_control_15days_exp1"])

def test_bystander_preprocessing_and_overlap(tmp_path: Path) -> None:
    t12 = tmp_path / "GSE269758_T_cell_hyper_beta_detp.csv.gz"
    t3 = tmp_path / "GSE269758_T_cell_R3_Hyperlibrary_beta_detp.csv.gz"
    msc = tmp_path / "GSE285251_MSC_Hyperlibrary_beta_detp.csv.gz"
    _write_gzip_text(
        t12,
        '"scrambleCTRL_D21_R1.beta","scrambleCTRL_D21_R1.detp","DNMT3A_HyperCpGs_D21_R1.beta","DNMT3A_HyperCpGs_D21_R1.detp","probeID"\n'
        '0.2,0.0,0.4,0.0,"cg1"\n0.3,0.0,0.5,0.0,"cg2"\n',
    )
    _write_gzip_text(
        t3,
        '"scrambleCTRL_D21_R3.beta","scrambleCTRL_D21_R3.detp","DNMT3A_HyperCpGs_D21_R3.beta","DNMT3A_HyperCpGs_D21_R3.detp","probeID"\n'
        '0.21,0.0,0.41,0.0,"cg1"\n0.31,0.0,0.51,0.0,"cg2"\n',
    )
    _write_gzip_text(
        msc,
        '"ID_REF","MSC170_scrDNMT.beta","MSC170_scrDNMT.detP","MSC170_HyperLDNMT.beta","MSC170_HyperLDNMT.detP","MSC179_scrDNMT.beta","MSC179_scrDNMT.detP","MSC179_HyperLDNMT.beta","MSC179_HyperLDNMT.detP"\n'
        '"cg1_BC21",0.2,0,0.4,0,0.3,0,0.5,0\n'
        '"cg2_TC21",0.3,0,0.5,0,0.4,0,0.6,0\n',
    )
    summary = preprocess_bystander_matrices([t12, t3], msc, tmp_path / "processed")
    assert summary["tcell"]["n_samples"] == 4
    assert summary["msc"]["n_samples"] == 4
    tmeta = pd.read_csv(tmp_path / "processed" / "tcell_metadata.csv")
    assert set(tmeta["treatment"]) == {"control", "targeted"}

    expr_meta = pd.DataFrame(
        {
            "sample_id": ["O1_negative_control_15days_exp1", "O1_transiently_reprogrammed_15days_exp1", "O2_negative_control_15days_exp1", "O2_transiently_reprogrammed_15days_exp1"],
            "donor_id": ["O1", "O1", "O2", "O2"],
            "condition": ["baseline", "reprogrammed", "baseline", "reprogrammed"],
        }
    )
    meth_meta = expr_meta.copy()
    expr_path = tmp_path / "expr_meta.csv"
    meth_path = tmp_path / "meth_meta.csv"
    expr_meta.to_csv(expr_path, index=False)
    meth_meta.to_csv(meth_path, index=False)
    report = build_overlap_report(expr_path, meth_path, tmp_path / "overlap.json")
    assert report["n_eligible_paired_donors"] == 2
    assert json.loads((tmp_path / "overlap.json").read_text())["n_overlapping_samples"] == 4
