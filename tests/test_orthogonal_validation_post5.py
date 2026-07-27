from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np

from ccid.orthogonal_validation import (
    _multidonor_lifespan_statistics,
    _parse_multidonor_lifespan_title,
    read_gse179847_processed_matrix,
)


def _write_metadata(path: Path) -> None:
    titles = [
        "Methylation_HC1_sp1_Normal_Control_ox21_p1",
        "Methylation_HC1_sp1_Normal_Control_ox21_p3",
        "Methylation_HC1_sp1_Normal_Control_ox21_p5",
        "Methylation_HC2_sp2_Normal_Control_ox21_p2",
        "Methylation_HC2_sp2_Normal_Control_ox21_p4",
        "Methylation_HC2_sp2_Normal_Control_ox21_p6",
    ]
    accessions = [f"GSM{i}" for i in range(1, 7)]
    descriptions = [f"ARR{i}" for i in range(1, 7)]
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("!Sample_title\t" + "\t".join(f'\"{x}\"' for x in titles) + "\n")
        handle.write("!Sample_geo_accession\t" + "\t".join(f'\"{x}\"' for x in accessions) + "\n")
        handle.write("!Sample_description\t" + "\t".join(f'\"{x}\"' for x in descriptions) + "\n")
        handle.write('!series_matrix_table_begin\n')
        handle.write('"ID_REF"\t' + "\t".join(f'\"{x}\"' for x in accessions) + "\n")
        handle.write('!series_matrix_table_end\n')


def _write_processed(path: Path) -> None:
    columns = []
    for i in range(1, 7):
        columns.extend([f"ARR{i} beta", f"ARR{i} Detection Pval"])
    beta = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60]
    values = []
    for value in beta:
        values.extend([value, 0.001])
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(','.join(['""'] + [f'"{x}"' for x in columns]) + "\n")
        handle.write(','.join(['"cg00000001"'] + [str(x) for x in values]) + "\n")


def test_post5_reader_maps_description_and_drops_detection_p(tmp_path: Path) -> None:
    metadata = tmp_path / "meta.txt.gz"
    processed = tmp_path / "processed.csv.gz"
    _write_metadata(metadata)
    _write_processed(processed)
    matrix = read_gse179847_processed_matrix(
        processed,
        target_ids=["cg00000001"],
        metadata_matrix_file=metadata,
    )
    assert matrix.values.shape == (1, 6)
    assert matrix.values.columns.tolist() == [f"ARR{i}" for i in range(1, 7)]
    assert matrix.metadata["sample_title"].str.startswith("Methylation_HC").all()
    assert np.isclose(matrix.values.loc["cg00000001", "ARR6"], 0.60)


def test_post5_parser_and_series_aware_statistics(tmp_path: Path) -> None:
    parsed = _parse_multidonor_lifespan_title(
        "Methylation_HC5_sp1_Normal_Control_ox21_p25"
    )
    assert parsed == ("HC5", "sp1", 21.0, 25.0, True)

    metadata = tmp_path / "meta.txt.gz"
    processed = tmp_path / "processed.csv.gz"
    _write_metadata(metadata)
    _write_processed(processed)
    matrix = read_gse179847_processed_matrix(
        processed,
        target_ids=["cg00000001"],
        metadata_matrix_file=metadata,
    )
    result = _multidonor_lifespan_statistics(matrix, ["cg00000001"]).iloc[0]
    assert bool(result["multidonor_available"])
    assert int(result["multidonor_n_donors"]) == 2
    assert int(result["multidonor_n_series"]) == 2
    assert result["multidonor_median_theil_sen_slope"] > 0
