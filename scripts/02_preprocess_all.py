#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

REQUIRED = [
    "GSE165177_Log2_RPM_Transient_reprogramming.txt.gz",
    "GSE165177_Log2_RPM_Transient_reprogramming_part2_170621.txt.gz",
    "GSE165179_Matrix_processed_transient.txt.gz",
    "GSE269758_T_cell_hyper_beta_detp.csv.gz",
    "GSE269758_T_cell_R3_Hyperlibrary_beta_detp.csv.gz",
    "GSE285251_MSC_Hyperlibrary_beta_detp.csv.gz",
]


def run(script: str, *args: str) -> None:
    command = [sys.executable, str(PROJECT_ROOT / "scripts" / script), *args]
    print("\nRUN", " ".join(command), flush=True)
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the complete public-data preprocessing stage")
    parser.add_argument("--raw-dir", default="data/raw")
    parser.add_argument("--processed-dir", default="data/processed")
    parser.add_argument("--chunksize", type=int, default=5000)
    parser.add_argument("--detection-p-threshold", type=float, default=0.01)
    parser.add_argument("--max-locus-missing", type=float, default=0.20)
    args = parser.parse_args()

    raw = PROJECT_ROOT / args.raw_dir
    missing = [name for name in REQUIRED if not (raw / name).exists()]
    if missing:
        print(json.dumps({"missing_required_files": missing}, indent=2))
        raise SystemExit(
            "Download the missing files with:\n"
            "python scripts/01_download_public_data.py "
            "--datasets gse165177_rna_part2,gse269758_tcell_bystander_r3 "
            "--output data/raw"
        )

    processed = args.processed_dir
    run(
        "02a_preprocess_fibroblast_expression.py",
        "--inputs",
        ",".join(
            [
                str(raw / REQUIRED[0]),
                str(raw / REQUIRED[1]),
            ]
        ),
        "--output-expression",
        f"{processed}/fibroblast_expression.csv.gz",
        "--output-metadata",
        f"{processed}/fibroblast_expression_metadata.csv",
        "--output-qc",
        f"{processed}/qc/fibroblast_expression_qc.json",
    )
    run(
        "02b_preprocess_fibroblast_methylation.py",
        "--input",
        str(raw / REQUIRED[2]),
        "--series",
        "mptr",
        "--output-dir",
        processed,
        "--chunksize",
        str(args.chunksize),
        "--detection-p-threshold",
        str(args.detection_p_threshold),
        "--max-locus-missing",
        str(args.max_locus_missing),
    )
    run(
        "02c_preprocess_bystander_methylation.py",
        "--tcell-inputs",
        ",".join([str(raw / REQUIRED[3]), str(raw / REQUIRED[4])]),
        "--msc-input",
        str(raw / REQUIRED[5]),
        "--output-dir",
        processed,
        "--max-locus-missing",
        str(args.max_locus_missing),
    )
    run(
        "02d_build_sample_overlap_report.py",
        "--expression-metadata",
        f"{processed}/fibroblast_expression_metadata.csv",
        "--methylation-metadata",
        f"{processed}/fibroblast_methylation_metadata.csv",
        "--output",
        f"{processed}/qc/fibroblast_sample_overlap.json",
    )
    print("\nPreprocessing completed successfully.")


if __name__ == "__main__":
    main()
