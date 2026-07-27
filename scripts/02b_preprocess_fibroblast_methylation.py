#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from ccid.preprocessing import preprocess_paired_beta_matrix


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stream-preprocess paired beta/detection-p fibroblast methylation matrices"
    )
    parser.add_argument(
        "--input",
        default="data/raw/GSE165179_Matrix_processed_transient.txt.gz",
    )
    parser.add_argument(
        "--series",
        choices=["mptr", "sendai"],
        default="mptr",
    )
    parser.add_argument("--output-dir", default="data/processed")
    parser.add_argument("--detection-p-threshold", type=float, default=0.01)
    parser.add_argument("--max-locus-missing", type=float, default=0.20)
    parser.add_argument("--chunksize", type=int, default=5000)
    args = parser.parse_args()

    output = Path(args.output_dir)
    if args.series == "mptr":
        prefix = "fibroblast"
        source_series = "GSE165179"
    else:
        prefix = "sendai_fibroblast"
        source_series = "GSE165178"

    result = preprocess_paired_beta_matrix(
        args.input,
        output / f"{prefix}_beta.csv.gz",
        output / f"{prefix}_methylation_metadata.csv",
        output / "qc" / f"{prefix}_methylation_sample_qc.csv",
        output / "qc" / f"{prefix}_methylation_summary.json",
        source_series=source_series,
        cell_type="dermal_fibroblast",
        detection_p_threshold=args.detection_p_threshold,
        max_locus_missing=args.max_locus_missing,
        chunksize=args.chunksize,
    )
    print(
        json.dumps(
            {
                "matrix": str(result.matrix_path),
                "metadata": str(result.metadata_path),
                "n_input_loci": result.n_input_loci,
                "n_output_loci": result.n_output_loci,
                "n_samples": result.n_samples,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
