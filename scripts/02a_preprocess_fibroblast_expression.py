#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from ccid.preprocessing import preprocess_fibroblast_expression


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge and preprocess GSE165177 MPTR fibroblast expression files"
    )
    parser.add_argument(
        "--inputs",
        default=(
            "data/raw/GSE165177_Log2_RPM_Transient_reprogramming.txt.gz,"
            "data/raw/GSE165177_Log2_RPM_Transient_reprogramming_part2_170621.txt.gz"
        ),
        help="Comma-separated GSE165177 processed expression files",
    )
    parser.add_argument(
        "--output-expression",
        default="data/processed/fibroblast_expression.csv.gz",
    )
    parser.add_argument(
        "--output-metadata",
        default="data/processed/fibroblast_expression_metadata.csv",
    )
    parser.add_argument(
        "--output-qc",
        default="data/processed/qc/fibroblast_expression_qc.json",
    )
    args = parser.parse_args()

    inputs = [x for x in args.inputs.split(",") if x]
    expression, metadata, qc = preprocess_fibroblast_expression(
        inputs,
        args.output_expression,
        args.output_metadata,
        args.output_qc,
    )
    print(
        json.dumps(
            {
                "n_genes": int(expression.shape[0]),
                "n_samples": int(expression.shape[1]),
                "conditions": metadata["condition"].value_counts().to_dict(),
                "qc": args.output_qc,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
