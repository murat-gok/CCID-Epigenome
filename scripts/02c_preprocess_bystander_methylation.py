#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from ccid.preprocessing import preprocess_bystander_matrices


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Preprocess T-cell and MSC epigenome-editing bystander matrices"
    )
    parser.add_argument(
        "--tcell-inputs",
        default=(
            "data/raw/GSE269758_T_cell_hyper_beta_detp.csv.gz,"
            "data/raw/GSE269758_T_cell_R3_Hyperlibrary_beta_detp.csv.gz"
        ),
    )
    parser.add_argument(
        "--msc-input",
        default="data/raw/GSE285251_MSC_Hyperlibrary_beta_detp.csv.gz",
    )
    parser.add_argument("--output-dir", default="data/processed")
    parser.add_argument("--max-locus-missing", type=float, default=0.20)
    args = parser.parse_args()

    summary = preprocess_bystander_matrices(
        [x for x in args.tcell_inputs.split(",") if x],
        args.msc_input,
        args.output_dir,
        max_locus_missing=args.max_locus_missing,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
