#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from ccid.preprocessing import build_overlap_report


def main() -> None:
    parser = argparse.ArgumentParser(description="Check expression/methylation sample overlap")
    parser.add_argument(
        "--expression-metadata",
        default="data/processed/fibroblast_expression_metadata.csv",
    )
    parser.add_argument(
        "--methylation-metadata",
        default="data/processed/fibroblast_methylation_metadata.csv",
    )
    parser.add_argument(
        "--output",
        default="data/processed/qc/fibroblast_sample_overlap.json",
    )
    args = parser.parse_args()
    report = build_overlap_report(
        args.expression_metadata,
        args.methylation_metadata,
        args.output,
    )
    print(json.dumps(report, indent=2))
    if report["n_eligible_paired_donors"] < 2:
        raise SystemExit(
            "ERROR: fewer than two donors have overlapping baseline and reprogrammed samples"
        )


if __name__ == "__main__":
    main()
