#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import pandas as pd

from ccid.bystander import analyze_bystander_effects
from ccid.dataio import read_matrix


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze reproducible bystander methylation effects")
    parser.add_argument("--beta", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--age-loci", help="One-column CSV/TXT of age-associated CpG IDs")
    parser.add_argument("--targeted-label", default="targeted")
    parser.add_argument("--control-label", default="control")
    parser.add_argument("--permutations", type=int, default=1000)
    args = parser.parse_args()

    beta = read_matrix(args.beta)
    metadata = pd.read_csv(args.metadata)
    age_loci: set[str] | None = None
    if args.age_loci:
        age_frame = pd.read_csv(args.age_loci, header=None)
        age_loci = set(age_frame.iloc[:, 0].astype(str))

    analysis = analyze_bystander_effects(
        beta,
        metadata,
        targeted_label=args.targeted_label,
        control_label=args.control_label,
        age_loci=age_loci,
        n_permutations=args.permutations,
    )
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    analysis.locus_results.to_csv(output / "bystander_locus_results.csv", index=False)
    analysis.replicate_correlations.to_csv(output / "bystander_replicate_correlations.csv", index=False)
    (output / "bystander_summary.json").write_text(json.dumps(analysis.summary, indent=2))
    print(json.dumps(analysis.summary, indent=2))


if __name__ == "__main__":
    main()
