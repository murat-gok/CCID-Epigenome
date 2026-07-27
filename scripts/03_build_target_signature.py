#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import pandas as pd

from ccid.signature import build_matched_methylation_target_signature_file


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build an exact donor-day-batch-stage matched safe-OSK methylation target signature"
        )
    )
    parser.add_argument("--beta", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--safe-window", required=True)
    parser.add_argument("--output", required=True, help="Selected/ranked target signature CSV")
    parser.add_argument("--all-loci-output")
    parser.add_argument("--pairs-output")
    parser.add_argument("--donor-effects-output")
    parser.add_argument("--qc-output")
    parser.add_argument("--baseline-condition", default="baseline")
    parser.add_argument("--discovery-donors", default="")
    parser.add_argument("--match-columns", default="donor_id,day,batch,stage")
    parser.add_argument("--min-abs-delta", type=float, default=0.03)
    parser.add_argument("--fdr", type=float, default=0.10)
    parser.add_argument("--min-sign-consistency", type=float, default=0.67)
    parser.add_argument("--min-matched-donors", type=int, default=2)
    parser.add_argument(
        "--min-inferential-donors",
        type=int,
        default=3,
        help=(
            "Minimum independent donors required for confirmatory t-test/FDR selection. "
            "With fewer donors, the script writes a ranked exploratory fallback."
        ),
    )
    parser.add_argument("--fallback-top-n", type=int, default=500)
    parser.add_argument("--chunksize", type=int, default=5000)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata)
    safe = pd.read_csv(args.safe_window)
    if "selected_safe_window" not in safe.columns:
        raise ValueError("Safe-window file lacks selected_safe_window column")
    safe_samples = safe.loc[safe["selected_safe_window"].astype(bool), "sample_id"].astype(str).tolist()
    if not safe_samples:
        raise ValueError("Safe-window file contains no selected samples")

    donors = [x.strip() for x in args.discovery_donors.split(",") if x.strip()] or None
    match_columns = [x.strip() for x in args.match_columns.split(",") if x.strip()]

    result = build_matched_methylation_target_signature_file(
        beta_path=args.beta,
        metadata=metadata,
        safe_samples=safe_samples,
        baseline_condition=args.baseline_condition,
        match_columns=match_columns,
        discovery_donors=donors,
        min_abs_delta=args.min_abs_delta,
        fdr=args.fdr,
        min_sign_consistency=args.min_sign_consistency,
        min_matched_donors=args.min_matched_donors,
        min_inferential_donors=args.min_inferential_donors,
        fallback_top_n=args.fallback_top_n,
        chunksize=args.chunksize,
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    stem = output.name[:-7] if output.name.endswith(".csv.gz") else output.stem
    all_loci_output = Path(args.all_loci_output) if args.all_loci_output else output.parent / f"{stem}_all_loci.csv.gz"
    pairs_output = Path(args.pairs_output) if args.pairs_output else output.parent / f"{stem}_matched_pairs.csv"
    donor_effects_output = (
        Path(args.donor_effects_output)
        if args.donor_effects_output
        else output.parent / f"{stem}_donor_effects.csv.gz"
    )
    qc_output = Path(args.qc_output) if args.qc_output else output.parent / f"{stem}_qc.json"

    result.signature.to_csv(output, index=False)
    result.all_loci.to_csv(all_loci_output, index=False, compression="infer")
    result.pair_table.to_csv(pairs_output, index=False)
    result.donor_effects.to_csv(donor_effects_output, index=True, index_label="cpg_id", compression="infer")
    qc_output.write_text(json.dumps(result.qc, indent=2), encoding="utf-8")

    print(f"Matched safe samples: {result.qc['n_matched_safe_samples']}")
    print(f"Matched donors: {result.qc['matched_donors']}")
    print(f"Strictly selected loci: {result.qc['n_strictly_selected_loci']}")
    print(f"Inferential status: {result.qc['inferential_status']}")
    print(f"Used exploratory fallback: {result.qc['used_underpowered_fallback']}")
    print(f"Signature: {output}")
    print(f"Pair table: {pairs_output}")
    print(f"QC: {qc_output}")


if __name__ == "__main__":
    main()
