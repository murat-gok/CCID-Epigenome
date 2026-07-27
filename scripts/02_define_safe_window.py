#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import pandas as pd

from ccid.dataio import read_matrix
from ccid.signature import define_safe_window

DEFAULT_IDENTITY = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "THY1"]
DEFAULT_PLURIPOTENCY = ["POU5F1", "SOX2", "NANOG", "LIN28A", "DPPA4", "ZFP42", "SALL4"]
DEFAULT_AGEING = ["CDKN2A", "CDKN1A", "SERPINE1", "IL6", "CXCL8", "MMP3"]


def main() -> None:
    parser = argparse.ArgumentParser(description="Define a leakage-free safe reprogramming window")
    parser.add_argument("--expression", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary-output")
    parser.add_argument("--discovery-donors", default="")
    parser.add_argument("--gene-sets-json")
    parser.add_argument("--identity-floor", type=float, default=-0.50)
    parser.add_argument("--pluripotency-ceiling", type=float, default=0.50)
    parser.add_argument("--min-rejuvenation", type=float, default=0.0)
    parser.add_argument("--min-donors-per-day", type=int, default=2)
    parser.add_argument("--utility-tolerance", type=float, default=0.10)
    parser.add_argument(
        "--allowed-conditions",
        default="reprogrammed,reprogrammed_intermediate",
        help="Comma-separated exact condition labels eligible for safe-window selection",
    )
    args = parser.parse_args()

    expression = read_matrix(args.expression)
    metadata = pd.read_csv(args.metadata)
    gene_sets = {
        "identity": DEFAULT_IDENTITY,
        "pluripotency": DEFAULT_PLURIPOTENCY,
        "ageing": DEFAULT_AGEING,
    }
    if args.gene_sets_json:
        gene_sets.update(json.loads(Path(args.gene_sets_json).read_text()))
    donors = [x.strip() for x in args.discovery_donors.split(",") if x.strip()] or None
    allowed_conditions = [x.strip() for x in args.allowed_conditions.split(",") if x.strip()]

    result = define_safe_window(
        expression=expression,
        metadata=metadata,
        identity_genes=gene_sets["identity"],
        pluripotency_genes=gene_sets["pluripotency"],
        rejuvenation_genes=gene_sets.get("ageing"),
        identity_floor=args.identity_floor,
        pluripotency_ceiling=args.pluripotency_ceiling,
        min_rejuvenation=args.min_rejuvenation,
        discovery_donors=donors,
        allowed_conditions=allowed_conditions,
        min_donors_per_day=args.min_donors_per_day,
        utility_tolerance=args.utility_tolerance,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.sample_scores.to_csv(output, index=False)

    selected = result.sample_scores[result.sample_scores["selected_safe_window"]].copy()
    summary = {
        "discovery_donors": donors,
        "allowed_conditions": allowed_conditions,
        "identity_floor": args.identity_floor,
        "pluripotency_ceiling": args.pluripotency_ceiling,
        "min_rejuvenation": args.min_rejuvenation,
        "min_donors_per_day": args.min_donors_per_day,
        "utility_tolerance": args.utility_tolerance,
        "selected_days": list(result.selected_days),
        "n_selected_samples": int(len(selected)),
        "selected_donors": sorted(selected["donor_id"].astype(str).unique().tolist()),
        "selected_conditions": sorted(selected["condition"].astype(str).unique().tolist()),
        "selected_samples": result.selected_samples,
    }
    summary_path = Path(args.summary_output) if args.summary_output else output.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Selected days: {result.selected_days}")
    print(f"Selected samples: {len(result.selected_samples)}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
