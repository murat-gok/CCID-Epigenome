#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import pandas as pd

from ccid.bystander_validation import (
    load_candidate_sets,
    load_probe_annotation,
    run_bystander_validation,
    write_bystander_validation_outputs,
)
from ccid.dataio import read_matrix


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate frozen safe-reprogramming CpGs in T-cell and MSC dCas9-DNMT3A bystander datasets"
    )
    parser.add_argument("--tcell-beta", default="data/processed/tcell_beta.csv.gz")
    parser.add_argument("--tcell-metadata", default="data/processed/tcell_metadata.csv")
    parser.add_argument("--msc-beta", default="data/processed/msc_beta.csv.gz")
    parser.add_argument("--msc-metadata", default="data/processed/msc_metadata.csv")
    parser.add_argument("--stable-core", default="results/real/o3_robustness/lodo_consensus_3of3_ids.csv")
    parser.add_argument("--consensus-2of3", default="results/real/o3_robustness/lodo_consensus_at_least_2of3.csv")
    parser.add_argument("--tier1", default="results/real/o3_validation/o3_tier1_candidates.csv.gz")
    parser.add_argument("--frozen-signature", default="results/real/safe_osk_target_signature_O1O2.csv.gz")
    parser.add_argument("--discovery-all-loci", default="results/real/safe_osk_target_signature_O1O2_all_loci.csv.gz")
    parser.add_argument("--probe-annotation", help="Optional EPIC/EPICv2 annotation CSV(.gz) with cpg_id, chromosome, position")
    parser.add_argument("--output-dir", default="results/real/bystander_validation")
    parser.add_argument("--min-abs-delta", type=float, default=0.02)
    parser.add_argument("--fdr", type=float, default=0.10)
    parser.add_argument("--null-iterations", type=int, default=1000)
    parser.add_argument("--null-bins", type=int, default=10)
    parser.add_argument("--cluster-gap-bp", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260721)
    args = parser.parse_args()

    candidate_sets, target_map = load_candidate_sets(
        stable_core_path=args.stable_core,
        consensus_2of3_path=args.consensus_2of3,
        tier1_path=args.tier1,
        frozen_signature_path=args.frozen_signature,
    )
    discovery_all = pd.read_csv(args.discovery_all_loci, compression="infer") if args.discovery_all_loci else None
    annotation = load_probe_annotation(args.probe_annotation) if args.probe_annotation else None

    result = run_bystander_validation(
        tcell_beta=read_matrix(args.tcell_beta),
        tcell_metadata=pd.read_csv(args.tcell_metadata),
        msc_beta=read_matrix(args.msc_beta),
        msc_metadata=pd.read_csv(args.msc_metadata),
        candidate_sets=candidate_sets,
        target_map=target_map,
        discovery_all_loci=discovery_all,
        probe_annotation=annotation,
        min_abs_delta=args.min_abs_delta,
        fdr=args.fdr,
        null_iterations=args.null_iterations,
        null_bins=args.null_bins,
        seed=args.seed,
        cluster_gap_bp=args.cluster_gap_bp,
    )
    write_bystander_validation_outputs(result, args.output_dir)

    print("T-cell pairs:", result.summary["n_tcell_pairs"])
    print("MSC pairs:", result.summary["n_msc_pairs"])
    print("Probe annotation supplied:", result.summary["probe_annotation_supplied"])
    print("Regional scores:", result.summary["n_regional_scores"])
    print("Results:", args.output_dir)
    print(json.dumps(result.summary["cross_celltype_summary"], indent=2))


if __name__ == "__main__":
    main()
