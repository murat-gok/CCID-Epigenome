#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import pandas as pd

from ccid.bystander_validation import load_probe_annotation
from ccid.bystander_validation_v071 import (
    build_candidate_metadata,
    run_bystander_validation_v071,
    write_bystander_validation_v071_outputs,
)
from ccid.dataio import read_matrix


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Positive-control, platform-aware T-cell, unique-region, coverage-bias and "
            "multiplicity-corrected validation of dCas9-DNMT3A bystander responsiveness"
        )
    )
    parser.add_argument(
        "--tcell-raw-inputs",
        default=(
            "data/raw/GSE269758_T_cell_hyper_beta_detp.csv.gz,"
            "data/raw/GSE269758_T_cell_R3_Hyperlibrary_beta_detp.csv.gz"
        ),
    )
    parser.add_argument("--msc-beta", default="data/processed/msc_beta.csv.gz")
    parser.add_argument("--msc-metadata", default="data/processed/msc_metadata.csv")
    parser.add_argument("--stable-core", default="results/real/o3_robustness/lodo_consensus_3of3_ids.csv")
    parser.add_argument("--consensus-2of3", default="results/real/o3_robustness/lodo_consensus_at_least_2of3.csv")
    parser.add_argument("--tier1", default="results/real/o3_validation/o3_tier1_candidates.csv.gz")
    parser.add_argument("--frozen-signature", default="results/real/safe_osk_target_signature_O1O2.csv.gz")
    parser.add_argument("--discovery-all-loci", default="results/real/safe_osk_target_signature_O1O2_all_loci.csv.gz")
    parser.add_argument("--probe-annotation", required=True)
    parser.add_argument("--output-dir", default="results/real/bystander_validation_v071")
    parser.add_argument("--min-abs-delta", type=float, default=0.02)
    parser.add_argument("--fdr", type=float, default=0.10)
    parser.add_argument("--null-iterations", type=int, default=1000)
    parser.add_argument("--null-bins", type=int, default=10)
    parser.add_argument("--cluster-gap-bp", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260721)
    args = parser.parse_args()

    candidate_sets, target_map, candidate_metadata = build_candidate_metadata(
        stable_core_path=args.stable_core,
        consensus_2of3_path=args.consensus_2of3,
        tier1_path=args.tier1,
        frozen_signature_path=args.frozen_signature,
    )
    discovery_all = pd.read_csv(args.discovery_all_loci, compression="infer")
    annotation = load_probe_annotation(args.probe_annotation)
    result = run_bystander_validation_v071(
        tcell_raw_paths=[x for x in args.tcell_raw_inputs.split(",") if x],
        msc_beta=read_matrix(args.msc_beta),
        msc_metadata=pd.read_csv(args.msc_metadata),
        candidate_sets=candidate_sets,
        target_map=target_map,
        candidate_metadata=candidate_metadata,
        discovery_all_loci=discovery_all,
        probe_annotation=annotation,
        min_abs_delta=args.min_abs_delta,
        fdr=args.fdr,
        null_iterations=args.null_iterations,
        null_bins=args.null_bins,
        cluster_gap_bp=args.cluster_gap_bp,
        seed=args.seed,
    )
    write_bystander_validation_v071_outputs(result, args.output_dir)
    print("T-cell native platforms:", result.summary["n_tcell_platforms"])
    print("T-cell pairs:", result.summary["n_tcell_pairs"])
    print("MSC pairs:", result.summary["n_msc_pairs"])
    print("Unique regions:", result.summary["n_unique_regions_all"])
    print("Unique multi-CpG regions:", result.summary["n_unique_multicpg_regions"])
    print("Results:", args.output_dir)
    print(json.dumps(result.summary["significant_matched_null_results_q_le_0_10"], indent=2))


if __name__ == "__main__":
    main()
