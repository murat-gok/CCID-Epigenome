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
from ccid.genomic_causality import (
    DEFAULT_WEIGHTS,
    load_gene_evidence,
    read_gmt,
    run_genomic_causality_scoring,
    write_genomic_causality_outputs,
)


def _optional_csv(path: str | None) -> pd.DataFrame | None:
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        print(f"Optional evidence file not found; continuing without it: {p}")
        return None
    return pd.read_csv(p, compression="infer", low_memory=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Integrate genomic context, biological-relevance panels, O3/temporal evidence and "
            "cell-context editor responsiveness into transparent CpG/region prioritisation scores"
        )
    )
    parser.add_argument("--v071-dir", default="results/real/bystander_validation_v071")
    parser.add_argument("--probe-annotation", default="data/reference/epicv2_probe_annotation_hg38.csv.gz")
    parser.add_argument("--gene-evidence-seed", default="data/reference/biological_evidence_seed_v080.csv")
    parser.add_argument("--external-gene-evidence", default="")
    parser.add_argument("--gmt", default="")
    parser.add_argument(
        "--successful-vs-failed",
        default="results/real/o3_robustness/o3_successful_vs_failed.csv.gz",
    )
    parser.add_argument(
        "--time-series-consistency",
        default="results/real/o3_robustness/time_series_cpg_consistency.csv.gz",
    )
    parser.add_argument("--weights-json", default="configs/causality_scoring_v080.json")
    parser.add_argument("--shortlist-size", type=int, default=100)
    parser.add_argument("--output-dir", default="results/real/genomic_causality_v080")
    args = parser.parse_args()

    v071 = Path(args.v071_dir)
    membership = pd.read_csv(v071 / "unique_region_membership.csv.gz", low_memory=False)
    regions = pd.read_csv(v071 / "unique_region_scores_all.csv", low_memory=False)
    annotation = load_probe_annotation(args.probe_annotation)
    # Preserve additional EPIC columns such as gene_group/strand that the generic loader leaves intact.
    full_annotation = pd.read_csv(args.probe_annotation, compression="infer", low_memory=False)
    if "cpg_id" not in full_annotation.columns:
        full_annotation = annotation
    else:
        annotation = annotation.merge(
            full_annotation[[c for c in full_annotation.columns if c not in annotation.columns or c == "cpg_id"]],
            on="cpg_id",
            how="left",
        )

    external = [x for x in args.external_gene_evidence.split(",") if x]
    gene_evidence = load_gene_evidence(args.gene_evidence_seed, external)
    gmt_paths = [x for x in args.gmt.split(",") if x]
    gmt = read_gmt(gmt_paths) if gmt_paths else None
    successful = _optional_csv(args.successful_vs_failed)
    time_series = _optional_csv(args.time_series_consistency)

    weights = dict(DEFAULT_WEIGHTS)
    weight_path = Path(args.weights_json)
    if weight_path.exists():
        payload = json.loads(weight_path.read_text(encoding="utf-8"))
        supplied = payload.get("weights", payload)
        weights = {str(k): float(v) for k, v in supplied.items()}
        total = sum(weights.values())
        if total <= 0:
            raise ValueError("Causality score weights must sum to a positive value")
        weights = {k: v / total for k, v in weights.items()}

    result = run_genomic_causality_scoring(
        membership=membership,
        region_scores=regions,
        annotation=annotation,
        gene_evidence=gene_evidence,
        successful_vs_failed=successful,
        time_series=time_series,
        gmt=gmt,
        weights=weights,
        shortlist_size=args.shortlist_size,
    )
    write_genomic_causality_outputs(result, args.output_dir)
    print("Candidate CpGs:", result.summary["n_candidate_cpgs"])
    print("Gene-annotated CpGs:", result.summary["n_annotated_to_gene"])
    print("Unique regions:", result.summary["n_unique_regions"])
    print("P1 candidates:", result.summary["n_priority_P1"])
    print("Safety-flagged CpGs:", result.summary["n_safety_flagged_cpgs"])
    print("Results:", args.output_dir)
    print(json.dumps(result.summary["top_cpgs"][:10], indent=2))


if __name__ == "__main__":
    main()
