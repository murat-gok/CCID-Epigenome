from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from ccid.manuscript_validation import (
    DEFAULT_V082_CONFIG,
    run_manuscript_validation,
    write_manuscript_validation_outputs,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Validate manuscript-readiness of candidate rankings with weight sensitivity, component uncertainty, "
            "leave-one-component-out ablation, permutation nulls, redundancy diagnostics and set enrichment."
        )
    )
    parser.add_argument("--v081-dir", default="results/real/regulatory_safety_v081")
    parser.add_argument("--candidate-file", default="")
    parser.add_argument("--config-json", default="configs/manuscript_validation_v082.json")
    parser.add_argument("--top-k", type=int, default=None)
    parser.add_argument("--iterations", type=int, default=None)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--output-dir", default="results/real/manuscript_validation_v082")
    args = parser.parse_args()

    if args.candidate_file:
        candidate_path = Path(args.candidate_file)
    else:
        filename = "candidate_cpg_regulatory_safety_v081.csv.gz"
        reviewed = Path("results/real/regulatory_safety_v081_reviewed") / filename
        base = Path(args.v081_dir) / filename
        candidate_path = reviewed if reviewed.exists() else base
    if not candidate_path.exists():
        raise FileNotFoundError(
            f"Candidate input not found: {candidate_path}. Run v0.8.1 first or pass --candidate-file."
        )

    config = dict(DEFAULT_V082_CONFIG)
    config_path = Path(args.config_json)
    if config_path.exists():
        config.update(json.loads(config_path.read_text(encoding="utf-8")))
    if args.top_k is not None:
        config["top_k"] = args.top_k
    if args.iterations is not None:
        config["weight_iterations"] = args.iterations
        config["uncertainty_iterations"] = args.iterations
        config["null_iterations"] = args.iterations
    if args.seed is not None:
        config["seed"] = args.seed

    cpg = pd.read_csv(candidate_path, low_memory=False)
    result = run_manuscript_validation(cpg, config=config)
    result.manifest["candidate_input"] = str(candidate_path)
    result.manifest["candidate_input_size_bytes"] = candidate_path.stat().st_size
    write_manuscript_validation_outputs(result, args.output_dir)

    print("Candidate CpGs:", result.summary["n_candidates"])
    print("Components:", len(result.summary["components"]))
    print("Stable primary shortlist:", result.summary["n_stable_primary_shortlist"])
    print("Median weight Spearman:", result.summary["median_weight_spearman"])
    print("Median uncertainty Spearman:", result.summary["median_uncertainty_spearman"])
    print("Results:", args.output_dir)
    print(json.dumps(result.summary["top_robust_candidates"][:10], indent=2))


if __name__ == "__main__":
    main()
