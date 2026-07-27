#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import numpy as np
import pandas as pd

from ccid.objectives import no_intervention_summary, relative_target_improvement
from ccid.optimizers import greedy_search, nsga2, random_search, select_target_priority
from ccid.reporting import records_to_frame
from ccid.schema import OptimizationContext, validate_candidate_table


def main() -> None:
    parser = argparse.ArgumentParser(description="Optimize minimal safe epigenome-editing protocols")
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--effect-mean", required=True)
    parser.add_argument("--effect-std", required=True)
    parser.add_argument("--target-signature", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--max-edits", type=int, default=4)
    parser.add_argument("--min-edits", type=int, default=1)
    parser.add_argument("--min-relative-improvement", type=float, default=0.10)
    parser.add_argument("--population", type=int, default=160)
    parser.add_argument("--generations", type=int, default=250)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    candidates = pd.read_csv(args.candidates)
    validate_candidate_table(candidates)
    mean = np.load(args.effect_mean)
    std = np.load(args.effect_std)
    target = pd.read_csv(args.target_signature)
    required_target = {"target_delta_beta", "weight"}
    if not required_target.issubset(target.columns):
        raise ValueError(f"Target signature requires columns {sorted(required_target)}")
    target_effect = target["target_delta_beta"].to_numpy(dtype=float)
    weights = target["weight"].to_numpy(dtype=float)
    if mean.shape != std.shape:
        raise ValueError("Effect mean and std matrices must have identical shapes")
    if mean.shape[0] != len(candidates) or mean.shape[1] != len(target):
        raise ValueError("Candidate/effect/target dimensions do not match")

    context = OptimizationContext(
        predicted_effects=mean,
        predicted_std=std,
        target_effect=target_effect,
        locus_weights=weights,
        candidate_table=candidates,
        max_edits=args.max_edits,
        min_edits=args.min_edits,
        min_relative_improvement=args.min_relative_improvement,
    )
    fronts = {
        "nsga2": nsga2(
            context,
            population_size=args.population,
            generations=args.generations,
            seed=args.seed,
        ),
        "random": random_search(context, n_samples=10000, seed=args.seed + 1),
        "greedy": greedy_search(context, seed=args.seed + 2),
    }

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    for method, front in fronts.items():
        records_to_frame(front).to_csv(output / f"{method}_pareto_front.csv", index=False)

    selected = select_target_priority(fronts["nsga2"])
    protocol = candidates.iloc[list(selected.selected_indices)].copy()
    protocol.insert(0, "candidate_index", list(selected.selected_indices))
    protocol["target_mismatch"] = selected.objectives[0]
    protocol["relative_target_improvement"] = relative_target_improvement(
        selected.objectives[0], context
    )
    protocol.to_csv(output / "selected_protocol.csv", index=False)
    pd.DataFrame([no_intervention_summary(context)]).to_csv(
        output / "no_intervention_baseline.csv", index=False
    )
    print(protocol[["candidate_index", "candidate_id", "region_id", "editor"]].to_string(index=False))
    print("No-intervention comparator written separately; it is not a protocol candidate.")


if __name__ == "__main__":
    main()
