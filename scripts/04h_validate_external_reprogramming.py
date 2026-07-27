#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ccid.orthogonal_validation import validate_external_reprogramming


def _write_json(payload: dict, path: Path) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate the frozen CCID CpG panel in independent GSE142439 skin fibroblasts."
    )
    parser.add_argument(
        "--candidate-file",
        default="results/real/regulatory_safety_v081/candidate_cpg_regulatory_safety_v081.csv.gz",
    )
    parser.add_argument(
        "--frozen-shortlist",
        default="data/reference/frozen_primary_shortlist_v090.csv",
    )
    parser.add_argument(
        "--series-matrix",
        default="data/raw/orthogonal_validation/GSE142439_series_matrix.txt.gz",
    )
    parser.add_argument("--bootstrap-iterations", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260724)
    parser.add_argument(
        "--output-dir", default="results/real/orthogonal_validation_v090/external_reprogramming"
    )
    args = parser.parse_args()

    for value in [args.candidate_file, args.frozen_shortlist, args.series_matrix]:
        if not Path(value).exists():
            raise FileNotFoundError(value)
    output = Path(args.output_dir)
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    result, summary = validate_external_reprogramming(
        args.candidate_file,
        args.frozen_shortlist,
        args.series_matrix,
        bootstrap_iterations=args.bootstrap_iterations,
        seed=args.seed,
    )
    result_path = output / "external_reprogramming_validation_v090.csv.gz"
    result.to_csv(result_path, index=False, compression="gzip")
    _write_json(summary, output / "external_reprogramming_summary_v090.json")

    plot = result[result["external_available"].fillna(False)].copy()
    plot = plot.sort_values("candidate_rank", ascending=False)
    if not plot.empty:
        y = np.arange(len(plot))
        means = pd.to_numeric(plot["external_mean_delta_beta"], errors="coerce").to_numpy()
        low = pd.to_numeric(plot["external_bootstrap_ci_low"], errors="coerce").to_numpy()
        high = pd.to_numeric(plot["external_bootstrap_ci_high"], errors="coerce").to_numpy()
        xerr = np.vstack([means - low, high - means])
        fig, ax = plt.subplots(figsize=(8.0, 7.0))
        ax.errorbar(means, y, xerr=xerr, fmt="o", capsize=3)
        ax.axvline(0, linewidth=1)
        ax.set_yticks(y, plot["cpg_id"].astype(str))
        ax.set_xlabel("Mean treated − normal DNA methylation (β)")
        ax.set_ylabel("Frozen candidate")
        ax.set_title("Independent transient-reprogramming effects in GSE142439")
        fig.tight_layout()
        for suffix in ["png", "pdf"]:
            fig.savefig(figures / f"fig_external_reprogramming_forest_v090.{suffix}", dpi=400)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6.5, 6.0))
        ax.scatter(plot["discovery_delta_beta"], plot["external_mean_delta_beta"])
        ax.axhline(0, linewidth=1)
        ax.axvline(0, linewidth=1)
        for row in plot.itertuples(index=False):
            ax.annotate(row.cpg_id, (row.discovery_delta_beta, row.external_mean_delta_beta), fontsize=7)
        ax.set_xlabel("Original MPTR discovery Δβ")
        ax.set_ylabel("Independent GSE142439 mean Δβ")
        ax.set_title("Direction and effect concordance without model refitting")
        fig.tight_layout()
        for suffix in ["png", "pdf"]:
            fig.savefig(figures / f"fig_external_effect_concordance_v090.{suffix}", dpi=400)
        plt.close(fig)

    print("Frozen candidates:", summary["n_frozen_candidates"])
    print("Candidates available:", summary["n_candidates_available"])
    print("Paired donors:", summary["n_paired_donors"])
    print("Direction concordant:", summary["direction_concordant_candidates"])
    print("Binomial p:", summary["direction_concordance_binomial_p"])
    print("Effect Spearman:", summary["effect_spearman"])
    print("Results:", output)


if __name__ == "__main__":
    main()
