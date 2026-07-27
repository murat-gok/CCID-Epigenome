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

from ccid.orthogonal_validation import validate_natural_aging


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Test whether frozen partial-reprogramming CpGs oppose natural fibroblast aging trajectories."
    )
    parser.add_argument(
        "--candidate-file",
        default="results/real/regulatory_safety_v081/candidate_cpg_regulatory_safety_v081.csv.gz",
    )
    parser.add_argument(
        "--frozen-shortlist", default="data/reference/frozen_primary_shortlist_v090.csv"
    )
    parser.add_argument(
        "--replicative-matrix",
        default="data/raw/orthogonal_validation/GSE131280_series_matrix.txt.gz",
    )
    parser.add_argument(
        "--chronological-matrix",
        default="data/raw/orthogonal_validation/GSE22595_series_matrix.txt.gz",
    )
    parser.add_argument(
        "--multidonor-matrix",
        default="",
        help="Optional GSE179847 processed beta CSV for six-donor within-donor aging trajectories.",
    )
    parser.add_argument(
        "--multidonor-metadata",
        default="",
        help="Optional small GSE179847 Series Matrix used to map GSM accessions to sample titles.",
    )
    parser.add_argument(
        "--output-dir", default="results/real/orthogonal_validation_v090/natural_aging"
    )
    args = parser.parse_args()

    for value in [args.candidate_file, args.frozen_shortlist]:
        if not Path(value).exists():
            raise FileNotFoundError(value)
    rep = args.replicative_matrix if Path(args.replicative_matrix).exists() else None
    age = args.chronological_matrix if Path(args.chronological_matrix).exists() else None
    multi = args.multidonor_matrix if args.multidonor_matrix and Path(args.multidonor_matrix).exists() else None
    multi_meta = (
        args.multidonor_metadata
        if args.multidonor_metadata and Path(args.multidonor_metadata).exists()
        else None
    )
    if rep is None and age is None and multi is None:
        raise FileNotFoundError("No natural-aging input matrix is available")

    output = Path(args.output_dir)
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    result, summary = validate_natural_aging(
        args.candidate_file,
        args.frozen_shortlist,
        replicative_matrix_file=rep,
        chronological_matrix_file=age,
        multidonor_matrix_file=multi,
        multidonor_metadata_file=multi_meta,
    )
    result.to_csv(output / "natural_aging_direction_validation_v090.csv.gz", index=False, compression="gzip")
    (output / "natural_aging_summary_v090.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    support_cols = [
        c
        for c in [
            "replicative_reversal_concordant",
            "chronological_reversal_concordant",
            "multidonor_reversal_concordant",
        ]
        if c in result.columns
    ]
    if support_cols:
        plot = result.sort_values("candidate_rank", ascending=True).copy()
        matrix = plot[support_cols].fillna(False).astype(int).to_numpy()
        fig, ax = plt.subplots(figsize=(7.2, 7.0))
        image = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1)
        ax.set_yticks(np.arange(len(plot)), plot["cpg_id"])
        ax.set_xticks(
            np.arange(len(support_cols)),
            [c.replace("_reversal_concordant", "").replace("_", " ") for c in support_cols],
            rotation=30,
            ha="right",
        )
        ax.set_title("Opposition between partial-reprogramming and aging directions")
        fig.colorbar(image, ax=ax, label="Directional reversal support")
        fig.tight_layout()
        for suffix in ["png", "pdf"]:
            fig.savefig(figures / f"fig_natural_aging_direction_matrix_v090.{suffix}", dpi=400)
        plt.close(fig)

    slope_specs = [
        ("replicative_theil_sen_slope", "GSE131280 replicative-aging slope"),
        ("chronological_theil_sen_slope", "GSE22595 chronological-aging slope"),
        ("multidonor_median_theil_sen_slope", "GSE179847 median donor slope"),
    ]
    available_specs = [(c, label) for c, label in slope_specs if c in result.columns]
    if available_specs:
        fig, axes = plt.subplots(1, len(available_specs), figsize=(6.2 * len(available_specs), 5.5), squeeze=False)
        for ax, (column, label) in zip(axes[0], available_specs):
            x = pd.to_numeric(result["discovery_delta_beta"], errors="coerce")
            y = pd.to_numeric(result[column], errors="coerce")
            ax.scatter(x, y)
            ax.axhline(0, linewidth=1)
            ax.axvline(0, linewidth=1)
            for row, xv, yv in zip(result.itertuples(index=False), x, y):
                if np.isfinite(xv) and np.isfinite(yv):
                    ax.annotate(row.cpg_id, (xv, yv), fontsize=6)
            ax.set_xlabel("Partial-reprogramming Δβ")
            ax.set_ylabel(label)
        fig.suptitle("Aging-reversal quadrants for the frozen CpG panel")
        fig.tight_layout()
        for suffix in ["png", "pdf"]:
            fig.savefig(figures / f"fig_aging_reversal_quadrants_v090.{suffix}", dpi=400)
        plt.close(fig)

    print(json.dumps(summary, indent=2))
    print("Results:", output)


if __name__ == "__main__":
    main()
