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

from ccid.orthogonal_validation import (
    design_candidate_guides,
    overlap_candidates_with_re2g,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Assess optional fibroblast rE2G links and SpCas9 sequence feasibility for frozen candidates."
    )
    parser.add_argument(
        "--candidate-file",
        default="results/real/regulatory_safety_v081/candidate_cpg_regulatory_safety_v081.csv.gz",
    )
    parser.add_argument("--frozen-shortlist", default="data/reference/frozen_primary_shortlist_v090.csv")
    parser.add_argument("--re2g-files", nargs="*", default=[])
    parser.add_argument("--biosample-names", nargs="*", default=[])
    parser.add_argument("--skip-guides", action="store_true")
    parser.add_argument("--flank", type=int, default=250)
    parser.add_argument("--max-distance", type=int, default=150)
    parser.add_argument("--genome", default="hg38")
    parser.add_argument("--output-dir", default="results/real/orthogonal_validation_v090/regulatory_editability")
    args = parser.parse_args()

    for value in [args.candidate_file, args.frozen_shortlist]:
        if not Path(value).exists():
            raise FileNotFoundError(value)
    output = Path(args.output_dir)
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    summary: dict[str, object] = {
        "scope": "Sequence feasibility and predicted regulatory linkage; neither establishes editing specificity or causality."
    }

    if args.re2g_files:
        for value in args.re2g_files:
            if not Path(value).exists():
                raise FileNotFoundError(value)
        names = args.biosample_names or [Path(p).stem for p in args.re2g_files]
        if len(names) != len(args.re2g_files):
            raise ValueError("--biosample-names must match --re2g-files")
        re2g = overlap_candidates_with_re2g(
            args.candidate_file,
            args.frozen_shortlist,
            args.re2g_files,
            biosample_names=names,
        )
        re2g.to_csv(output / "fibroblast_re2g_candidate_links_v090.csv", index=False)
        summary["re2g_candidate_rows"] = int(len(re2g))
        summary["re2g_linked_candidates"] = int(re2g.groupby("cpg_id")["re2g_linked"].max().sum())

    if not args.skip_guides:
        guides, guide_summary = design_candidate_guides(
            args.candidate_file,
            args.frozen_shortlist,
            flank=args.flank,
            max_distance=args.max_distance,
            genome=args.genome,
        )
        guides.to_csv(output / "candidate_spcas9_guides_v090.csv.gz", index=False, compression="gzip")
        guide_summary.to_csv(output / "candidate_guide_feasibility_summary_v090.csv", index=False)
        if not guides.empty:
            batch = guides[guides["basic_sequence_pass"].fillna(False)].copy()
            batch["guide_name"] = batch["cpg_id"].astype(str) + "_" + batch.groupby("cpg_id").cumcount().add(1).astype(str)
            batch[["guide_name", "cpg_id", "chromosome", "cpg_position", "protospacer", "pam", "strand", "distance_to_cpg"]].to_csv(
                output / "guide_sequences_for_genomewide_offtarget_tools_v090.tsv", sep="\t", index=False
            )
        summary["guide_feasible_candidates"] = int(guide_summary["guide_feasible"].fillna(False).sum())
        summary["sequence_fetch_failures"] = int(guide_summary["sequence_fetch_error"].fillna("").ne("").sum())

        plot = guide_summary.sort_values("position", ascending=False).copy()
        fig, ax = plt.subplots(figsize=(8.0, 6.5))
        y = np.arange(len(plot))
        ax.barh(y, plot["n_basic_sequence_pass"])
        ax.set_yticks(y, plot["cpg_id"])
        ax.set_xlabel("SpCas9 NGG guides passing basic sequence filters")
        ax.set_title("Sequence-level editability of the frozen CpG panel")
        fig.tight_layout()
        for suffix in ["png", "pdf"]:
            fig.savefig(figures / f"fig_candidate_guide_feasibility_v090.{suffix}", dpi=400)
        plt.close(fig)

    (output / "regulatory_editability_summary_v090.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    print("Results:", output)


if __name__ == "__main__":
    main()
