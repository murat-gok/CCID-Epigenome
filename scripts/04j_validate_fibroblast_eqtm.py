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

from ccid.orthogonal_validation import validate_fibroblast_eqtm


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Test CpG–gene expression associations in 62 primary fibroblast lines."
    )
    parser.add_argument(
        "--candidate-file",
        default="results/real/regulatory_safety_v081/candidate_cpg_regulatory_safety_v081.csv.gz",
    )
    parser.add_argument("--frozen-shortlist", default="data/reference/frozen_primary_shortlist_v090.csv")
    parser.add_argument(
        "--methylation-matrix",
        default="data/raw/orthogonal_validation/GSE52025_series_matrix.txt.gz",
    )
    parser.add_argument(
        "--expression-matrix",
        default="data/raw/orthogonal_validation/GSE53243_series_matrix.txt.gz",
    )
    parser.add_argument(
        "--expression-annotation",
        default="data/raw/orthogonal_validation/GPL6883_data.txt",
    )
    parser.add_argument("--gene-map", default="data/reference/frozen_primary_shortlist_v090.csv")
    parser.add_argument("--output-dir", default="results/real/orthogonal_validation_v090/eqtm")
    args = parser.parse_args()

    for value in [
        args.candidate_file,
        args.frozen_shortlist,
        args.methylation_matrix,
        args.expression_matrix,
        args.expression_annotation,
        args.gene_map,
    ]:
        if not Path(value).exists():
            raise FileNotFoundError(value)
    output = Path(args.output_dir)
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    result, summary = validate_fibroblast_eqtm(
        args.candidate_file,
        args.frozen_shortlist,
        args.methylation_matrix,
        args.expression_matrix,
        args.expression_annotation,
        args.gene_map,
    )
    result.to_csv(output / "fibroblast_eqtm_validation_v090.csv", index=False)
    (output / "fibroblast_eqtm_summary_v090.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    if not result.empty and "eqtm_spearman_rho" in result.columns:
        plot = result[result["eqtm_spearman_rho"].notna()].copy()
        plot["label"] = plot["cpg_id"].astype(str) + " → " + plot["gene"].astype(str)
        plot = plot.sort_values("eqtm_spearman_rho")
        fig, ax = plt.subplots(figsize=(8.0, max(4.0, 0.35 * len(plot) + 1.5)))
        y = np.arange(len(plot))
        ax.barh(y, plot["eqtm_spearman_rho"])
        ax.axvline(0, linewidth=1)
        ax.set_yticks(y, plot["label"])
        ax.set_xlabel("Spearman correlation between CpG methylation and gene expression")
        ax.set_title("Cross-sectional fibroblast eQTM evidence")
        for yi, row in enumerate(plot.itertuples(index=False)):
            if np.isfinite(row.eqtm_fdr):
                ax.text(row.eqtm_spearman_rho, yi, f"  q={row.eqtm_fdr:.3g}", va="center", fontsize=7)
        fig.tight_layout()
        for suffix in ["png", "pdf"]:
            fig.savefig(figures / f"fig_fibroblast_eqtm_v090.{suffix}", dpi=400)
        plt.close(fig)

    print(json.dumps(summary, indent=2))
    print("Results:", output)


if __name__ == "__main__":
    main()
