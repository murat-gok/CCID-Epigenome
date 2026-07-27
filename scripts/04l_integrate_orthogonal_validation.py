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

from ccid.orthogonal_validation import integrate_orthogonal_evidence


def existing(path: str) -> str | None:
    return path if path and Path(path).exists() else None


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Integrate frozen-panel external replication, aging direction, eQTM, rE2G and editability evidence."
    )
    parser.add_argument("--frozen-shortlist", default="data/reference/frozen_primary_shortlist_v090.csv")
    parser.add_argument(
        "--external-file",
        default="results/real/orthogonal_validation_v090/external_reprogramming/external_reprogramming_validation_v090.csv.gz",
    )
    parser.add_argument(
        "--aging-file",
        default="results/real/orthogonal_validation_v090/natural_aging/natural_aging_direction_validation_v090.csv.gz",
    )
    parser.add_argument(
        "--eqtm-file",
        default="results/real/orthogonal_validation_v090/eqtm/fibroblast_eqtm_validation_v090.csv",
    )
    parser.add_argument(
        "--re2g-file",
        default="results/real/orthogonal_validation_v090/regulatory_editability/fibroblast_re2g_candidate_links_v090.csv",
    )
    parser.add_argument(
        "--guide-summary-file",
        default="results/real/orthogonal_validation_v090/regulatory_editability/candidate_guide_feasibility_summary_v090.csv",
    )
    parser.add_argument("--output-dir", default="results/real/orthogonal_validation_v090/integrated")
    args = parser.parse_args()

    for value in [args.frozen_shortlist, args.external_file, args.aging_file]:
        if not Path(value).exists():
            raise FileNotFoundError(value)
    output = Path(args.output_dir)
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    result, summary = integrate_orthogonal_evidence(
        args.frozen_shortlist,
        args.external_file,
        args.aging_file,
        eqtm_file=existing(args.eqtm_file),
        re2g_file=existing(args.re2g_file),
        guide_summary_file=existing(args.guide_summary_file),
    )
    result.to_csv(output / "orthogonally_validated_candidate_panel_v090.csv", index=False)
    (output / "orthogonal_validation_summary_v090.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    evidence_cols = [
        c
        for c in [
            "external_mean_direction_concordant",
            "replicative_reversal_concordant",
            "chronological_reversal_concordant",
            "multidonor_reversal_concordant",
            "eqtm_supported",
            "re2g_linked",
            "guide_feasible",
        ]
        if c in result.columns
    ]
    if evidence_cols:
        plot = result.sort_values("robust_rank", ascending=True)
        matrix = plot[evidence_cols].fillna(False).astype(int).to_numpy()
        fig, ax = plt.subplots(figsize=(9.0, 7.0))
        image = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1)
        ax.set_yticks(np.arange(len(plot)), plot["cpg_id"])
        ax.set_xticks(
            np.arange(len(evidence_cols)),
            [c.replace("_mean_direction_concordant", "").replace("_reversal_concordant", "").replace("_", " ") for c in evidence_cols],
            rotation=35,
            ha="right",
        )
        ax.set_title("Orthogonal evidence matrix for the frozen CpG panel")
        fig.colorbar(image, ax=ax, label="Evidence present")
        fig.tight_layout()
        for suffix in ["png", "pdf"]:
            fig.savefig(figures / f"fig_integrated_orthogonal_evidence_v090.{suffix}", dpi=400)
        plt.close(fig)

    print(json.dumps(summary, indent=2))
    cols = [c for c in ["robust_rank", "cpg_id", "candidate_identity", "orthogonal_support_count", "orthogonal_validation_tier"] if c in result.columns]
    print(result[cols].to_string(index=False))
    print("Results:", output)


if __name__ == "__main__":
    main()
