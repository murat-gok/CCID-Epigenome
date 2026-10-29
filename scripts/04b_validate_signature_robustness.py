#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import matplotlib
# Figures are written to disk and never displayed, so force the
# non-interactive backend. Without this, matplotlib may select an
# interactive backend (e.g. TkAgg) and fail on machines with a
# broken or absent Tcl/Tk, and in any headless environment.
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ccid.robust_validation import run_robust_validation


def _json_clean(value: object) -> object:
    if isinstance(value, dict):
        return {str(k): _json_clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return [_json_clean(v) for v in value.tolist()]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_figures(result, output_dir: Path) -> None:
    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    sensitivity = result.baseline_sensitivity_metrics.copy()
    if not sensitivity.empty:
        labels = sensitivity["comparison"].str.replace("successful_vs_", "", regex=False)
        fig, ax = plt.subplots(figsize=(9, 5))
        positions = np.arange(len(sensitivity))
        ax.bar(positions - 0.18, sensitivity["raw_sign_concordance"], width=0.36, label="Raw")
        ax.bar(
            positions + 0.18,
            sensitivity["centered_sign_concordance"],
            width=0.36,
            label="Global-median centered",
        )
        ax.set_xticks(positions, labels, rotation=25, ha="right")
        ax.set_ylim(0, 1)
        ax.set_ylabel("Sign concordance")
        ax.set_title("O3 baseline and successful-vs-failed sensitivity")
        ax.legend()
        fig.tight_layout()
        fig.savefig(figures / "baseline_sensitivity.png", dpi=220)
        plt.close(fig)

    nulls = result.direction_balanced_null
    if not nulls.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        scales = list(nulls["analysis_scale"].drop_duplicates())
        values = [nulls.loc[nulls["analysis_scale"] == scale, "sign_concordance"] for scale in scales]
        ax.boxplot(values, tick_labels=scales, showfliers=False)
        for index, scale in enumerate(scales, start=1):
            observed = float(
                nulls.loc[nulls["analysis_scale"] == scale, "observed_sign_concordance"].iloc[0]
            )
            ax.scatter([index], [observed], marker="D", s=55, label="Observed" if index == 1 else None)
        ax.set_ylim(0, 1)
        ax.set_ylabel("Sign concordance")
        ax.set_title("Direction- and effect-matched empirical null")
        ax.legend()
        fig.tight_layout()
        fig.savefig(figures / "direction_balanced_null.png", dpi=220)
        plt.close(fig)

    lodo = result.lodo_summary
    if not lodo.empty:
        fig, ax = plt.subplots(figsize=(7, 5))
        x = np.arange(len(lodo))
        ax.bar(x - 0.18, lodo["raw_sign_concordance"], width=0.36, label="Raw")
        ax.bar(x + 0.18, lodo["centered_sign_concordance"], width=0.36, label="Centered")
        ax.set_xticks(x, lodo["holdout_donor"])
        ax.set_ylim(0, 1)
        ax.set_xlabel("Held-out donor")
        ax.set_ylabel("Sign concordance")
        ax.set_title("Leave-one-donor-out robustness")
        ax.legend()
        fig.tight_layout()
        fig.savefig(figures / "lodo_sign_concordance.png", dpi=220)
        plt.close(fig)

    time = result.time_series_metrics.copy()
    if not time.empty:
        fig, ax = plt.subplots(figsize=(9, 5))
        for (donor, condition, stage), group in time.groupby(
            ["donor_id", "condition", "stage"], sort=True
        ):
            group = group.sort_values("day")
            ax.plot(
                group["day"],
                group["sign_concordance"],
                marker="o",
                label=f"{donor} | {condition} | {stage}",
            )
        ax.set_ylim(0, 1)
        ax.set_xlabel("Day")
        ax.set_ylabel("Sign concordance")
        ax.set_title("Frozen-signature consistency across the reprogramming time series")
        ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(figures / "time_series_sign_concordance.png", dpi=220)
        plt.close(fig)

    if not result.successful_vs_failed.empty:
        frame = result.successful_vs_failed
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(
            frame["discovery_delta_beta"],
            frame["validation_delta_beta"],
            s=18,
            alpha=0.65,
        )
        lower = float(
            min(frame["discovery_delta_beta"].min(), frame["validation_delta_beta"].min(), 0)
        )
        upper = float(
            max(frame["discovery_delta_beta"].max(), frame["validation_delta_beta"].max(), 0)
        )
        ax.plot([lower, upper], [lower, upper], linestyle="--", linewidth=1)
        ax.axhline(0, linewidth=0.8)
        ax.axvline(0, linewidth=0.8)
        ax.set_xlabel("O1–O2 discovery Δβ")
        ax.set_ylabel("O3 successful − failed Δβ")
        ax.set_title("Direct successful-versus-failed comparison")
        fig.tight_layout()
        fig.savefig(figures / "successful_vs_failed.png", dpi=220)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run prespecified O3 robustness tests: successful-vs-failed, multi-baseline "
            "sensitivity, direction-balanced nulls, leave-one-donor-out, and time series"
        )
    )
    parser.add_argument("--beta", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--signature", required=True)
    parser.add_argument("--discovery-pairs", required=True)
    parser.add_argument("--all-loci", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--validation-donor", default="O3")
    parser.add_argument("--donors", default="O1,O2,O3")
    parser.add_argument("--success-condition", default="reprogrammed")
    parser.add_argument("--baseline-condition", default="baseline")
    parser.add_argument("--failed-condition", default="failed")
    parser.add_argument("--lodo-top-n", type=int, default=500)
    parser.add_argument("--lodo-min-abs-delta", type=float, default=0.03)
    parser.add_argument("--null-iterations", type=int, default=1000)
    parser.add_argument("--null-bins", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260721)
    parser.add_argument("--chunksize", type=int, default=10000)
    args = parser.parse_args()

    metadata_path = Path(args.metadata)
    signature_path = Path(args.signature)
    pairs_path = Path(args.discovery_pairs)
    all_loci_path = Path(args.all_loci)
    metadata = pd.read_csv(metadata_path)
    signature = pd.read_csv(signature_path)
    pairs = pd.read_csv(pairs_path)
    all_loci = pd.read_csv(
        all_loci_path,
        usecols=lambda column: column in {"cpg_id", "target_delta_beta"},
    )
    donors = [x.strip() for x in args.donors.split(",") if x.strip()]

    result = run_robust_validation(
        beta_path=args.beta,
        metadata=metadata,
        signature=signature,
        discovery_pair_table=pairs,
        all_loci=all_loci,
        validation_donor=args.validation_donor,
        donors=donors,
        success_condition=args.success_condition,
        baseline_condition=args.baseline_condition,
        failed_condition=args.failed_condition,
        lodo_top_n=args.lodo_top_n,
        lodo_min_abs_delta=args.lodo_min_abs_delta,
        null_iterations=args.null_iterations,
        null_bins=args.null_bins,
        seed=args.seed,
        chunksize=args.chunksize,
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    result.successful_vs_failed.to_csv(
        output_dir / "o3_successful_vs_failed.csv.gz", index=False
    )
    result.baseline_sensitivity_metrics.to_csv(
        output_dir / "o3_baseline_sensitivity_metrics.csv", index=False
    )
    result.baseline_sensitivity_cpg.to_csv(
        output_dir / "o3_baseline_sensitivity_cpg.csv.gz", index=False
    )
    result.direction_balanced_null.to_csv(
        output_dir / "o3_direction_balanced_null.csv.gz", index=False
    )
    result.lodo_summary.to_csv(output_dir / "lodo_summary.csv", index=False)
    result.lodo_candidates.to_csv(output_dir / "lodo_candidates.csv.gz", index=False)
    result.lodo_overlap.to_csv(output_dir / "lodo_overlap.csv", index=False)
    result.time_series_metrics.to_csv(output_dir / "time_series_metrics.csv", index=False)
    result.time_series_cpg_consistency.to_csv(
        output_dir / "time_series_cpg_consistency.csv.gz", index=False
    )
    result.sample_manifest.to_csv(output_dir / "time_series_sample_manifest.csv", index=False)
    (output_dir / "robustness_summary.json").write_text(
        json.dumps(_json_clean(result.summary), indent=2, allow_nan=False),
        encoding="utf-8",
    )
    manifest = {
        "package_version": "0.6.0",
        "analysis_prespecified": True,
        "candidate_list_frozen": True,
        "beta_path": str(Path(args.beta)),
        "metadata_path": str(metadata_path),
        "metadata_sha256": _sha256(metadata_path),
        "signature_path": str(signature_path),
        "signature_sha256": _sha256(signature_path),
        "discovery_pairs_path": str(pairs_path),
        "discovery_pairs_sha256": _sha256(pairs_path),
        "all_loci_path": str(all_loci_path),
        "all_loci_sha256": _sha256(all_loci_path),
        "validation_donor": args.validation_donor,
        "donors": donors,
        "lodo_top_n": args.lodo_top_n,
        "lodo_min_abs_delta": args.lodo_min_abs_delta,
        "null_iterations": args.null_iterations,
        "null_bins": args.null_bins,
        "seed": args.seed,
        "chunksize": args.chunksize,
    }
    (output_dir / "robustness_analysis_manifest.json").write_text(
        json.dumps(_json_clean(manifest), indent=2, allow_nan=False),
        encoding="utf-8",
    )
    _write_figures(result, output_dir)

    print(f"Frozen stratum: {result.summary['frozen_stratum']}")
    if not result.successful_vs_failed.empty:
        metric = result.baseline_sensitivity_metrics
        row = metric[metric["comparison"] == "successful_vs_failed_same_stratum"].iloc[0]
        print(f"Successful-vs-failed sign concordance: {row['raw_sign_concordance']:.4f}")
        print(f"Successful-vs-failed centered concordance: {row['centered_sign_concordance']:.4f}")
    print("LODO holdouts:")
    for row in result.lodo_summary.itertuples(index=False):
        print(
            f"  {row.holdout_donor}: raw={row.raw_sign_concordance:.4f}, "
            f"centered={row.centered_sign_concordance:.4f}"
        )
    print(f"Time-series comparisons: {len(result.time_series_metrics)}")
    print(f"Results: {output_dir}")


if __name__ == "__main__":
    main()
