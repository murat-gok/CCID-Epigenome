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

from ccid.o3_validation import (
    build_o3_pair_table,
    infer_validation_stratum,
    validate_o3_signature_file,
)


def _json_clean(value: object) -> object:
    if isinstance(value, dict):
        return {str(k): _json_clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return [_json_clean(v) for v in value.tolist()]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_figures(result, output_dir: Path, seed: int) -> None:
    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    frame = result.cpg_validation
    valid = frame[frame["validation_delta_beta"].notna()].copy()

    fig, ax = plt.subplots(figsize=(7, 6))
    for tier, group in valid.groupby("tier", sort=True):
        ax.scatter(
            group["discovery_delta_beta"],
            group["validation_delta_beta"],
            alpha=0.60,
            s=18,
            label=tier,
        )
    if len(valid):
        lower = float(
            min(
                valid["discovery_delta_beta"].min(),
                valid["validation_delta_beta"].min(),
                0.0,
            )
        )
        upper = float(
            max(
                valid["discovery_delta_beta"].max(),
                valid["validation_delta_beta"].max(),
                0.0,
            )
        )
        ax.plot([lower, upper], [lower, upper], linestyle="--", linewidth=1, label="y=x")
    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)
    ax.set_xlabel("O1–O2 discovery Δβ")
    ax.set_ylabel("O3 validation Δβ")
    ax.set_title("Held-out O3 validation of the frozen CpG signature")
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures / "delta_correlation.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(frame["discovery_delta_beta"].dropna(), bins=40, alpha=0.55, label="Discovery")
    ax.hist(frame["validation_delta_beta"].dropna(), bins=40, alpha=0.55, label="O3")
    ax.set_xlabel("Δβ")
    ax.set_ylabel("CpG count")
    ax.set_title("Discovery and O3 effect distributions")
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures / "delta_distribution.png", dpi=220)
    plt.close(fig)

    topk = result.topk_metrics
    fig, ax = plt.subplots(figsize=(7, 5))
    if not topk.empty:
        ax.plot(topk["top_k"], topk["sign_concordance"], marker="o", label="Observed")
    ax.axhline(0.5, linestyle="--", linewidth=1, label="Chance direction")
    ax.set_ylim(0, 1)
    ax.set_xlabel("Top-k frozen discovery CpGs")
    ax.set_ylabel("O3 sign concordance")
    ax.set_title("Direction concordance across candidate-list depth")
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures / "sign_concordance_topk.png", dpi=220)
    plt.close(fig)

    rng = np.random.default_rng(seed)
    sampled_values: list[np.ndarray] = []
    labels: list[str] = []
    for sample, values in result.sample_values.items():
        finite = values[np.isfinite(values)]
        if len(finite) > 20000:
            finite = rng.choice(finite, size=20000, replace=False)
        sampled_values.append(finite)
        labels.append(sample)
    fig, ax = plt.subplots(figsize=(max(7, 1.4 * len(labels)), 5))
    if sampled_values:
        ax.boxplot(sampled_values, tick_labels=labels, showfliers=False)
    ax.set_ylabel("Beta value")
    ax.set_title("Genome-wide beta-distribution QC for O3 and diagnostic controls")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(figures / "sample_qc.png", dpi=220)
    plt.close(fig)

    summary_rows: list[tuple[str, float, float]] = []
    for null_type, label in (
        ("cpg_label_permutation", "CpG-label permutation"),
        ("effect_matched_random_cpg_set", "Effect-matched random sets"),
    ):
        subset = result.null_tests[result.null_tests["null_type"] == null_type]
        if subset.empty or "sign_concordance" not in subset:
            continue
        summary_rows.append(
            (
                label,
                float(subset["sign_concordance"].mean()),
                float(subset["sign_concordance"].std(ddof=1)) if len(subset) > 1 else 0.0,
            )
        )
    if summary_rows:
        fig, ax = plt.subplots(figsize=(7, 5))
        positions = np.arange(len(summary_rows))
        ax.bar(
            positions,
            [row[1] for row in summary_rows],
            yerr=[row[2] for row in summary_rows],
            capsize=4,
        )
        observed = float(result.summary["overall_metrics"]["sign_concordance"])
        ax.axhline(observed, linestyle="--", linewidth=2, label="Observed O3")
        ax.set_xticks(positions, [row[0] for row in summary_rows], rotation=20, ha="right")
        ax.set_ylim(0, 1)
        ax.set_ylabel("Sign concordance")
        ax.set_title("Observed O3 concordance versus empirical nulls")
        ax.legend()
        fig.tight_layout()
        fig.savefig(figures / "null_sign_concordance.png", dpi=220)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate the frozen O1/O2 CpG signature in held-out donor O3"
    )
    parser.add_argument("--beta", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--signature", required=True)
    parser.add_argument(
        "--discovery-pairs",
        required=True,
        help="Frozen O1/O2 matched-pair table used to infer day/batch/stage",
    )
    parser.add_argument(
        "--all-loci",
        help="Optional all-loci discovery file for effect-matched random CpG-set nulls",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--validation-donor", default="O3")
    parser.add_argument("--discovery-donors", default="O1,O2")
    parser.add_argument("--validation-condition", default="reprogrammed")
    parser.add_argument("--baseline-condition", default="baseline")
    parser.add_argument("--top-k", default="50,100,250,500")
    parser.add_argument("--tier1-min-abs-delta", type=float, default=0.03)
    parser.add_argument("--tier1-min-retention", type=float, default=0.25)
    parser.add_argument("--tier1-max-retention", type=float, default=2.0)
    parser.add_argument("--tier2-min-abs-delta", type=float, default=0.01)
    parser.add_argument("--tier2-min-retention", type=float, default=0.10)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--random-sets", type=int, default=1000)
    parser.add_argument("--random-effect-bins", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260721)
    parser.add_argument("--chunksize", type=int, default=10000)
    args = parser.parse_args()

    metadata_path = Path(args.metadata)
    signature_path = Path(args.signature)
    discovery_pairs_path = Path(args.discovery_pairs)
    metadata = pd.read_csv(metadata_path)
    signature = pd.read_csv(signature_path)
    discovery_pairs = pd.read_csv(discovery_pairs_path)
    frozen_stratum = infer_validation_stratum(discovery_pairs)

    all_loci_path: Path | None = Path(args.all_loci) if args.all_loci else None
    if all_loci_path is None:
        stem = (
            signature_path.name[:-7]
            if signature_path.name.endswith(".csv.gz")
            else signature_path.stem
        )
        inferred = signature_path.parent / f"{stem}_all_loci.csv.gz"
        if inferred.exists():
            all_loci_path = inferred
    all_loci = (
        pd.read_csv(
            all_loci_path,
            usecols=lambda column: column in {"cpg_id", "target_delta_beta"},
        )
        if all_loci_path is not None
        else None
    )
    donors = [x.strip() for x in args.discovery_donors.split(",") if x.strip()]
    top_k = [int(x.strip()) for x in args.top_k.split(",") if x.strip()]

    result = validate_o3_signature_file(
        beta_path=args.beta,
        metadata=metadata,
        discovery_signature=signature,
        discovery_pair_table=discovery_pairs,
        discovery_all_loci=all_loci,
        validation_donor=args.validation_donor,
        validation_condition=args.validation_condition,
        baseline_condition=args.baseline_condition,
        discovery_donors=donors,
        top_k_values=top_k,
        tier1_min_abs_delta=args.tier1_min_abs_delta,
        tier1_min_effect_ratio=args.tier1_min_retention,
        tier1_max_effect_ratio=args.tier1_max_retention,
        tier2_min_abs_delta=args.tier2_min_abs_delta,
        tier2_min_effect_ratio=args.tier2_min_retention,
        n_permutations=args.permutations,
        random_seed=args.seed,
        n_random_sets=args.random_sets,
        random_effect_bins=args.random_effect_bins,
        chunksize=args.chunksize,
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    primary_pair, auxiliary_controls = build_o3_pair_table(
        metadata,
        discovery_pairs,
        validation_donor=args.validation_donor,
        validation_condition=args.validation_condition,
        baseline_condition=args.baseline_condition,
        discovery_donors=donors,
    )
    primary_pair.to_csv(output_dir / "o3_primary_pair.csv", index=False)
    auxiliary_controls.to_csv(output_dir / "o3_auxiliary_controls.csv", index=False)
    result.cpg_validation.to_csv(output_dir / "o3_cpg_validation.csv.gz", index=False)
    for tier in ("Tier1", "Tier2", "Tier3"):
        result.cpg_validation[result.cpg_validation["tier"] == tier].to_csv(
            output_dir / f"o3_{tier.lower()}_candidates.csv.gz", index=False
        )
    result.null_tests.to_csv(output_dir / "o3_null_tests.csv", index=False)
    result.sample_qc.to_csv(output_dir / "o3_sample_qc.csv", index=False)
    result.topk_metrics.to_csv(output_dir / "o3_topk_metrics.csv", index=False)
    (output_dir / "o3_validation_summary.json").write_text(
        json.dumps(_json_clean(result.summary), indent=2, allow_nan=False),
        encoding="utf-8",
    )

    manifest = {
        "package_version": "0.5.0",
        "candidate_list_frozen_before_o3": True,
        "signature_path": str(signature_path),
        "signature_sha256": _sha256(signature_path),
        "discovery_pairs_path": str(discovery_pairs_path),
        "discovery_pairs_sha256": _sha256(discovery_pairs_path),
        "metadata_path": str(metadata_path),
        "metadata_sha256": _sha256(metadata_path),
        "all_loci_path": str(all_loci_path) if all_loci_path else None,
        "frozen_validation_stratum": frozen_stratum,
        "validation_donor": args.validation_donor,
        "random_seed": args.seed,
        "tier_thresholds": {
            "tier1_min_abs_delta": args.tier1_min_abs_delta,
            "tier1_min_retention": args.tier1_min_retention,
            "tier1_max_retention": args.tier1_max_retention,
            "tier2_min_abs_delta": args.tier2_min_abs_delta,
            "tier2_min_retention": args.tier2_min_retention,
        },
        "n_permutations": args.permutations,
        "n_effect_matched_random_sets": args.random_sets,
    }
    (output_dir / "o3_frozen_analysis_manifest.json").write_text(
        json.dumps(_json_clean(manifest), indent=2, allow_nan=False),
        encoding="utf-8",
    )
    _write_figures(result, output_dir, args.seed)

    metrics = result.summary["overall_metrics"]
    print(
        "Frozen stratum: "
        f"day={result.summary['day']}, batch={result.summary['batch']}, "
        f"stage={result.summary['stage']}"
    )
    print(f"Validation samples: {result.summary['validation_sample_ids']}")
    print(f"Baseline samples: {result.summary['baseline_sample_ids']}")
    print(f"Valid frozen CpGs: {result.summary['n_o3_valid_cpgs']}")
    print(f"Sign concordance: {metrics['sign_concordance']:.4f}")
    print(f"Spearman rho: {metrics['spearman_rho']:.4f}")
    print(f"Tier counts: {result.summary['tier_counts']}")
    print(f"Results: {output_dir}")


if __name__ == "__main__":
    main()
