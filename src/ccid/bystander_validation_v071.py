from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import json
import math
import re

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .bystander_validation import (
    TARGETED_GENES,
    _empirical_p,
    _normalise_cpg_column,
    _safe_pearson,
    _safe_spearman,
    build_coverage_and_candidate_effects,
    build_cross_celltype_results,
    build_paired_celltype_effects,
    load_candidate_sets,
    load_probe_annotation,
    run_matched_null_tests,
)
from .preprocessing import _read_bystander_wide


@dataclass(frozen=True)
class PlatformAwareEffects:
    locus_effects: pd.DataFrame
    pair_effects: pd.DataFrame
    pair_manifest: pd.DataFrame
    platform_qc: pd.DataFrame


@dataclass(frozen=True)
class BystanderValidationV071Result:
    candidate_coverage: pd.DataFrame
    candidate_set_summary: pd.DataFrame
    candidate_effects: pd.DataFrame
    cross_celltype_effects: pd.DataFrame
    cross_celltype_summary: pd.DataFrame
    tcell_platform_meta_effects: pd.DataFrame
    tcell_platform_pair_effects: pd.DataFrame
    tcell_platform_manifest: pd.DataFrame
    tcell_platform_qc: pd.DataFrame
    msc_effects: pd.DataFrame
    positive_control_probe_effects: pd.DataFrame
    positive_control_gene_summary: pd.DataFrame
    positive_control_distance_decay: pd.DataFrame
    unique_region_membership: pd.DataFrame
    unique_region_scores_all: pd.DataFrame
    unique_region_scores_multicpg: pd.DataFrame
    coverage_bias_per_candidate: pd.DataFrame
    coverage_bias_summary: pd.DataFrame
    matched_null_tests: pd.DataFrame
    matched_null_summary: pd.DataFrame
    cross_celltype_null_tests: pd.DataFrame
    cross_celltype_null_summary: pd.DataFrame
    summary: dict[str, object]
    manifest: dict[str, object]


def _jsonable(value: object) -> object:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def _mad_variance(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 10:
        return 1.0
    med = np.median(values)
    mad = np.median(np.abs(values - med))
    scale = 1.4826 * mad
    variance = float(scale * scale)
    if not np.isfinite(variance) or variance <= 1e-10:
        variance = float(np.nanvar(values))
    return max(variance, 1e-8)


def _pair_deltas_from_platform(
    beta: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    platform_id: str,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    beta = beta.copy()
    beta.index = beta.index.astype(str).str.replace(r"_(BC|TC)\d+$", "", regex=True)
    if beta.index.duplicated().any():
        beta = beta.groupby(level=0, sort=False).mean()
    beta = beta.apply(pd.to_numeric, errors="coerce")

    meta = metadata.copy()
    required = {"sample_id", "treatment", "replicate"}
    missing = required.difference(meta.columns)
    if missing:
        raise ValueError(f"Platform metadata missing columns: {sorted(missing)}")
    meta["sample_id"] = meta["sample_id"].astype(str)
    meta = meta[meta["sample_id"].isin(beta.columns)].copy()

    effects: dict[str, pd.Series] = {}
    baselines: dict[str, pd.Series] = {}
    rows: list[dict[str, object]] = []
    for replicate, group in meta.groupby("replicate", sort=True):
        targeted = group.loc[group["treatment"].astype(str) == "targeted", "sample_id"].tolist()
        controls = group.loc[group["treatment"].astype(str) == "control", "sample_id"].tolist()
        if not targeted or not controls:
            continue
        pair_id = str(replicate)
        col = f"{platform_id}:{pair_id}"
        effects[col] = beta[targeted].mean(axis=1) - beta[controls].mean(axis=1)
        baselines[col] = beta[controls].mean(axis=1)
        rows.append(
            {
                "platform_id": platform_id,
                "pair_id": pair_id,
                "pair_column": col,
                "targeted_sample_ids": ";".join(targeted),
                "control_sample_ids": ";".join(controls),
                "n_targeted": len(targeted),
                "n_controls": len(controls),
                "n_platform_loci": int(beta.shape[0]),
            }
        )
    if not effects:
        raise ValueError(f"No complete targeted-control pairs found for {platform_id}")
    return pd.DataFrame(effects), pd.DataFrame(baselines), rows


def build_platform_aware_tcell_effects(
    tcell_raw_paths: Sequence[str | Path],
    *,
    min_abs_delta: float = 0.02,
    fdr: float = 0.10,
) -> PlatformAwareEffects:
    effect_frames: list[pd.DataFrame] = []
    baseline_frames: list[pd.DataFrame] = []
    manifest_rows: list[dict[str, object]] = []
    platform_qc_rows: list[dict[str, object]] = []

    for index, path_like in enumerate(tcell_raw_paths, start=1):
        path = Path(path_like)
        if not path.exists():
            raise FileNotFoundError(path)
        beta, metadata, qc = _read_bystander_wide(path)
        platform_id = f"platform_{index}_{re.sub(r'[^A-Za-z0-9]+', '_', path.stem)[:48]}"
        effects, baselines, manifest = _pair_deltas_from_platform(
            beta, metadata, platform_id=platform_id
        )
        effect_frames.append(effects)
        baseline_frames.append(baselines)
        manifest_rows.extend(manifest)
        platform_qc_rows.append(
            {
                "platform_id": platform_id,
                "source_file": str(path),
                "n_loci": int(beta.shape[0]),
                "n_samples": int(beta.shape[1]),
                "n_complete_pairs": len(manifest),
                **{k: v for k, v in qc.items() if k not in {"source_file", "n_loci", "n_samples"}},
            }
        )

    pair_effects = pd.concat(effect_frames, axis=1, join="outer")
    pair_baselines = pd.concat(baseline_frames, axis=1, join="outer")
    pair_manifest = pd.DataFrame(manifest_rows)

    pair_variance: dict[str, float] = {}
    for col in pair_effects.columns:
        pair_variance[col] = _mad_variance(pair_effects[col].to_numpy(dtype=float))
    pair_manifest["empirical_pair_variance"] = pair_manifest["pair_column"].map(pair_variance)

    values = pair_effects.to_numpy(dtype=float, copy=True)
    baselines = pair_baselines.reindex(index=pair_effects.index, columns=pair_effects.columns).to_numpy(dtype=float)
    variances = np.array([pair_variance[c] for c in pair_effects.columns], dtype=float)

    rows: list[dict[str, object]] = []
    for idx, cpg_id in enumerate(pair_effects.index.astype(str)):
        y = values[idx]
        b = baselines[idx]
        valid = np.isfinite(y)
        k = int(valid.sum())
        if k == 0:
            rows.append({"cpg_id": cpg_id, "cell_type": "t_cell", "n_pairs": 0})
            continue
        yv = y[valid]
        vv = variances[valid]
        w = 1.0 / vv
        fixed = float(np.sum(w * yv) / np.sum(w))
        fixed_se = float(math.sqrt(1.0 / np.sum(w)))
        q = float(np.sum(w * (yv - fixed) ** 2)) if k > 1 else 0.0
        c = float(np.sum(w) - np.sum(w**2) / np.sum(w)) if k > 1 else 0.0
        tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
        wr = 1.0 / (vv + tau2)
        random_mean = float(np.sum(wr * yv) / np.sum(wr))
        random_se = float(math.sqrt(1.0 / np.sum(wr)))
        random_p = float(2 * stats.norm.sf(abs(random_mean / random_se))) if random_se > 0 else 1.0
        direction = np.sign(random_mean)
        sign_consistency = float(np.mean(np.sign(yv) == direction)) if direction != 0 else 0.0
        baseline = float(np.nanmean(b[valid])) if np.isfinite(b[valid]).any() else float("nan")
        rows.append(
            {
                "cpg_id": cpg_id,
                "cell_type": "t_cell",
                "baseline_beta": baseline,
                "targeted_beta": baseline + random_mean if np.isfinite(baseline) else float("nan"),
                "editor_delta_beta": random_mean,
                "available_case_mean_delta": float(np.mean(yv)),
                "median_editor_delta_beta": float(np.median(yv)),
                "fixed_effect_delta": fixed,
                "fixed_effect_se": fixed_se,
                "random_effect_delta": random_mean,
                "random_effect_se": random_se,
                "tau2": tau2,
                "cochran_q": q,
                "n_pairs": k,
                "sign_consistency": sign_consistency,
                "p_value": random_p,
            }
        )

    locus = pd.DataFrame(rows)
    valid_p = locus["p_value"].fillna(1.0).to_numpy(dtype=float)
    locus["q_value"] = multipletests(valid_p, method="fdr_bh")[1]
    locus["significant_editor_effect"] = (
        (locus["q_value"] <= fdr)
        & (locus["editor_delta_beta"].abs() >= min_abs_delta)
        & (locus["n_pairs"] >= 2)
    )
    locus["editor_direction"] = np.select(
        [locus["editor_delta_beta"] > 0, locus["editor_delta_beta"] < 0],
        ["hypermethylating", "hypomethylating"],
        default="neutral",
    )

    pair_long = (
        pair_effects.rename_axis("cpg_id")
        .reset_index()
        .melt(id_vars="cpg_id", var_name="pair_column", value_name="pair_delta_beta")
        .merge(pair_manifest[["pair_column", "platform_id", "pair_id"]], on="pair_column", how="left")
    )
    pair_long.insert(1, "cell_type", "t_cell")
    return PlatformAwareEffects(
        locus_effects=locus,
        pair_effects=pair_long,
        pair_manifest=pair_manifest,
        platform_qc=pd.DataFrame(platform_qc_rows),
    )


def build_candidate_metadata(
    *,
    stable_core_path: str | Path,
    consensus_2of3_path: str | Path,
    tier1_path: str | Path,
    frozen_signature_path: str | Path,
) -> tuple[dict[str, set[str]], pd.DataFrame, pd.DataFrame]:
    candidate_sets, target_map = load_candidate_sets(
        stable_core_path=stable_core_path,
        consensus_2of3_path=consensus_2of3_path,
        tier1_path=tier1_path,
        frozen_signature_path=frozen_signature_path,
    )
    frozen = _normalise_cpg_column(pd.read_csv(frozen_signature_path, compression="infer"))
    keep = [c for c in ["cpg_id", "target_delta_beta", "rank_score", "weight"] if c in frozen.columns]
    metadata = frozen[keep].copy()

    consensus = _normalise_cpg_column(pd.read_csv(consensus_2of3_path, compression="infer"))
    if "n_lodo_rotations" in consensus.columns:
        metadata = metadata.merge(
            consensus[["cpg_id", "n_lodo_rotations"]], on="cpg_id", how="outer"
        )
    else:
        metadata["n_lodo_rotations"] = np.nan

    membership = pd.DataFrame({"cpg_id": sorted(set().union(*candidate_sets.values()))})
    for name, cpgs in candidate_sets.items():
        membership[f"in_{name}"] = membership["cpg_id"].isin(cpgs)
    metadata = membership.merge(metadata, on="cpg_id", how="left")
    return candidate_sets, target_map, metadata


def build_positive_control_analysis(
    annotation: pd.DataFrame,
    tcell_effects: pd.DataFrame,
    msc_effects: pd.DataFrame,
    *,
    max_distance_bp: int = 5000,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if annotation.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    ann = annotation.copy()
    ann["gene_tokens"] = ann["gene"].fillna("").astype(str).map(
        lambda x: {g.strip().upper() for g in re.split(r"[;,]", x) if g.strip()}
    )
    anchors: list[dict[str, object]] = []
    probe_rows: list[pd.DataFrame] = []
    for gene in sorted(TARGETED_GENES):
        gene_probes = ann[ann["gene_tokens"].map(lambda s: gene in s)].copy()
        if gene_probes.empty:
            continue
        for chrom, group in gene_probes.groupby("chromosome"):
            anchor = int(np.median(group["position"].to_numpy(dtype=int)))
            anchors.append(
                {
                    "gene": gene,
                    "chromosome": str(chrom),
                    "annotation_derived_anchor": anchor,
                    "n_gene_annotated_probes": int(group["cpg_id"].nunique()),
                }
            )
            local = ann[
                (ann["chromosome"].astype(str) == str(chrom))
                & ((ann["position"] - anchor).abs() <= max_distance_bp)
            ].copy()
            local["positive_control_gene"] = gene
            local["annotation_derived_anchor"] = anchor
            local["distance_to_anchor_bp"] = (local["position"] - anchor).abs().astype(int)
            local["is_gene_annotated_probe"] = local["gene_tokens"].map(lambda s: gene in s)
            probe_rows.append(local)
    if not probe_rows:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(anchors)
    probes = pd.concat(probe_rows, ignore_index=True).drop(columns="gene_tokens")

    t = tcell_effects[["cpg_id", "editor_delta_beta", "sign_consistency", "n_pairs"]].rename(
        columns={
            "editor_delta_beta": "tcell_delta",
            "sign_consistency": "tcell_sign_consistency",
            "n_pairs": "tcell_n_pairs",
        }
    )
    m = msc_effects[["cpg_id", "editor_delta_beta", "sign_consistency", "n_pairs"]].rename(
        columns={
            "editor_delta_beta": "msc_delta",
            "sign_consistency": "msc_sign_consistency",
            "n_pairs": "msc_n_pairs",
        }
    )
    probes = probes.merge(t, on="cpg_id", how="left").merge(m, on="cpg_id", how="left")
    probes["mean_abs_editor_delta"] = np.nanmean(
        np.abs(probes[["tcell_delta", "msc_delta"]].to_numpy(dtype=float)), axis=1
    )

    gene_rows: list[dict[str, object]] = []
    for gene, group in probes.groupby("positive_control_gene", sort=True):
        annotated = group[group["is_gene_annotated_probe"]]
        for label, subset in [("gene_annotated", annotated), ("within_5kb_anchor", group)]:
            gene_rows.append(
                {
                    "positive_control_gene": gene,
                    "scope": label,
                    "n_probes": int(subset["cpg_id"].nunique()),
                    "n_tcell_measured": int(subset["tcell_delta"].notna().sum()),
                    "n_msc_measured": int(subset["msc_delta"].notna().sum()),
                    "median_abs_tcell_delta": float(subset["tcell_delta"].abs().median()),
                    "median_abs_msc_delta": float(subset["msc_delta"].abs().median()),
                    "mean_tcell_sign_consistency": float(subset["tcell_sign_consistency"].mean()),
                    "mean_msc_sign_consistency": float(subset["msc_sign_consistency"].mean()),
                    "max_mean_abs_editor_delta": float(subset["mean_abs_editor_delta"].max()),
                }
            )

    bins = [-1, 250, 500, 1000, 2000, 5000]
    labels = ["0-250", "251-500", "501-1000", "1001-2000", "2001-5000"]
    probes["distance_bin_bp"] = pd.cut(
        probes["distance_to_anchor_bp"], bins=bins, labels=labels
    )
    decay_rows: list[dict[str, object]] = []
    for (gene, distance_bin), group in probes.groupby(
        ["positive_control_gene", "distance_bin_bp"], observed=True, sort=True
    ):
        decay_rows.append(
            {
                "positive_control_gene": gene,
                "distance_bin_bp": str(distance_bin),
                "n_probes": int(group["cpg_id"].nunique()),
                "median_abs_tcell_delta": float(group["tcell_delta"].abs().median()),
                "median_abs_msc_delta": float(group["msc_delta"].abs().median()),
                "mean_abs_editor_delta": float(group["mean_abs_editor_delta"].mean()),
            }
        )
    return probes, pd.DataFrame(gene_rows), pd.DataFrame(decay_rows)


def build_unique_regions(
    candidate_sets: Mapping[str, set[str]],
    candidate_metadata: pd.DataFrame,
    annotation: pd.DataFrame,
    tcell_effects: pd.DataFrame,
    msc_effects: pd.DataFrame,
    *,
    cluster_gap_bp: int = 2000,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if annotation.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    union = set().union(*candidate_sets.values())
    loci = pd.DataFrame({"cpg_id": sorted(union)}).merge(annotation, on="cpg_id", how="left")
    loci = loci.dropna(subset=["chromosome", "position"]).copy()
    for name, cpgs in candidate_sets.items():
        loci[f"in_{name}"] = loci["cpg_id"].isin(cpgs)
    loci = loci.merge(
        candidate_metadata[[c for c in candidate_metadata.columns if c not in loci.columns or c == "cpg_id"]],
        on="cpg_id",
        how="left",
        suffixes=("", "_meta"),
    )
    region_parts: list[pd.DataFrame] = []
    for chrom, group in loci.groupby("chromosome", sort=True):
        group = group.sort_values(["position", "cpg_id"]).copy()
        cluster = []
        current = 0
        previous: int | None = None
        for pos in group["position"].astype(int):
            if previous is None or pos - previous > cluster_gap_bp:
                current += 1
            cluster.append(current)
            previous = pos
        group["cluster_index"] = cluster
        region_parts.append(group)
    membership = pd.concat(region_parts, ignore_index=True)
    membership["unique_region_id"] = membership.apply(
        lambda r: f"chr{r['chromosome']}:{int(r['position'])}:{int(r['cluster_index']):04d}", axis=1
    )
    # Replace provisional IDs with coordinate-based stable IDs.
    coords = membership.groupby(["chromosome", "cluster_index"])["position"].agg(["min", "max"]).reset_index()
    coords["unique_region_id"] = coords.apply(
        lambda r: f"chr{r['chromosome']}:{int(r['min'])}-{int(r['max'])}", axis=1
    )
    membership = membership.drop(columns="unique_region_id").merge(
        coords[["chromosome", "cluster_index", "unique_region_id"]],
        on=["chromosome", "cluster_index"],
        how="left",
    )

    t = tcell_effects[["cpg_id", "editor_delta_beta", "sign_consistency"]].rename(
        columns={"editor_delta_beta": "tcell_delta", "sign_consistency": "tcell_consistency"}
    )
    m = msc_effects[["cpg_id", "editor_delta_beta", "sign_consistency"]].rename(
        columns={"editor_delta_beta": "msc_delta", "sign_consistency": "msc_consistency"}
    )
    membership = membership.merge(t, on="cpg_id", how="left").merge(m, on="cpg_id", how="left")

    rows: list[dict[str, object]] = []
    membership_cols = [f"in_{name}" for name in candidate_sets]
    for region_id, group in membership.groupby("unique_region_id", sort=True):
        td = group["tcell_delta"].to_numpy(dtype=float)
        md = group["msc_delta"].to_numpy(dtype=float)
        paired = np.isfinite(td) & np.isfinite(md)
        sign_conc = float(np.mean(np.sign(td[paired]) == np.sign(md[paired]))) if paired.any() else float("nan")
        matrix = np.column_stack([td, md])
        finite_abs = np.abs(matrix[np.isfinite(matrix)])
        mean_abs = float(np.mean(finite_abs)) if len(finite_abs) else float("nan")
        cons_matrix = group[["tcell_consistency", "msc_consistency"]].to_numpy(dtype=float)
        finite_cons = cons_matrix[np.isfinite(cons_matrix)]
        mean_cons = float(np.mean(finite_cons)) if len(finite_cons) else float("nan")
        measured_fraction = float(np.mean(np.isfinite(td) | np.isfinite(md)))
        factor = sign_conc if np.isfinite(sign_conc) else 0.5
        score = mean_abs * mean_cons * measured_fraction * factor if np.isfinite(mean_abs) and np.isfinite(mean_cons) else float("nan")
        row: dict[str, object] = {
            "unique_region_id": region_id,
            "chromosome": str(group["chromosome"].iloc[0]),
            "start": int(group["position"].min()),
            "end": int(group["position"].max()),
            "width_bp": int(group["position"].max() - group["position"].min() + 1),
            "n_candidate_cpgs": int(group["cpg_id"].nunique()),
            "region_class": "multi_cpg_region" if group["cpg_id"].nunique() >= 2 else "singleton_locus",
            "candidate_cpg_ids": ";".join(sorted(group["cpg_id"].unique())),
            "genes": ";".join(sorted({g for x in group["gene"].fillna("") for g in re.split(r"[;,]", str(x)) if g})),
            "measured_fraction": measured_fraction,
            "cross_cell_sign_concordance": sign_conc,
            "mean_abs_editor_delta": mean_abs,
            "mean_sign_consistency": mean_cons,
            "regional_responsiveness_score": score,
            "interpretation": "unique_region_bystander_responsiveness_not_direct_editability",
        }
        for col in membership_cols:
            row[col] = bool(group[col].fillna(False).any())
        rows.append(row)
    scores = pd.DataFrame(rows).sort_values("regional_responsiveness_score", ascending=False)
    return membership, scores, scores[scores["n_candidate_cpgs"] >= 2].copy()


def _robust_standardised_difference(measured: pd.Series, missing: pd.Series) -> float:
    a = pd.to_numeric(measured, errors="coerce").dropna().to_numpy(dtype=float)
    b = pd.to_numeric(missing, errors="coerce").dropna().to_numpy(dtype=float)
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    pooled = np.median(np.abs(np.concatenate([a, b]) - np.median(np.concatenate([a, b]))))
    return float((np.median(a) - np.median(b)) / (1.4826 * pooled)) if pooled > 0 else float("nan")


def build_coverage_bias_analysis(
    candidate_sets: Mapping[str, set[str]],
    candidate_metadata: pd.DataFrame,
    cell_effects: Mapping[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[pd.DataFrame] = []
    summary: list[dict[str, object]] = []
    metrics = ["abs_target_delta", "rank_score", "n_lodo_rotations"]
    meta = candidate_metadata.copy()
    if "target_delta_beta" in meta.columns:
        meta["abs_target_delta"] = pd.to_numeric(meta["target_delta_beta"], errors="coerce").abs()
    else:
        meta["abs_target_delta"] = np.nan
    for set_name, cpgs in candidate_sets.items():
        subset = meta[meta["cpg_id"].isin(cpgs)].copy()
        for cell_type, effects in cell_effects.items():
            measured_set = set(effects.loc[effects["editor_delta_beta"].notna(), "cpg_id"].astype(str))
            frame = subset.copy()
            frame["candidate_set"] = set_name
            frame["cell_type"] = cell_type
            frame["measured"] = frame["cpg_id"].isin(measured_set)
            rows.append(frame)
            for metric in metrics:
                if metric not in frame.columns:
                    continue
                a = pd.to_numeric(frame.loc[frame["measured"], metric], errors="coerce").dropna()
                b = pd.to_numeric(frame.loc[~frame["measured"], metric], errors="coerce").dropna()
                p = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue) if len(a) and len(b) else float("nan")
                summary.append(
                    {
                        "candidate_set": set_name,
                        "cell_type": cell_type,
                        "metric": metric,
                        "n_measured_with_metric": int(len(a)),
                        "n_missing_with_metric": int(len(b)),
                        "measured_median": float(a.median()) if len(a) else float("nan"),
                        "missing_median": float(b.median()) if len(b) else float("nan"),
                        "median_difference": float(a.median() - b.median()) if len(a) and len(b) else float("nan"),
                        "robust_standardised_difference": _robust_standardised_difference(a, b),
                        "mannwhitney_p": p,
                    }
                )
    per_candidate = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    summary_df = pd.DataFrame(summary)
    if not summary_df.empty:
        valid = summary_df["mannwhitney_p"].fillna(1.0).to_numpy(dtype=float)
        summary_df["bh_q_global"] = multipletests(valid, method="fdr_bh")[1]
    return per_candidate, summary_df


def summarise_matched_nulls(null_tests: pd.DataFrame) -> pd.DataFrame:
    if null_tests.empty:
        return pd.DataFrame()
    summary = (
        null_tests.groupby(["candidate_set", "cell_type", "metric"], as_index=False)
        .agg(
            observed=("observed", "first"),
            null_mean=("null_value", "mean"),
            null_median=("null_value", "median"),
            null_q95=("null_value", lambda x: x.quantile(0.95)),
            empirical_p=("empirical_p", "first"),
            n_candidates_measured=("n_candidates_measured", "first"),
        )
    )
    summary["bh_q_global"] = multipletests(summary["empirical_p"].fillna(1.0), method="fdr_bh")[1]
    summary["bh_q_within_celltype"] = np.nan
    for cell, idx in summary.groupby("cell_type").groups.items():
        summary.loc[idx, "bh_q_within_celltype"] = multipletests(
            summary.loc[idx, "empirical_p"].fillna(1.0), method="fdr_bh"
        )[1]
    return summary


def _quantile_codes(series: pd.Series, bins: int) -> pd.Series:
    valid = pd.to_numeric(series, errors="coerce")
    try:
        return pd.qcut(valid.rank(method="first"), q=min(bins, valid.notna().sum()), labels=False, duplicates="drop")
    except ValueError:
        return pd.Series(0, index=series.index, dtype=float)


def run_cross_celltype_empirical_null(
    candidate_sets: Mapping[str, set[str]],
    target_map: pd.DataFrame,
    discovery_all_loci: pd.DataFrame | None,
    tcell_effects: pd.DataFrame,
    msc_effects: pd.DataFrame,
    *,
    n_iterations: int = 1000,
    n_bins: int = 10,
    seed: int = 20260721,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    t = tcell_effects[["cpg_id", "baseline_beta", "editor_delta_beta"]].rename(
        columns={"baseline_beta": "t_baseline", "editor_delta_beta": "t_delta"}
    )
    m = msc_effects[["cpg_id", "baseline_beta", "editor_delta_beta"]].rename(
        columns={"baseline_beta": "m_baseline", "editor_delta_beta": "m_delta"}
    )
    pool = t.merge(m, on="cpg_id", how="inner")
    pool = pool.dropna(subset=["t_delta", "m_delta"]).copy()
    pool = pool.merge(target_map[["cpg_id", "target_delta_beta"]], on="cpg_id", how="left")
    if discovery_all_loci is not None and not discovery_all_loci.empty:
        disc = _normalise_cpg_column(discovery_all_loci)
        delta_col = next((c for c in ["target_delta_beta", "mean_delta", "discovery_delta_beta"] if c in disc.columns), None)
        if delta_col:
            disc = disc[["cpg_id", delta_col]].rename(columns={delta_col: "background_discovery_delta"})
            pool = pool.merge(disc, on="cpg_id", how="left")
    if "background_discovery_delta" not in pool.columns:
        pool["background_discovery_delta"] = pool["target_delta_beta"]
    pool["mean_baseline"] = pool[["t_baseline", "m_baseline"]].mean(axis=1)
    pool["abs_discovery_delta"] = pd.to_numeric(pool["background_discovery_delta"], errors="coerce").abs()
    pool["baseline_bin"] = _quantile_codes(pool["mean_baseline"], n_bins)
    pool["effect_bin"] = _quantile_codes(pool["abs_discovery_delta"], n_bins)
    pool["stratum"] = pool["baseline_bin"].astype(str) + "|" + pool["effect_bin"].astype(str)

    rows: list[dict[str, object]] = []
    metrics = ["sign_concordance", "abs_spearman", "mean_abs_tcell", "mean_abs_msc"]
    for set_name, cpgs in candidate_sets.items():
        observed = pool[pool["cpg_id"].isin(cpgs)].copy()
        if len(observed) < 3:
            continue
        rho, _ = _safe_spearman(observed["t_delta"].to_numpy(float), observed["m_delta"].to_numpy(float))
        obs = {
            "sign_concordance": float(np.mean(np.sign(observed["t_delta"]) == np.sign(observed["m_delta"]))),
            "abs_spearman": abs(rho) if np.isfinite(rho) else float("nan"),
            "mean_abs_tcell": float(observed["t_delta"].abs().mean()),
            "mean_abs_msc": float(observed["m_delta"].abs().mean()),
        }
        stratum_counts = observed["stratum"].value_counts().to_dict()
        null_values = {metric: np.full(n_iterations, np.nan) for metric in metrics}
        eligible = pool[~pool["cpg_id"].isin(cpgs)].copy()
        for i in range(n_iterations):
            parts: list[pd.DataFrame] = []
            for stratum, count in stratum_counts.items():
                choices = eligible[eligible["stratum"] == stratum]
                if len(choices) < count:
                    choices = eligible
                if choices.empty:
                    continue
                take = choices.iloc[rng.choice(len(choices), size=count, replace=len(choices) < count)]
                parts.append(take)
            if not parts:
                continue
            sampled = pd.concat(parts, ignore_index=True).head(len(observed))
            srho, _ = _safe_spearman(sampled["t_delta"].to_numpy(float), sampled["m_delta"].to_numpy(float))
            null_values["sign_concordance"][i] = np.mean(np.sign(sampled["t_delta"]) == np.sign(sampled["m_delta"]))
            null_values["abs_spearman"][i] = abs(srho) if np.isfinite(srho) else np.nan
            null_values["mean_abs_tcell"][i] = sampled["t_delta"].abs().mean()
            null_values["mean_abs_msc"][i] = sampled["m_delta"].abs().mean()
        for metric in metrics:
            p = _empirical_p(obs[metric], null_values[metric], alternative="greater")
            for i, value in enumerate(null_values[metric], start=1):
                rows.append(
                    {
                        "candidate_set": set_name,
                        "metric": metric,
                        "iteration": i,
                        "observed": obs[metric],
                        "null_value": value,
                        "empirical_p": p,
                        "matching": "common_measured_baseline_and_discovery_effect_bins",
                        "n_common_measured": len(observed),
                    }
                )
    tests = pd.DataFrame(rows)
    if tests.empty:
        return tests, pd.DataFrame()
    summary = (
        tests.groupby(["candidate_set", "metric"], as_index=False)
        .agg(
            observed=("observed", "first"),
            null_mean=("null_value", "mean"),
            null_q95=("null_value", lambda x: x.quantile(0.95)),
            empirical_p=("empirical_p", "first"),
            n_common_measured=("n_common_measured", "first"),
        )
    )
    summary["bh_q_global"] = multipletests(summary["empirical_p"].fillna(1.0), method="fdr_bh")[1]
    return tests, summary


def run_bystander_validation_v071(
    *,
    tcell_raw_paths: Sequence[str | Path],
    msc_beta: pd.DataFrame,
    msc_metadata: pd.DataFrame,
    candidate_sets: Mapping[str, set[str]],
    target_map: pd.DataFrame,
    candidate_metadata: pd.DataFrame,
    discovery_all_loci: pd.DataFrame | None,
    probe_annotation: pd.DataFrame,
    min_abs_delta: float = 0.02,
    fdr: float = 0.10,
    null_iterations: int = 1000,
    null_bins: int = 10,
    cluster_gap_bp: int = 2000,
    seed: int = 20260721,
) -> BystanderValidationV071Result:
    tcell = build_platform_aware_tcell_effects(
        tcell_raw_paths, min_abs_delta=min_abs_delta, fdr=fdr
    )
    msc = build_paired_celltype_effects(
        msc_beta, msc_metadata, cell_type="mesenchymal_stromal_cell",
        min_abs_delta=min_abs_delta, fdr=fdr
    )
    cell_effects = {
        "t_cell": tcell.locus_effects,
        "mesenchymal_stromal_cell": msc.locus_effects,
    }
    coverage, candidate_effects, candidate_summary = build_coverage_and_candidate_effects(
        candidate_sets, target_map, cell_effects
    )
    cross_effects, cross_summary, _ = build_cross_celltype_results(
        candidate_sets, target_map, tcell.locus_effects, msc.locus_effects
    )

    matched_nulls = run_matched_null_tests(
        candidate_sets, cell_effects, discovery_all_loci,
        n_iterations=null_iterations, n_bins=null_bins, seed=seed
    )
    matched_summary = summarise_matched_nulls(matched_nulls)
    cross_nulls, cross_null_summary = run_cross_celltype_empirical_null(
        candidate_sets, target_map, discovery_all_loci,
        tcell.locus_effects, msc.locus_effects,
        n_iterations=null_iterations, n_bins=null_bins, seed=seed + 101
    )
    positive_probes, positive_genes, positive_decay = build_positive_control_analysis(
        probe_annotation, tcell.locus_effects, msc.locus_effects
    )
    region_membership, regions_all, regions_multi = build_unique_regions(
        candidate_sets, candidate_metadata, probe_annotation,
        tcell.locus_effects, msc.locus_effects,
        cluster_gap_bp=cluster_gap_bp
    )
    coverage_bias, coverage_bias_summary = build_coverage_bias_analysis(
        candidate_sets, candidate_metadata, cell_effects
    )

    significant_nulls = matched_summary[
        matched_summary["bh_q_global"] <= 0.10
    ].to_dict("records") if not matched_summary.empty else []
    summary = {
        "analysis_status": "positive_control_platform_aware_unique_region_validation",
        "interpretation_scope": (
            "T-cell R1/R2 and R3 are analysed on their native probe spaces and combined by "
            "available-pair empirical-variance meta-analysis. Positive-control distances use an "
            "annotation-derived gene anchor, not the unpublished guide cut site. Unique-region scores "
            "are bystander-responsiveness priorities, not causal or direct-editability proof."
        ),
        "n_tcell_platforms": int(len(tcell.platform_qc)),
        "n_tcell_pairs": int(len(tcell.pair_manifest)),
        "n_msc_pairs": int(len(msc.pair_manifest)),
        "candidate_set_sizes": {k: len(v) for k, v in candidate_sets.items()},
        "coverage": coverage.to_dict("records"),
        "candidate_set_summary": candidate_summary.to_dict("records"),
        "cross_celltype_summary": cross_summary.to_dict("records"),
        "n_positive_control_probes": int(len(positive_probes)),
        "positive_control_genes_found": sorted(positive_genes["positive_control_gene"].unique().tolist()) if not positive_genes.empty else [],
        "n_unique_regions_all": int(len(regions_all)),
        "n_unique_multicpg_regions": int(len(regions_multi)),
        "n_significant_matched_null_results_q_le_0_10": int(len(significant_nulls)),
        "significant_matched_null_results_q_le_0_10": significant_nulls,
        "min_abs_delta": min_abs_delta,
        "fdr": fdr,
        "null_iterations": null_iterations,
        "null_bins": null_bins,
        "cluster_gap_bp": cluster_gap_bp,
        "seed": seed,
    }
    manifest = {
        "candidate_lists_frozen_before_bystander_analysis": True,
        "tcell_platform_policy": "native_probe_space_available_pair_meta_analysis",
        "tcell_meta_weighting": "inverse_empirical_pair_variance_with_derSimonian_laird_sensitivity",
        "positive_control_anchor": "median_coordinate_of_gene_annotated_EPICv2_probes",
        "positive_control_guardrail": "distance_is_not_distance_to_guide_cut_site",
        "unique_region_policy": "cluster_union_once_then_attach_candidate_set_memberships",
        "main_region_policy": "multi_cpg_regions_n_at_least_2",
        "multiple_testing": ["BH_within_celltype", "BH_global"],
        "cross_celltype_null": "common_measured_baseline_and_discovery_effect_matched_random_sets",
        "editor": "dCas9-DNMT3A/3L",
    }
    return BystanderValidationV071Result(
        candidate_coverage=coverage,
        candidate_set_summary=candidate_summary,
        candidate_effects=candidate_effects,
        cross_celltype_effects=cross_effects,
        cross_celltype_summary=cross_summary,
        tcell_platform_meta_effects=tcell.locus_effects,
        tcell_platform_pair_effects=tcell.pair_effects,
        tcell_platform_manifest=tcell.pair_manifest,
        tcell_platform_qc=tcell.platform_qc,
        msc_effects=msc.locus_effects,
        positive_control_probe_effects=positive_probes,
        positive_control_gene_summary=positive_genes,
        positive_control_distance_decay=positive_decay,
        unique_region_membership=region_membership,
        unique_region_scores_all=regions_all,
        unique_region_scores_multicpg=regions_multi,
        coverage_bias_per_candidate=coverage_bias,
        coverage_bias_summary=coverage_bias_summary,
        matched_null_tests=matched_nulls,
        matched_null_summary=matched_summary,
        cross_celltype_null_tests=cross_nulls,
        cross_celltype_null_summary=cross_null_summary,
        summary=_jsonable(summary),
        manifest=_jsonable(manifest),
    )


def write_bystander_validation_v071_outputs(
    result: BystanderValidationV071Result,
    output_dir: str | Path,
) -> None:
    out = Path(output_dir)
    figures = out / "figures"
    out.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    csvs = {
        "candidate_probe_coverage_platform_aware.csv": result.candidate_coverage,
        "candidate_set_summary_platform_aware.csv": result.candidate_set_summary,
        "cross_celltype_replication_platform_aware.csv": result.cross_celltype_summary,
        "tcell_platform_pair_manifest.csv": result.tcell_platform_manifest,
        "tcell_platform_qc.csv": result.tcell_platform_qc,
        "positive_control_gene_summary.csv": result.positive_control_gene_summary,
        "positive_control_distance_decay.csv": result.positive_control_distance_decay,
        "unique_region_scores_all.csv": result.unique_region_scores_all,
        "unique_region_scores_multicpg.csv": result.unique_region_scores_multicpg,
        "coverage_bias_summary.csv": result.coverage_bias_summary,
        "matched_null_summary_bh.csv": result.matched_null_summary,
        "cross_celltype_null_summary.csv": result.cross_celltype_null_summary,
    }
    for name, frame in csvs.items():
        frame.to_csv(out / name, index=False)

    gz = {
        "candidate_bystander_effects_platform_aware.csv.gz": result.candidate_effects,
        "cross_celltype_candidate_effects_platform_aware.csv.gz": result.cross_celltype_effects,
        "tcell_platform_meta_effects.csv.gz": result.tcell_platform_meta_effects,
        "tcell_platform_pair_effects.csv.gz": result.tcell_platform_pair_effects,
        "msc_effects_v071.csv.gz": result.msc_effects,
        "positive_control_probe_effects.csv.gz": result.positive_control_probe_effects,
        "unique_region_membership.csv.gz": result.unique_region_membership,
        "coverage_bias_per_candidate.csv.gz": result.coverage_bias_per_candidate,
        "matched_null_tests_v071.csv.gz": result.matched_null_tests,
        "cross_celltype_null_tests.csv.gz": result.cross_celltype_null_tests,
    }
    for name, frame in gz.items():
        frame.to_csv(out / name, index=False, compression="gzip")

    (out / "bystander_validation_v071_summary.json").write_text(
        json.dumps(result.summary, indent=2), encoding="utf-8"
    )
    (out / "bystander_validation_v071_manifest.json").write_text(
        json.dumps(result.manifest, indent=2), encoding="utf-8"
    )
    _write_figures(result, figures)


def _write_figures(result: BystanderValidationV071Result, figures: Path) -> None:
    import matplotlib.pyplot as plt

    if not result.positive_control_distance_decay.empty:
        for cell_col, label in [
            ("median_abs_tcell_delta", "T-cell"),
            ("median_abs_msc_delta", "MSC"),
        ]:
            pivot = result.positive_control_distance_decay.pivot(
                index="distance_bin_bp", columns="positive_control_gene", values=cell_col
            )
            ax = pivot.plot(marker="o", figsize=(9, 5))
            ax.set_ylabel(f"Median |Δβ| ({label})")
            ax.set_xlabel("Distance from annotation-derived gene anchor (bp bin)")
            ax.set_title(f"Published-target positive controls: {label}")
            plt.tight_layout()
            plt.savefig(figures / f"positive_control_decay_{label.lower().replace('-', '')}.png", dpi=200)
            plt.close()

    if not result.matched_null_summary.empty:
        subset = result.matched_null_summary[
            result.matched_null_summary["metric"].isin(["mean_abs_delta", "median_abs_delta", "mean_sign_consistency"])
        ].copy()
        subset["label"] = subset["candidate_set"] + ":" + subset["cell_type"] + ":" + subset["metric"]
        subset = subset.sort_values("bh_q_global").head(20)
        plt.figure(figsize=(10, 7))
        plt.barh(subset["label"][::-1], -np.log10(subset["bh_q_global"].clip(lower=1e-12))[::-1])
        plt.axvline(-math.log10(0.10), linestyle="--")
        plt.xlabel("-log10(BH q, global)")
        plt.title("Matched-null evidence after multiplicity correction")
        plt.tight_layout()
        plt.savefig(figures / "matched_null_bh_summary.png", dpi=200)
        plt.close()

    if not result.unique_region_scores_multicpg.empty:
        top = result.unique_region_scores_multicpg.head(20)
        plt.figure(figsize=(9, 6))
        plt.barh(top["unique_region_id"][::-1], top["regional_responsiveness_score"][::-1])
        plt.xlabel("Regional responsiveness score")
        plt.title("Unique multi-CpG candidate regions")
        plt.tight_layout()
        plt.savefig(figures / "unique_multicpg_regions.png", dpi=200)
        plt.close()

    if not result.coverage_bias_summary.empty:
        plot = result.coverage_bias_summary.dropna(subset=["robust_standardised_difference"]).copy()
        plot["label"] = plot["candidate_set"] + ":" + plot["cell_type"] + ":" + plot["metric"]
        plot = plot.head(30)
        plt.figure(figsize=(10, 7))
        plt.barh(plot["label"][::-1], plot["robust_standardised_difference"][::-1])
        plt.axvline(0, linewidth=1)
        plt.xlabel("Robust standardised difference: measured - missing")
        plt.title("Candidate coverage-bias diagnostics")
        plt.tight_layout()
        plt.savefig(figures / "coverage_bias.png", dpi=200)
        plt.close()

    if not result.cross_celltype_null_summary.empty:
        subset = result.cross_celltype_null_summary[
            result.cross_celltype_null_summary["metric"] == "sign_concordance"
        ]
        x = np.arange(len(subset))
        plt.figure(figsize=(9, 5))
        plt.bar(x - 0.18, subset["observed"], width=0.36, label="Observed")
        plt.bar(x + 0.18, subset["null_mean"], width=0.36, label="Matched null mean")
        plt.xticks(x, subset["candidate_set"], rotation=30, ha="right")
        plt.ylabel("Cross-cell sign concordance")
        plt.title("Cross-cell-type empirical null")
        plt.legend()
        plt.tight_layout()
        plt.savefig(figures / "cross_celltype_empirical_null.png", dpi=200)
        plt.close()
