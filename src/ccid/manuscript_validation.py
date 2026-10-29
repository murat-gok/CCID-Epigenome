from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import json
import math

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


DEFAULT_COMPONENT_WEIGHTS: dict[str, float] = {
    "discovery_strength_score": 0.10,
    "donor_stability_score": 0.18,
    "o3_support_score": 0.14,
    "successful_failed_specificity_score": 0.12,
    "temporal_consistency_score": 0.12,
    "editor_responsiveness_score": 0.10,
    "regional_support_score": 0.08,
    "genomic_context_score": 0.06,
    "regulatory_evidence_score": 0.10,
}

DEFAULT_V082_CONFIG: dict[str, object] = {
    "component_weights": DEFAULT_COMPONENT_WEIGHTS,
    "top_k": 20,
    "secondary_top_k": 50,
    "weight_iterations": 1000,
    "weight_log_sigma": 0.30,
    "uncertainty_iterations": 1000,
    "component_noise_sd": 0.05,
    "null_iterations": 1000,
    "seed": 20260722,
    "stable_topk_frequency": 0.70,
    "shortlist_size": 50,
}


@dataclass(frozen=True)
class ManuscriptValidationResult:
    candidate_stability: pd.DataFrame
    weight_iteration_summary: pd.DataFrame
    uncertainty_iteration_summary: pd.DataFrame
    ablation_summary: pd.DataFrame
    ablation_candidate_shifts: pd.DataFrame
    null_summary: pd.DataFrame
    null_iterations: pd.DataFrame
    component_correlation: pd.DataFrame
    redundancy_summary: pd.DataFrame
    candidate_set_enrichment: pd.DataFrame
    robust_shortlist: pd.DataFrame
    summary: dict[str, object]
    manifest: dict[str, object]


def _jsonable(value: object) -> object:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def _normalise_weights(weights: Mapping[str, float], available: Sequence[str]) -> dict[str, float]:
    clean = {c: max(float(weights.get(c, 0.0)), 0.0) for c in available}
    total = sum(clean.values())
    if total <= 0:
        return {c: 1.0 / max(len(available), 1) for c in available}
    return {c: w / total for c, w in clean.items()}


def _score_available(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    valid = np.isfinite(values)
    numerator = np.nansum(np.where(valid, values * weights[None, :], 0.0), axis=1)
    denominator = np.sum(valid * weights[None, :], axis=1)
    return np.divide(
        numerator,
        denominator,
        out=np.full(values.shape[0], np.nan, dtype=float),
        where=denominator > 0,
    )


def _rank_desc(scores: np.ndarray) -> np.ndarray:
    return pd.Series(scores).rank(method="average", ascending=False, na_option="bottom").to_numpy(dtype=float)


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    valid = np.isfinite(a) & np.isfinite(b)
    if valid.sum() < 3:
        return float("nan")
    result = stats.spearmanr(a[valid], b[valid])
    return float(result.statistic if hasattr(result, "statistic") else result[0])


def _top_indices(scores: np.ndarray, k: int) -> np.ndarray:
    finite = np.where(np.isfinite(scores))[0]
    if len(finite) == 0:
        return np.array([], dtype=int)
    order = finite[np.argsort(-scores[finite], kind="mergesort")]
    return order[: min(k, len(order))]


def _jaccard(a: Sequence[int], b: Sequence[int]) -> float:
    sa, sb = set(map(int, a)), set(map(int, b))
    if not sa and not sb:
        return 1.0
    return len(sa & sb) / max(len(sa | sb), 1)


def _prepare_candidates(cpg: pd.DataFrame, component_weights: Mapping[str, float]) -> tuple[pd.DataFrame, list[str], np.ndarray, np.ndarray]:
    if "cpg_id" not in cpg.columns:
        raise ValueError("Candidate table must contain cpg_id")
    frame = cpg.copy().drop_duplicates("cpg_id", keep="first").reset_index(drop=True)
    components = [c for c in component_weights if c in frame.columns and pd.to_numeric(frame[c], errors="coerce").notna().any()]
    if len(components) < 2:
        raise ValueError("At least two evidence component columns are required")
    for col in components:
        frame[col] = pd.to_numeric(frame[col], errors="coerce").clip(0, 1)
    values = frame[components].to_numpy(dtype=float)
    weights_map = _normalise_weights(component_weights, components)
    weights = np.array([weights_map[c] for c in components], dtype=float)
    baseline = 100.0 * _score_available(values, weights)
    frame["v082_baseline_score"] = baseline
    frame["v082_baseline_rank"] = _rank_desc(baseline)
    frame["v082_component_availability"] = np.isfinite(values).mean(axis=1)
    return frame, components, values, weights


def _perturb_weights(
    values: np.ndarray,
    base_weights: np.ndarray,
    baseline_scores: np.ndarray,
    *,
    iterations: int,
    sigma: float,
    top_k: int,
    secondary_top_k: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    n = values.shape[0]
    scores = np.empty((iterations, n), dtype=np.float32)
    ranks = np.empty((iterations, n), dtype=np.float32)
    baseline_top = _top_indices(baseline_scores, top_k)
    baseline_top2 = _top_indices(baseline_scores, secondary_top_k)
    rows: list[dict[str, float | int]] = []
    safe = np.clip(base_weights, 1e-12, None)
    for i in range(iterations):
        sampled = rng.lognormal(mean=np.log(safe), sigma=sigma)
        sampled = sampled / sampled.sum()
        s = 100.0 * _score_available(values, sampled)
        r = _rank_desc(s)
        scores[i] = s
        ranks[i] = r
        top = _top_indices(s, top_k)
        top2 = _top_indices(s, secondary_top_k)
        rows.append(
            {
                "iteration": i + 1,
                "spearman_vs_baseline": _spearman(baseline_scores, s),
                "top_k_jaccard": _jaccard(baseline_top, top),
                "secondary_top_k_jaccard": _jaccard(baseline_top2, top2),
                "maximum_score": float(np.nanmax(s)),
                "mean_top_k_score": float(np.nanmean(s[top])) if len(top) else float("nan"),
            }
        )
    return scores, ranks, pd.DataFrame(rows)


def _perturb_components(
    values: np.ndarray,
    weights: np.ndarray,
    baseline_scores: np.ndarray,
    *,
    iterations: int,
    noise_sd: float,
    top_k: int,
    secondary_top_k: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    n = values.shape[0]
    scores = np.empty((iterations, n), dtype=np.float32)
    ranks = np.empty((iterations, n), dtype=np.float32)
    baseline_top = _top_indices(baseline_scores, top_k)
    baseline_top2 = _top_indices(baseline_scores, secondary_top_k)
    rows: list[dict[str, float | int]] = []
    valid = np.isfinite(values)
    for i in range(iterations):
        jitter = rng.normal(0.0, noise_sd, size=values.shape)
        perturbed = np.where(valid, np.clip(values + jitter, 0.0, 1.0), np.nan)
        s = 100.0 * _score_available(perturbed, weights)
        r = _rank_desc(s)
        scores[i] = s
        ranks[i] = r
        top = _top_indices(s, top_k)
        top2 = _top_indices(s, secondary_top_k)
        rows.append(
            {
                "iteration": i + 1,
                "spearman_vs_baseline": _spearman(baseline_scores, s),
                "top_k_jaccard": _jaccard(baseline_top, top),
                "secondary_top_k_jaccard": _jaccard(baseline_top2, top2),
                "maximum_score": float(np.nanmax(s)),
                "mean_top_k_score": float(np.nanmean(s[top])) if len(top) else float("nan"),
            }
        )
    return scores, ranks, pd.DataFrame(rows)


def _summarise_stability(
    frame: pd.DataFrame,
    weight_scores: np.ndarray,
    weight_ranks: np.ndarray,
    uncertainty_scores: np.ndarray,
    uncertainty_ranks: np.ndarray,
    *,
    top_k: int,
    secondary_top_k: int,
) -> pd.DataFrame:
    out = frame.copy()
    for prefix, score_mat, rank_mat in [
        ("weight", weight_scores, weight_ranks),
        ("uncertainty", uncertainty_scores, uncertainty_ranks),
    ]:
        out[f"{prefix}_score_mean"] = np.nanmean(score_mat, axis=0)
        out[f"{prefix}_score_sd"] = np.nanstd(score_mat, axis=0, ddof=1)
        out[f"{prefix}_rank_median"] = np.nanmedian(rank_mat, axis=0)
        out[f"{prefix}_rank_q25"] = np.nanquantile(rank_mat, 0.25, axis=0)
        out[f"{prefix}_rank_q75"] = np.nanquantile(rank_mat, 0.75, axis=0)
        out[f"{prefix}_rank_iqr"] = out[f"{prefix}_rank_q75"] - out[f"{prefix}_rank_q25"]
        out[f"{prefix}_top_k_frequency"] = np.mean(rank_mat <= top_k, axis=0)
        out[f"{prefix}_secondary_top_k_frequency"] = np.mean(rank_mat <= secondary_top_k, axis=0)
    n = max(len(out), 1)
    baseline_percentile = 1.0 - (out["v082_baseline_rank"] - 1.0) / max(n - 1, 1)
    rank_iqr_scaled = np.clip(
        (out["weight_rank_iqr"] + out["uncertainty_rank_iqr"]) / (2.0 * max(n - 1, 1)),
        0,
        1,
    )
    out["robustness_composite_score"] = 100.0 * (
        0.35 * baseline_percentile
        + 0.25 * out["weight_top_k_frequency"]
        + 0.25 * out["uncertainty_top_k_frequency"]
        + 0.15 * (1.0 - rank_iqr_scaled)
    )
    out["robustness_rank"] = out["robustness_composite_score"].rank(method="min", ascending=False).astype(int)
    return out


def _ablation(
    frame: pd.DataFrame,
    components: Sequence[str],
    values: np.ndarray,
    weights: np.ndarray,
    baseline_scores: np.ndarray,
    *,
    top_k: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    baseline_ranks = _rank_desc(baseline_scores)
    baseline_top = _top_indices(baseline_scores, top_k)
    summary_rows: list[dict[str, object]] = []
    shift_rows: list[pd.DataFrame] = []
    for j, component in enumerate(components):
        keep = np.ones(len(components), dtype=bool)
        keep[j] = False
        new_weights = weights[keep]
        new_weights = new_weights / new_weights.sum()
        scores = 100.0 * _score_available(values[:, keep], new_weights)
        ranks = _rank_desc(scores)
        top = _top_indices(scores, top_k)
        shift = ranks - baseline_ranks
        summary_rows.append(
            {
                "removed_component": component,
                "spearman_vs_baseline": _spearman(baseline_scores, scores),
                "top_k_jaccard": _jaccard(baseline_top, top),
                "median_absolute_rank_change": float(np.nanmedian(np.abs(shift))),
                "mean_absolute_rank_change": float(np.nanmean(np.abs(shift))),
                "maximum_absolute_rank_change": float(np.nanmax(np.abs(shift))),
            }
        )
        part = frame[["cpg_id", "v082_baseline_rank"]].copy()
        part["removed_component"] = component
        part["ablated_score"] = scores
        part["ablated_rank"] = ranks
        part["rank_change"] = shift
        shift_rows.append(part)
    return pd.DataFrame(summary_rows), pd.concat(shift_rows, ignore_index=True)


def _ranking_null(
    values: np.ndarray,
    weights: np.ndarray,
    baseline_scores: np.ndarray,
    *,
    iterations: int,
    top_k: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    observed_top = _top_indices(baseline_scores, top_k)
    observed_mean_top = float(np.nanmean(baseline_scores[observed_top]))
    observed_max = float(np.nanmax(baseline_scores))
    rows: list[dict[str, float | int]] = []
    for i in range(iterations):
        permuted = values.copy()
        for j in range(values.shape[1]):
            idx = np.where(np.isfinite(permuted[:, j]))[0]
            if len(idx) > 1:
                permuted[idx, j] = rng.permutation(permuted[idx, j])
        scores = 100.0 * _score_available(permuted, weights)
        top = _top_indices(scores, top_k)
        rows.append(
            {
                "iteration": i + 1,
                "null_mean_top_k_score": float(np.nanmean(scores[top])),
                "null_maximum_score": float(np.nanmax(scores)),
            }
        )
    it = pd.DataFrame(rows)
    p_top = (1.0 + float((it["null_mean_top_k_score"] >= observed_mean_top).sum())) / (iterations + 1.0)
    p_max = (1.0 + float((it["null_maximum_score"] >= observed_max).sum())) / (iterations + 1.0)
    summary = pd.DataFrame(
        [
            {
                "statistic": "mean_top_k_score",
                "observed": observed_mean_top,
                "null_mean": float(it["null_mean_top_k_score"].mean()),
                "null_q95": float(it["null_mean_top_k_score"].quantile(0.95)),
                "empirical_p": p_top,
            },
            {
                "statistic": "maximum_score",
                "observed": observed_max,
                "null_mean": float(it["null_maximum_score"].mean()),
                "null_q95": float(it["null_maximum_score"].quantile(0.95)),
                "empirical_p": p_max,
            },
        ]
    )
    return summary, it


def _component_redundancy(frame: pd.DataFrame, components: Sequence[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    corr = frame[list(components)].corr(method="spearman", min_periods=5)
    matrix = corr.to_numpy(dtype=float)
    safe = np.nan_to_num(matrix, nan=0.0)
    np.fill_diagonal(safe, 1.0)
    eig = np.linalg.eigvalsh((safe + safe.T) / 2.0)
    eig = np.clip(eig, 0.0, None)
    total = eig.sum()
    if total > 0:
        p = eig[eig > 0] / total
        effective_rank = float(np.exp(-np.sum(p * np.log(p))))
    else:
        effective_rank = float("nan")
    off_diag = np.abs(safe[np.triu_indices_from(safe, k=1)])
    high_pairs: list[tuple[str, str, float]] = []
    for i in range(len(components)):
        for j in range(i + 1, len(components)):
            value = corr.iloc[i, j]
            if pd.notna(value) and abs(value) >= 0.80:
                high_pairs.append((components[i], components[j], float(value)))
    summary = pd.DataFrame(
        [
            {
                "n_components": len(components),
                "effective_rank": effective_rank,
                "maximum_absolute_pairwise_spearman": float(np.nanmax(off_diag)) if len(off_diag) else float("nan"),
                "n_pairs_abs_spearman_ge_0_80": len(high_pairs),
                "high_correlation_pairs": ";".join(f"{a}|{b}|{r:.3f}" for a, b, r in high_pairs),
            }
        ]
    )
    corr.index.name = "component"
    return corr.reset_index(), summary


def _fisher_enrichment(frame: pd.DataFrame, *, rank_col: str, top_values: Sequence[int]) -> pd.DataFrame:
    flags: dict[str, pd.Series] = {}
    for col in [
        "in_stable_core_3of3",
        "in_consensus_at_least_2of3",
        "in_o3_tier1",
        "in_frozen_500",
    ]:
        if col in frame.columns:
            flags[col] = frame[col].fillna(False).astype(bool)
    if "n_candidate_cpgs" in frame.columns:
        flags["in_multicpg_region"] = pd.to_numeric(frame["n_candidate_cpgs"], errors="coerce").fillna(1) >= 2
    rows: list[dict[str, object]] = []
    n = len(frame)
    for k in top_values:
        top = pd.to_numeric(frame[rank_col], errors="coerce") <= min(int(k), n)
        for name, flag in flags.items():
            a = int((top & flag).sum())
            b = int((top & ~flag).sum())
            c = int((~top & flag).sum())
            d = int((~top & ~flag).sum())
            odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append(
                {
                    "top_k": int(k),
                    "feature": name,
                    "top_with_feature": a,
                    "top_without_feature": b,
                    "background_with_feature": c,
                    "background_without_feature": d,
                    "top_fraction": a / max(a + b, 1),
                    "background_fraction": c / max(c + d, 1),
                    "odds_ratio": float(odds),
                    "p_value": float(p),
                }
            )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["bh_q_global"] = multipletests(result["p_value"], method="fdr_bh")[1]
    return result


def run_manuscript_validation(cpg: pd.DataFrame, config: Mapping[str, object] | None = None) -> ManuscriptValidationResult:
    cfg = dict(DEFAULT_V082_CONFIG)
    if config:
        cfg.update(config)
    weights_cfg = dict(DEFAULT_COMPONENT_WEIGHTS)
    weights_cfg.update(dict(cfg.get("component_weights", {})))
    frame, components, values, weights = _prepare_candidates(cpg, weights_cfg)
    baseline = frame["v082_baseline_score"].to_numpy(dtype=float)
    top_k = int(cfg["top_k"])
    secondary_top_k = int(cfg["secondary_top_k"])
    seed = int(cfg["seed"])
    root_rng = np.random.default_rng(seed)
    child_seeds = root_rng.integers(0, np.iinfo(np.uint32).max, size=3, dtype=np.uint32)

    w_scores, w_ranks, w_iter = _perturb_weights(
        values,
        weights,
        baseline,
        iterations=int(cfg["weight_iterations"]),
        sigma=float(cfg["weight_log_sigma"]),
        top_k=top_k,
        secondary_top_k=secondary_top_k,
        rng=np.random.default_rng(int(child_seeds[0])),
    )
    u_scores, u_ranks, u_iter = _perturb_components(
        values,
        weights,
        baseline,
        iterations=int(cfg["uncertainty_iterations"]),
        noise_sd=float(cfg["component_noise_sd"]),
        top_k=top_k,
        secondary_top_k=secondary_top_k,
        rng=np.random.default_rng(int(child_seeds[1])),
    )
    stability = _summarise_stability(
        frame,
        w_scores,
        w_ranks,
        u_scores,
        u_ranks,
        top_k=top_k,
        secondary_top_k=secondary_top_k,
    )
    ablation_summary, ablation_shifts = _ablation(
        frame,
        components,
        values,
        weights,
        baseline,
        top_k=top_k,
    )
    null_summary, null_iterations = _ranking_null(
        values,
        weights,
        baseline,
        iterations=int(cfg["null_iterations"]),
        top_k=top_k,
        rng=np.random.default_rng(int(child_seeds[2])),
    )
    corr, redundancy = _component_redundancy(frame, components)
    enrichment = _fisher_enrichment(stability, rank_col="robustness_rank", top_values=[top_k, secondary_top_k])

    worst_drop = (
        ablation_shifts.groupby("cpg_id", as_index=False)["rank_change"]
        .agg(worst_ablation_rank_drop="max", best_ablation_rank_gain="min", mean_absolute_ablation_rank_change=lambda x: float(np.mean(np.abs(x))))
    )
    stability = stability.merge(worst_drop, on="cpg_id", how="left")
    threshold = float(cfg["stable_topk_frequency"])
    stability["manuscript_stability_status"] = np.select(
        [
            (stability["weight_top_k_frequency"] >= threshold)
            & (stability["uncertainty_top_k_frequency"] >= threshold),
            (stability["weight_secondary_top_k_frequency"] >= threshold)
            & (stability["uncertainty_secondary_top_k_frequency"] >= threshold),
        ],
        ["stable_primary_shortlist", "stable_extended_shortlist"],
        default="sensitivity_dependent",
    )
    shortlist = stability.sort_values(
        ["robustness_composite_score", "v082_baseline_score"], ascending=False
    ).head(int(cfg["shortlist_size"])).copy()
    shortlist["manuscript_shortlist_rank"] = np.arange(1, len(shortlist) + 1)

    stable_primary = stability["manuscript_stability_status"].eq("stable_primary_shortlist")
    summary = {
        "analysis_status": "manuscript_readiness_sensitivity_ablation_and_negative_control_validation",
        "interpretation_scope": (
            "The suite evaluates ranking robustness to plausible weight and component uncertainty, "
            "leave-one-component-out ablation, marginal-preserving permutation nulls, evidence redundancy, "
            "and enrichment of pre-specified robust candidate sets. It does not replace biological safety review or wet-lab validation."
        ),
        "n_candidates": len(frame),
        "components": components,
        "component_weights": {c: float(w) for c, w in zip(components, weights)},
        "top_k": top_k,
        "secondary_top_k": secondary_top_k,
        "n_stable_primary_shortlist": int(stable_primary.sum()),
        "median_weight_spearman": float(w_iter["spearman_vs_baseline"].median()),
        "median_weight_top_k_jaccard": float(w_iter["top_k_jaccard"].median()),
        "median_uncertainty_spearman": float(u_iter["spearman_vs_baseline"].median()),
        "median_uncertainty_top_k_jaccard": float(u_iter["top_k_jaccard"].median()),
        "ranking_null": null_summary.to_dict(orient="records"),
        "most_influential_ablation": (
            ablation_summary.sort_values("mean_absolute_rank_change", ascending=False).head(3).to_dict(orient="records")
        ),
        "redundancy": redundancy.to_dict(orient="records")[0],
        "top_robust_candidates": shortlist[
            [c for c in [
                "cpg_id", "gene", "nearest_tss_gene", "unique_region_id", "v082_baseline_score",
                "robustness_composite_score", "weight_top_k_frequency", "uncertainty_top_k_frequency",
                "worst_ablation_rank_drop", "manuscript_stability_status"
            ] if c in shortlist.columns]
        ].head(20).to_dict(orient="records"),
        "seed": seed,
    }
    manifest = {
        "version": "0.8.2",
        "config": cfg,
        "component_columns_used": components,
        "input_columns": list(cpg.columns),
        "guardrails": [
            "Missing evidence components remain unavailable and are not imputed as zero.",
            "Permutation nulls preserve each evidence component's marginal distribution.",
            "Stable ranking is prioritisation evidence, not causal or clinical-safety proof.",
            "Literature and safety evidence completion remains mandatory before design-ready claims.",
        ],
    }
    return ManuscriptValidationResult(
        candidate_stability=stability,
        weight_iteration_summary=w_iter,
        uncertainty_iteration_summary=u_iter,
        ablation_summary=ablation_summary,
        ablation_candidate_shifts=ablation_shifts,
        null_summary=null_summary,
        null_iterations=null_iterations,
        component_correlation=corr,
        redundancy_summary=redundancy,
        candidate_set_enrichment=enrichment,
        robust_shortlist=shortlist,
        summary=_jsonable(summary),
        manifest=_jsonable(manifest),
    )


def _plot_outputs(result: ManuscriptValidationResult, output_dir: Path) -> None:
    import matplotlib
    # Figures are written to disk and never displayed, so force the
    # non-interactive backend. Without this, matplotlib may select an
    # interactive backend (e.g. TkAgg) and fail on machines with a
    # broken or absent Tcl/Tk, and in any headless environment.
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    top = result.robust_shortlist.head(30).sort_values("robustness_composite_score")
    plt.figure(figsize=(9, 8))
    plt.barh(top["cpg_id"].astype(str), top["robustness_composite_score"])
    plt.xlabel("Robustness composite score")
    plt.ylabel("CpG")
    plt.tight_layout()
    plt.savefig(fig_dir / "top_robust_candidates.png", dpi=220)
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.scatter(
        result.candidate_stability["v082_baseline_rank"],
        result.candidate_stability["weight_rank_iqr"],
        alpha=0.65,
    )
    plt.xlabel("Baseline rank")
    plt.ylabel("Rank IQR under weight perturbation")
    plt.tight_layout()
    plt.savefig(fig_dir / "weight_rank_stability.png", dpi=220)
    plt.close()

    abl = result.ablation_summary.sort_values("mean_absolute_rank_change")
    plt.figure(figsize=(9, 6))
    plt.barh(abl["removed_component"], abl["mean_absolute_rank_change"])
    plt.xlabel("Mean absolute rank change")
    plt.ylabel("Removed evidence component")
    plt.tight_layout()
    plt.savefig(fig_dir / "ablation_rank_impact.png", dpi=220)
    plt.close()

    null = result.null_iterations
    observed = float(result.null_summary.loc[result.null_summary["statistic"].eq("mean_top_k_score"), "observed"].iloc[0])
    plt.figure(figsize=(8, 5))
    plt.hist(null["null_mean_top_k_score"], bins=35)
    plt.axvline(observed, linestyle="--", linewidth=2)
    plt.xlabel("Mean top-k score under permutation null")
    plt.ylabel("Iterations")
    plt.tight_layout()
    plt.savefig(fig_dir / "ranking_specificity_null.png", dpi=220)
    plt.close()

    corr = result.component_correlation.set_index("component")
    plt.figure(figsize=(9, 8))
    image = plt.imshow(corr.to_numpy(dtype=float), vmin=-1, vmax=1)
    plt.xticks(range(len(corr.columns)), corr.columns, rotation=90)
    plt.yticks(range(len(corr.index)), corr.index)
    plt.colorbar(image, label="Spearman correlation")
    plt.tight_layout()
    plt.savefig(fig_dir / "component_correlation.png", dpi=220)
    plt.close()


def write_manuscript_validation_outputs(result: ManuscriptValidationResult, output_dir: str | Path) -> None:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    result.candidate_stability.to_csv(out / "candidate_rank_stability_v082.csv.gz", index=False, compression="gzip")
    result.weight_iteration_summary.to_csv(out / "weight_perturbation_iterations_v082.csv.gz", index=False, compression="gzip")
    result.uncertainty_iteration_summary.to_csv(out / "component_uncertainty_iterations_v082.csv.gz", index=False, compression="gzip")
    result.ablation_summary.to_csv(out / "ablation_summary_v082.csv", index=False)
    result.ablation_candidate_shifts.to_csv(out / "ablation_candidate_rank_shifts_v082.csv.gz", index=False, compression="gzip")
    result.null_summary.to_csv(out / "ranking_specificity_null_summary_v082.csv", index=False)
    result.null_iterations.to_csv(out / "ranking_specificity_null_iterations_v082.csv.gz", index=False, compression="gzip")
    result.component_correlation.to_csv(out / "component_correlation_v082.csv", index=False)
    result.redundancy_summary.to_csv(out / "component_redundancy_summary_v082.csv", index=False)
    result.candidate_set_enrichment.to_csv(out / "candidate_set_enrichment_v082.csv", index=False)
    result.robust_shortlist.to_csv(out / "manuscript_ready_robust_shortlist_v082.csv", index=False)
    (out / "manuscript_validation_v082_summary.json").write_text(
        json.dumps(result.summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (out / "manuscript_validation_v082_manifest.json").write_text(
        json.dumps(result.manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    _plot_outputs(result, out)
