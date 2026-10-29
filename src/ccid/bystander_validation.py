from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import json
import math

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


TARGETED_GENES = {"PDE4C", "FHL2", "ELOVL2", "KLF14", "TEAD1"}


@dataclass(frozen=True)
class CellTypeEffectResult:
    locus_effects: pd.DataFrame
    pair_effects: pd.DataFrame
    pair_manifest: pd.DataFrame
    sample_qc: pd.DataFrame


@dataclass(frozen=True)
class BystanderValidationResult:
    coverage: pd.DataFrame
    candidate_set_summary: pd.DataFrame
    candidate_effects: pd.DataFrame
    direction_compatibility: pd.DataFrame
    cross_celltype_effects: pd.DataFrame
    cross_celltype_summary: pd.DataFrame
    null_tests: pd.DataFrame
    regional_clusters: pd.DataFrame
    regional_scores: pd.DataFrame
    pair_manifest: pd.DataFrame
    sample_qc: pd.DataFrame
    tcell_locus_effects: pd.DataFrame
    msc_locus_effects: pd.DataFrame
    tcell_pair_effects: pd.DataFrame
    msc_pair_effects: pd.DataFrame
    summary: dict[str, object]
    manifest: dict[str, object]


def _safe_spearman(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    valid = np.isfinite(a) & np.isfinite(b)
    if valid.sum() < 3 or np.nanstd(a[valid]) == 0 or np.nanstd(b[valid]) == 0:
        return float("nan"), float("nan")
    rho, p = stats.spearmanr(a[valid], b[valid])
    return float(rho), float(p)


def _safe_pearson(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    valid = np.isfinite(a) & np.isfinite(b)
    if valid.sum() < 3 or np.nanstd(a[valid]) == 0 or np.nanstd(b[valid]) == 0:
        return float("nan"), float("nan")
    r, p = stats.pearsonr(a[valid], b[valid])
    return float(r), float(p)


def _empirical_p(observed: float, null: np.ndarray, *, alternative: str = "greater") -> float:
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if not np.isfinite(observed) or len(null) == 0:
        return float("nan")
    if alternative == "greater":
        extreme = np.sum(null >= observed)
    elif alternative == "less":
        extreme = np.sum(null <= observed)
    elif alternative == "two-sided":
        extreme = np.sum(np.abs(null) >= abs(observed))
    else:
        raise ValueError(f"Unknown alternative: {alternative}")
    return float((1 + extreme) / (1 + len(null)))


def _normalise_cpg_column(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    candidates = [c for c in frame.columns if str(c).lower() in {"cpg_id", "locus_id", "probe_id", "id_ref", "name"}]
    if candidates:
        frame = frame.rename(columns={candidates[0]: "cpg_id"})
    elif frame.index.name and str(frame.index.name).lower() in {"cpg_id", "locus_id", "probe_id", "id_ref", "name"}:
        frame = frame.reset_index().rename(columns={frame.index.name: "cpg_id"})
    elif frame.shape[1] == 1:
        frame = frame.rename(columns={frame.columns[0]: "cpg_id"})
    else:
        raise ValueError("Could not identify a CpG ID column")
    frame["cpg_id"] = frame["cpg_id"].astype(str).str.replace(r"_(BC|TC)\d+$", "", regex=True)
    return frame.drop_duplicates("cpg_id")


def load_candidate_sets(
    *,
    stable_core_path: str | Path,
    consensus_2of3_path: str | Path,
    tier1_path: str | Path,
    frozen_signature_path: str | Path,
) -> tuple[dict[str, set[str]], pd.DataFrame]:
    frozen = _normalise_cpg_column(pd.read_csv(frozen_signature_path))
    if "target_delta_beta" not in frozen.columns:
        raise ValueError("Frozen signature must contain target_delta_beta")
    frozen["target_delta_beta"] = pd.to_numeric(frozen["target_delta_beta"], errors="coerce")
    target_map = frozen[["cpg_id", "target_delta_beta"]].drop_duplicates("cpg_id")

    stable = _normalise_cpg_column(pd.read_csv(stable_core_path))
    consensus = _normalise_cpg_column(pd.read_csv(consensus_2of3_path))
    tier1 = _normalise_cpg_column(pd.read_csv(tier1_path))

    sets = {
        "stable_core_3of3": set(stable["cpg_id"]),
        "consensus_at_least_2of3": set(consensus["cpg_id"]),
        "o3_tier1": set(tier1["cpg_id"]),
        "frozen_500": set(frozen["cpg_id"]),
    }
    if not sets["stable_core_3of3"].issubset(sets["consensus_at_least_2of3"]):
        raise ValueError("The 3/3 stable core is not a subset of the >=2/3 consensus")
    if not sets["consensus_at_least_2of3"].issubset(sets["frozen_500"]):
        # LODO rotations may contain loci not present in the original O1/O2 frozen set.
        # Preserve them but target direction will be missing where unavailable.
        pass
    return sets, target_map


def build_paired_celltype_effects(
    beta: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    cell_type: str,
    targeted_label: str = "targeted",
    control_label: str = "control",
    treatment_column: str = "treatment",
    pair_column: str = "replicate",
    min_abs_delta: float = 0.02,
    fdr: float = 0.10,
) -> CellTypeEffectResult:
    required = {"sample_id", treatment_column, pair_column}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing columns: {sorted(missing)}")

    beta = beta.copy()
    beta.index = beta.index.astype(str).str.replace(r"_(BC|TC)\d+$", "", regex=True)
    if beta.index.duplicated().any():
        beta = beta.groupby(level=0, sort=False).mean()
    beta = beta.apply(pd.to_numeric, errors="coerce")

    meta = metadata.copy()
    meta["sample_id"] = meta["sample_id"].astype(str)
    meta = meta[meta["sample_id"].isin(beta.columns)].copy()
    if meta.empty:
        raise ValueError(f"No metadata samples overlap {cell_type} beta matrix")

    pair_rows: list[dict[str, object]] = []
    pair_effects: dict[str, pd.Series] = {}
    for pair, group in meta.groupby(pair_column, sort=True, dropna=False):
        targeted = group.loc[group[treatment_column].astype(str) == targeted_label, "sample_id"].tolist()
        controls = group.loc[group[treatment_column].astype(str) == control_label, "sample_id"].tolist()
        if not targeted or not controls:
            continue
        pair_id = str(pair)
        pair_effects[pair_id] = beta[targeted].mean(axis=1) - beta[controls].mean(axis=1)
        pair_rows.append(
            {
                "cell_type": cell_type,
                "pair_id": pair_id,
                "targeted_sample_ids": ";".join(targeted),
                "control_sample_ids": ";".join(controls),
                "n_targeted": len(targeted),
                "n_controls": len(controls),
            }
        )
    if len(pair_effects) < 2:
        raise ValueError(f"At least two complete targeted-control pairs are required for {cell_type}")

    pair_frame = pd.DataFrame(pair_effects, index=beta.index)
    values = pair_frame.to_numpy(dtype=float, copy=True)
    valid_n = np.isfinite(values).sum(axis=1)
    mean_delta = np.nanmean(values, axis=1)
    median_delta = np.nanmedian(values, axis=1)
    sd_delta = np.nanstd(values, axis=1, ddof=1)
    se_delta = sd_delta / np.sqrt(np.maximum(valid_n, 1))
    _, p_values = stats.ttest_1samp(values, 0.0, axis=1, nan_policy="omit")
    p_values = np.nan_to_num(p_values, nan=1.0, posinf=1.0, neginf=1.0)
    p_values = np.where(valid_n >= 2, p_values, 1.0)
    q_values = multipletests(p_values, method="fdr_bh")[1]
    direction = np.sign(mean_delta)
    same = np.where(
        direction[:, None] == 0,
        0.0,
        (np.sign(values) == direction[:, None]).astype(float),
    )
    same[~np.isfinite(values)] = np.nan
    sign_consistency = np.nanmean(same, axis=1)

    control_ids = meta.loc[meta[treatment_column].astype(str) == control_label, "sample_id"].tolist()
    targeted_ids = meta.loc[meta[treatment_column].astype(str) == targeted_label, "sample_id"].tolist()
    baseline_beta = beta[control_ids].mean(axis=1)
    targeted_beta = beta[targeted_ids].mean(axis=1)

    locus = pd.DataFrame(
        {
            "cpg_id": beta.index.astype(str),
            "cell_type": cell_type,
            "baseline_beta": baseline_beta.to_numpy(dtype=float),
            "targeted_beta": targeted_beta.to_numpy(dtype=float),
            "editor_delta_beta": mean_delta,
            "median_editor_delta_beta": median_delta,
            "editor_sd": sd_delta,
            "editor_se": se_delta,
            "n_pairs": valid_n,
            "sign_consistency": sign_consistency,
            "p_value": p_values,
            "q_value": q_values,
        }
    )
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
        pair_frame.rename_axis("cpg_id")
        .reset_index()
        .melt(id_vars="cpg_id", var_name="pair_id", value_name="pair_delta_beta")
    )
    pair_long.insert(1, "cell_type", cell_type)

    sample_qc_rows: list[dict[str, object]] = []
    meta_index = meta.set_index("sample_id")
    for sample_id in meta["sample_id"]:
        values_s = pd.to_numeric(beta[sample_id], errors="coerce")
        row = {
            "cell_type": cell_type,
            "sample_id": sample_id,
            "treatment": meta_index.loc[sample_id, treatment_column],
            "pair_id": meta_index.loc[sample_id, pair_column],
            "n_loci": int(len(values_s)),
            "valid_loci": int(values_s.notna().sum()),
            "valid_fraction": float(values_s.notna().mean()),
            "mean_beta": float(values_s.mean()),
            "median_beta": float(values_s.median()),
            "sd_beta": float(values_s.std()),
        }
        sample_qc_rows.append(row)

    return CellTypeEffectResult(
        locus_effects=locus,
        pair_effects=pair_long,
        pair_manifest=pd.DataFrame(pair_rows),
        sample_qc=pd.DataFrame(sample_qc_rows),
    )


def _set_membership_frame(candidate_sets: Mapping[str, set[str]]) -> pd.DataFrame:
    rows = [
        {"candidate_set": set_name, "cpg_id": cpg}
        for set_name, cpgs in candidate_sets.items()
        for cpg in sorted(cpgs)
    ]
    return pd.DataFrame(rows)


def build_coverage_and_candidate_effects(
    candidate_sets: Mapping[str, set[str]],
    target_map: pd.DataFrame,
    cell_effects: Mapping[str, pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    membership = _set_membership_frame(candidate_sets)
    membership = membership.merge(target_map, on="cpg_id", how="left")
    coverage_rows: list[dict[str, object]] = []
    effect_frames: list[pd.DataFrame] = []

    for cell_type, effects in cell_effects.items():
        measured = set(effects.loc[effects["editor_delta_beta"].notna(), "cpg_id"].astype(str))
        for set_name, cpgs in candidate_sets.items():
            observed = sorted(set(cpgs).intersection(measured))
            coverage_rows.append(
                {
                    "candidate_set": set_name,
                    "cell_type": cell_type,
                    "n_candidates": len(cpgs),
                    "n_measured": len(observed),
                    "coverage_fraction": len(observed) / max(len(cpgs), 1),
                    "n_missing": len(cpgs) - len(observed),
                }
            )
        merged = membership.merge(effects, on="cpg_id", how="left")
        merged["cell_type"] = cell_type
        effect_frames.append(merged)

    candidate_effects = pd.concat(effect_frames, ignore_index=True)
    candidate_effects["target_direction"] = np.select(
        [candidate_effects["target_delta_beta"] > 0, candidate_effects["target_delta_beta"] < 0],
        ["hypermethylation", "hypomethylation"],
        default="unknown",
    )
    candidate_effects["same_as_target_direction"] = (
        np.sign(candidate_effects["editor_delta_beta"])
        == np.sign(candidate_effects["target_delta_beta"])
    ) & candidate_effects["editor_delta_beta"].notna() & candidate_effects["target_delta_beta"].notna()
    candidate_effects["opposite_to_target_direction"] = (
        np.sign(candidate_effects["editor_delta_beta"])
        == -np.sign(candidate_effects["target_delta_beta"])
    ) & candidate_effects["editor_delta_beta"].notna() & candidate_effects["target_delta_beta"].notna()

    summary_rows: list[dict[str, object]] = []
    for (set_name, cell_type), group in candidate_effects.groupby(["candidate_set", "cell_type"], sort=True):
        valid = group[group["editor_delta_beta"].notna()]
        summary_rows.append(
            {
                "candidate_set": set_name,
                "cell_type": cell_type,
                "n_measured": int(len(valid)),
                "median_editor_delta_beta": float(valid["editor_delta_beta"].median()) if len(valid) else float("nan"),
                "median_abs_editor_delta_beta": float(valid["editor_delta_beta"].abs().median()) if len(valid) else float("nan"),
                "mean_sign_consistency": float(valid["sign_consistency"].mean()) if len(valid) else float("nan"),
                "significant_fraction": float(valid["significant_editor_effect"].mean()) if len(valid) else float("nan"),
                "same_as_target_fraction": float(valid["same_as_target_direction"].mean()) if len(valid) else float("nan"),
                "opposite_to_target_fraction": float(valid["opposite_to_target_direction"].mean()) if len(valid) else float("nan"),
            }
        )
    return pd.DataFrame(coverage_rows), candidate_effects, pd.DataFrame(summary_rows)


def build_cross_celltype_results(
    candidate_sets: Mapping[str, set[str]],
    target_map: pd.DataFrame,
    tcell_effects: pd.DataFrame,
    msc_effects: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    t = tcell_effects.rename(
        columns={c: f"tcell_{c}" for c in tcell_effects.columns if c not in {"cpg_id", "cell_type"}}
    ).drop(columns=["cell_type"])
    m = msc_effects.rename(
        columns={c: f"msc_{c}" for c in msc_effects.columns if c not in {"cpg_id", "cell_type"}}
    ).drop(columns=["cell_type"])
    merged = target_map.merge(t, on="cpg_id", how="outer").merge(m, on="cpg_id", how="outer")

    membership = _set_membership_frame(candidate_sets)
    cross = membership.merge(merged, on="cpg_id", how="left")
    tdelta = cross["tcell_editor_delta_beta"].to_numpy(dtype=float)
    mdelta = cross["msc_editor_delta_beta"].to_numpy(dtype=float)
    target = cross["target_delta_beta"].to_numpy(dtype=float)
    valid_both = np.isfinite(tdelta) & np.isfinite(mdelta)
    cross["cross_cell_same_direction"] = valid_both & (np.sign(tdelta) == np.sign(mdelta))
    cross["both_significant"] = (
        cross["tcell_significant_editor_effect"].fillna(False).astype(bool)
        & cross["msc_significant_editor_effect"].fillna(False).astype(bool)
    )
    cross["mean_abs_editor_delta"] = np.nanmean(np.abs(np.column_stack([tdelta, mdelta])), axis=1)
    cross["mean_sign_consistency"] = np.nanmean(
        np.column_stack([
            cross["tcell_sign_consistency"].to_numpy(dtype=float),
            cross["msc_sign_consistency"].to_numpy(dtype=float),
        ]),
        axis=1,
    )
    cross["target_same_both"] = valid_both & np.isfinite(target) & (
        np.sign(tdelta) == np.sign(target)
    ) & (np.sign(mdelta) == np.sign(target))
    cross["target_opposite_both"] = valid_both & np.isfinite(target) & (
        np.sign(tdelta) == -np.sign(target)
    ) & (np.sign(mdelta) == -np.sign(target))

    summaries: list[dict[str, object]] = []
    direction_rows: list[dict[str, object]] = []
    for set_name, group in cross.groupby("candidate_set", sort=True):
        valid = group[group["tcell_editor_delta_beta"].notna() & group["msc_editor_delta_beta"].notna()]
        a = valid["tcell_editor_delta_beta"].to_numpy(dtype=float)
        b = valid["msc_editor_delta_beta"].to_numpy(dtype=float)
        pearson_r, pearson_p = _safe_pearson(a, b)
        spearman_rho, spearman_p = _safe_spearman(a, b)
        summaries.append(
            {
                "candidate_set": set_name,
                "n_common_measured": int(len(valid)),
                "cross_cell_sign_concordance": float(valid["cross_cell_same_direction"].mean()) if len(valid) else float("nan"),
                "pearson_r": pearson_r,
                "pearson_p": pearson_p,
                "spearman_rho": spearman_rho,
                "spearman_p": spearman_p,
                "median_abs_tcell_delta": float(valid["tcell_editor_delta_beta"].abs().median()) if len(valid) else float("nan"),
                "median_abs_msc_delta": float(valid["msc_editor_delta_beta"].abs().median()) if len(valid) else float("nan"),
                "both_significant_fraction": float(valid["both_significant"].mean()) if len(valid) else float("nan"),
                "target_same_both_fraction": float(valid["target_same_both"].mean()) if len(valid) else float("nan"),
                "target_opposite_both_fraction": float(valid["target_opposite_both"].mean()) if len(valid) else float("nan"),
            }
        )

    unique = cross.drop_duplicates("cpg_id").copy()
    for row in unique.itertuples(index=False):
        tval = getattr(row, "tcell_editor_delta_beta")
        mval = getattr(row, "msc_editor_delta_beta")
        target_val = getattr(row, "target_delta_beta")
        if not np.isfinite(tval) and not np.isfinite(mval):
            category = "not_measured"
        elif np.isfinite(tval) and np.isfinite(mval) and np.sign(tval) == np.sign(mval):
            if np.isfinite(target_val) and np.sign(tval) == np.sign(target_val):
                category = "direct_direction_support"
            elif np.isfinite(target_val) and np.sign(tval) == -np.sign(target_val):
                category = "opposite_direction_responsiveness"
            else:
                category = "cross_cell_consistent_target_unknown"
        else:
            category = "cell_type_specific_or_mixed"
        score = np.nanmean([abs(tval) if np.isfinite(tval) else np.nan, abs(mval) if np.isfinite(mval) else np.nan])
        consistency = np.nanmean([
            getattr(row, "tcell_sign_consistency") if np.isfinite(getattr(row, "tcell_sign_consistency")) else np.nan,
            getattr(row, "msc_sign_consistency") if np.isfinite(getattr(row, "msc_sign_consistency")) else np.nan,
        ])
        direction_rows.append(
            {
                "cpg_id": getattr(row, "cpg_id"),
                "target_delta_beta": target_val,
                "tcell_editor_delta_beta": tval,
                "msc_editor_delta_beta": mval,
                "direction_category": category,
                "bystander_responsiveness_score": float(score * consistency) if np.isfinite(score) and np.isfinite(consistency) else float("nan"),
            }
        )
    return cross, pd.DataFrame(summaries), pd.DataFrame(direction_rows)


def _quantile_bin(series: pd.Series, n_bins: int) -> pd.Series:
    valid = pd.to_numeric(series, errors="coerce")
    if valid.notna().sum() == 0:
        return pd.Series(pd.NA, index=series.index, dtype="Int64")
    ranks = valid.rank(method="average", pct=True)
    bins = np.minimum((ranks * n_bins).fillna(0).astype(int), n_bins - 1)
    bins[valid.isna()] = -1
    return bins.astype(int)


def run_matched_null_tests(
    candidate_sets: Mapping[str, set[str]],
    cell_effects: Mapping[str, pd.DataFrame],
    discovery_all_loci: pd.DataFrame | None,
    *,
    n_iterations: int = 1000,
    n_bins: int = 10,
    seed: int = 0,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    discovery_map: pd.DataFrame | None = None
    if discovery_all_loci is not None:
        discovery_map = _normalise_cpg_column(discovery_all_loci)
        if "target_delta_beta" in discovery_map.columns:
            discovery_map = discovery_map[["cpg_id", "target_delta_beta"]]
        else:
            discovery_map = None

    rows: list[dict[str, object]] = []
    for cell_type, effects in cell_effects.items():
        universe = effects.loc[effects["editor_delta_beta"].notna(), [
            "cpg_id", "baseline_beta", "editor_delta_beta", "sign_consistency", "significant_editor_effect"
        ]].copy()
        if discovery_map is not None:
            universe = universe.merge(discovery_map, on="cpg_id", how="left")
        else:
            universe["target_delta_beta"] = np.nan
        universe["baseline_bin"] = _quantile_bin(universe["baseline_beta"], n_bins)
        universe["target_abs_bin"] = _quantile_bin(universe["target_delta_beta"].abs(), n_bins)
        universe["match_bin"] = (
            universe["baseline_bin"].astype(str)
            + ":"
            + universe["target_abs_bin"].astype(str)
        )

        for set_name, cpgs in candidate_sets.items():
            candidate = universe[universe["cpg_id"].isin(cpgs)].copy()
            if candidate.empty:
                continue
            observed = {
                "median_abs_delta": float(candidate["editor_delta_beta"].abs().median()),
                "mean_abs_delta": float(candidate["editor_delta_beta"].abs().mean()),
                "mean_sign_consistency": float(candidate["sign_consistency"].mean()),
                "significant_fraction": float(candidate["significant_editor_effect"].mean()),
            }
            bin_counts = candidate["match_bin"].value_counts().to_dict()
            pools = {b: universe[universe["match_bin"] == b] for b in bin_counts}
            null_metrics = {key: np.full(n_iterations, np.nan) for key in observed}
            for i in range(n_iterations):
                sampled_parts: list[pd.DataFrame] = []
                for b, count in bin_counts.items():
                    pool = pools[b]
                    if len(pool) == 0:
                        continue
                    replace = len(pool) < count
                    take = rng.choice(len(pool), size=count, replace=replace)
                    sampled_parts.append(pool.iloc[take])
                if not sampled_parts:
                    continue
                sampled = pd.concat(sampled_parts, ignore_index=True)
                null_metrics["median_abs_delta"][i] = sampled["editor_delta_beta"].abs().median()
                null_metrics["mean_abs_delta"][i] = sampled["editor_delta_beta"].abs().mean()
                null_metrics["mean_sign_consistency"][i] = sampled["sign_consistency"].mean()
                null_metrics["significant_fraction"][i] = sampled["significant_editor_effect"].mean()
            for metric, observed_value in observed.items():
                null = null_metrics[metric]
                for i, null_value in enumerate(null, start=1):
                    rows.append(
                        {
                            "candidate_set": set_name,
                            "cell_type": cell_type,
                            "metric": metric,
                            "iteration": i,
                            "observed": observed_value,
                            "null_value": null_value,
                            "empirical_p": _empirical_p(observed_value, null, alternative="greater"),
                            "matching": "baseline_beta_and_discovery_abs_delta",
                            "n_candidates_measured": len(candidate),
                        }
                    )
    return pd.DataFrame(rows)


def load_probe_annotation(path: str | Path | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame(columns=["cpg_id", "chromosome", "position"])
    frame = pd.read_csv(path, compression="infer", low_memory=False)
    frame = _normalise_cpg_column(frame)
    rename: dict[str, str] = {}
    for col in frame.columns:
        low = str(col).lower()
        if low in {"chr", "chrom", "chromosome"}:
            rename[col] = "chromosome"
        elif low in {"pos", "position", "mapinfo", "coordinate"}:
            rename[col] = "position"
        elif low in {"gene", "gene_name", "ucsc_refgene_name"}:
            rename[col] = "gene"
        elif low in {"relation_to_island", "relation_to_ucsc_cpg_island", "island_relation"}:
            rename[col] = "island_relation"
    frame = frame.rename(columns=rename)
    required = {"cpg_id", "chromosome", "position"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Probe annotation is missing columns: {sorted(missing)}")
    frame["chromosome"] = frame["chromosome"].astype(str).str.replace("^chr", "", regex=True)
    frame["position"] = pd.to_numeric(frame["position"], errors="coerce")
    frame = frame.dropna(subset=["chromosome", "position"]).drop_duplicates("cpg_id")
    frame["position"] = frame["position"].astype(int)
    if "gene" not in frame.columns:
        frame["gene"] = ""
    if "island_relation" not in frame.columns:
        frame["island_relation"] = ""
    return frame


def build_regional_scores(
    candidate_sets: Mapping[str, set[str]],
    target_map: pd.DataFrame,
    annotation: pd.DataFrame,
    tcell_effects: pd.DataFrame,
    msc_effects: pd.DataFrame,
    *,
    cluster_gap_bp: int = 2000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if annotation.empty:
        return (
            pd.DataFrame(columns=["region_id", "chromosome", "start", "end", "n_candidate_cpgs"]),
            pd.DataFrame(columns=["region_id", "candidate_set", "regional_responsiveness_score"]),
        )

    membership = _set_membership_frame(candidate_sets)
    annotated = membership.merge(annotation, on="cpg_id", how="left").dropna(subset=["chromosome", "position"])
    annotated = annotated.sort_values(["candidate_set", "chromosome", "position", "cpg_id"])
    cluster_rows: list[dict[str, object]] = []
    for (set_name, chromosome), group in annotated.groupby(["candidate_set", "chromosome"], sort=True):
        group = group.sort_values("position")
        cluster_index = 0
        previous: int | None = None
        region_ids: list[str] = []
        for position in group["position"].astype(int):
            if previous is None or position - previous > cluster_gap_bp:
                cluster_index += 1
            region_ids.append(f"{set_name}:{chromosome}:{cluster_index:04d}")
            previous = position
        group = group.copy()
        group["region_id"] = region_ids
        cluster_rows.append(group)
    clusters_long = pd.concat(cluster_rows, ignore_index=True) if cluster_rows else pd.DataFrame()
    if clusters_long.empty:
        return clusters_long, pd.DataFrame()

    t = tcell_effects[["cpg_id", "editor_delta_beta", "sign_consistency", "significant_editor_effect"]].rename(
        columns={
            "editor_delta_beta": "tcell_delta",
            "sign_consistency": "tcell_consistency",
            "significant_editor_effect": "tcell_significant",
        }
    )
    m = msc_effects[["cpg_id", "editor_delta_beta", "sign_consistency", "significant_editor_effect"]].rename(
        columns={
            "editor_delta_beta": "msc_delta",
            "sign_consistency": "msc_consistency",
            "significant_editor_effect": "msc_significant",
        }
    )
    clusters_long = clusters_long.merge(target_map, on="cpg_id", how="left").merge(t, on="cpg_id", how="left").merge(m, on="cpg_id", how="left")
    clusters_long["direct_target_gene_flag"] = clusters_long["gene"].fillna("").astype(str).apply(
        lambda value: any(gene in {x.strip() for x in value.replace(";", ",").split(",") if x.strip()} for gene in TARGETED_GENES)
    )

    score_rows: list[dict[str, object]] = []
    for (set_name, region_id), group in clusters_long.groupby(["candidate_set", "region_id"], sort=True):
        tdelta = group["tcell_delta"].to_numpy(dtype=float)
        mdelta = group["msc_delta"].to_numpy(dtype=float)
        paired = np.isfinite(tdelta) & np.isfinite(mdelta)
        cross_concordance = float(np.mean(np.sign(tdelta[paired]) == np.sign(mdelta[paired]))) if paired.any() else float("nan")
        mean_abs = float(np.nanmean(np.abs(np.column_stack([tdelta, mdelta]))))
        mean_consistency = float(np.nanmean(np.column_stack([
            group["tcell_consistency"].to_numpy(dtype=float),
            group["msc_consistency"].to_numpy(dtype=float),
        ])))
        measured_fraction = float(np.mean(np.isfinite(tdelta) | np.isfinite(mdelta)))
        concordance_factor = cross_concordance if np.isfinite(cross_concordance) else 0.5
        regional_score = mean_abs * mean_consistency * measured_fraction * concordance_factor
        score_rows.append(
            {
                "region_id": region_id,
                "candidate_set": set_name,
                "chromosome": str(group["chromosome"].iloc[0]),
                "start": int(group["position"].min()),
                "end": int(group["position"].max()),
                "width_bp": int(group["position"].max() - group["position"].min() + 1),
                "n_candidate_cpgs": int(group["cpg_id"].nunique()),
                "candidate_cpg_ids": ";".join(sorted(group["cpg_id"].unique())),
                "genes": ";".join(sorted({g for value in group["gene"].fillna("") for g in str(value).replace(",", ";").split(";") if g})),
                "contains_published_direct_target_gene": bool(group["direct_target_gene_flag"].any()),
                "measured_fraction": measured_fraction,
                "cross_cell_sign_concordance": cross_concordance,
                "mean_abs_editor_delta": mean_abs,
                "mean_sign_consistency": mean_consistency,
                "regional_responsiveness_score": regional_score,
                "interpretation": "regional_bystander_responsiveness_not_direct_editability",
            }
        )
    return clusters_long, pd.DataFrame(score_rows).sort_values("regional_responsiveness_score", ascending=False)


def _jsonable(value: object) -> object:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def run_bystander_validation(
    *,
    tcell_beta: pd.DataFrame,
    tcell_metadata: pd.DataFrame,
    msc_beta: pd.DataFrame,
    msc_metadata: pd.DataFrame,
    candidate_sets: Mapping[str, set[str]],
    target_map: pd.DataFrame,
    discovery_all_loci: pd.DataFrame | None = None,
    probe_annotation: pd.DataFrame | None = None,
    min_abs_delta: float = 0.02,
    fdr: float = 0.10,
    null_iterations: int = 1000,
    null_bins: int = 10,
    seed: int = 20260721,
    cluster_gap_bp: int = 2000,
) -> BystanderValidationResult:
    tcell = build_paired_celltype_effects(
        tcell_beta,
        tcell_metadata,
        cell_type="t_cell",
        min_abs_delta=min_abs_delta,
        fdr=fdr,
    )
    msc = build_paired_celltype_effects(
        msc_beta,
        msc_metadata,
        cell_type="mesenchymal_stromal_cell",
        min_abs_delta=min_abs_delta,
        fdr=fdr,
    )
    cell_effects = {
        "t_cell": tcell.locus_effects,
        "mesenchymal_stromal_cell": msc.locus_effects,
    }
    coverage, candidate_effects, candidate_summary = build_coverage_and_candidate_effects(
        candidate_sets, target_map, cell_effects
    )
    cross_effects, cross_summary, direction = build_cross_celltype_results(
        candidate_sets, target_map, tcell.locus_effects, msc.locus_effects
    )
    null_tests = run_matched_null_tests(
        candidate_sets,
        cell_effects,
        discovery_all_loci,
        n_iterations=null_iterations,
        n_bins=null_bins,
        seed=seed,
    )
    annotation = probe_annotation if probe_annotation is not None else pd.DataFrame()
    clusters, regional = build_regional_scores(
        candidate_sets,
        target_map,
        annotation,
        tcell.locus_effects,
        msc.locus_effects,
        cluster_gap_bp=cluster_gap_bp,
    )

    null_summary_rows: list[dict[str, object]] = []
    if not null_tests.empty:
        for keys, group in null_tests.groupby(["candidate_set", "cell_type", "metric"], sort=True):
            set_name, cell_type, metric = keys
            null_summary_rows.append(
                {
                    "candidate_set": set_name,
                    "cell_type": cell_type,
                    "metric": metric,
                    "observed": float(group["observed"].iloc[0]),
                    "null_mean": float(group["null_value"].mean()),
                    "null_q95": float(group["null_value"].quantile(0.95)),
                    "empirical_p": float(group["empirical_p"].iloc[0]),
                }
            )

    summary: dict[str, object] = {
        "analysis_status": "paired_cross_celltype_bystander_responsiveness_validation",
        "interpretation_scope": (
            "The dCas9-DNMT3A datasets support locus and region responsiveness to a hypermethylating "
            "editor. They do not directly prove that hypomethylating safe-reprogramming targets can be "
            "reproduced with dCas9-TET1. Regional scores are prioritisation evidence, not causal proof."
        ),
        "published_editor_targets": sorted(TARGETED_GENES),
        "n_tcell_pairs": int(len(tcell.pair_manifest)),
        "n_msc_pairs": int(len(msc.pair_manifest)),
        "candidate_set_sizes": {name: len(values) for name, values in candidate_sets.items()},
        "coverage": coverage.to_dict("records"),
        "candidate_set_summary": candidate_summary.to_dict("records"),
        "cross_celltype_summary": cross_summary.to_dict("records"),
        "null_summary": null_summary_rows,
        "probe_annotation_supplied": bool(probe_annotation is not None and not probe_annotation.empty),
        "n_regional_scores": int(len(regional)),
        "min_abs_delta": min_abs_delta,
        "fdr": fdr,
        "null_iterations": null_iterations,
        "null_bins": null_bins,
        "seed": seed,
    }
    manifest: dict[str, object] = {
        "candidate_sets_frozen_before_bystander_analysis": True,
        "candidate_set_names": list(candidate_sets),
        "targeted_label": "targeted",
        "control_label": "control",
        "pair_column": "replicate",
        "cell_types": ["t_cell", "mesenchymal_stromal_cell"],
        "editor": "dCas9-DNMT3A/3L",
        "null_matching": "baseline beta plus discovery absolute target effect quantile bins when available",
        "regional_cluster_gap_bp": cluster_gap_bp,
        "regional_claim_guardrail": "regional_bystander_responsiveness_not_direct_editability",
    }

    return BystanderValidationResult(
        coverage=coverage,
        candidate_set_summary=candidate_summary,
        candidate_effects=candidate_effects,
        direction_compatibility=direction,
        cross_celltype_effects=cross_effects,
        cross_celltype_summary=cross_summary,
        null_tests=null_tests,
        regional_clusters=clusters,
        regional_scores=regional,
        pair_manifest=pd.concat([tcell.pair_manifest, msc.pair_manifest], ignore_index=True),
        sample_qc=pd.concat([tcell.sample_qc, msc.sample_qc], ignore_index=True),
        tcell_locus_effects=tcell.locus_effects,
        msc_locus_effects=msc.locus_effects,
        tcell_pair_effects=tcell.pair_effects,
        msc_pair_effects=msc.pair_effects,
        summary=_jsonable(summary),
        manifest=_jsonable(manifest),
    )


def write_bystander_validation_outputs(
    result: BystanderValidationResult,
    output_dir: str | Path,
) -> None:
    output = Path(output_dir)
    figures = output / "figures"
    output.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    result.coverage.to_csv(output / "candidate_probe_coverage.csv", index=False)
    result.candidate_set_summary.to_csv(output / "candidate_set_summary.csv", index=False)
    result.candidate_effects.to_csv(output / "candidate_bystander_effects.csv.gz", index=False, compression="gzip")
    result.direction_compatibility.to_csv(output / "direction_compatibility.csv", index=False)
    result.cross_celltype_effects.to_csv(output / "cross_celltype_candidate_effects.csv.gz", index=False, compression="gzip")
    result.cross_celltype_summary.to_csv(output / "cross_celltype_replication.csv", index=False)
    result.null_tests.to_csv(output / "empirical_null_tests.csv.gz", index=False, compression="gzip")
    result.regional_clusters.to_csv(output / "regional_candidate_clusters.csv.gz", index=False, compression="gzip")
    result.regional_scores.to_csv(output / "regional_editability_scores.csv", index=False)
    result.pair_manifest.to_csv(output / "sample_pair_manifest.csv", index=False)
    result.sample_qc.to_csv(output / "sample_qc.csv", index=False)
    result.tcell_locus_effects.to_csv(output / "tcell_cpg_effects.csv.gz", index=False, compression="gzip")
    result.msc_locus_effects.to_csv(output / "msc_cpg_effects.csv.gz", index=False, compression="gzip")
    result.tcell_pair_effects.to_csv(output / "tcell_pair_effects.csv.gz", index=False, compression="gzip")
    result.msc_pair_effects.to_csv(output / "msc_pair_effects.csv.gz", index=False, compression="gzip")
    (output / "bystander_validation_summary.json").write_text(json.dumps(result.summary, indent=2))
    (output / "bystander_analysis_manifest.json").write_text(json.dumps(result.manifest, indent=2))

    _write_figures(result, figures)


def _write_figures(result: BystanderValidationResult, figures: Path) -> None:
    import matplotlib
    # Figures are written to disk and never displayed, so force the
    # non-interactive backend. Without this, matplotlib may select an
    # interactive backend (e.g. TkAgg) and fail on machines with a
    # broken or absent Tcl/Tk, and in any headless environment.
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    coverage = result.coverage.copy()
    if not coverage.empty:
        pivot = coverage.pivot(index="candidate_set", columns="cell_type", values="coverage_fraction")
        ax = pivot.plot(kind="bar", figsize=(9, 5))
        ax.set_ylabel("Measured fraction")
        ax.set_xlabel("Candidate set")
        ax.set_ylim(0, 1.05)
        ax.set_title("Candidate CpG coverage")
        plt.tight_layout()
        plt.savefig(figures / "candidate_coverage.png", dpi=200)
        plt.close()

    cross = result.cross_celltype_effects
    stable = cross[cross["candidate_set"] == "stable_core_3of3"].dropna(
        subset=["tcell_editor_delta_beta", "msc_editor_delta_beta"]
    )
    if stable.empty:
        stable = cross.dropna(subset=["tcell_editor_delta_beta", "msc_editor_delta_beta"])
    if not stable.empty:
        plt.figure(figsize=(6, 6))
        plt.scatter(stable["tcell_editor_delta_beta"], stable["msc_editor_delta_beta"], alpha=0.7)
        low = float(min(stable["tcell_editor_delta_beta"].min(), stable["msc_editor_delta_beta"].min()))
        high = float(max(stable["tcell_editor_delta_beta"].max(), stable["msc_editor_delta_beta"].max()))
        plt.plot([low, high], [low, high], linestyle="--")
        plt.xlabel("T-cell editor Δβ")
        plt.ylabel("MSC editor Δβ")
        plt.title("Cross-cell-type bystander response")
        plt.tight_layout()
        plt.savefig(figures / "tcell_msc_effect_correlation.png", dpi=200)
        plt.close()

    candidate_summary = result.candidate_set_summary
    if not candidate_summary.empty:
        pivot = candidate_summary.pivot(index="candidate_set", columns="cell_type", values="median_abs_editor_delta_beta")
        ax = pivot.plot(kind="bar", figsize=(9, 5))
        ax.set_ylabel("Median |editor Δβ|")
        ax.set_xlabel("Candidate set")
        ax.set_title("Candidate-set bystander responsiveness")
        plt.tight_layout()
        plt.savefig(figures / "candidate_set_comparison.png", dpi=200)
        plt.close()

    pair_effects = pd.concat([result.tcell_pair_effects, result.msc_pair_effects], ignore_index=True)
    candidates = set(result.direction_compatibility["cpg_id"])
    pair_effects = pair_effects[pair_effects["cpg_id"].isin(candidates)]
    if not pair_effects.empty:
        grouped = pair_effects.groupby(["cell_type", "pair_id"])["pair_delta_beta"].apply(lambda x: x.abs().median()).reset_index()
        labels = grouped["cell_type"] + ":" + grouped["pair_id"].astype(str)
        plt.figure(figsize=(9, 5))
        plt.bar(labels, grouped["pair_delta_beta"])
        plt.xticks(rotation=45, ha="right")
        plt.ylabel("Median candidate |Δβ|")
        plt.title("Pair/donor consistency")
        plt.tight_layout()
        plt.savefig(figures / "donor_consistency.png", dpi=200)
        plt.close()

    regional = result.regional_scores.head(20)
    if not regional.empty:
        plt.figure(figsize=(9, 6))
        labels = regional["region_id"].astype(str)
        plt.barh(labels[::-1], regional["regional_responsiveness_score"][::-1])
        plt.xlabel("Regional responsiveness score")
        plt.title("Top candidate regions")
        plt.tight_layout()
        plt.savefig(figures / "regional_editability.png", dpi=200)
        plt.close()

    nulls = result.null_tests
    subset = nulls[
        (nulls["candidate_set"] == "stable_core_3of3")
        & (nulls["metric"] == "median_abs_delta")
    ]
    if not subset.empty:
        cells = subset["cell_type"].drop_duplicates().tolist()
        plt.figure(figsize=(8, 5))
        data = [subset.loc[subset["cell_type"] == cell, "null_value"].dropna() for cell in cells]
        plt.boxplot(data, tick_labels=cells)
        for idx, cell in enumerate(cells, start=1):
            observed = subset.loc[subset["cell_type"] == cell, "observed"].iloc[0]
            plt.scatter([idx], [observed], marker="x", s=80)
        plt.ylabel("Median |editor Δβ|")
        plt.title("Stable-core matched-null comparison")
        plt.tight_layout()
        plt.savefig(figures / "matched_null_comparison.png", dpi=200)
        plt.close()
