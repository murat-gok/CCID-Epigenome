from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from ccid.o3_validation import (
    _comparison_metrics,
    _match_column,
    _normalise_scalar,
    infer_validation_stratum,
    read_validation_beta,
)


@dataclass(frozen=True)
class RobustValidationResult:
    successful_vs_failed: pd.DataFrame
    baseline_sensitivity_metrics: pd.DataFrame
    baseline_sensitivity_cpg: pd.DataFrame
    direction_balanced_null: pd.DataFrame
    lodo_summary: pd.DataFrame
    lodo_candidates: pd.DataFrame
    lodo_overlap: pd.DataFrame
    time_series_metrics: pd.DataFrame
    time_series_cpg_consistency: pd.DataFrame
    sample_manifest: pd.DataFrame
    summary: dict[str, object]


def _unique(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    output: list[str] = []
    for value in values:
        text = str(value)
        if text not in seen:
            seen.add(text)
            output.append(text)
    return output


def _require_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    required = {"sample_id", "donor_id", "day", "condition", "batch", "stage"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing: {sorted(missing)}")
    meta = metadata.copy()
    meta["sample_id"] = meta["sample_id"].astype(str)
    if meta["sample_id"].duplicated().any():
        raise ValueError("Metadata contains duplicate sample_id values")
    return meta


def _mask_stratum(
    metadata: pd.DataFrame,
    *,
    donor: str | None = None,
    day: object | None = None,
    batch: object | None = None,
    stage: object | None = None,
    condition: str | None = None,
) -> pd.Series:
    mask = pd.Series(True, index=metadata.index)
    if donor is not None:
        mask &= _match_column(metadata, "donor_id", donor)
    if day is not None:
        mask &= _match_column(metadata, "day", day)
    if batch is not None:
        mask &= _match_column(metadata, "batch", batch)
    if stage is not None:
        mask &= _match_column(metadata, "stage", stage)
    if condition is not None:
        mask &= _match_column(metadata, "condition", condition)
    return mask


def _ids(metadata: pd.DataFrame, mask: pd.Series) -> tuple[str, ...]:
    return tuple(metadata.loc[mask, "sample_id"].astype(str).tolist())


def _frozen_donor_samples(
    metadata: pd.DataFrame,
    *,
    donor: str,
    day: object,
    batch: object,
    stage: object,
    success_condition: str = "reprogrammed",
    baseline_condition: str = "baseline",
    failed_condition: str = "failed",
) -> dict[str, tuple[str, ...]]:
    success = _ids(
        metadata,
        _mask_stratum(
            metadata,
            donor=donor,
            day=day,
            batch=batch,
            stage=stage,
            condition=success_condition,
        ),
    )
    baseline = _ids(
        metadata,
        _mask_stratum(
            metadata,
            donor=donor,
            day=day,
            batch=batch,
            stage=stage,
            condition=baseline_condition,
        ),
    )
    failed = _ids(
        metadata,
        _mask_stratum(
            metadata,
            donor=donor,
            day=day,
            batch=batch,
            stage=stage,
            condition=failed_condition,
        ),
    )
    if not success:
        raise ValueError(f"No successful sample for donor={donor} in frozen stratum")
    if not baseline:
        raise ValueError(f"No baseline sample for donor={donor} in frozen stratum")
    return {"success": success, "baseline": baseline, "failed": failed}


def _mean(beta: pd.DataFrame, sample_ids: Sequence[str]) -> pd.Series:
    if not sample_ids:
        return pd.Series(np.nan, index=beta.index, dtype=float)
    return beta.loc[:, list(sample_ids)].mean(axis=1)


def _normalise_signature(signature: pd.DataFrame) -> pd.DataFrame:
    required = {"cpg_id", "target_delta_beta"}
    missing = required.difference(signature.columns)
    if missing:
        raise ValueError(f"Signature is missing: {sorted(missing)}")
    frame = signature.copy()
    frame["cpg_id"] = frame["cpg_id"].astype(str)
    frame = frame.drop_duplicates("cpg_id", keep="first")
    if "rank_score" in frame.columns:
        frame = frame.sort_values("rank_score", ascending=False, kind="stable")
    frame = frame.reset_index(drop=True)
    frame["candidate_rank"] = np.arange(1, len(frame) + 1)
    frame = frame.rename(columns={"target_delta_beta": "discovery_delta_beta"})
    if "weight" not in frame.columns:
        frame["weight"] = 1.0
    return frame


def _candidate_comparison(
    signature: pd.DataFrame,
    delta: pd.Series,
    *,
    comparison: str,
    genomewide_center: float,
    control_ids: Sequence[str],
) -> tuple[pd.DataFrame, dict[str, object]]:
    frame = signature.copy()
    frame["comparison"] = comparison
    frame["validation_delta_beta"] = frame["cpg_id"].map(delta)
    frame["genomewide_median_delta"] = float(genomewide_center)
    frame["centered_discovery_delta_beta"] = (
        frame["discovery_delta_beta"] - float(frame.attrs.get("discovery_center", 0.0))
    )
    frame["centered_validation_delta_beta"] = (
        frame["validation_delta_beta"] - float(genomewide_center)
    )
    frame["same_direction"] = (
        np.sign(frame["discovery_delta_beta"])
        == np.sign(frame["validation_delta_beta"])
    ) & frame["validation_delta_beta"].notna() & (frame["validation_delta_beta"] != 0)
    frame["centered_same_direction"] = (
        np.sign(frame["centered_discovery_delta_beta"])
        == np.sign(frame["centered_validation_delta_beta"])
    ) & frame["centered_validation_delta_beta"].notna() & (
        frame["centered_validation_delta_beta"] != 0
    )
    frame["effect_retention"] = frame["validation_delta_beta"].abs() / frame[
        "discovery_delta_beta"
    ].abs().clip(lower=1e-12)
    frame["control_sample_ids"] = ";".join(str(x) for x in control_ids)

    raw = _comparison_metrics(
        frame["discovery_delta_beta"].to_numpy(dtype=float),
        frame["validation_delta_beta"].to_numpy(dtype=float),
        frame["weight"].to_numpy(dtype=float),
    )
    centered = _comparison_metrics(
        frame["centered_discovery_delta_beta"].to_numpy(dtype=float),
        frame["centered_validation_delta_beta"].to_numpy(dtype=float),
        frame["weight"].to_numpy(dtype=float),
    )
    metrics: dict[str, object] = {
        "comparison": comparison,
        "control_sample_ids": ";".join(str(x) for x in control_ids),
        "n_control_samples": int(len(control_ids)),
        "genomewide_median_delta": float(genomewide_center),
        **{f"raw_{key}": value for key, value in raw.items()},
        **{f"centered_{key}": value for key, value in centered.items()},
    }
    return frame, metrics


def _direction_effect_bins(
    values: pd.Series,
    reference: pd.Series,
    n_bins: int,
) -> tuple[pd.Series, pd.Series]:
    combined = pd.concat([values.abs(), reference.abs()], ignore_index=True)
    finite = combined[np.isfinite(combined)]
    if finite.empty:
        raise ValueError("No finite effects available for null stratification")
    edges = np.unique(np.quantile(finite, np.linspace(0, 1, max(int(n_bins), 2) + 1)))
    if len(edges) < 3:
        edges = np.asarray([-np.inf, np.inf], dtype=float)
    else:
        edges[0] = -np.inf
        edges[-1] = np.inf
    return (
        pd.cut(values.abs(), bins=edges, labels=False, include_lowest=True),
        pd.cut(reference.abs(), bins=edges, labels=False, include_lowest=True),
    )


def direction_balanced_effect_matched_null(
    candidates: pd.DataFrame,
    all_loci: pd.DataFrame,
    validation_delta: pd.Series,
    *,
    n_iterations: int = 1000,
    n_bins: int = 10,
    seed: int = 20260721,
    global_centered: bool = False,
) -> pd.DataFrame:
    """Sample null CpGs matched on discovery direction and absolute-effect bin."""

    if n_iterations < 1:
        return pd.DataFrame()
    required = {"cpg_id", "target_delta_beta"}
    if not required.issubset(all_loci.columns):
        raise ValueError("all_loci must contain cpg_id and target_delta_beta")
    cand = candidates.copy()
    if "discovery_delta_beta" not in cand.columns:
        if "target_delta_beta" not in cand.columns:
            raise ValueError("candidates must contain discovery/target delta")
        cand = cand.rename(columns={"target_delta_beta": "discovery_delta_beta"})
    cand["cpg_id"] = cand["cpg_id"].astype(str)

    pool = all_loci.loc[:, ["cpg_id", "target_delta_beta"]].copy()
    pool["cpg_id"] = pool["cpg_id"].astype(str)
    pool["validation_delta_beta"] = pool["cpg_id"].map(validation_delta)
    pool = pool[~pool["cpg_id"].isin(set(cand["cpg_id"]))].copy()
    cand["validation_delta_beta"] = cand["cpg_id"].map(validation_delta)

    discovery_center = float(
        np.nanmedian(all_loci["target_delta_beta"].to_numpy(dtype=float))
    )
    validation_center = float(np.nanmedian(validation_delta.to_numpy(dtype=float)))
    if global_centered:
        cand["analysis_discovery"] = cand["discovery_delta_beta"] - discovery_center
        cand["analysis_validation"] = cand["validation_delta_beta"] - validation_center
        pool["analysis_discovery"] = pool["target_delta_beta"] - discovery_center
        pool["analysis_validation"] = pool["validation_delta_beta"] - validation_center
        analysis_scale = "global_centered"
    else:
        cand["analysis_discovery"] = cand["discovery_delta_beta"]
        cand["analysis_validation"] = cand["validation_delta_beta"]
        pool["analysis_discovery"] = pool["target_delta_beta"]
        pool["analysis_validation"] = pool["validation_delta_beta"]
        analysis_scale = "raw"

    cand = cand[
        np.isfinite(cand["analysis_discovery"])
        & np.isfinite(cand["analysis_validation"])
        & (cand["analysis_discovery"] != 0)
    ].copy()
    pool = pool[
        np.isfinite(pool["analysis_discovery"])
        & np.isfinite(pool["analysis_validation"])
        & (pool["analysis_discovery"] != 0)
    ].copy()
    if cand.empty or pool.empty:
        raise ValueError("Candidate or null pool is empty after finite-value filtering")

    cand["direction"] = np.where(cand["analysis_discovery"] > 0, "positive", "negative")
    pool["direction"] = np.where(pool["analysis_discovery"] > 0, "positive", "negative")
    cand_bins, pool_bins = _direction_effect_bins(
        cand["analysis_discovery"], pool["analysis_discovery"], n_bins
    )
    cand["effect_bin"] = cand_bins
    pool["effect_bin"] = pool_bins
    strata_counts = cand.groupby(["direction", "effect_bin"], dropna=False).size()
    grouped = {
        key: group for key, group in pool.groupby(["direction", "effect_bin"], dropna=False)
    }
    direction_groups = {key: group for key, group in pool.groupby("direction")}

    observed = _comparison_metrics(
        cand["analysis_discovery"].to_numpy(dtype=float),
        cand["analysis_validation"].to_numpy(dtype=float),
    )
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for iteration in range(int(n_iterations)):
        parts: list[pd.DataFrame] = []
        replacement_count = 0
        fallback_count = 0
        for (direction, effect_bin), count in strata_counts.items():
            group = grouped.get((direction, effect_bin))
            if group is None or group.empty:
                group = direction_groups.get(direction, pool)
                fallback_count += int(count)
            replace = len(group) < int(count)
            replacement_count += int(count) if replace else 0
            chosen = rng.choice(group.index.to_numpy(), size=int(count), replace=replace)
            parts.append(group.loc[chosen])
        sampled = pd.concat(parts, ignore_index=True)
        metrics = _comparison_metrics(
            sampled["analysis_discovery"].to_numpy(dtype=float),
            sampled["analysis_validation"].to_numpy(dtype=float),
        )
        rows.append(
            {
                "null_type": "direction_and_effect_matched_random_cpg_set",
                "analysis_scale": analysis_scale,
                "iteration": iteration + 1,
                "pool_size": int(len(pool)),
                "replacement_count": replacement_count,
                "fallback_count": fallback_count,
                "sign_concordance": metrics["sign_concordance"],
                "abs_spearman": abs(float(metrics["spearman_rho"]))
                if np.isfinite(metrics["spearman_rho"])
                else np.nan,
                "observed_sign_concordance": observed["sign_concordance"],
                "observed_abs_spearman": abs(float(observed["spearman_rho"]))
                if np.isfinite(observed["spearman_rho"])
                else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _rotation_signature(
    donor_deltas: pd.DataFrame,
    discovery_donors: Sequence[str],
    *,
    top_n: int,
    min_abs_delta: float,
) -> pd.DataFrame:
    effects = donor_deltas.loc[:, [str(x) for x in discovery_donors]].copy()
    values = effects.to_numpy(dtype=float, copy=True)
    valid_n = np.isfinite(values).sum(axis=1)
    mean_delta = np.nanmean(values, axis=1)
    sd = np.nanstd(values, axis=1, ddof=1)
    se = sd / np.sqrt(np.maximum(valid_n, 1))
    positive = np.nansum(values > 0, axis=1) / np.maximum(valid_n, 1)
    negative = np.nansum(values < 0, axis=1) / np.maximum(valid_n, 1)
    sign_consistency = np.maximum(positive, negative)
    table = pd.DataFrame(
        {
            "cpg_id": effects.index.astype(str),
            "discovery_delta_beta": mean_delta,
            "discovery_sd": sd,
            "discovery_se": se,
            "sign_consistency": sign_consistency,
            "n_discovery_donors": valid_n,
        }
    )
    table["passes_primary_filter"] = (
        (table["n_discovery_donors"] == len(discovery_donors))
        & (table["sign_consistency"] >= 1.0)
        & (table["discovery_delta_beta"].abs() >= min_abs_delta)
    )
    table["rank_score"] = (
        table["discovery_delta_beta"].abs()
        * table["sign_consistency"]
        * np.sqrt(table["n_discovery_donors"].clip(lower=1))
        / np.maximum(table["discovery_se"], 0.01)
    )
    primary = table[table["passes_primary_filter"]].sort_values(
        ["rank_score", "cpg_id"], ascending=[False, True], kind="stable"
    )
    if len(primary) >= top_n:
        selected = primary.head(top_n).copy()
        selected["rotation_fallback"] = False
    else:
        remaining = table[~table["cpg_id"].isin(set(primary["cpg_id"]))].sort_values(
            ["rank_score", "cpg_id"], ascending=[False, True], kind="stable"
        )
        selected = pd.concat([primary, remaining.head(max(top_n - len(primary), 0))])
        selected = selected.head(top_n).copy()
        selected["rotation_fallback"] = ~selected["passes_primary_filter"]
    selected["candidate_rank"] = np.arange(1, len(selected) + 1)
    return selected.reset_index(drop=True)


def leave_one_donor_out_analysis(
    donor_deltas: pd.DataFrame,
    *,
    donors: Sequence[str] = ("O1", "O2", "O3"),
    top_n: int = 500,
    min_abs_delta: float = 0.03,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    donors = [str(x) for x in donors]
    missing = [x for x in donors if x not in donor_deltas.columns]
    if missing:
        raise ValueError(f"Donor deltas missing: {missing}")
    summaries: list[dict[str, object]] = []
    candidate_frames: list[pd.DataFrame] = []
    candidate_sets: dict[str, set[str]] = {}
    for holdout in donors:
        discovery = [x for x in donors if x != holdout]
        selected = _rotation_signature(
            donor_deltas, discovery, top_n=top_n, min_abs_delta=min_abs_delta
        )
        selected["validation_delta_beta"] = selected["cpg_id"].map(donor_deltas[holdout])
        selected["same_direction"] = (
            np.sign(selected["discovery_delta_beta"])
            == np.sign(selected["validation_delta_beta"])
        ) & selected["validation_delta_beta"].notna() & (
            selected["validation_delta_beta"] != 0
        )
        selected["effect_retention"] = selected["validation_delta_beta"].abs() / selected[
            "discovery_delta_beta"
        ].abs().clip(lower=1e-12)
        selected["holdout_donor"] = holdout
        selected["discovery_donors"] = ";".join(discovery)
        raw = _comparison_metrics(
            selected["discovery_delta_beta"].to_numpy(dtype=float),
            selected["validation_delta_beta"].to_numpy(dtype=float),
        )
        discovery_center = float(
            np.nanmedian(donor_deltas.loc[:, discovery].mean(axis=1).to_numpy(dtype=float))
        )
        validation_center = float(np.nanmedian(donor_deltas[holdout].to_numpy(dtype=float)))
        centered = _comparison_metrics(
            (selected["discovery_delta_beta"] - discovery_center).to_numpy(dtype=float),
            (selected["validation_delta_beta"] - validation_center).to_numpy(dtype=float),
        )
        summaries.append(
            {
                "holdout_donor": holdout,
                "discovery_donors": ";".join(discovery),
                "n_candidates": int(len(selected)),
                "n_primary_filter": int(selected["passes_primary_filter"].sum()),
                "n_rotation_fallback": int(selected["rotation_fallback"].sum()),
                "discovery_global_median_delta": discovery_center,
                "validation_global_median_delta": validation_center,
                **{f"raw_{key}": value for key, value in raw.items()},
                **{f"centered_{key}": value for key, value in centered.items()},
            }
        )
        candidate_frames.append(selected)
        candidate_sets[holdout] = set(selected["cpg_id"])

    overlap_rows: list[dict[str, object]] = []
    for i, left in enumerate(donors):
        for right in donors[i + 1 :]:
            intersection = candidate_sets[left] & candidate_sets[right]
            union = candidate_sets[left] | candidate_sets[right]
            overlap_rows.append(
                {
                    "holdout_a": left,
                    "holdout_b": right,
                    "intersection": int(len(intersection)),
                    "union": int(len(union)),
                    "jaccard": float(len(intersection) / len(union)) if union else np.nan,
                }
            )
    membership: dict[str, int] = {}
    for values in candidate_sets.values():
        for cpg in values:
            membership[cpg] = membership.get(cpg, 0) + 1
    overlap_rows.append(
        {
            "holdout_a": "consensus",
            "holdout_b": "all_three_rotations",
            "intersection": int(sum(count == len(donors) for count in membership.values())),
            "union": int(len(membership)),
            "jaccard": np.nan,
        }
    )
    return (
        pd.DataFrame(summaries),
        pd.concat(candidate_frames, ignore_index=True),
        pd.DataFrame(overlap_rows),
    )


def _read_candidate_beta(
    beta_path: str | Path,
    sample_ids: Sequence[str],
    candidate_ids: Sequence[str],
    *,
    chunksize: int,
) -> pd.DataFrame:
    beta_path = Path(beta_path)
    sample_order = _unique(sample_ids)
    candidate_set = {str(x) for x in candidate_ids}
    header = pd.read_csv(beta_path, compression="infer", nrows=0)
    index_name = str(header.columns[0])
    available = set(str(x) for x in header.columns[1:])
    absent = [x for x in sample_order if x not in available]
    if absent:
        raise ValueError(f"Time-series samples absent from beta matrix: {absent[:5]}")
    pieces: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        beta_path,
        compression="infer",
        usecols=[index_name, *sample_order],
        index_col=0,
        chunksize=chunksize,
        low_memory=False,
    ):
        chunk.index = chunk.index.astype(str)
        subset = chunk.loc[chunk.index.isin(candidate_set)]
        if not subset.empty:
            pieces.append(subset.apply(pd.to_numeric, errors="coerce").astype(np.float32))
    if not pieces:
        raise ValueError("No frozen candidate CpGs were found in the beta matrix")
    output = pd.concat(pieces, axis=0)
    if output.index.duplicated().any():
        output = output.groupby(level=0, sort=False).mean()
    return output


def time_series_consistency_analysis(
    beta_path: str | Path,
    metadata: pd.DataFrame,
    signature: pd.DataFrame,
    *,
    donors: Sequence[str],
    batch: object,
    success_conditions: Sequence[str] = ("reprogrammed", "reprogrammed_intermediate"),
    failed_condition: str = "failed",
    baseline_condition: str = "baseline",
    chunksize: int = 10000,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    meta = _require_metadata(metadata)
    allowed_donors = {str(x).upper() for x in donors}
    conditions = {str(x).lower() for x in [*success_conditions, failed_condition]}
    sample_mask = (
        meta["donor_id"].astype(str).str.upper().isin(allowed_donors)
        & meta["condition"].astype(str).str.lower().isin(conditions)
        & _match_column(meta, "batch", batch)
        & meta["day"].notna()
    )
    candidates = meta.loc[sample_mask].copy()
    rows: list[dict[str, object]] = []
    required_samples: list[str] = []
    for sample in candidates.itertuples(index=False):
        baseline_mask = _mask_stratum(
            meta,
            donor=str(sample.donor_id),
            day=sample.day,
            batch=sample.batch,
            stage=sample.stage,
            condition=baseline_condition,
        )
        baseline_ids = _ids(meta, baseline_mask)
        if not baseline_ids:
            continue
        required_samples.append(str(sample.sample_id))
        required_samples.extend(baseline_ids)
        rows.append(
            {
                "sample_id": str(sample.sample_id),
                "baseline_sample_ids": ";".join(baseline_ids),
                "donor_id": str(sample.donor_id),
                "day": sample.day,
                "batch": sample.batch,
                "stage": sample.stage,
                "condition": str(sample.condition),
            }
        )
    if not rows:
        raise ValueError("No time-series sample had an exact matched baseline")
    manifest = pd.DataFrame(rows)
    frozen = _normalise_signature(signature)
    beta = _read_candidate_beta(
        beta_path,
        required_samples,
        frozen["cpg_id"].tolist(),
        chunksize=chunksize,
    )
    metric_rows: list[dict[str, object]] = []
    cpg_rows: list[pd.DataFrame] = []
    discovery = frozen.set_index("cpg_id")["discovery_delta_beta"]
    weights = frozen.set_index("cpg_id")["weight"]
    for row in manifest.itertuples(index=False):
        baseline_ids = str(row.baseline_sample_ids).split(";")
        delta = beta[str(row.sample_id)] - beta[baseline_ids].mean(axis=1)
        aligned_discovery = discovery.reindex(beta.index)
        aligned_weights = weights.reindex(beta.index).fillna(1.0)
        metrics = _comparison_metrics(
            aligned_discovery.to_numpy(dtype=float),
            delta.to_numpy(dtype=float),
            aligned_weights.to_numpy(dtype=float),
        )
        metric_rows.append(
            {
                **row._asdict(),
                **metrics,
                "median_delta_beta": float(np.nanmedian(delta.to_numpy(dtype=float))),
            }
        )
        detail = pd.DataFrame(
            {
                "cpg_id": beta.index.astype(str),
                "sample_id": str(row.sample_id),
                "donor_id": str(row.donor_id),
                "day": row.day,
                "stage": row.stage,
                "condition": str(row.condition),
                "discovery_delta_beta": aligned_discovery.to_numpy(dtype=float),
                "timepoint_delta_beta": delta.to_numpy(dtype=float),
            }
        )
        detail["same_direction"] = (
            np.sign(detail["discovery_delta_beta"])
            == np.sign(detail["timepoint_delta_beta"])
        ) & detail["timepoint_delta_beta"].notna() & (detail["timepoint_delta_beta"] != 0)
        cpg_rows.append(detail)
    detail_all = pd.concat(cpg_rows, ignore_index=True)

    consistency_rows: list[dict[str, object]] = []
    success_names = {str(x).lower() for x in success_conditions}
    for cpg, group in detail_all.groupby("cpg_id", sort=False):
        success = group[group["condition"].str.lower().isin(success_names)]
        failed = group[group["condition"].str.lower() == str(failed_condition).lower()]
        o3_success = success[success["donor_id"].str.upper() == "O3"]
        supported_days = sorted(
            pd.to_numeric(
                o3_success.loc[o3_success["same_direction"], "day"], errors="coerce"
            ).dropna().unique().tolist()
        )
        consistency_rows.append(
            {
                "cpg_id": cpg,
                "n_success_timepoints": int(len(success)),
                "success_sign_consistency": float(success["same_direction"].mean())
                if len(success)
                else np.nan,
                "n_failed_timepoints": int(len(failed)),
                "failed_sign_consistency": float(failed["same_direction"].mean())
                if len(failed)
                else np.nan,
                "n_o3_success_timepoints": int(len(o3_success)),
                "o3_success_sign_consistency": float(o3_success["same_direction"].mean())
                if len(o3_success)
                else np.nan,
                "earliest_o3_supported_day": supported_days[0] if supported_days else np.nan,
                "latest_o3_supported_day": supported_days[-1] if supported_days else np.nan,
            }
        )
    return pd.DataFrame(metric_rows), pd.DataFrame(consistency_rows), manifest


def run_robust_validation(
    *,
    beta_path: str | Path,
    metadata: pd.DataFrame,
    signature: pd.DataFrame,
    discovery_pair_table: pd.DataFrame,
    all_loci: pd.DataFrame,
    validation_donor: str = "O3",
    donors: Sequence[str] = ("O1", "O2", "O3"),
    success_condition: str = "reprogrammed",
    baseline_condition: str = "baseline",
    failed_condition: str = "failed",
    lodo_top_n: int = 500,
    lodo_min_abs_delta: float = 0.03,
    null_iterations: int = 1000,
    null_bins: int = 10,
    seed: int = 20260721,
    chunksize: int = 10000,
) -> RobustValidationResult:
    meta = _require_metadata(metadata)
    frozen_stratum = infer_validation_stratum(discovery_pair_table)
    day = frozen_stratum["day"]
    batch = frozen_stratum["batch"]
    stage = frozen_stratum["stage"]
    frozen = _normalise_signature(signature)
    if not {"cpg_id", "target_delta_beta"}.issubset(all_loci.columns):
        raise ValueError("all_loci must contain cpg_id and target_delta_beta")
    all_loci = all_loci.copy()
    all_loci["cpg_id"] = all_loci["cpg_id"].astype(str)
    discovery_center = float(np.nanmedian(all_loci["target_delta_beta"].to_numpy(dtype=float)))
    frozen.attrs["discovery_center"] = discovery_center

    donor_samples: dict[str, dict[str, tuple[str, ...]]] = {}
    required_core: list[str] = []
    for donor in donors:
        donor_samples[str(donor)] = _frozen_donor_samples(
            meta,
            donor=str(donor),
            day=day,
            batch=batch,
            stage=stage,
            success_condition=success_condition,
            baseline_condition=baseline_condition,
            failed_condition=failed_condition,
        )
        required_core.extend(donor_samples[str(donor)]["success"])
        required_core.extend(donor_samples[str(donor)]["baseline"])
        required_core.extend(donor_samples[str(donor)]["failed"])

    validation_baselines = _ids(
        meta,
        _mask_stratum(
            meta,
            donor=validation_donor,
            batch=batch,
            stage=stage,
            condition=baseline_condition,
        ),
    )
    required_core.extend(validation_baselines)
    beta, _, _ = read_validation_beta(beta_path, _unique(required_core), chunksize=chunksize)

    donor_delta_dict: dict[str, pd.Series] = {}
    for donor, sample_map in donor_samples.items():
        donor_delta_dict[donor] = _mean(beta, sample_map["success"]) - _mean(
            beta, sample_map["baseline"]
        )
    donor_deltas = pd.DataFrame(donor_delta_dict, index=beta.index)

    validation_map = donor_samples[str(validation_donor)]
    validation_mean = _mean(beta, validation_map["success"])
    exact_baseline = _mean(beta, validation_map["baseline"])
    failed_mean = _mean(beta, validation_map["failed"])
    other_baselines = [x for x in validation_baselines if x not in validation_map["baseline"]]
    all_baselines = list(validation_baselines)

    comparisons: list[tuple[str, pd.Series, Sequence[str]]] = [
        ("successful_vs_exact_baseline", validation_mean - exact_baseline, validation_map["baseline"]),
    ]
    if other_baselines:
        comparisons.append(
            (
                "successful_vs_other_day_baseline_mean",
                validation_mean - _mean(beta, other_baselines),
                other_baselines,
            )
        )
    if all_baselines:
        comparisons.append(
            (
                "successful_vs_all_same_donor_baseline_mean",
                validation_mean - _mean(beta, all_baselines),
                all_baselines,
            )
        )
    if validation_map["failed"]:
        comparisons.append(
            ("successful_vs_failed_same_stratum", validation_mean - failed_mean, validation_map["failed"])
        )

    sensitivity_frames: list[pd.DataFrame] = []
    sensitivity_metrics: list[dict[str, object]] = []
    comparison_deltas: dict[str, pd.Series] = {}
    for name, delta, controls in comparisons:
        comparison_deltas[name] = delta
        center = float(np.nanmedian(delta.to_numpy(dtype=float)))
        frame, metrics = _candidate_comparison(
            frozen,
            delta,
            comparison=name,
            genomewide_center=center,
            control_ids=controls,
        )
        sensitivity_frames.append(frame)
        sensitivity_metrics.append(metrics)
    sensitivity_cpg = pd.concat(sensitivity_frames, ignore_index=True)
    sensitivity_metric_frame = pd.DataFrame(sensitivity_metrics)
    successful_vs_failed = sensitivity_cpg[
        sensitivity_cpg["comparison"] == "successful_vs_failed_same_stratum"
    ].copy()

    exact_delta = comparison_deltas["successful_vs_exact_baseline"]
    null_raw = direction_balanced_effect_matched_null(
        frozen,
        all_loci,
        exact_delta,
        n_iterations=null_iterations,
        n_bins=null_bins,
        seed=seed,
        global_centered=False,
    )
    null_centered = direction_balanced_effect_matched_null(
        frozen,
        all_loci,
        exact_delta,
        n_iterations=null_iterations,
        n_bins=null_bins,
        seed=seed + 1,
        global_centered=True,
    )
    null_frame = pd.concat([null_raw, null_centered], ignore_index=True)

    lodo_summary, lodo_candidates, lodo_overlap = leave_one_donor_out_analysis(
        donor_deltas,
        donors=donors,
        top_n=lodo_top_n,
        min_abs_delta=lodo_min_abs_delta,
    )
    time_metrics, time_consistency, time_manifest = time_series_consistency_analysis(
        beta_path,
        meta,
        frozen.rename(columns={"discovery_delta_beta": "target_delta_beta"}),
        donors=donors,
        batch=batch,
        success_conditions=("reprogrammed", "reprogrammed_intermediate"),
        failed_condition=failed_condition,
        baseline_condition=baseline_condition,
        chunksize=chunksize,
    )

    empirical_null_summary: dict[str, object] = {}
    for scale, group in null_frame.groupby("analysis_scale"):
        observed_sign = float(group["observed_sign_concordance"].iloc[0])
        observed_rho = float(group["observed_abs_spearman"].iloc[0])
        empirical_null_summary[str(scale)] = {
            "n_iterations": int(len(group)),
            "observed_sign_concordance": observed_sign,
            "null_mean_sign_concordance": float(group["sign_concordance"].mean()),
            "empirical_p_sign": float(
                (1 + (group["sign_concordance"] >= observed_sign).sum()) / (len(group) + 1)
            ),
            "observed_abs_spearman": observed_rho,
            "null_mean_abs_spearman": float(group["abs_spearman"].mean()),
            "empirical_p_spearman": float(
                (1 + (group["abs_spearman"] >= observed_rho).sum()) / (len(group) + 1)
            ),
        }

    summary: dict[str, object] = {
        "analysis_status": "prespecified_robustness_and_sensitivity_analysis",
        "candidate_list_frozen_before_validation": True,
        "frozen_stratum": {
            "day": _normalise_scalar(day),
            "batch": _normalise_scalar(batch),
            "stage": _normalise_scalar(stage),
        },
        "validation_donor": str(validation_donor),
        "donors": [str(x) for x in donors],
        "n_frozen_candidates": int(len(frozen)),
        "discovery_global_median_delta": discovery_center,
        "successful_vs_failed_available": bool(len(successful_vs_failed)),
        "baseline_comparisons": sensitivity_metric_frame.to_dict(orient="records"),
        "direction_balanced_null_summary": empirical_null_summary,
        "lodo_summary": lodo_summary.to_dict(orient="records"),
        "lodo_overlap": lodo_overlap.to_dict(orient="records"),
        "n_time_series_comparisons": int(len(time_metrics)),
        "time_series_conditions": sorted(time_metrics["condition"].astype(str).unique().tolist()),
        "interpretation_note": (
            "The O1/O2 frozen candidate list remains unchanged. LODO rotations are "
            "separate sensitivity analyses that rebuild candidates using only the two "
            "non-held-out donors. Time-series and baseline analyses are descriptive." 
        ),
    }
    return RobustValidationResult(
        successful_vs_failed=successful_vs_failed,
        baseline_sensitivity_metrics=sensitivity_metric_frame,
        baseline_sensitivity_cpg=sensitivity_cpg,
        direction_balanced_null=null_frame,
        lodo_summary=lodo_summary,
        lodo_candidates=lodo_candidates,
        lodo_overlap=lodo_overlap,
        time_series_metrics=time_metrics,
        time_series_cpg_consistency=time_consistency,
        sample_manifest=time_manifest,
        summary=summary,
    )
