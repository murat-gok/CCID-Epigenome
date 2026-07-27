from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy import stats


@dataclass(frozen=True)
class ValidationPair:
    validation_sample_ids: tuple[str, ...]
    baseline_sample_ids: tuple[str, ...]
    wrong_day_baseline_ids: tuple[str, ...]
    cross_donor_baseline_ids: tuple[str, ...]
    failed_sample_ids: tuple[str, ...]
    donor_id: str
    day: int | float
    batch: str
    stage: str


@dataclass(frozen=True)
class O3ValidationResult:
    cpg_validation: pd.DataFrame
    summary: dict[str, object]
    null_tests: pd.DataFrame
    sample_qc: pd.DataFrame
    topk_metrics: pd.DataFrame
    sample_values: dict[str, np.ndarray]


def _normalise_scalar(value: object) -> object:
    if pd.isna(value):
        return None
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return int(value)
    if isinstance(value, str):
        return value.strip().lower()
    return value


def _match_column(frame: pd.DataFrame, column: str, value: object) -> pd.Series:
    return frame[column].map(_normalise_scalar) == _normalise_scalar(value)


def infer_validation_stratum(
    discovery_pair_table: pd.DataFrame,
    *,
    match_columns: Sequence[str] = ("day", "batch", "stage"),
) -> dict[str, object]:
    """Infer the frozen validation design stratum from discovery matched pairs.

    A single held-out validation analysis must reproduce the exact discovery
    day/batch/stage stratum. Multiple strata are rejected rather than silently
    pooled because pooling would change the prespecified validation target.
    """

    missing = set(match_columns).difference(discovery_pair_table.columns)
    if missing:
        raise ValueError(
            f"Discovery pair table lacks frozen design columns: {sorted(missing)}"
        )
    strata = discovery_pair_table.loc[:, list(match_columns)].drop_duplicates()
    if len(strata) != 1:
        raise ValueError(
            "O3 validation currently requires exactly one frozen discovery stratum; "
            f"found {len(strata)}"
        )
    row = strata.iloc[0]
    return {column: row[column] for column in match_columns}


def select_validation_pair(
    metadata: pd.DataFrame,
    *,
    validation_donor: str = "O3",
    day: int | float = 10,
    batch: str = "exp2",
    stage: str = "endpoint",
    validation_condition: str = "reprogrammed",
    baseline_condition: str = "baseline",
    discovery_donors: Sequence[str] = ("O1", "O2"),
) -> ValidationPair:
    """Select a held-out donor pair and diagnostic negative-control samples.

    The primary pair is exact on donor, day, batch, and stage.  Wrong-day,
    cross-donor, and failed-reprogramming samples are returned only for diagnostic
    null analyses and never influence candidate selection or tier thresholds.
    """

    required = {"sample_id", "donor_id", "day", "condition", "batch", "stage"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing: {sorted(missing)}")

    meta = metadata.copy()
    meta["sample_id"] = meta["sample_id"].astype(str)
    donor_mask = _match_column(meta, "donor_id", validation_donor)
    day_mask = _match_column(meta, "day", day)
    batch_mask = _match_column(meta, "batch", batch)
    stage_mask = _match_column(meta, "stage", stage)

    validation_mask = (
        donor_mask
        & day_mask
        & batch_mask
        & stage_mask
        & _match_column(meta, "condition", validation_condition)
    )
    baseline_mask = (
        donor_mask
        & day_mask
        & batch_mask
        & stage_mask
        & _match_column(meta, "condition", baseline_condition)
    )
    validation_ids = tuple(meta.loc[validation_mask, "sample_id"].tolist())
    baseline_ids = tuple(meta.loc[baseline_mask, "sample_id"].tolist())
    if not validation_ids:
        raise ValueError(
            "No held-out validation sample matched donor/day/batch/stage/condition"
        )
    if not baseline_ids:
        raise ValueError(
            "No held-out baseline sample matched donor/day/batch/stage/condition"
        )

    wrong_day_mask = (
        donor_mask
        & ~day_mask
        & batch_mask
        & stage_mask
        & _match_column(meta, "condition", baseline_condition)
    )
    wrong_day_ids = tuple(meta.loc[wrong_day_mask, "sample_id"].tolist())

    allowed_discovery = {str(x).upper() for x in discovery_donors}
    cross_donor_mask = (
        meta["donor_id"].astype(str).str.upper().isin(allowed_discovery)
        & day_mask
        & batch_mask
        & stage_mask
        & _match_column(meta, "condition", baseline_condition)
    )
    cross_donor_ids = tuple(meta.loc[cross_donor_mask, "sample_id"].tolist())

    failed_mask = (
        donor_mask
        & day_mask
        & batch_mask
        & stage_mask
        & _match_column(meta, "condition", "failed")
    )
    failed_ids = tuple(meta.loc[failed_mask, "sample_id"].tolist())

    return ValidationPair(
        validation_sample_ids=validation_ids,
        baseline_sample_ids=baseline_ids,
        wrong_day_baseline_ids=wrong_day_ids,
        cross_donor_baseline_ids=cross_donor_ids,
        failed_sample_ids=failed_ids,
        donor_id=str(validation_donor),
        day=day,
        batch=str(batch),
        stage=str(stage),
    )


def _unique_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    output: list[str] = []
    for value in values:
        value = str(value)
        if value not in seen:
            seen.add(value)
            output.append(value)
    return output


def read_validation_beta(
    beta_path: str | Path,
    sample_ids: Sequence[str],
    *,
    chunksize: int = 10000,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, np.ndarray]]:
    """Read only required samples from an EPIC-scale processed beta matrix.

    Returns a locus-by-sample beta frame, per-sample QC, and finite value arrays
    used for distribution plots.  The full 96-sample matrix is never loaded.
    """

    beta_path = Path(beta_path)
    if not beta_path.exists():
        raise FileNotFoundError(beta_path)
    if chunksize < 1:
        raise ValueError("chunksize must be at least 1")

    sample_order = _unique_order(sample_ids)
    header = pd.read_csv(beta_path, compression="infer", nrows=0)
    index_name = str(header.columns[0])
    available = {str(x) for x in header.columns[1:]}
    absent = [x for x in sample_order if x not in available]
    if absent:
        raise ValueError(f"Required samples absent from beta matrix: {absent[:5]}")

    chunks: list[pd.DataFrame] = []
    values_by_sample: dict[str, list[np.ndarray]] = {sample: [] for sample in sample_order}
    for chunk in pd.read_csv(
        beta_path,
        compression="infer",
        usecols=[index_name, *sample_order],
        index_col=0,
        chunksize=chunksize,
        low_memory=False,
    ):
        chunk.index = chunk.index.astype(str)
        chunk = chunk.apply(pd.to_numeric, errors="coerce")
        chunks.append(chunk.astype(np.float32))
        for sample in sample_order:
            finite = chunk[sample].to_numpy(dtype=float, copy=False)
            values_by_sample[sample].append(finite[np.isfinite(finite)])

    if not chunks:
        raise ValueError("Beta matrix contained no rows")
    beta = pd.concat(chunks, axis=0)
    if beta.index.duplicated().any():
        beta = beta.groupby(level=0, sort=False).mean()

    sample_values = {
        sample: np.concatenate(parts) if parts else np.asarray([], dtype=float)
        for sample, parts in values_by_sample.items()
    }
    qc_rows: list[dict[str, object]] = []
    total = int(len(beta))
    for sample in sample_order:
        values = beta[sample].to_numpy(dtype=float, copy=False)
        finite = values[np.isfinite(values)]
        qc_rows.append(
            {
                "sample_id": sample,
                "total_loci": total,
                "valid_loci": int(len(finite)),
                "valid_fraction": float(len(finite) / max(total, 1)),
                "mean_beta": float(np.mean(finite)) if len(finite) else np.nan,
                "median_beta": float(np.median(finite)) if len(finite) else np.nan,
                "sd_beta": float(np.std(finite, ddof=1)) if len(finite) > 1 else np.nan,
                "q05_beta": float(np.quantile(finite, 0.05)) if len(finite) else np.nan,
                "q95_beta": float(np.quantile(finite, 0.95)) if len(finite) else np.nan,
            }
        )
    return beta, pd.DataFrame(qc_rows), sample_values


def _safe_corr(x: np.ndarray, y: np.ndarray, method: str) -> tuple[float, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    x_valid = x[valid]
    y_valid = y[valid]
    if len(x_valid) < 3 or np.nanstd(x_valid) == 0 or np.nanstd(y_valid) == 0:
        return np.nan, np.nan
    if method == "pearson":
        result = stats.pearsonr(x_valid, y_valid)
    elif method == "spearman":
        result = stats.spearmanr(x_valid, y_valid)
    else:
        raise ValueError(method)
    return float(result.statistic), float(result.pvalue)


def _comparison_metrics(
    discovery: np.ndarray,
    validation: np.ndarray,
    weights: np.ndarray | None = None,
) -> dict[str, float | int]:
    valid = np.isfinite(discovery) & np.isfinite(validation) & (discovery != 0)
    d = discovery[valid]
    v = validation[valid]
    n = int(len(d))
    if n == 0:
        return {
            "n_valid": 0,
            "sign_concordance": np.nan,
            "weighted_sign_concordance": np.nan,
            "sign_binomial_p": np.nan,
            "pearson_r": np.nan,
            "pearson_p": np.nan,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "mae": np.nan,
            "rmse": np.nan,
            "median_abs_validation_delta": np.nan,
            "median_effect_retention": np.nan,
            "slope_through_origin": np.nan,
        }

    same = np.sign(d) == np.sign(v)
    sign_count = int(same.sum())
    binom_p = float(stats.binomtest(sign_count, n, p=0.5, alternative="greater").pvalue)
    if weights is None:
        weighted = float(np.mean(same))
    else:
        w = np.asarray(weights, dtype=float)[valid]
        w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
        weighted = float(np.average(same.astype(float), weights=w)) if w.sum() > 0 else float(np.mean(same))
    pearson_r, pearson_p = _safe_corr(d, v, "pearson")
    spearman_r, spearman_p = _safe_corr(d, v, "spearman")
    retention = np.abs(v) / np.maximum(np.abs(d), 1e-12)
    denominator = float(np.dot(d, d))
    slope = float(np.dot(d, v) / denominator) if denominator > 0 else np.nan
    return {
        "n_valid": n,
        "sign_concordance": float(np.mean(same)),
        "weighted_sign_concordance": weighted,
        "sign_binomial_p": binom_p,
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_rho": spearman_r,
        "spearman_p": spearman_p,
        "mae": float(np.mean(np.abs(v - d))),
        "rmse": float(np.sqrt(np.mean((v - d) ** 2))),
        "median_abs_validation_delta": float(np.median(np.abs(v))),
        "median_effect_retention": float(np.median(retention)),
        "slope_through_origin": slope,
    }


def _assign_tiers(
    frame: pd.DataFrame,
    *,
    tier1_min_abs_delta: float,
    tier1_min_retention: float,
    tier1_max_retention: float,
    tier2_min_abs_delta: float,
    tier2_min_retention: float,
) -> pd.Series:
    same = frame["same_direction"].fillna(False)
    abs_delta = frame["abs_validation_delta"]
    retention = frame["effect_retention"]
    tier1 = (
        same
        & (abs_delta >= tier1_min_abs_delta)
        & (retention >= tier1_min_retention)
        & (retention <= tier1_max_retention)
    )
    tier2 = (
        same
        & ~tier1
        & (abs_delta >= tier2_min_abs_delta)
        & (retention >= tier2_min_retention)
    )
    output = pd.Series("Tier3", index=frame.index, dtype="object")
    output.loc[tier2] = "Tier2"
    output.loc[tier1] = "Tier1"
    return output


def _topk_metrics(frame: pd.DataFrame, top_k: Sequence[int]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for k in sorted({int(x) for x in top_k if int(x) > 0}):
        subset = frame.head(min(k, len(frame)))
        metrics = _comparison_metrics(
            subset["discovery_delta_beta"].to_numpy(dtype=float),
            subset["validation_delta_beta"].to_numpy(dtype=float),
            subset["weight"].to_numpy(dtype=float) if "weight" in subset else None,
        )
        rows.append({"top_k": min(k, len(frame)), **metrics})
    return pd.DataFrame(rows).drop_duplicates("top_k", keep="first")


def _permutation_null(
    discovery: np.ndarray,
    validation: np.ndarray,
    *,
    n_permutations: int,
    seed: int,
) -> dict[str, object]:
    valid = np.isfinite(discovery) & np.isfinite(validation) & (discovery != 0)
    d = discovery[valid]
    v = validation[valid]
    observed = _comparison_metrics(d, v)
    if len(d) < 3 or n_permutations < 1:
        return {
            "null_type": "cpg_label_permutation",
            "n_iterations": int(n_permutations),
            "observed_sign_concordance": observed["sign_concordance"],
            "null_mean_sign_concordance": np.nan,
            "empirical_p_sign": np.nan,
            "observed_abs_spearman": abs(float(observed["spearman_rho"])) if np.isfinite(observed["spearman_rho"]) else np.nan,
            "null_mean_abs_spearman": np.nan,
            "empirical_p_spearman": np.nan,
        }
    rng = np.random.default_rng(seed)
    sign_values = np.empty(n_permutations, dtype=float)
    rho_values = np.empty(n_permutations, dtype=float)
    for i in range(n_permutations):
        shuffled = rng.permutation(v)
        sign_values[i] = np.mean(np.sign(d) == np.sign(shuffled))
        rho, _ = _safe_corr(d, shuffled, "spearman")
        rho_values[i] = abs(rho) if np.isfinite(rho) else 0.0
    observed_sign = float(observed["sign_concordance"])
    observed_rho = abs(float(observed["spearman_rho"])) if np.isfinite(observed["spearman_rho"]) else 0.0
    return {
        "null_type": "cpg_label_permutation",
        "n_iterations": int(n_permutations),
        "observed_sign_concordance": observed_sign,
        "null_mean_sign_concordance": float(sign_values.mean()),
        "null_sd_sign_concordance": float(sign_values.std(ddof=1)) if n_permutations > 1 else 0.0,
        "empirical_p_sign": float((1 + np.sum(sign_values >= observed_sign)) / (n_permutations + 1)),
        "observed_abs_spearman": observed_rho,
        "null_mean_abs_spearman": float(rho_values.mean()),
        "null_sd_abs_spearman": float(rho_values.std(ddof=1)) if n_permutations > 1 else 0.0,
        "empirical_p_spearman": float((1 + np.sum(rho_values >= observed_rho)) / (n_permutations + 1)),
    }


def _effect_matched_random_null(
    candidates: pd.DataFrame,
    all_loci: pd.DataFrame,
    o3_delta: pd.Series,
    *,
    n_iterations: int,
    n_bins: int,
    seed: int,
) -> dict[str, object]:
    required = {"cpg_id", "target_delta_beta"}
    if not required.issubset(all_loci.columns):
        raise ValueError("All-loci file must contain cpg_id and target_delta_beta")
    pool = all_loci.loc[:, ["cpg_id", "target_delta_beta"]].copy()
    pool["cpg_id"] = pool["cpg_id"].astype(str)
    pool["validation_delta_beta"] = pool["cpg_id"].map(o3_delta)
    candidate_ids = set(candidates["cpg_id"].astype(str))
    pool = pool[~pool["cpg_id"].isin(candidate_ids)]
    pool = pool[
        np.isfinite(pool["target_delta_beta"])
        & np.isfinite(pool["validation_delta_beta"])
        & (pool["target_delta_beta"] != 0)
    ].copy()
    if len(pool) < len(candidates) or n_iterations < 1:
        return {
            "null_type": "effect_matched_random_cpg_sets",
            "n_iterations": int(n_iterations),
            "pool_size": int(len(pool)),
            "observed_sign_concordance": np.nan,
            "null_mean_sign_concordance": np.nan,
            "empirical_p_sign": np.nan,
            "observed_abs_spearman": np.nan,
            "null_mean_abs_spearman": np.nan,
            "empirical_p_spearman": np.nan,
        }

    candidate_valid = candidates[
        np.isfinite(candidates["discovery_delta_beta"])
        & np.isfinite(candidates["validation_delta_beta"])
    ].copy()
    observed = _comparison_metrics(
        candidate_valid["discovery_delta_beta"].to_numpy(dtype=float),
        candidate_valid["validation_delta_beta"].to_numpy(dtype=float),
    )

    combined_abs = pd.concat(
        [pool["target_delta_beta"].abs(), candidate_valid["discovery_delta_beta"].abs()],
        ignore_index=True,
    )
    quantiles = np.unique(np.quantile(combined_abs, np.linspace(0, 1, max(n_bins, 2) + 1)))
    if len(quantiles) < 3:
        quantiles = np.array([-np.inf, np.inf])
    else:
        quantiles[0] = -np.inf
        quantiles[-1] = np.inf
    pool["effect_bin"] = pd.cut(pool["target_delta_beta"].abs(), bins=quantiles, labels=False, include_lowest=True)
    candidate_valid["effect_bin"] = pd.cut(
        candidate_valid["discovery_delta_beta"].abs(), bins=quantiles, labels=False, include_lowest=True
    )
    bin_counts = candidate_valid["effect_bin"].value_counts(dropna=False).to_dict()
    rng = np.random.default_rng(seed)
    sign_values = np.empty(n_iterations, dtype=float)
    rho_values = np.empty(n_iterations, dtype=float)
    grouped = {key: group for key, group in pool.groupby("effect_bin", dropna=False)}
    for i in range(n_iterations):
        samples: list[pd.DataFrame] = []
        for key, count in bin_counts.items():
            group = grouped.get(key)
            if group is None or group.empty:
                group = pool
            indices = rng.choice(group.index.to_numpy(), size=int(count), replace=len(group) < int(count))
            samples.append(group.loc[indices])
        sampled = pd.concat(samples, ignore_index=True)
        metrics = _comparison_metrics(
            sampled["target_delta_beta"].to_numpy(dtype=float),
            sampled["validation_delta_beta"].to_numpy(dtype=float),
        )
        sign_values[i] = float(metrics["sign_concordance"])
        rho = float(metrics["spearman_rho"])
        rho_values[i] = abs(rho) if np.isfinite(rho) else 0.0

    observed_sign = float(observed["sign_concordance"])
    observed_rho = abs(float(observed["spearman_rho"])) if np.isfinite(observed["spearman_rho"]) else 0.0
    return {
        "null_type": "effect_matched_random_cpg_sets",
        "n_iterations": int(n_iterations),
        "pool_size": int(len(pool)),
        "observed_sign_concordance": observed_sign,
        "null_mean_sign_concordance": float(sign_values.mean()),
        "null_sd_sign_concordance": float(sign_values.std(ddof=1)) if n_iterations > 1 else 0.0,
        "empirical_p_sign": float((1 + np.sum(sign_values >= observed_sign)) / (n_iterations + 1)),
        "observed_abs_spearman": observed_rho,
        "null_mean_abs_spearman": float(rho_values.mean()),
        "null_sd_abs_spearman": float(rho_values.std(ddof=1)) if n_iterations > 1 else 0.0,
        "empirical_p_spearman": float((1 + np.sum(rho_values >= observed_rho)) / (n_iterations + 1)),
    }


def validate_o3_signature(
    beta_path: str | Path,
    metadata: pd.DataFrame,
    signature: pd.DataFrame,
    *,
    validation_donor: str = "O3",
    day: int | float = 10,
    batch: str = "exp2",
    stage: str = "endpoint",
    validation_condition: str = "reprogrammed",
    baseline_condition: str = "baseline",
    discovery_donors: Sequence[str] = ("O1", "O2"),
    top_k: Sequence[int] = (50, 100, 250, 500),
    tier1_min_abs_delta: float = 0.03,
    tier1_min_retention: float = 0.25,
    tier1_max_retention: float = 2.0,
    tier2_min_abs_delta: float = 0.01,
    tier2_min_retention: float = 0.10,
    n_permutations: int = 1000,
    random_seed: int = 2026,
    discovery_all_loci: pd.DataFrame | None = None,
    n_random_sets: int = 1000,
    random_effect_bins: int = 10,
    chunksize: int = 10000,
) -> O3ValidationResult:
    """Validate a frozen O1/O2 signature in the held-out O3 donor.

    Candidate selection and thresholds are never re-estimated from O3. Tier labels
    are descriptive prioritisation categories because O3 contributes one independent
    donor pair, not a confirmatory cohort.
    """

    required_signature = {"cpg_id", "target_delta_beta"}
    missing = required_signature.difference(signature.columns)
    if missing:
        raise ValueError(f"Signature is missing: {sorted(missing)}")
    if not (0 <= tier2_min_abs_delta <= tier1_min_abs_delta):
        raise ValueError("Tier absolute-delta thresholds are inconsistent")
    if not (0 <= tier2_min_retention <= tier1_min_retention):
        raise ValueError("Tier retention thresholds are inconsistent")
    if tier1_max_retention < tier1_min_retention:
        raise ValueError("tier1_max_retention must be >= tier1_min_retention")

    pair = select_validation_pair(
        metadata,
        validation_donor=validation_donor,
        day=day,
        batch=batch,
        stage=stage,
        validation_condition=validation_condition,
        baseline_condition=baseline_condition,
        discovery_donors=discovery_donors,
    )
    required_samples = _unique_order(
        [
            *pair.validation_sample_ids,
            *pair.baseline_sample_ids,
            *pair.wrong_day_baseline_ids,
            *pair.cross_donor_baseline_ids,
            *pair.failed_sample_ids,
        ]
    )
    beta, sample_qc, sample_values = read_validation_beta(
        beta_path, required_samples, chunksize=chunksize
    )
    validation_mean = beta.loc[:, list(pair.validation_sample_ids)].mean(axis=1)
    baseline_mean = beta.loc[:, list(pair.baseline_sample_ids)].mean(axis=1)
    o3_delta = validation_mean - baseline_mean

    frozen = signature.copy()
    frozen["cpg_id"] = frozen["cpg_id"].astype(str)
    frozen = frozen.drop_duplicates("cpg_id", keep="first")
    order_column = "rank_score" if "rank_score" in frozen.columns else None
    if order_column is not None:
        frozen = frozen.sort_values(order_column, ascending=False, kind="stable")
    frozen = frozen.reset_index(drop=True)
    frozen["candidate_rank"] = np.arange(1, len(frozen) + 1)
    frozen = frozen.rename(columns={"target_delta_beta": "discovery_delta_beta"})
    frozen["validation_delta_beta"] = frozen["cpg_id"].map(o3_delta)
    frozen["abs_discovery_delta"] = frozen["discovery_delta_beta"].abs()
    frozen["abs_validation_delta"] = frozen["validation_delta_beta"].abs()
    frozen["same_direction"] = (
        np.sign(frozen["discovery_delta_beta"])
        == np.sign(frozen["validation_delta_beta"])
    ) & frozen["validation_delta_beta"].notna() & (frozen["validation_delta_beta"] != 0)
    frozen["effect_retention"] = (
        frozen["abs_validation_delta"] / frozen["abs_discovery_delta"].clip(lower=1e-12)
    )
    frozen["effect_shrinkage"] = 1.0 - frozen["effect_retention"]
    frozen["tier"] = _assign_tiers(
        frozen,
        tier1_min_abs_delta=tier1_min_abs_delta,
        tier1_min_retention=tier1_min_retention,
        tier1_max_retention=tier1_max_retention,
        tier2_min_abs_delta=tier2_min_abs_delta,
        tier2_min_retention=tier2_min_retention,
    )
    if "weight" not in frozen.columns:
        frozen["weight"] = 1.0

    overall = _comparison_metrics(
        frozen["discovery_delta_beta"].to_numpy(dtype=float),
        frozen["validation_delta_beta"].to_numpy(dtype=float),
        frozen["weight"].to_numpy(dtype=float),
    )
    topk = _topk_metrics(frozen, top_k)

    pair_valid = np.isfinite(validation_mean) & np.isfinite(baseline_mean)
    sample_pair_pearson, sample_pair_pearson_p = _safe_corr(
        baseline_mean.to_numpy(dtype=float)[pair_valid],
        validation_mean.to_numpy(dtype=float)[pair_valid],
        "pearson",
    )
    sample_pair_spearman, sample_pair_spearman_p = _safe_corr(
        baseline_mean.to_numpy(dtype=float)[pair_valid],
        validation_mean.to_numpy(dtype=float)[pair_valid],
        "spearman",
    )

    null_rows: list[dict[str, object]] = []
    null_rows.append(
        _permutation_null(
            frozen["discovery_delta_beta"].to_numpy(dtype=float),
            frozen["validation_delta_beta"].to_numpy(dtype=float),
            n_permutations=n_permutations,
            seed=random_seed,
        )
    )
    if discovery_all_loci is not None:
        null_rows.append(
            _effect_matched_random_null(
                frozen,
                discovery_all_loci,
                o3_delta,
                n_iterations=n_random_sets,
                n_bins=random_effect_bins,
                seed=random_seed + 1,
            )
        )

    discovery_delta = frozen["discovery_delta_beta"].to_numpy(dtype=float)
    weights = frozen["weight"].to_numpy(dtype=float)
    for control_id in pair.wrong_day_baseline_ids:
        diagnostic_delta = validation_mean - beta[control_id]
        aligned = frozen["cpg_id"].map(diagnostic_delta).to_numpy(dtype=float)
        metrics = _comparison_metrics(discovery_delta, aligned, weights)
        null_rows.append(
            {
                "null_type": "wrong_day_baseline",
                "sample_id": ";".join(pair.validation_sample_ids),
                "control_id": control_id,
                **metrics,
            }
        )
    for control_id in pair.cross_donor_baseline_ids:
        diagnostic_delta = validation_mean - beta[control_id]
        aligned = frozen["cpg_id"].map(diagnostic_delta).to_numpy(dtype=float)
        metrics = _comparison_metrics(discovery_delta, aligned, weights)
        null_rows.append(
            {
                "null_type": "cross_donor_baseline",
                "sample_id": ";".join(pair.validation_sample_ids),
                "control_id": control_id,
                **metrics,
            }
        )
    for failed_id in pair.failed_sample_ids:
        diagnostic_delta = beta[failed_id] - baseline_mean
        aligned = frozen["cpg_id"].map(diagnostic_delta).to_numpy(dtype=float)
        metrics = _comparison_metrics(discovery_delta, aligned, weights)
        null_rows.append(
            {
                "null_type": "failed_reprogramming_vs_baseline",
                "sample_id": failed_id,
                "control_id": ";".join(pair.baseline_sample_ids),
                **metrics,
            }
        )

    tier_counts = frozen["tier"].value_counts().reindex(["Tier1", "Tier2", "Tier3"], fill_value=0)
    summary: dict[str, object] = {
        "analysis_status": "held_out_single_donor_descriptive_validation",
        "candidate_list_frozen_before_o3": True,
        "validation_donor": pair.donor_id,
        "day": pair.day,
        "batch": pair.batch,
        "stage": pair.stage,
        "validation_sample_ids": list(pair.validation_sample_ids),
        "baseline_sample_ids": list(pair.baseline_sample_ids),
        "wrong_day_baseline_ids": list(pair.wrong_day_baseline_ids),
        "cross_donor_baseline_ids": list(pair.cross_donor_baseline_ids),
        "failed_sample_ids": list(pair.failed_sample_ids),
        "n_signature_cpgs": int(len(frozen)),
        "n_o3_valid_cpgs": int(overall["n_valid"]),
        "tier_counts": {str(k): int(v) for k, v in tier_counts.items()},
        "tier_thresholds": {
            "tier1_min_abs_delta": tier1_min_abs_delta,
            "tier1_min_retention": tier1_min_retention,
            "tier1_max_retention": tier1_max_retention,
            "tier2_min_abs_delta": tier2_min_abs_delta,
            "tier2_min_retention": tier2_min_retention,
        },
        "overall_metrics": overall,
        "sample_pair_genomewide_qc": {
            "n_joint_valid_loci": int(pair_valid.sum()),
            "pearson_r": sample_pair_pearson,
            "pearson_p": sample_pair_pearson_p,
            "spearman_rho": sample_pair_spearman,
            "spearman_p": sample_pair_spearman_p,
        },
        "top_k_metrics": topk.to_dict(orient="records"),
        "interpretation_note": (
            "Tier labels are descriptive prioritisation categories. O3 is one independent "
            "held-out donor and does not provide confirmatory population-level inference."
        ),
        "correlated_cpg_warning": (
            "CpGs are correlated; the exact binomial and parametric correlation p-values are "
            "descriptive. The permutation and effect-matched random-set tests are the preferred "
            "empirical null comparisons."
        ),
    }
    return O3ValidationResult(
        cpg_validation=frozen,
        summary=summary,
        null_tests=pd.DataFrame(null_rows),
        sample_qc=sample_qc,
        topk_metrics=topk,
        sample_values=sample_values,
    )


def build_o3_pair_table(
    metadata: pd.DataFrame,
    discovery_pair_table: pd.DataFrame,
    *,
    validation_donor: str = "O3",
    validation_condition: str = "reprogrammed",
    baseline_condition: str = "baseline",
    discovery_donors: Sequence[str] = ("O1", "O2"),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build the O3 primary pair from the frozen O1/O2 discovery stratum.

    The discovery matched-pair table is the source of truth for day, batch, and
    stage.  This prevents post-hoc selection of a more favourable O3 time point.
    Auxiliary rows are diagnostic null comparisons and do not alter the primary
    candidate list or tier thresholds.
    """

    stratum = infer_validation_stratum(discovery_pair_table)
    try:
        pair = select_validation_pair(
            metadata,
            validation_donor=validation_donor,
            day=stratum["day"],
            batch=stratum["batch"],
            stage=stratum["stage"],
            validation_condition=validation_condition,
            baseline_condition=baseline_condition,
            discovery_donors=discovery_donors,
        )
    except ValueError as exc:
        raise ValueError(
            "No exact O3 validation/baseline pair was found for the frozen "
            f"discovery stratum {stratum}"
        ) from exc

    primary = pd.DataFrame(
        [
            {
                "validation_donor": pair.donor_id,
                "validation_sample_ids": ";".join(pair.validation_sample_ids),
                "baseline_sample_ids": ";".join(pair.baseline_sample_ids),
                "n_validation_samples": len(pair.validation_sample_ids),
                "n_baseline_samples": len(pair.baseline_sample_ids),
                "day": pair.day,
                "batch": pair.batch,
                "stage": pair.stage,
                "validation_condition": validation_condition,
                "baseline_condition": baseline_condition,
            }
        ]
    )

    auxiliary_rows: list[dict[str, object]] = []
    for sample_id in pair.failed_sample_ids:
        auxiliary_rows.append(
            {
                "null_type": "failed_reprogramming_same_stratum",
                "sample_id": sample_id,
                "control_id": ";".join(pair.baseline_sample_ids),
                "day": pair.day,
                "batch": pair.batch,
                "stage": pair.stage,
            }
        )
    meta = metadata.set_index("sample_id", drop=False)
    for control_id in pair.wrong_day_baseline_ids:
        auxiliary_rows.append(
            {
                "null_type": "wrong_day_baseline",
                "sample_id": ";".join(pair.validation_sample_ids),
                "control_id": control_id,
                "day": meta.loc[control_id, "day"] if control_id in meta.index else np.nan,
                "batch": pair.batch,
                "stage": pair.stage,
            }
        )
    for control_id in pair.cross_donor_baseline_ids:
        auxiliary_rows.append(
            {
                "null_type": "cross_donor_baseline",
                "sample_id": ";".join(pair.validation_sample_ids),
                "control_id": control_id,
                "day": pair.day,
                "batch": pair.batch,
                "stage": pair.stage,
            }
        )
    return primary, pd.DataFrame(auxiliary_rows)


def _permutation_distribution_rows(
    discovery: np.ndarray,
    validation: np.ndarray,
    *,
    n_permutations: int,
    seed: int,
) -> pd.DataFrame:
    valid = np.isfinite(discovery) & np.isfinite(validation) & (discovery != 0)
    d = discovery[valid]
    v = validation[valid]
    if len(d) < 3 or n_permutations <= 0:
        return pd.DataFrame()
    observed = _comparison_metrics(d, v)
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for iteration in range(n_permutations):
        shuffled = rng.permutation(v)
        metrics = _comparison_metrics(d, shuffled)
        rows.append(
            {
                "null_type": "cpg_label_permutation",
                "iteration": iteration + 1,
                "n_valid": metrics["n_valid"],
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


def _effect_matched_random_distribution_rows(
    candidates: pd.DataFrame,
    all_loci: pd.DataFrame,
    o3_delta: pd.Series,
    *,
    n_iterations: int,
    n_bins: int,
    seed: int,
) -> pd.DataFrame:
    if n_iterations <= 0:
        return pd.DataFrame()
    required = {"cpg_id", "target_delta_beta"}
    if not required.issubset(all_loci.columns):
        raise ValueError("All-loci file must contain cpg_id and target_delta_beta")

    candidate_valid = candidates[
        np.isfinite(candidates["discovery_delta_beta"])
        & np.isfinite(candidates["validation_delta_beta"])
    ].copy()
    pool = all_loci.loc[:, ["cpg_id", "target_delta_beta"]].copy()
    pool["cpg_id"] = pool["cpg_id"].astype(str)
    pool["validation_delta_beta"] = pool["cpg_id"].map(o3_delta)
    pool = pool[~pool["cpg_id"].isin(set(candidates["cpg_id"].astype(str)))]
    pool = pool[
        np.isfinite(pool["target_delta_beta"])
        & np.isfinite(pool["validation_delta_beta"])
        & (pool["target_delta_beta"] != 0)
    ].copy()
    if candidate_valid.empty or pool.empty:
        return pd.DataFrame()

    combined_abs = pd.concat(
        [pool["target_delta_beta"].abs(), candidate_valid["discovery_delta_beta"].abs()],
        ignore_index=True,
    )
    quantiles = np.unique(
        np.quantile(combined_abs, np.linspace(0, 1, max(int(n_bins), 2) + 1))
    )
    if len(quantiles) < 3:
        quantiles = np.asarray([-np.inf, np.inf])
    else:
        quantiles[0], quantiles[-1] = -np.inf, np.inf
    pool["effect_bin"] = pd.cut(
        pool["target_delta_beta"].abs(), bins=quantiles, labels=False, include_lowest=True
    )
    candidate_valid["effect_bin"] = pd.cut(
        candidate_valid["discovery_delta_beta"].abs(),
        bins=quantiles,
        labels=False,
        include_lowest=True,
    )
    counts = candidate_valid["effect_bin"].value_counts(dropna=False).to_dict()
    grouped = {key: group for key, group in pool.groupby("effect_bin", dropna=False)}
    observed = _comparison_metrics(
        candidate_valid["discovery_delta_beta"].to_numpy(dtype=float),
        candidate_valid["validation_delta_beta"].to_numpy(dtype=float),
    )
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    for iteration in range(n_iterations):
        sampled_parts: list[pd.DataFrame] = []
        for key, count in counts.items():
            group = grouped.get(key)
            if group is None or group.empty:
                group = pool
            chosen = rng.choice(
                group.index.to_numpy(), size=int(count), replace=len(group) < int(count)
            )
            sampled_parts.append(group.loc[chosen])
        sampled = pd.concat(sampled_parts, ignore_index=True)
        metrics = _comparison_metrics(
            sampled["target_delta_beta"].to_numpy(dtype=float),
            sampled["validation_delta_beta"].to_numpy(dtype=float),
        )
        rows.append(
            {
                "null_type": "effect_matched_random_cpg_set",
                "iteration": iteration + 1,
                "pool_size": int(len(pool)),
                "n_valid": metrics["n_valid"],
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


def validate_o3_signature_file(
    *,
    beta_path: str | Path,
    metadata: pd.DataFrame,
    discovery_signature: pd.DataFrame,
    discovery_pair_table: pd.DataFrame,
    discovery_all_loci: pd.DataFrame | None = None,
    validation_donor: str = "O3",
    validation_condition: str = "reprogrammed",
    baseline_condition: str = "baseline",
    discovery_donors: Sequence[str] = ("O1", "O2"),
    tier1_min_abs_delta: float = 0.03,
    tier1_min_effect_ratio: float = 0.25,
    tier1_max_effect_ratio: float = 2.0,
    tier2_min_abs_delta: float = 0.01,
    tier2_min_effect_ratio: float = 0.0,
    top_k_values: Sequence[int] = (50, 100, 250, 500),
    n_permutations: int = 1000,
    n_random_sets: int = 1000,
    random_effect_bins: int = 10,
    random_seed: int = 2026,
    chunksize: int = 10000,
) -> O3ValidationResult:
    """File-backed O3 validation using the frozen O1/O2 matched-pair stratum."""

    primary, _ = build_o3_pair_table(
        metadata,
        discovery_pair_table,
        validation_donor=validation_donor,
        validation_condition=validation_condition,
        baseline_condition=baseline_condition,
        discovery_donors=discovery_donors,
    )
    row = primary.iloc[0]
    base_result = validate_o3_signature(
        beta_path,
        metadata,
        discovery_signature,
        validation_donor=validation_donor,
        day=row["day"],
        batch=row["batch"],
        stage=row["stage"],
        validation_condition=validation_condition,
        baseline_condition=baseline_condition,
        discovery_donors=discovery_donors,
        top_k=top_k_values,
        tier1_min_abs_delta=tier1_min_abs_delta,
        tier1_min_retention=tier1_min_effect_ratio,
        tier1_max_retention=tier1_max_effect_ratio,
        tier2_min_abs_delta=tier2_min_abs_delta,
        tier2_min_retention=tier2_min_effect_ratio,
        n_permutations=0,
        random_seed=random_seed,
        discovery_all_loci=None,
        n_random_sets=0,
        chunksize=chunksize,
    )

    frame = base_result.cpg_validation
    permutation_rows = _permutation_distribution_rows(
        frame["discovery_delta_beta"].to_numpy(dtype=float),
        frame["validation_delta_beta"].to_numpy(dtype=float),
        n_permutations=n_permutations,
        seed=random_seed,
    )

    pair = select_validation_pair(
        metadata,
        validation_donor=validation_donor,
        day=row["day"],
        batch=row["batch"],
        stage=row["stage"],
        validation_condition=validation_condition,
        baseline_condition=baseline_condition,
        discovery_donors=discovery_donors,
    )
    required_samples = _unique_order(
        [*pair.validation_sample_ids, *pair.baseline_sample_ids]
    )
    beta_primary, _, _ = read_validation_beta(
        beta_path, required_samples, chunksize=chunksize
    )
    o3_delta = (
        beta_primary.loc[:, list(pair.validation_sample_ids)].mean(axis=1)
        - beta_primary.loc[:, list(pair.baseline_sample_ids)].mean(axis=1)
    )
    random_rows = (
        _effect_matched_random_distribution_rows(
            frame,
            discovery_all_loci,
            o3_delta,
            n_iterations=n_random_sets,
            n_bins=random_effect_bins,
            seed=random_seed + 1,
        )
        if discovery_all_loci is not None
        else pd.DataFrame()
    )
    diagnostics = base_result.null_tests[
        ~base_result.null_tests["null_type"].isin(
            ["cpg_label_permutation", "effect_matched_random_cpg_sets"]
        )
    ].copy() if not base_result.null_tests.empty else pd.DataFrame()
    null_tests = pd.concat(
        [x for x in [permutation_rows, random_rows, diagnostics] if not x.empty],
        ignore_index=True,
        sort=False,
    ) if any(not x.empty for x in [permutation_rows, random_rows, diagnostics]) else pd.DataFrame()

    summary = dict(base_result.summary)
    counts = frame["tier"].value_counts()
    empirical_null_summary: dict[str, object] = {}
    for null_name in ["cpg_label_permutation", "effect_matched_random_cpg_set"]:
        subset = null_tests[null_tests["null_type"] == null_name] if not null_tests.empty else pd.DataFrame()
        if subset.empty:
            continue
        observed_sign = float(subset["observed_sign_concordance"].iloc[0])
        observed_rho = float(subset["observed_abs_spearman"].iloc[0])
        empirical_null_summary[null_name] = {
            "n_iterations": int(len(subset)),
            "observed_sign_concordance": observed_sign,
            "null_mean_sign_concordance": float(subset["sign_concordance"].mean()),
            "empirical_p_sign": float(
                (1 + (subset["sign_concordance"] >= observed_sign).sum()) / (len(subset) + 1)
            ),
            "observed_abs_spearman": observed_rho,
            "null_mean_abs_spearman": float(subset["abs_spearman"].mean()),
            "empirical_p_spearman": float(
                (1 + (subset["abs_spearman"] >= observed_rho).sum()) / (len(subset) + 1)
            ),
        }

    summary.update(
        {
            "independent_validation_donors": 1,
            "frozen_discovery_stratum": {
                "day": row["day"],
                "batch": row["batch"],
                "stage": row["stage"],
            },
            "n_tier1": int(counts.get("Tier1", 0)),
            "n_tier2": int(counts.get("Tier2", 0)),
            "n_tier3": int(counts.get("Tier3", 0)),
            "n_permutations": int(n_permutations),
            "n_effect_matched_random_sets": int(n_random_sets),
            "empirical_null_summary": empirical_null_summary,
        }
    )
    return O3ValidationResult(
        cpg_validation=frame,
        summary=summary,
        null_tests=null_tests,
        sample_qc=base_result.sample_qc,
        topk_metrics=base_result.topk_metrics,
        sample_values=base_result.sample_values,
    )


# Backward/forward-compatible semantic alias used by the command-line layer.
O3ValidationResult.validation_table = property(lambda self: self.cpg_validation)  # type: ignore[attr-defined]
