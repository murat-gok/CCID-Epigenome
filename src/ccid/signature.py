from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


@dataclass(frozen=True)
class SafeWindowResult:
    sample_scores: pd.DataFrame
    selected_samples: tuple[str, ...]
    selected_days: tuple[float, ...]


@dataclass(frozen=True)
class MatchedSignatureResult:
    signature: pd.DataFrame
    all_loci: pd.DataFrame
    pair_table: pd.DataFrame
    donor_effects: pd.DataFrame
    qc: dict[str, object]


def zscore_rows(
    matrix: pd.DataFrame,
    reference_samples: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Row-wise z scores fitted on a prespecified reference sample set.

    When ``reference_samples`` is provided, row means and standard deviations are
    estimated only from those samples and then applied to every column. This keeps
    held-out donors from influencing discovery-set score calibration.
    """

    if reference_samples is None:
        reference = matrix
    else:
        reference_order = [str(x) for x in reference_samples]
        missing = [x for x in reference_order if x not in matrix.columns]
        if missing:
            raise ValueError(f"Reference samples absent from matrix: {missing[:5]}")
        reference = matrix.loc[:, reference_order]
    mean = reference.mean(axis=1)
    std = reference.std(axis=1, ddof=1).replace(0, np.nan)
    return matrix.sub(mean, axis=0).div(std, axis=0).fillna(0.0)


def gene_set_score(
    expression: pd.DataFrame,
    genes: list[str],
    reference_samples: Sequence[str] | None = None,
) -> pd.Series:
    available = expression.index.intersection(genes)
    if len(available) == 0:
        raise ValueError("None of the requested genes are present in the expression matrix")
    return zscore_rows(
        expression.loc[available], reference_samples=reference_samples
    ).mean(axis=0)


def define_safe_window(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    identity_genes: list[str],
    pluripotency_genes: list[str],
    rejuvenation_genes: list[str] | None = None,
    identity_floor: float = -0.50,
    pluripotency_ceiling: float = 0.50,
    min_rejuvenation: float = 0.0,
    discovery_donors: list[str] | None = None,
    allowed_conditions: Sequence[str] = ("reprogrammed", "reprogrammed_intermediate"),
    min_donors_per_day: int = 2,
    utility_tolerance: float = 0.10,
) -> SafeWindowResult:
    """Select a leakage-free safe reprogramming window using discovery donors only.

    Required metadata columns are ``sample_id``, ``donor_id``, ``day`` and
    ``condition``. Conditions are matched exactly, preventing failed samples from
    entering merely because their raw names contain the word "reprogram".
    Candidate days must be represented by at least ``min_donors_per_day`` donors.
    """

    required = {"sample_id", "donor_id", "day", "condition"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing: {sorted(missing)}")
    if min_donors_per_day < 1:
        raise ValueError("min_donors_per_day must be at least 1")

    metadata = metadata.copy().set_index("sample_id", drop=False)
    common = metadata.index.intersection(expression.columns)
    if len(common) == 0:
        raise ValueError("No sample IDs overlap between expression and metadata")
    metadata = metadata.loc[common]
    expression = expression.loc[:, common]

    if discovery_donors is None:
        calibration_samples = common.astype(str).tolist()
    else:
        allowed_donors = {str(x) for x in discovery_donors}
        calibration_samples = metadata.index[
            metadata["donor_id"].astype(str).isin(allowed_donors)
        ].astype(str).tolist()
        if len(calibration_samples) < 2:
            raise ValueError(
                "At least two discovery-set samples are required for score calibration"
            )

    identity = gene_set_score(
        expression, identity_genes, reference_samples=calibration_samples
    )
    pluripotency = gene_set_score(
        expression, pluripotency_genes, reference_samples=calibration_samples
    )
    if rejuvenation_genes:
        rejuvenation = -gene_set_score(
            expression, rejuvenation_genes, reference_samples=calibration_samples
        )
    else:
        rejuvenation = pd.Series(0.0, index=common)

    score_columns = ["sample_id", "donor_id", "day", "condition"]
    for optional in ("batch", "stage"):
        if optional in metadata.columns:
            score_columns.append(optional)
    scores = metadata[score_columns].copy()
    scores["identity_score"] = identity
    scores["pluripotency_score"] = pluripotency
    scores["rejuvenation_score"] = rejuvenation

    discovery_mask = pd.Series(True, index=scores.index)
    if discovery_donors is not None:
        discovery_mask = scores["donor_id"].astype(str).isin([str(x) for x in discovery_donors])
    scores["discovery_set"] = discovery_mask
    scores["score_calibration_set"] = scores.index.isin(calibration_samples)

    allowed = {str(x).lower() for x in allowed_conditions}
    candidate_mask = (
        discovery_mask
        & scores["condition"].astype(str).str.lower().isin(allowed)
        & scores["day"].notna()
        & (scores["identity_score"] >= identity_floor)
        & (scores["pluripotency_score"] <= pluripotency_ceiling)
        & (scores["rejuvenation_score"] >= min_rejuvenation)
    )
    candidate = scores.loc[candidate_mask].copy()
    if candidate.empty:
        raise ValueError("No samples satisfy the prespecified safe-window rule")

    day_summary = candidate.groupby("day", dropna=False).agg(
        identity_score=("identity_score", "mean"),
        pluripotency_score=("pluripotency_score", "mean"),
        rejuvenation_score=("rejuvenation_score", "mean"),
        n_samples=("sample_id", "size"),
        n_donors=("donor_id", "nunique"),
    )
    day_summary = day_summary[day_summary["n_donors"] >= min_donors_per_day].copy()
    if day_summary.empty:
        raise ValueError(
            "No candidate day is represented by the required number of discovery donors "
            f"(min_donors_per_day={min_donors_per_day})"
        )

    day_summary["safe_utility"] = (
        day_summary["rejuvenation_score"]
        + 0.5 * day_summary["identity_score"]
        - 0.8 * np.maximum(day_summary["pluripotency_score"], 0)
    )
    best_utility = float(day_summary["safe_utility"].max())
    selected_day_index = day_summary.index[
        day_summary["safe_utility"] >= best_utility - utility_tolerance
    ]
    selected_days = tuple(float(x) for x in selected_day_index)
    selected_samples = tuple(
        candidate.index[candidate["day"].astype(float).isin(selected_days)].tolist()
    )

    scores["safe_day_eligible"] = scores["day"].isin(day_summary.index)
    scores["safe_day_utility"] = scores["day"].map(day_summary["safe_utility"])
    scores["selected_safe_window"] = scores.index.isin(selected_samples)
    return SafeWindowResult(scores.reset_index(drop=True), selected_samples, selected_days)


def _normalise_match_value(value: object) -> object:
    if pd.isna(value):
        return None
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return int(value)
    return str(value).strip().lower() if isinstance(value, str) else value


def build_exact_match_table(
    metadata: pd.DataFrame,
    safe_samples: Sequence[str],
    *,
    baseline_condition: str = "baseline",
    match_columns: Sequence[str] = ("donor_id", "day", "batch", "stage"),
    discovery_donors: Sequence[str] | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    """Match each safe sample to baseline controls from the same design stratum.

    Matching is exact on donor, day, experimental batch and stage by default.
    Multiple controls in the same stratum are averaged. No cross-day or cross-batch
    fallback is performed because that would reintroduce time/batch confounding.
    """

    required = {"sample_id", "condition", *match_columns}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing matching columns: {sorted(missing)}")

    meta = metadata.copy()
    meta["sample_id"] = meta["sample_id"].astype(str)
    meta = meta.drop_duplicates("sample_id", keep="first").set_index("sample_id", drop=False)

    safe_order = [str(x) for x in safe_samples]
    absent = [x for x in safe_order if x not in meta.index]
    if absent:
        raise ValueError(f"Safe-window samples absent from metadata: {absent[:5]}")

    safe_meta = meta.loc[safe_order].copy()
    if discovery_donors is not None:
        allowed_donors = {str(x) for x in discovery_donors}
        safe_meta = safe_meta[safe_meta["donor_id"].astype(str).isin(allowed_donors)]

    baseline = meta[
        meta["condition"].astype(str).str.lower() == str(baseline_condition).lower()
    ].copy()
    rows: list[dict[str, object]] = []
    unmatched: list[str] = []

    for safe_id, safe_row in safe_meta.iterrows():
        candidates = baseline
        key_values: dict[str, object] = {}
        for column in match_columns:
            target = _normalise_match_value(safe_row[column])
            key_values[column] = safe_row[column]
            norm_series = candidates[column].map(_normalise_match_value)
            candidates = candidates[norm_series == target]
        baseline_ids = candidates["sample_id"].astype(str).tolist()
        if not baseline_ids:
            unmatched.append(str(safe_id))
            continue
        row: dict[str, object] = {
            "pair_id": f"pair_{len(rows) + 1:03d}",
            "safe_sample_id": str(safe_id),
            "baseline_sample_ids": ";".join(baseline_ids),
            "n_baseline_samples": len(baseline_ids),
            "donor_id": str(safe_row["donor_id"]),
            "condition": str(safe_row["condition"]),
        }
        row.update(key_values)
        rows.append(row)

    if not rows:
        raise ValueError("No safe-window sample had an exact matched baseline control")
    return pd.DataFrame(rows), unmatched


def _summarise_donor_effects(
    donor_effects: pd.DataFrame,
    *,
    min_abs_delta: float,
    fdr: float,
    min_sign_consistency: float,
    min_matched_donors: int,
    min_inferential_donors: int,
) -> pd.DataFrame:
    values = donor_effects.to_numpy(dtype=float, copy=True)
    valid_n = np.isfinite(values).sum(axis=1)
    mean_delta = np.nanmean(values, axis=1)
    donor_sd = np.nanstd(values, axis=1, ddof=1)
    se = donor_sd / np.sqrt(np.maximum(valid_n, 1))
    _, p_value = stats.ttest_1samp(values, popmean=0.0, axis=1, nan_policy="omit")
    p_value = np.nan_to_num(p_value, nan=1.0, posinf=1.0, neginf=1.0)
    inferentially_supported = valid_n >= min_inferential_donors
    p_value = np.where(inferentially_supported, p_value, 1.0)
    q_value = multipletests(p_value, method="fdr_bh")[1]
    positive = np.nansum(values > 0, axis=1) / np.maximum(valid_n, 1)
    negative = np.nansum(values < 0, axis=1) / np.maximum(valid_n, 1)
    sign_consistency = np.maximum(positive, negative)

    result = pd.DataFrame(
        {
            "cpg_id": donor_effects.index.astype(str),
            "target_delta_beta": mean_delta,
            "donor_sd": donor_sd,
            "standard_error": se,
            "p_value": p_value,
            "q_value": q_value,
            "sign_consistency": sign_consistency,
            "n_donors": valid_n,
            "inferentially_supported": inferentially_supported,
        }
    )
    result["selected"] = (
        (result["n_donors"] >= min_matched_donors)
        & result["inferentially_supported"]
        & (result["target_delta_beta"].abs() >= min_abs_delta)
        & (result["q_value"] <= fdr)
        & (result["sign_consistency"] >= min_sign_consistency)
    )
    result["rank_score"] = (
        result["target_delta_beta"].abs()
        * result["sign_consistency"]
        * np.sqrt(result["n_donors"].clip(lower=1))
        / np.maximum(result["standard_error"], 0.01)
    )
    return result


def build_matched_methylation_target_signature(
    beta: pd.DataFrame,
    metadata: pd.DataFrame,
    safe_samples: Sequence[str],
    *,
    baseline_condition: str = "baseline",
    match_columns: Sequence[str] = ("donor_id", "day", "batch", "stage"),
    discovery_donors: Sequence[str] | None = None,
    min_abs_delta: float = 0.03,
    fdr: float = 0.10,
    min_sign_consistency: float = 0.67,
    min_matched_donors: int = 2,
    min_inferential_donors: int = 3,
    fallback_top_n: int = 500,
) -> MatchedSignatureResult:
    """Build an exact-matched, donor-independent safe-OSK methylation signature.

    Each selected reprogramming sample is compared with negative controls from the
    same donor, day, batch and stage. Pair-level deltas are averaged within donor;
    inferential statistics then use donors—not technical samples—as independent
    units.
    """

    pair_table, unmatched = build_exact_match_table(
        metadata,
        safe_samples,
        baseline_condition=baseline_condition,
        match_columns=match_columns,
        discovery_donors=discovery_donors,
    )

    needed = set(pair_table["safe_sample_id"])
    for item in pair_table["baseline_sample_ids"]:
        needed.update(str(item).split(";"))
    missing_columns = sorted(needed.difference(beta.columns.astype(str)))
    if missing_columns:
        raise ValueError(f"Matched samples absent from beta matrix: {missing_columns[:5]}")

    pair_effects: dict[str, pd.Series] = {}
    for row in pair_table.itertuples(index=False):
        baseline_ids = str(row.baseline_sample_ids).split(";")
        pair_effects[str(row.pair_id)] = (
            beta[str(row.safe_sample_id)] - beta[baseline_ids].mean(axis=1)
        )
    pair_effect_frame = pd.DataFrame(pair_effects, index=beta.index)

    donor_effects: dict[str, pd.Series] = {}
    for donor_id, donor_pairs in pair_table.groupby("donor_id", sort=True):
        donor_effects[str(donor_id)] = pair_effect_frame[donor_pairs["pair_id"].tolist()].mean(axis=1)
    donor_effect_frame = pd.DataFrame(donor_effects, index=beta.index)
    if donor_effect_frame.shape[1] < min_matched_donors:
        raise ValueError(
            f"Only {donor_effect_frame.shape[1]} matched donors; "
            f"at least {min_matched_donors} are required"
        )

    all_loci = _summarise_donor_effects(
        donor_effect_frame,
        min_abs_delta=min_abs_delta,
        fdr=fdr,
        min_sign_consistency=min_sign_consistency,
        min_matched_donors=min_matched_donors,
        min_inferential_donors=min_inferential_donors,
    )
    selected = all_loci[all_loci["selected"]].copy()
    underpowered = selected.empty
    if underpowered:
        selected = all_loci.sort_values(
            ["rank_score", "q_value"], ascending=[False, True]
        ).head(min(fallback_top_n, len(all_loci))).copy()
    selected["underpowered_fallback"] = underpowered
    selected["weight"] = selected["rank_score"]
    mean_weight = float(selected["weight"].replace([np.inf, -np.inf], np.nan).mean())
    if not np.isfinite(mean_weight) or mean_weight <= 0:
        selected["weight"] = 1.0
    else:
        selected["weight"] = selected["weight"] / mean_weight

    donor_pair_counts = pair_table.groupby("donor_id").size().astype(int).to_dict()
    qc: dict[str, object] = {
        "n_requested_safe_samples": int(len(safe_samples)),
        "n_matched_safe_samples": int(len(pair_table)),
        "n_unmatched_safe_samples": int(len(unmatched)),
        "unmatched_safe_samples": unmatched,
        "n_matched_donors": int(pair_table["donor_id"].nunique()),
        "matched_donors": sorted(pair_table["donor_id"].astype(str).unique().tolist()),
        "donor_pair_counts": {str(k): int(v) for k, v in donor_pair_counts.items()},
        "match_columns": list(match_columns),
        "baseline_condition": baseline_condition,
        "n_all_loci": int(len(all_loci)),
        "n_strictly_selected_loci": int(all_loci["selected"].sum()),
        "min_inferential_donors": int(min_inferential_donors),
        "inferential_status": (
            "confirmatory"
            if pair_table["donor_id"].nunique() >= min_inferential_donors
            else "exploratory_underpowered"
        ),
        "used_underpowered_fallback": bool(underpowered),
    }
    return MatchedSignatureResult(
        signature=selected.reset_index(drop=True),
        all_loci=all_loci.reset_index(drop=True),
        pair_table=pair_table.reset_index(drop=True),
        donor_effects=donor_effect_frame,
        qc=qc,
    )



def build_matched_methylation_target_signature_file(
    beta_path: str | Path,
    metadata: pd.DataFrame,
    safe_samples: Sequence[str],
    *,
    baseline_condition: str = "baseline",
    match_columns: Sequence[str] = ("donor_id", "day", "batch", "stage"),
    discovery_donors: Sequence[str] | None = None,
    min_abs_delta: float = 0.03,
    fdr: float = 0.10,
    min_sign_consistency: float = 0.67,
    min_matched_donors: int = 2,
    min_inferential_donors: int = 3,
    fallback_top_n: int = 500,
    chunksize: int = 5000,
) -> MatchedSignatureResult:
    """Chunked file-backed version of the exact-matched signature builder.

    The full CpG-by-sample beta matrix is never loaded at once. Only the samples
    needed by the exact matched pairs are read, and chunk-level pair effects are
    reduced to donor-level effects before concatenation. This substantially lowers
    memory use for EPIC-scale matrices.
    """

    beta_path = Path(beta_path)
    if not beta_path.exists():
        raise FileNotFoundError(beta_path)
    if chunksize < 1:
        raise ValueError("chunksize must be at least 1")

    pair_table, unmatched = build_exact_match_table(
        metadata,
        safe_samples,
        baseline_condition=baseline_condition,
        match_columns=match_columns,
        discovery_donors=discovery_donors,
    )
    needed_samples: list[str] = []
    for sample in pair_table["safe_sample_id"].astype(str):
        if sample not in needed_samples:
            needed_samples.append(sample)
    for item in pair_table["baseline_sample_ids"].astype(str):
        for sample in item.split(";"):
            if sample not in needed_samples:
                needed_samples.append(sample)

    header = pd.read_csv(beta_path, nrows=0, compression="infer")
    index_name = str(header.columns[0])
    available = {str(x) for x in header.columns[1:]}
    missing_columns = [x for x in needed_samples if x not in available]
    if missing_columns:
        raise ValueError(f"Matched samples absent from beta matrix: {missing_columns[:5]}")

    usecols = [index_name, *needed_samples]
    donor_chunks: list[pd.DataFrame] = []
    reader = pd.read_csv(
        beta_path,
        compression="infer",
        usecols=usecols,
        index_col=0,
        chunksize=chunksize,
        low_memory=False,
    )
    for chunk in reader:
        chunk.index = chunk.index.astype(str)
        chunk = chunk.apply(pd.to_numeric, errors="coerce")
        pair_effects: dict[str, pd.Series] = {}
        for row in pair_table.itertuples(index=False):
            baseline_ids = str(row.baseline_sample_ids).split(";")
            pair_effects[str(row.pair_id)] = (
                chunk[str(row.safe_sample_id)] - chunk[baseline_ids].mean(axis=1)
            )
        pair_frame = pd.DataFrame(pair_effects, index=chunk.index)
        donor_effects: dict[str, pd.Series] = {}
        for donor_id, donor_pairs in pair_table.groupby("donor_id", sort=True):
            donor_effects[str(donor_id)] = pair_frame[
                donor_pairs["pair_id"].tolist()
            ].mean(axis=1)
        donor_chunks.append(pd.DataFrame(donor_effects, index=chunk.index).astype(np.float32))

    if not donor_chunks:
        raise ValueError("Beta matrix contained no data rows")
    donor_effect_frame = pd.concat(donor_chunks, axis=0)
    if donor_effect_frame.index.duplicated().any():
        donor_effect_frame = donor_effect_frame.groupby(level=0, sort=False).mean()
    if donor_effect_frame.shape[1] < min_matched_donors:
        raise ValueError(
            f"Only {donor_effect_frame.shape[1]} matched donors; "
            f"at least {min_matched_donors} are required"
        )

    all_loci = _summarise_donor_effects(
        donor_effect_frame,
        min_abs_delta=min_abs_delta,
        fdr=fdr,
        min_sign_consistency=min_sign_consistency,
        min_matched_donors=min_matched_donors,
        min_inferential_donors=min_inferential_donors,
    )
    selected = all_loci[all_loci["selected"]].copy()
    underpowered = selected.empty
    if underpowered:
        selected = all_loci.sort_values(
            ["rank_score", "q_value"], ascending=[False, True]
        ).head(min(fallback_top_n, len(all_loci))).copy()
    selected["underpowered_fallback"] = underpowered
    selected["weight"] = selected["rank_score"]
    mean_weight = float(selected["weight"].replace([np.inf, -np.inf], np.nan).mean())
    if not np.isfinite(mean_weight) or mean_weight <= 0:
        selected["weight"] = 1.0
    else:
        selected["weight"] = selected["weight"] / mean_weight

    donor_pair_counts = pair_table.groupby("donor_id").size().astype(int).to_dict()
    qc: dict[str, object] = {
        "n_requested_safe_samples": int(len(safe_samples)),
        "n_matched_safe_samples": int(len(pair_table)),
        "n_unmatched_safe_samples": int(len(unmatched)),
        "unmatched_safe_samples": unmatched,
        "n_matched_donors": int(pair_table["donor_id"].nunique()),
        "matched_donors": sorted(pair_table["donor_id"].astype(str).unique().tolist()),
        "donor_pair_counts": {str(k): int(v) for k, v in donor_pair_counts.items()},
        "match_columns": list(match_columns),
        "baseline_condition": baseline_condition,
        "chunksize": int(chunksize),
        "n_beta_samples_read": int(len(needed_samples)),
        "n_all_loci": int(len(all_loci)),
        "n_strictly_selected_loci": int(all_loci["selected"].sum()),
        "min_inferential_donors": int(min_inferential_donors),
        "inferential_status": (
            "confirmatory"
            if pair_table["donor_id"].nunique() >= min_inferential_donors
            else "exploratory_underpowered"
        ),
        "used_underpowered_fallback": bool(underpowered),
    }
    return MatchedSignatureResult(
        signature=selected.reset_index(drop=True),
        all_loci=all_loci.reset_index(drop=True),
        pair_table=pair_table.reset_index(drop=True),
        donor_effects=donor_effect_frame,
        qc=qc,
    )

def build_methylation_target_signature(
    beta: pd.DataFrame,
    metadata: pd.DataFrame,
    safe_samples: list[str] | tuple[str, ...],
    baseline_condition: str = "baseline",
    min_abs_delta: float = 0.03,
    fdr: float = 0.10,
    min_sign_consistency: float = 0.67,
) -> pd.DataFrame:
    """Legacy donor-pooled signature builder retained for backward compatibility.

    New real-data analyses should use :func:`build_matched_methylation_target_signature`.
    """

    required = {"sample_id", "donor_id", "condition"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing: {sorted(missing)}")

    meta = metadata.copy().set_index("sample_id", drop=False)
    common = meta.index.intersection(beta.columns)
    meta = meta.loc[common]
    beta = beta.loc[:, common]
    safe_set = set(safe_samples)

    donor_effects: list[pd.Series] = []
    donor_ids: list[str] = []
    for donor, donor_meta in meta.groupby("donor_id"):
        safe_ids = [s for s in donor_meta.index if s in safe_set]
        base_ids = donor_meta.index[
            donor_meta["condition"].str.lower() == baseline_condition.lower()
        ].tolist()
        if safe_ids and base_ids:
            effect = beta[safe_ids].mean(axis=1) - beta[base_ids].mean(axis=1)
            donor_effects.append(effect)
            donor_ids.append(str(donor))
    if len(donor_effects) < 2:
        raise ValueError("At least two donors with safe-window and baseline samples are required")

    effects = pd.concat(donor_effects, axis=1)
    effects.columns = donor_ids
    all_loci = _summarise_donor_effects(
        effects,
        min_abs_delta=min_abs_delta,
        fdr=fdr,
        min_sign_consistency=min_sign_consistency,
        min_matched_donors=2,
        min_inferential_donors=2,
    )
    selected = all_loci[all_loci["selected"]].copy()
    if selected.empty:
        selected = all_loci.sort_values(["rank_score", "q_value"], ascending=[False, True]).head(
            min(500, len(all_loci))
        ).copy()
        selected["underpowered_fallback"] = True
    else:
        selected["underpowered_fallback"] = False
    selected["weight"] = selected["rank_score"]
    selected["weight"] /= max(float(selected["weight"].mean()), 1e-12)
    return selected.reset_index(drop=True)
