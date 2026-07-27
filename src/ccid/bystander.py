from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


@dataclass(frozen=True)
class BystanderAnalysis:
    locus_results: pd.DataFrame
    replicate_correlations: pd.DataFrame
    summary: dict[str, float]


def analyze_bystander_effects(
    beta: pd.DataFrame,
    metadata: pd.DataFrame,
    targeted_label: str = "targeted",
    control_label: str = "control",
    treatment_column: str = "treatment",
    pair_column: str = "replicate",
    age_loci: set[str] | None = None,
    min_abs_delta: float = 0.02,
    fdr: float = 0.10,
    seed: int = 0,
    n_permutations: int = 500,
) -> BystanderAnalysis:
    """Quantify genome-wide, reproducible bystander methylation changes."""

    required = {"sample_id", treatment_column, pair_column}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Metadata is missing columns: {sorted(missing)}")

    meta = metadata.copy().set_index("sample_id", drop=False)
    common = meta.index.intersection(beta.columns)
    meta = meta.loc[common]
    beta = beta.loc[:, common].apply(pd.to_numeric, errors="coerce")

    targeted_ids = meta.index[meta[treatment_column].astype(str) == targeted_label].tolist()
    control_ids = meta.index[meta[treatment_column].astype(str) == control_label].tolist()
    if len(targeted_ids) < 2 or len(control_ids) < 2:
        raise ValueError("At least two targeted and two control samples are required")

    target_values = beta[targeted_ids].to_numpy(dtype=float)
    control_values = beta[control_ids].to_numpy(dtype=float)
    delta = np.nanmean(target_values, axis=1) - np.nanmean(control_values, axis=1)
    t_stat, p_value = stats.ttest_ind(
        target_values,
        control_values,
        axis=1,
        equal_var=False,
        nan_policy="omit",
    )
    p_value = np.nan_to_num(p_value, nan=1.0)
    q_value = multipletests(p_value, method="fdr_bh")[1]

    loci = pd.DataFrame(
        {
            "locus_id": beta.index.astype(str),
            "delta_beta": delta,
            "t_stat": np.nan_to_num(t_stat),
            "p_value": p_value,
            "q_value": q_value,
        }
    )
    loci["significant"] = (loci["q_value"] <= fdr) & (loci["delta_beta"].abs() >= min_abs_delta)
    loci["age_associated"] = loci["locus_id"].isin(age_loci or set())

    paired_effects: dict[str, pd.Series] = {}
    for pair, pair_meta in meta.groupby(pair_column):
        target_pair = pair_meta.index[pair_meta[treatment_column].astype(str) == targeted_label].tolist()
        control_pair = pair_meta.index[pair_meta[treatment_column].astype(str) == control_label].tolist()
        if target_pair and control_pair:
            paired_effects[str(pair)] = beta[target_pair].mean(axis=1) - beta[control_pair].mean(axis=1)

    correlation_rows: list[dict[str, float | str]] = []
    pair_names = list(paired_effects)
    for i in range(len(pair_names)):
        for j in range(i + 1, len(pair_names)):
            a = paired_effects[pair_names[i]]
            b = paired_effects[pair_names[j]]
            valid = a.notna() & b.notna()
            rho, p = stats.spearmanr(a[valid], b[valid])
            correlation_rows.append(
                {
                    "replicate_a": pair_names[i],
                    "replicate_b": pair_names[j],
                    "spearman_rho": float(rho),
                    "p_value": float(p),
                }
            )
    correlations = pd.DataFrame(correlation_rows)

    enrichment_odds = float("nan")
    enrichment_p = float("nan")
    permutation_p = float("nan")
    if age_loci:
        selected = loci["significant"].to_numpy(dtype=bool)
        age = loci["age_associated"].to_numpy(dtype=bool)
        contingency = np.array(
            [
                [np.sum(selected & age), np.sum(selected & ~age)],
                [np.sum(~selected & age), np.sum(~selected & ~age)],
            ]
        )
        enrichment_odds, enrichment_p = stats.fisher_exact(contingency, alternative="greater")

        rng = np.random.default_rng(seed)
        observed = np.sum(selected & age)
        null = np.empty(n_permutations, dtype=int)
        for i in range(n_permutations):
            null[i] = np.sum(selected & rng.permutation(age))
        permutation_p = float((1 + np.sum(null >= observed)) / (n_permutations + 1))

    summary = {
        "n_loci": float(len(loci)),
        "n_significant": float(loci["significant"].sum()),
        "median_abs_delta": float(loci["delta_beta"].abs().median()),
        "mean_replicate_spearman": float(correlations["spearman_rho"].mean()) if not correlations.empty else float("nan"),
        "age_enrichment_odds_ratio": float(enrichment_odds),
        "age_enrichment_fisher_p": float(enrichment_p),
        "age_enrichment_permutation_p": float(permutation_p),
    }
    return BystanderAnalysis(loci, correlations, summary)
