from __future__ import annotations

import numpy as np
import pandas as pd


EVIDENCE_WEIGHTS = {
    "mr_support": 0.22,
    "coloc_support": 0.16,
    "mediation_support": 0.10,
    "perturbation_support": 0.22,
    "damage_clock_support": 0.12,
    "safe_reprogramming_support": 0.18,
}


def build_evidence_score(table: pd.DataFrame) -> pd.DataFrame:
    """Combine heterogeneous evidence without treating missingness as absence."""

    output = table.copy()
    weighted_sum = np.zeros(len(output), dtype=float)
    observed_weight = np.zeros(len(output), dtype=float)

    for column, weight in EVIDENCE_WEIGHTS.items():
        if column not in output:
            output[column] = np.nan
        values = output[column].astype(float).to_numpy()
        observed = np.isfinite(values)
        weighted_sum[observed] += weight * np.clip(values[observed], 0, 1)
        observed_weight[observed] += weight

    output["causal_score"] = np.divide(
        weighted_sum,
        observed_weight,
        out=np.zeros_like(weighted_sum),
        where=observed_weight > 0,
    )
    output["evidence_completeness"] = observed_weight / sum(EVIDENCE_WEIGHTS.values())
    output["causal_score_shrunk"] = (
        output["causal_score"] * output["evidence_completeness"]
        + 0.25 * (1 - output["evidence_completeness"])
    )

    strong = (
        (output["mr_support"].fillna(0) >= 0.7)
        & (output["coloc_support"].fillna(0) >= 0.7)
        & (output["perturbation_support"].fillna(0) >= 0.5)
    )
    moderate = (
        (output["causal_score_shrunk"] >= 0.50)
        & (output["evidence_completeness"] >= 0.45)
    )
    output["evidence_tier"] = np.select(
        [strong, moderate],
        ["Tier-1", "Tier-2"],
        default="Tier-3/exploratory",
    )
    return output


def validate_guide_feasibility(table: pd.DataFrame) -> pd.DataFrame:
    required = {
        "candidate_id",
        "region_id",
        "editor",
        "guide_sequence",
        "pam_feasible",
        "direction_feasible",
        "chromatin_accessible",
        "amplicon_feasible",
        "editability",
        "sequence_offtarget",
        "bystander_risk",
        "oncogenic_risk",
        "identity_risk",
    }
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"Guide table is missing columns: {sorted(missing)}")

    output = table.copy()
    output["hard_feasible"] = (
        output["pam_feasible"].astype(bool)
        & output["direction_feasible"].astype(bool)
        & output["chromatin_accessible"].astype(bool)
        & output["amplicon_feasible"].astype(bool)
    )
    output["guide_quality"] = (
        0.45 * output["editability"]
        + 0.20 * (1 - output["sequence_offtarget"])
        + 0.15 * (1 - output["bystander_risk"])
        + 0.10 * (1 - output["oncogenic_risk"])
        + 0.10 * (1 - output["identity_risk"])
    )
    return output
