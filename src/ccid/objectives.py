from __future__ import annotations

import numpy as np

from .schema import OptimizationContext


def _bool_column(context: OptimizationContext, name: str, default: bool = True) -> np.ndarray:
    table = context.candidate_table
    if name not in table.columns:
        return np.full(len(table), default, dtype=bool)
    return table[name].to_numpy(dtype=bool)


def baseline_target_mismatch(context: OptimizationContext) -> float:
    """Weighted RMSE for the no-intervention state (zero predicted effect)."""
    return float(
        np.sqrt(
            np.average(
                np.square(context.target_effect),
                weights=context.locus_weights,
            )
        )
    )


def relative_target_improvement(target_mismatch: float, context: OptimizationContext) -> float:
    """Relative mismatch reduction versus the prespecified no-intervention baseline."""
    baseline = baseline_target_mismatch(context)
    if baseline <= 1e-12:
        return 0.0
    return float((baseline - target_mismatch) / baseline)


def no_intervention_summary(context: OptimizationContext) -> dict[str, float | bool | str]:
    """Return the zero-edit comparator without treating it as a candidate protocol."""
    mismatch = baseline_target_mismatch(context)
    return {
        "baseline": "No intervention",
        "target_mismatch": mismatch,
        "relative_target_improvement": 0.0,
        "safety_risk": 0.0,
        "edit_count": 0.0,
        "uncertainty": 0.0,
        "is_protocol_candidate": False,
    }


def repair_mask(mask: np.ndarray, context: OptimizationContext, rng: np.random.Generator) -> np.ndarray:
    """Repair hard-invalid, oversized and undersized protocol masks.

    The no-intervention state is deliberately excluded from protocol search when
    ``context.min_edits >= 1``. It remains available through
    :func:`no_intervention_summary` as a separate comparator.
    """

    if context.min_edits < 0:
        raise ValueError("min_edits must be non-negative")
    if context.max_edits < context.min_edits:
        raise ValueError("max_edits must be >= min_edits")

    repaired = mask.astype(bool, copy=True)
    table = context.candidate_table

    hard_valid = (
        _bool_column(context, "pam_feasible")
        & _bool_column(context, "direction_feasible")
        & _bool_column(context, "chromatin_accessible")
        & _bool_column(context, "amplicon_feasible")
    )
    repaired[~hard_valid] = False

    def desirability(indices: np.ndarray) -> np.ndarray:
        if len(indices) == 0:
            return np.empty(0, dtype=float)
        return (
            table.loc[indices, "causal_score"].to_numpy(dtype=float)
            + table.loc[indices, "editability"].to_numpy(dtype=float)
            - table.loc[indices, "sequence_offtarget"].to_numpy(dtype=float)
            - table.loc[indices, "oncogenic_risk"].to_numpy(dtype=float)
            - table.loc[indices, "identity_risk"].to_numpy(dtype=float)
            + rng.normal(0, 1e-9, len(indices))
        )

    selected = np.flatnonzero(repaired)
    if len(selected) > context.max_edits:
        keep = selected[np.argsort(desirability(selected))[-context.max_edits :]]
        repaired[:] = False
        repaired[keep] = True

    selected = np.flatnonzero(repaired)
    if len(selected) < context.min_edits:
        eligible = np.flatnonzero(hard_valid & ~repaired)
        if len(eligible):
            # Prefer candidates from regions not already represented.
            selected_regions = set(table.loc[selected, "region_id"].astype(str))
            unique_region = np.array(
                [
                    idx
                    for idx in eligible
                    if str(table.iloc[int(idx)]["region_id"]) not in selected_regions
                ],
                dtype=int,
            )
            ordered_groups = [unique_region, np.setdiff1d(eligible, unique_region)]
            needed = context.min_edits - len(selected)
            for group in ordered_groups:
                if needed <= 0 or len(group) == 0:
                    continue
                ordered = group[np.argsort(desirability(group))[::-1]]
                take = ordered[:needed]
                repaired[take] = True
                needed -= len(take)

    return repaired


def predict_combined_effect(mask: np.ndarray, context: OptimizationContext) -> tuple[np.ndarray, np.ndarray]:
    selected = np.flatnonzero(mask)
    if len(selected) == 0:
        zeros = np.zeros_like(context.target_effect)
        return zeros, zeros

    mean = context.predicted_effects[selected].sum(axis=0)
    # Conservative independent-error aggregation. Real analysis should also
    # test correlated residuals by study and editor.
    std = np.sqrt(np.square(context.predicted_std[selected]).sum(axis=0))

    # Small nonlinear interaction proxy for multiplexed edits.
    if len(selected) > 1:
        interaction_scale = 0.018 * (len(selected) - 1)
        mean = mean + interaction_scale * np.tanh(mean) * np.cos(np.arange(len(mean)) / 5)
    return mean, std


def evaluate_mask(mask: np.ndarray, context: OptimizationContext) -> tuple[np.ndarray, bool, float]:
    table = context.candidate_table
    selected = np.flatnonzero(mask)
    mean, std = predict_combined_effect(mask, context)

    robust_effect = mean - context.uncertainty_quantile * np.sign(mean) * std
    residual = robust_effect - context.target_effect
    target_mismatch = float(
        np.sqrt(np.average(np.square(residual), weights=context.locus_weights))
    )
    relative_improvement = relative_target_improvement(target_mismatch, context)

    if len(selected) == 0:
        safety_risk = 0.0
        identity_risk = 0.0
        causal_shortfall = 0.0
    else:
        sequence_risk = table.loc[selected, "sequence_offtarget"].sum()
        bystander_risk = table.loc[selected, "bystander_risk"].sum()
        oncogenic_risk = table.loc[selected, "oncogenic_risk"].sum()
        identity_risk = table.loc[selected, "identity_risk"].sum()
        low_editability = np.maximum(
            0.5 - table.loc[selected, "editability"].to_numpy(), 0
        ).sum()
        safety_risk = float(
            0.30 * sequence_risk
            + 0.25 * bystander_risk
            + 0.25 * oncogenic_risk
            + 0.20 * identity_risk
            + 0.15 * low_editability
        )
        causal_shortfall = float(
            np.maximum(
                context.causal_threshold
                - table.loc[selected, "causal_score"].to_numpy(),
                0,
            ).sum()
        )

    edit_count = float(len(selected))
    uncertainty = float(np.average(std, weights=context.locus_weights))

    violations = 0.0
    violations += max(0.0, context.min_edits - edit_count)
    violations += max(0.0, edit_count - context.max_edits)
    violations += max(0.0, context.min_relative_improvement - relative_improvement)
    violations += max(0.0, identity_risk - context.identity_risk_limit)
    violations += max(0.0, safety_risk - context.safety_risk_limit)
    violations += causal_shortfall

    if len(selected) > 0:
        duplicate_regions = len(selected) - table.loc[selected, "region_id"].nunique()
        violations += 0.5 * duplicate_regions
        hard_ok = (
            table.loc[selected, "pam_feasible"].all()
            and table.loc[selected, "direction_feasible"].all()
            and _bool_column(context, "chromatin_accessible")[selected].all()
            and _bool_column(context, "amplicon_feasible")[selected].all()
        )
        if not hard_ok:
            violations += 10.0

    objectives = np.array([target_mismatch, safety_risk, edit_count, uncertainty], dtype=float)
    return objectives, violations <= 1e-12, float(violations)


def scalarized_score(objectives: np.ndarray, violations: float = 0.0) -> float:
    scales = np.array([0.12, 1.2, 4.0, 0.05])
    weights = np.array([0.48, 0.24, 0.18, 0.10])
    return float(np.dot(weights, objectives / scales) + 100.0 * violations)
