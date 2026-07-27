from __future__ import annotations

import numpy as np
import pandas as pd

from .schema import SyntheticProblem, validate_candidate_table


def make_synthetic_problem(
    seed: int = 7,
    n_candidates: int = 28,
    n_loci: int = 48,
    n_features: int = 12,
    planted_k: int = 3,
) -> SyntheticProblem:
    """Create a benchmark with nonlinear, network-like bystander effects.

    The planted target is generated from a small feasible edit set, allowing
    recovery tests without claiming biological validity.
    """

    if planted_k >= n_candidates:
        raise ValueError("planted_k must be smaller than n_candidates")

    rng = np.random.default_rng(seed)
    editor_names = ("dCas9-TET1", "dCas9-DNMT3A", "CRISPRoff")
    editor_idx = rng.integers(0, len(editor_names), size=n_candidates)

    region_latent = rng.normal(size=(n_candidates, 5))
    sequence_features = rng.normal(size=(n_candidates, n_features - 5))
    features = np.concatenate([region_latent, sequence_features], axis=1)

    locus_groups = np.repeat(np.arange(8), repeats=int(np.ceil(n_loci / 8)))[:n_loci]
    group_loadings = rng.normal(0, 0.22, size=(8, n_features))
    local_loadings = rng.normal(0, 0.10, size=(n_loci, n_features))
    editor_loadings = rng.normal(0, 0.12, size=(len(editor_names), n_loci))

    linear = features @ (group_loadings[locus_groups] + local_loadings).T
    editor_term = editor_loadings[editor_idx]
    nonlinear = 0.10 * np.tanh(features[:, [0]] * features[:, [3]]) @ rng.normal(
        size=(1, n_loci)
    )
    true_effects = 0.12 * np.tanh(linear + editor_term + nonlinear)

    causal_score = np.clip(0.55 + 0.18 * features[:, 0] - 0.10 * features[:, 4], 0, 1)
    editability = np.clip(0.65 + 0.15 * features[:, 1] + rng.normal(0, 0.07, n_candidates), 0, 1)
    sequence_offtarget = np.clip(0.30 - 0.12 * features[:, 2] + rng.normal(0, 0.08, n_candidates), 0, 1)
    bystander_risk = np.clip(np.mean(np.abs(true_effects), axis=1) * 4.5, 0, 1)
    oncogenic_risk = np.clip(0.22 + 0.14 * np.abs(features[:, 3]) + rng.normal(0, 0.05, n_candidates), 0, 1)
    identity_risk = np.clip(0.18 + 0.15 * np.maximum(features[:, 4], 0) + rng.normal(0, 0.05, n_candidates), 0, 1)
    pam_feasible = rng.random(n_candidates) > 0.08
    direction_feasible = rng.random(n_candidates) > 0.05

    table = pd.DataFrame(
        {
            "candidate_id": [f"cand_{i:03d}" for i in range(n_candidates)],
            "region_id": [f"region_{i // 2:03d}" for i in range(n_candidates)],
            "editor": [editor_names[i] for i in editor_idx],
            "causal_score": causal_score,
            "editability": editability,
            "sequence_offtarget": sequence_offtarget,
            "bystander_risk": bystander_risk,
            "oncogenic_risk": oncogenic_risk,
            "identity_risk": identity_risk,
            "pam_feasible": pam_feasible,
            "direction_feasible": direction_feasible,
            "chromatin_accessible": rng.random(n_candidates) > 0.10,
            "amplicon_feasible": rng.random(n_candidates) > 0.06,
        }
    )

    valid = (
        table["pam_feasible"]
        & table["direction_feasible"]
        & table["chromatin_accessible"]
        & table["amplicon_feasible"]
        & (table["causal_score"] >= 0.35)
        & (table["identity_risk"] <= 0.55)
    )
    feasible_indices = np.flatnonzero(valid.to_numpy())
    if len(feasible_indices) < planted_k:
        # Deterministically repair the first candidates in the rare case that
        # random feasibility removes too many options.
        feasible_indices = np.arange(planted_k)
        table.loc[feasible_indices, [
            "pam_feasible",
            "direction_feasible",
            "chromatin_accessible",
            "amplicon_feasible",
        ]] = True
        table.loc[feasible_indices, "causal_score"] = 0.8
        table.loc[feasible_indices, "identity_risk"] = 0.15

    desirability = (
        table.loc[feasible_indices, "causal_score"].to_numpy()
        + table.loc[feasible_indices, "editability"].to_numpy()
        - table.loc[feasible_indices, "sequence_offtarget"].to_numpy()
        - table.loc[feasible_indices, "oncogenic_risk"].to_numpy()
        - table.loc[feasible_indices, "identity_risk"].to_numpy()
    )
    planted_indices = tuple(feasible_indices[np.argsort(desirability)[-planted_k:]].tolist())

    # A mild interaction term makes the inverse problem non-additive.
    target_effect = true_effects[list(planted_indices)].sum(axis=0)
    for a in range(len(planted_indices)):
        for b in range(a + 1, len(planted_indices)):
            ia, ib = planted_indices[a], planted_indices[b]
            interaction = 0.025 * np.tanh(features[ia, 0] * features[ib, 1])
            target_effect += interaction * np.cos(np.arange(n_loci) / 4)

    target_effect += rng.normal(0, 0.004, size=n_loci)
    locus_weights = np.ones(n_loci, dtype=float)
    locus_weights[locus_groups <= 2] = 1.6
    locus_weights /= locus_weights.mean()

    validate_candidate_table(table)
    return SyntheticProblem(
        candidate_table=table,
        candidate_features=features,
        true_effects=true_effects,
        target_effect=target_effect,
        locus_weights=locus_weights,
        planted_indices=planted_indices,
        locus_groups=locus_groups,
        editor_names=editor_names,
    )


def make_surrogate_training_data(
    problem: SyntheticProblem,
    seed: int = 11,
    n_replicates: int = 16,
    n_studies: int = 8,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Generate study-structured perturbation observations.

    Study offsets intentionally create a leakage trap: random CV is optimistic,
    while GroupKFold approximates external-study generalization.
    """

    rng = np.random.default_rng(seed)
    x_rows: list[np.ndarray] = []
    y_rows: list[np.ndarray] = []
    groups: list[int] = []
    candidate_indices: list[int] = []

    study_offsets = rng.normal(0, 0.025, size=(n_studies, problem.true_effects.shape[1]))
    study_feature_offsets = rng.normal(0, 0.22, size=(n_studies, problem.candidate_features.shape[1]))
    for candidate_idx, base_x in enumerate(problem.candidate_features):
        for rep in range(n_replicates):
            study = (candidate_idx + rep) % n_studies
            assay_noise = rng.normal(0, 0.05, size=base_x.shape)
            x_rows.append(base_x + study_feature_offsets[study] + assay_noise)
            y_rows.append(
                problem.true_effects[candidate_idx]
                + study_offsets[study]
                + rng.normal(0, 0.012, size=problem.true_effects.shape[1])
            )
            groups.append(study)
            candidate_indices.append(candidate_idx)

    return (
        np.vstack(x_rows),
        np.vstack(y_rows),
        np.asarray(groups, dtype=int),
        np.asarray(candidate_indices, dtype=int),
    )
