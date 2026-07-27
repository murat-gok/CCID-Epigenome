from __future__ import annotations

import pandas as pd

from ccid.candidates import build_evidence_score, validate_guide_feasibility


def test_candidate_scoring_handles_missing_evidence() -> None:
    table = pd.DataFrame(
        {
            "candidate_id": ["a"],
            "region_id": ["r1"],
            "editor": ["dCas9-TET1"],
            "guide_sequence": ["ACTG" * 5],
            "pam_feasible": [True],
            "direction_feasible": [True],
            "chromatin_accessible": [True],
            "amplicon_feasible": [True],
            "editability": [0.8],
            "sequence_offtarget": [0.1],
            "bystander_risk": [0.2],
            "oncogenic_risk": [0.1],
            "identity_risk": [0.1],
            "mr_support": [0.8],
            "coloc_support": [0.9],
            "mediation_support": [float("nan")],
            "perturbation_support": [0.7],
            "damage_clock_support": [0.6],
            "safe_reprogramming_support": [0.8],
        }
    )
    scored = build_evidence_score(table)
    scored["causal_score"] = scored["causal_score_shrunk"]
    checked = validate_guide_feasibility(scored)
    assert checked.loc[0, "hard_feasible"]
    assert checked.loc[0, "causal_score"] > 0.5
    assert checked.loc[0, "evidence_tier"] in {"Tier-1", "Tier-2"}
