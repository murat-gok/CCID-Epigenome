#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import pandas as pd

from ccid.candidates import build_evidence_score, validate_guide_feasibility


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge guide feasibility and causal evidence")
    parser.add_argument("--guides", required=True, help="Guide-level CSV")
    parser.add_argument("--evidence", required=True, help="Region/CpG evidence CSV")
    parser.add_argument("--merge-key", default="region_id")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    guides = pd.read_csv(args.guides)
    evidence = pd.read_csv(args.evidence)
    merged = guides.merge(evidence, on=args.merge_key, how="left", validate="many_to_one")
    merged = build_evidence_score(merged)
    merged["causal_score"] = merged["causal_score_shrunk"]
    merged = validate_guide_feasibility(merged)
    merged = merged.sort_values(
        ["hard_feasible", "evidence_tier", "guide_quality"],
        ascending=[False, True, False],
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output, index=False)
    print(f"Wrote {len(merged)} guide-editor candidates to {output}")


if __name__ == "__main__":
    main()
