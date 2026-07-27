#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import pandas as pd

from ccid.statistics import compare_algorithms


def main() -> None:
    parser = argparse.ArgumentParser(description="Run paper-level optimizer statistics")
    parser.add_argument("--results", required=True)
    parser.add_argument("--metric", default="hypervolume")
    parser.add_argument("--higher-is-better", action="store_true")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    results = pd.read_csv(args.results)
    pairwise, omnibus = compare_algorithms(
        results,
        metric=args.metric,
        lower_is_better=not args.higher_is_better,
    )
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    pairwise.to_csv(output / f"pairwise_{args.metric}.csv", index=False)
    (output / f"omnibus_{args.metric}.json").write_text(json.dumps(omnibus, indent=2))
    print(json.dumps(omnibus, indent=2))


if __name__ == "__main__":
    main()
