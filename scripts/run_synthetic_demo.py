#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from ccid.pipeline import run_demo


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the CCID synthetic validation suite")
    parser.add_argument("--output", default=str(PROJECT_ROOT / "results"))
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()
    summary = run_demo(args.output, seed=args.seed)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
