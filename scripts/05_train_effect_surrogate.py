#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import numpy as np
import pandas as pd

from ccid.surrogate import BootstrapEffectSurrogate, grouped_train_cal_test_split


def main() -> None:
    parser = argparse.ArgumentParser(description="Train an editor- and study-aware effect surrogate")
    parser.add_argument("--features", required=True, help="NPY matrix: observations x features")
    parser.add_argument("--effects", required=True, help="NPY matrix: observations x target loci")
    parser.add_argument("--groups", required=True, help="NPY vector of study/donor groups")
    parser.add_argument("--candidate-features", required=True, help="NPY matrix: candidates x features")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--test-groups", default="6,7")
    parser.add_argument("--calibration-groups", default="5")
    parser.add_argument("--seed", type=int, default=19)
    args = parser.parse_args()

    x = np.load(args.features)
    y = np.load(args.effects)
    groups = np.load(args.groups)
    candidate_features = np.load(args.candidate_features)
    test_groups = tuple(int(x) for x in args.test_groups.split(",") if x)
    calibration_groups = tuple(int(x) for x in args.calibration_groups.split(",") if x)
    train_idx, cal_idx, test_idx = grouped_train_cal_test_split(
        groups,
        test_groups=test_groups,
        calibration_groups=calibration_groups,
    )

    model = BootstrapEffectSurrogate(random_state=args.seed)
    model.fit(x[train_idx], y[train_idx])
    model.calibrate(x[cal_idx], y[cal_idx], level=0.90)
    evaluation = model.evaluate(x[test_idx], y[test_idx])
    mean, std = model.predict_raw(candidate_features)

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    np.save(output / "candidate_effect_mean.npy", mean)
    np.save(output / "candidate_effect_std.npy", std)
    metrics = {
        "test_rmse": evaluation.rmse,
        "empirical_90pct_coverage": evaluation.empirical_coverage,
        "conformal_q": evaluation.conformal_q,
        "n_train": int(len(train_idx)),
        "n_calibration": int(len(cal_idx)),
        "n_test": int(len(test_idx)),
    }
    (output / "surrogate_metrics.json").write_text(json.dumps(metrics, indent=2))
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
