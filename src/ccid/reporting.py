from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .schema import SolutionRecord


def records_to_frame(records: list[SolutionRecord]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for solution_id, record in enumerate(records):
        rows.append(
            {
                "solution_id": solution_id,
                "method": record.method,
                "seed": record.seed,
                "feasible": record.feasible,
                "violations": record.violations,
                "target_mismatch": record.objectives[0],
                "safety_risk": record.objectives[1],
                "edit_count": int(record.objectives[2]),
                "uncertainty": record.objectives[3],
                "selected_indices": ";".join(map(str, record.selected_indices)),
            }
        )
    return pd.DataFrame(rows)


def save_pareto_plot(fronts: dict[str, list[SolutionRecord]], output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    markers = ["o", "s", "^", "D", "P"]
    for marker, (method, records) in zip(markers, fronts.items()):
        if not records:
            continue
        matrix = np.vstack([r.objectives for r in records])
        ax.scatter(
            matrix[:, 0],
            matrix[:, 1],
            label=method,
            marker=marker,
            alpha=0.78,
            s=46,
        )
    ax.set_xlabel("Robust target mismatch (lower is better)")
    ax.set_ylabel("Safety risk (lower is better)")
    ax.set_title("Pareto fronts on the synthetic inverse-design benchmark")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def save_effect_plot(
    target: np.ndarray,
    predicted: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    output_path: Path,
) -> None:
    x = np.arange(len(target))
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.plot(x, target, label="Target signature", linewidth=2)
    ax.plot(x, predicted, label="Predicted protocol effect", linewidth=1.8)
    ax.fill_between(x, lower, upper, alpha=0.22, label="Conformal uncertainty band")
    ax.axhline(0, linewidth=0.8)
    ax.set_xlabel("Synthetic methylation locus")
    ax.set_ylabel("Expected methylation shift")
    ax.set_title("Target-signature matching for the selected protocol")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
