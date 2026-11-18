#!/usr/bin/env python3
"""Figure source for the dose-response figure: distance bands overall, and the
same gradient computed separately within narrow bands of control methylation."""
from pathlib import Path
import numpy as np, pandas as pd

SRC = Path("/mnt/user-data/outputs/dose/dose_response_per_locus.csv.gz")
OUT = Path("/mnt/user-data/outputs/dose")

df = pd.read_csv(SRC, index_col=0)
d = df["d_KLF4_or_OCT4"]

BANDS = [(0, 0, "Inside"), (1, 500, "1-500 bp"), (501, 1000, "0.5-1 kb"),
         (1001, 5000, "1-5 kb"), (5001, 20000, "5-20 kb"),
         (20001, 100000, "20-100 kb"), (100001, np.inf, ">100 kb")]
rows = []
for lo, hi, lab in BANDS:
    m = (d >= lo) & (d <= hi)
    v = df.loc[m, "trough"]
    rows.append(dict(band=lab, n=int(m.sum()), mean=v.mean(),
                     sem=v.std(ddof=1) / np.sqrt(len(v)),
                     mean_control_beta=df.loc[m, "control_beta"].mean()))
overall = pd.DataFrame(rows)
overall.to_csv(OUT / "fig7_bands_overall.csv", index=False)
print(overall.round(4).to_string(index=False))

BETA = [(0.6, 0.7), (0.7, 0.8), (0.8, 0.9)]
rows = []
for blo, bhi in BETA:
    s = df[(df["control_beta"] >= blo) & (df["control_beta"] < bhi)]
    for lo, hi, lab in BANDS:
        m = (s["d_KLF4_or_OCT4"] >= lo) & (s["d_KLF4_or_OCT4"] <= hi)
        if m.sum() < 30:
            continue
        v = s.loc[m, "trough"]
        rows.append(dict(beta_band=f"{blo:.1f}-{bhi:.1f}", band=lab, n=int(m.sum()),
                         mean=v.mean(), sem=v.std(ddof=1) / np.sqrt(len(v))))
strat = pd.DataFrame(rows)
strat.to_csv(OUT / "fig7_bands_by_control_beta.csv", index=False)
print()
print(strat.pivot(index="band", columns="beta_band", values="mean")
      .reindex([b[2] for b in BANDS]).round(4).to_string())
