#!/usr/bin/env python3
"""Cache mean control methylation for every evaluable locus.

One streaming pass over the beta matrix, taking only the negative-control
samples matched to the day-10 discovery contrast. Cached so the matched-
background sampling can be re-run cheaply.
"""
from __future__ import annotations
import gzip
from pathlib import Path
import pandas as pd

ROOT = Path(".")
OUTF = ROOT / "results/faz0/control_beta_all_loci.csv.gz"

pairs = pd.read_csv(ROOT / "results/real/safe_osk_target_signature_O1O2_matched_pairs.csv")
ctrl = sorted({s for v in pairs["baseline_sample_ids"] for s in str(v).split(";") if s.strip()})

vals = {}
with gzip.open(ROOT / "data/processed/fibroblast_beta.csv.gz", "rt") as fh:
    header = fh.readline().rstrip("\n").split(",")
    idx = [i for i, c in enumerate(header) if c in set(ctrl)]
    assert idx, "control samples not found in the beta matrix"
    for line in fh:
        p = line.rstrip("\n").split(",")
        xs = [float(p[i]) for i in idx if p[i] not in ("", "NA", "nan", "NaN")]
        if xs:
            vals[p[0]] = sum(xs) / len(xs)

s = pd.Series(vals, name="control_beta")
s.index.name = "cpg_id"
s.to_csv(OUTF, compression="gzip")
print(f"control samples: {ctrl}")
print(f"loci with control beta: {len(s):,}")
print(f"mean {s.mean():.4f}  median {s.median():.4f}")
print(f"wrote {OUTF} ({OUTF.stat().st_size/1024:.0f} KB)")
