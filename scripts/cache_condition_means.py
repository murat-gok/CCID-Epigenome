#!/usr/bin/env python3
"""Cache per-locus condition means for the whole array.

One streaming pass, three groups: baseline, the reprogramming intermediate, and
induced pluripotent stem cells. With these cached, the per-locus trough depth
used by the dose-response test is available for every evaluable locus rather
than only for the panel, which is what gives that test any power.
"""
from __future__ import annotations
import gzip
from pathlib import Path
import pandas as pd

ROOT = Path(".")
OUTF = ROOT / "results/faz0/condition_means_all_loci.csv.gz"
GROUPS = ("baseline", "reprogrammed_intermediate", "ipsc")

meta = pd.read_csv(ROOT / "data/processed/fibroblast_methylation_metadata.csv")
want = {g: set(meta.loc[meta["condition"] == g, "sample_id"]) for g in GROUPS}

with gzip.open(ROOT / "data/processed/fibroblast_beta.csv.gz", "rt") as fh:
    header = fh.readline().rstrip("\n").split(",")
    idx = {g: [i for i, c in enumerate(header) if c in want[g]] for g in GROUPS}
    print({g: len(v) for g, v in idx.items()})
    rows = []
    bad = ("", "NA", "nan", "NaN")
    for line in fh:
        p = line.rstrip("\n").split(",")
        out = [p[0]]
        for g in GROUPS:
            xs = [float(p[i]) for i in idx[g] if p[i] not in bad]
            out.append(sum(xs) / len(xs) if xs else float("nan"))
        rows.append(out)

df = pd.DataFrame(rows, columns=["cpg_id", *GROUPS]).set_index("cpg_id")
df.to_csv(OUTF, compression="gzip", float_format="%.6f")
print(f"loci: {len(df):,}  wrote {OUTF} ({OUTF.stat().st_size/1e6:.1f} MB)")
print(df.mean().round(4).to_dict())
