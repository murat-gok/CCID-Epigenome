#!/usr/bin/env python3
"""Is the 16-CpG panel enriched at OSKM binding sites?

The reviewer's strongest point: reprogramming factors demethylate their own
binding sites while they are expressed, and those sites remethylate once the
factors are withdrawn. That single mechanism would produce everything the panel
shows - a trough in the intermediate, a return in iPSCs, and partial
demethylation in failed reprogramming, where factors were also expressed. It
has to be tested before the panel can be called a marker of anything.

Data: Soufi, Donahue and Zaret (Cell 2012), OCT4/SOX2/KLF4/c-MYC ChIP-seq in
human fibroblasts 48 h after induction (GSE36570), lifted hg18 -> hg38.

Null: the same baseline-methylation-matched background used for the
state-specificity test, so the comparison is not confounded by the panel's high
control methylation.
"""

from __future__ import annotations

import gzip
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from pyliftover import LiftOver

UP = Path("/mnt/user-data/uploads/CCID-Epigenome")
OUT = Path("/mnt/user-data/outputs/oskm")
RNG = np.random.default_rng(20261118)
N_DRAWS = 1000
BIN = 0.02
WINDOWS = [0, 500, 1000, 5000]


def load_peaks() -> dict[str, list[tuple[str, int, int]]]:
    lo = LiftOver(str(UP / "data/reference/oskm_chip/hg18ToHg38.over.chain.gz"))
    peaks, dropped, total = defaultdict(list), 0, 0
    with gzip.open(UP / "data/reference/oskm_chip/GSE36570_All_48hrs_MTFBRs.bed.gz", "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or not p[1].isdigit():
                continue
            total += 1
            a = lo.convert_coordinate(p[0], int(p[1]))
            b = lo.convert_coordinate(p[0], int(p[2]))
            if not a or not b or a[0][0] != b[0][0]:
                dropped += 1
                continue
            s, e = sorted((a[0][1], b[0][1]))
            label = p[3].upper().strip('"')
            peaks[label].append((a[0][0], s, e))
    print(f"peaks: {total:,} read, {dropped:,} failed liftOver "
          f"({dropped/total:.1%}), factors: {sorted(peaks)}")
    for k in peaks:
        print(f"   {k:6s} {len(peaks[k]):>7,}")
    return peaks


def index_peaks(regions, pad: int):
    by_chr = defaultdict(list)
    for c, s, e in regions:
        by_chr[c].append((max(0, s - pad), e + pad))
    for c in by_chr:
        by_chr[c].sort()
    return {c: (np.array([x[0] for x in v]), np.array([x[1] for x in v]))
            for c, v in by_chr.items()}


def hits(idx, coords) -> np.ndarray:
    out = np.zeros(len(coords), dtype=bool)
    for i, (c, pos) in enumerate(coords):
        if c not in idx:
            continue
        starts, ends = idx[c]
        j = np.searchsorted(starts, pos, side="right") - 1
        # walk back a little: intervals are padded and may overlap
        lo_j = max(0, j - 50)
        if j >= 0 and np.any((starts[lo_j:j + 1] <= pos) & (ends[lo_j:j + 1] >= pos)):
            out[i] = True
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    panel = pd.read_csv(UP / "results/faz0/testA_panel_floor_anatomy.csv")["cpg_id"].tolist()
    frozen = set(pd.read_csv(UP / "results/real/safe_osk_target_signature_O1O2.csv.gz")["cpg_id"])
    ctrl = pd.read_csv(UP / "results/faz0/control_beta_all_loci.csv.gz",
                       index_col=0)["control_beta"]
    zh = pd.read_csv(UP / "data/reference/probe_qc/EPIC.hg38.zhou_annotation.tsv.gz",
                     sep="\t", compression="gzip", low_memory=False,
                     usecols=["probe_id", "chrm", "beg"]).dropna()
    zh = zh[zh["probe_id"].isin(ctrl.index)]
    coord = dict(zip(zh["probe_id"], zip(zh["chrm"].astype(str), zh["beg"].astype(int))))
    print(f"loci with hg38 coordinates: {len(coord):,}")

    raw = load_peaks()
    # a region labelled "KLF4,OCT4" is a KLF4 region and an OCT4 region alike,
    # so membership is by substring, not by exact label
    peaks = {f: [r for lab, v in raw.items() if f in lab.split(",") for r in v]
             for f in ("OCT4", "SOX2", "KLF4", "CMYC")}
    peaks["UNION"] = [r for v in raw.values() for r in v]
    print("\nregions per factor (multi-factor regions counted for each):")
    for f, v in peaks.items():
        print(f"   {f:6s} {len(v):>7,}")

    # baseline-matched background draws, same recipe as the state-specificity test
    excl = frozen | set(panel)
    pool = ctrl.drop(index=[i for i in excl if i in ctrl.index]).dropna()
    pool = pool[pool.index.isin(coord)]
    bins = np.floor(pool.values / BIN).astype(int)
    by_bin = defaultdict(list)
    for cpg, b in zip(pool.index.values, bins):
        by_bin[int(b)].append(cpg)
    draws = np.empty((N_DRAWS, len(panel)), dtype=object)
    for j, p in enumerate(panel):
        b = int(np.floor(ctrl.get(p, np.nan) / BIN))
        cand = by_bin.get(b, [])
        if len(cand) < 50:
            cand = sum((by_bin.get(k, []) for k in (b - 1, b, b + 1)), [])
        draws[:, j] = RNG.choice(np.asarray(cand, dtype=object), size=N_DRAWS, replace=True)

    panel_coords = [coord[p] for p in panel if p in coord]
    print(f"panel loci with coordinates: {len(panel_coords)}/{len(panel)}")
    flat = np.unique(draws)
    bg_coords = {c: coord[c] for c in flat if c in coord}

    rows = []
    for factor in ["OCT4", "SOX2", "KLF4", "CMYC", "UNION"]:
        if factor not in peaks:
            continue
        for pad in WINDOWS:
            idx = index_peaks(peaks[factor], pad)
            obs = int(hits(idx, panel_coords).sum())
            bg_keys = list(bg_coords)
            bg_hit = dict(zip(bg_keys, hits(idx, [bg_coords[k] for k in bg_keys])))
            null = np.array([sum(bg_hit.get(c, False) for c in draws[i])
                             for i in range(N_DRAWS)])
            p_emp = (1 + int((null >= obs).sum())) / (1 + N_DRAWS)
            rows.append({"factor": factor, "window_bp": pad,
                         "panel_hits": obs, "panel_n": len(panel_coords),
                         "null_mean": null.mean(), "null_sd": null.std(ddof=1),
                         "null_99pct": float(np.quantile(null, 0.99)),
                         "empirical_p": p_emp})
            print(f"  {factor:7s} ±{pad:>5d} bp   panel {obs:2d}/{len(panel_coords)}   "
                  f"null {null.mean():5.2f} ± {null.std(ddof=1):4.2f}   p = {p_emp:.3g}")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "oskm_binding_enrichment.csv", index=False)
    print(f"\nwrote {OUT/'oskm_binding_enrichment.csv'}")


if __name__ == "__main__":
    main()
