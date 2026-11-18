#!/usr/bin/env python3
"""Does the depth of the reprogramming trough scale with proximity to a
reprogramming-factor binding region?

The enrichment result is binary: a candidate is near a binding region or it is
not. A dose-response would be harder to explain away. For every evaluable locus
this measures trough depth, baseline methylation minus methylation in the
reprogramming intermediate, and the distance to the nearest KLF4, OCT4 or
combined binding region, and asks whether the two are related.

Distance is confounded with baseline methylation - a locus cannot fall far if
it starts low - so the relationship is also measured within narrow bands of
control methylation, and partial correlation controlling for it is reported.
"""

from __future__ import annotations

import gzip
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from pyliftover import LiftOver
from scipy import stats

UP = Path("/mnt/user-data/uploads/CCID-Epigenome")
OUT = Path("/mnt/user-data/outputs/dose")
MAXD = 100_000          # distances are capped; beyond this "near" is meaningless


def peaks_hg38() -> dict[str, list[tuple[str, int, int]]]:
    lo = LiftOver(str(UP / "data/reference/oskm_chip/hg18ToHg38.over.chain.gz"))
    raw = defaultdict(list)
    with gzip.open(UP / "data/reference/oskm_chip/GSE36570_All_48hrs_MTFBRs.bed.gz", "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or not p[1].isdigit():
                continue
            a = lo.convert_coordinate(p[0], int(p[1]))
            b = lo.convert_coordinate(p[0], int(p[2]))
            if not a or not b or a[0][0] != b[0][0]:
                continue
            s, e = sorted((a[0][1], b[0][1]))
            raw[p[3].upper().strip('"')].append((a[0][0], s, e))
    out = {f: [r for lab, v in raw.items() if f in lab.split(",") for r in v]
           for f in ("KLF4", "OCT4")}
    out["KLF4_or_OCT4"] = out["KLF4"] + out["OCT4"]
    return out


def distances(regions, coords: pd.DataFrame) -> np.ndarray:
    by_chr = defaultdict(list)
    for c, s, e in regions:
        by_chr[c].append((s, e))
    idx = {}
    for c, v in by_chr.items():
        v.sort()
        idx[c] = (np.array([x[0] for x in v]), np.array([x[1] for x in v]))
    d = np.full(len(coords), np.nan)
    for c, grp in coords.groupby("chrm", sort=False):
        if c not in idx:
            continue
        starts, ends = idx[c]
        pos = grp["beg"].to_numpy()
        j = np.searchsorted(starts, pos, side="right") - 1
        best = np.full(len(pos), np.inf)
        for off in range(-2, 3):                 # neighbours either side
            k = np.clip(j + off, 0, len(starts) - 1)
            before = np.maximum(0, starts[k] - pos)
            after = np.maximum(0, pos - ends[k])
            best = np.minimum(best, np.maximum(before, after) * (
                (pos < starts[k]) | (pos > ends[k])))
        d[coords.index.get_indexer(grp.index)] = best
    return d


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cm = pd.read_csv(UP / "results/faz0/condition_means_all_loci.csv.gz", index_col=0)
    ctrl = pd.read_csv(UP / "results/faz0/control_beta_all_loci.csv.gz",
                       index_col=0)["control_beta"]
    panel = pd.read_csv(UP / "results/faz0/testA_panel_floor_anatomy.csv")["cpg_id"].tolist()
    zh = pd.read_csv(UP / "data/reference/probe_qc/EPIC.hg38.zhou_annotation.tsv.gz",
                     sep="\t", compression="gzip", low_memory=False,
                     usecols=["probe_id", "chrm", "beg"]).dropna()
    zh["beg"] = zh["beg"].astype(int)
    zh = zh.set_index("probe_id")

    df = cm.join(ctrl, how="inner").join(zh, how="inner").dropna(
        subset=["baseline", "reprogrammed_intermediate", "control_beta", "chrm", "beg"])
    df["trough"] = df["baseline"] - df["reprogrammed_intermediate"]
    df["chrm"] = df["chrm"].astype(str)
    print(f"loci with methylation and coordinates: {len(df):,}")

    pk = peaks_hg38()
    for f, regions in pk.items():
        df[f"d_{f}"] = distances(regions, df[["chrm", "beg"]])
    df = df.dropna(subset=[f"d_{f}" for f in pk])
    print(f"loci with a distance on every factor: {len(df):,}")

    print("\nTrough depth versus distance to the nearest binding region")
    rows = []
    for f in pk:
        d = df[f"d_{f}"].clip(upper=MAXD)
        rho, p = stats.spearmanr(d, df["trough"])
        # the same, holding control methylation fixed, by residualising both
        lo_ = df["control_beta"]
        rx = df["trough"] - np.poly1d(np.polyfit(lo_, df["trough"], 2))(lo_)
        ry = d - np.poly1d(np.polyfit(lo_, d, 2))(lo_)
        rho_p, p_p = stats.spearmanr(ry, rx)
        rows.append({"factor": f, "spearman_rho": rho, "p": p,
                     "partial_rho_given_control_beta": rho_p, "partial_p": p_p})
        print(f"  {f:14s} rho = {rho:+.4f} (p = {p:.3g})   "
              f"partial rho = {rho_p:+.4f} (p = {p_p:.3g})")

    print("\nMean trough depth by distance band (KLF4 or OCT4)")
    bands = [(0, 0), (1, 500), (501, 1000), (1001, 5000), (5001, 20000),
             (20001, 100000), (100001, np.inf)]
    band_rows = []
    for lo_b, hi_b in bands:
        m = (df["d_KLF4_or_OCT4"] >= lo_b) & (df["d_KLF4_or_OCT4"] <= hi_b)
        if m.sum() == 0:
            continue
        lab = "inside" if hi_b == 0 else (f">{lo_b-1:,}" if hi_b == np.inf
                                          else f"{lo_b:,}-{hi_b:,}")
        band_rows.append({"band_bp": lab, "n": int(m.sum()),
                          "mean_trough": df.loc[m, "trough"].mean(),
                          "mean_control_beta": df.loc[m, "control_beta"].mean()})
        print(f"  {lab:>16s}  n={m.sum():>7,}  trough {df.loc[m,'trough'].mean():+.4f}  "
              f"control beta {df.loc[m,'control_beta'].mean():.3f}")

    print("\nWithin bands of control methylation (KLF4 or OCT4)")
    strat = []
    for lo_b, hi_b in [(0.6, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 1.0)]:
        s = df[(df["control_beta"] >= lo_b) & (df["control_beta"] < hi_b)]
        if len(s) < 500:
            continue
        near = s[s["d_KLF4_or_OCT4"] <= 1000]["trough"]
        far = s[s["d_KLF4_or_OCT4"] > 20000]["trough"]
        u, pu = stats.mannwhitneyu(near, far, alternative="greater")
        strat.append({"control_beta_band": f"{lo_b:.1f}-{hi_b:.1f}",
                      "n_near": len(near), "n_far": len(far),
                      "mean_trough_near": near.mean(), "mean_trough_far": far.mean(),
                      "mannwhitney_p": pu})
        print(f"  beta {lo_b:.1f}-{hi_b:.1f}: near (<=1 kb) n={len(near):>6,} "
              f"trough {near.mean():+.4f} | far (>20 kb) n={len(far):>6,} "
              f"trough {far.mean():+.4f} | p = {pu:.3g}")

    pan = df.loc[[c for c in panel if c in df.index]]
    rho_pan, p_pan = stats.spearmanr(pan["d_KLF4_or_OCT4"], pan["trough"])
    print(f"\nWithin the 16-CpG panel: rho = {rho_pan:+.3f} (p = {p_pan:.3f}, n = {len(pan)})")
    print(pan[["control_beta", "trough", "d_KLF4", "d_OCT4", "d_KLF4_or_OCT4"]]
          .round(3).to_string())

    keep = df[["control_beta", "baseline", "reprogrammed_intermediate", "trough",
               "d_KLF4", "d_OCT4", "d_KLF4_or_OCT4"]]
    keep.to_csv(OUT / "dose_response_per_locus.csv.gz", compression="gzip",
                float_format="%.6g")
    pd.DataFrame(rows).to_csv(OUT / "dose_response_correlations.csv", index=False)
    pd.DataFrame(band_rows).to_csv(OUT / "dose_response_distance_bands.csv", index=False)
    pd.DataFrame(strat).to_csv(OUT / "dose_response_within_beta_bands.csv", index=False)
    pan.to_csv(OUT / "panel_distances.csv")
    print(f"\nwrote four tables to {OUT}")


if __name__ == "__main__":
    main()
