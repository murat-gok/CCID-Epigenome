#!/usr/bin/env python3
"""Baseline-matched background, an operational state-specificity test, and
cluster-level independence for the 16-CpG panel.

Why each piece exists:

  The published background was a systematic 1% sample sitting at control
  methylation ~0.47, while the panel sits at 0.766. Ceiling effects and
  regression to the mean therefore confound any comparison of methylation
  change between them. Here every background locus is drawn from the same
  control-methylation bin as the panel locus it stands in for.

  The state-specificity challenge, which the manuscript proposes as its
  generalisable contribution, is given a statistic, a null and a threshold:

      T = beta(baseline) - beta(intermediate)           depth of the trough
      R = 1 - |beta(iPSC) - beta(baseline)| / T         return to baseline

  A locus set qualifies as a transient-state marker when T exceeds the 99th
  percentile of the matched-background null and R exceeds 0.8.

  The 16 CpGs are also collapsed into genomic clusters, so that counting
  arguments are not inflated by co-located probes.

Run from the repository root, after scripts/cache_control_beta.py.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(".")
OUT = ROOT / "results/faz0"
RNG = np.random.default_rng(20261118)
N_DRAWS = 1000
BIN = 0.02
CLUSTER_BP = 5000


def main() -> None:
    panel = pd.read_csv(OUT / "testA_panel_floor_anatomy.csv")["cpg_id"].tolist()
    frozen = set(pd.read_csv(
        ROOT / "results/real/safe_osk_target_signature_O1O2.csv.gz")["cpg_id"])
    meta = pd.read_csv(ROOT / "data/processed/fibroblast_methylation_metadata.csv")
    ctrl = pd.read_csv(OUT / "control_beta_all_loci.csv.gz",
                       index_col=0)["control_beta"]

    excl = frozen | set(panel)
    pool = ctrl.drop(index=[i for i in excl if i in ctrl.index]).dropna()
    bins = np.floor(pool.values / BIN).astype(int)
    by_bin: dict[int, np.ndarray] = {}
    order = np.argsort(bins, kind="stable")
    b_sorted, ids_sorted = bins[order], pool.index.values[order]
    edges = np.flatnonzero(np.diff(b_sorted)) + 1
    for chunk in np.split(np.arange(len(b_sorted)), edges):
        by_bin[int(b_sorted[chunk[0]])] = ids_sorted[chunk]

    draws = np.empty((N_DRAWS, len(panel)), dtype=object)
    for j, p in enumerate(panel):
        b = int(np.floor(ctrl.get(p, np.nan) / BIN))
        cand = by_bin.get(b)
        if cand is None or len(cand) < 50:          # widen if the bin is thin
            cand = np.concatenate([by_bin.get(k, np.empty(0, dtype=object))
                                   for k in (b - 1, b, b + 1)])
        draws[:, j] = RNG.choice(cand, size=N_DRAWS, replace=True)
    need = set(panel) | set(np.unique(draws))
    print(f"panel control beta {ctrl.reindex(panel).mean():.4f} | "
          f"pool {len(pool):,} | distinct background loci drawn {len(need)-len(panel):,}")

    rows, cols = {}, None
    with gzip.open(ROOT / "data/processed/fibroblast_beta.csv.gz", "rt") as fh:
        cols = fh.readline().rstrip("\n").split(",")[1:]
        for line in fh:
            p = line.rstrip("\n").split(",")
            if p[0] in need:
                rows[p[0]] = p[1:]
    beta = pd.DataFrame.from_dict(rows, orient="index", columns=cols).apply(
        pd.to_numeric, errors="coerce")
    print(f"read {len(beta):,} loci x {beta.shape[1]} samples")

    m = meta.set_index("sample_id").reindex(beta.columns)
    cond = m["condition"].values

    def stats(ids):
        sub = beta.reindex([i for i in ids if i in beta.index]).to_numpy(dtype=float)
        per_sample = np.nanmean(sub, axis=0)
        g = {c: np.nanmean(per_sample[cond == c]) for c in np.unique(cond.astype(str))}
        base, inter, ips = g.get("baseline"), g.get("reprogrammed_intermediate"), g.get("ipsc")
        T = base - inter
        R = 1 - abs(ips - base) / T if T else np.nan
        return base, inter, ips, T, R

    base, inter, ips, T_obs, R_obs = stats(panel)
    print(f"\npanel  baseline {base:.4f}  intermediate {inter:.4f}  iPSC {ips:.4f}"
          f"  T {T_obs:.4f}  R {R_obs:.4f}")

    nT, nR = np.empty(N_DRAWS), np.empty(N_DRAWS)
    for i in range(N_DRAWS):
        _, _, _, nT[i], nR[i] = stats(list(draws[i]))
    q99 = float(np.nanquantile(nT, 0.99))
    p_T = (1 + int(np.nansum(nT >= T_obs))) / (1 + N_DRAWS)
    ok = (T_obs > q99) and (R_obs > 0.8)
    print(f"matched null (n={N_DRAWS})  T mean {np.nanmean(nT):.4f} "
          f"sd {np.nanstd(nT, ddof=1):.4f}  99th pct {q99:.4f}  p = {p_T:.3g}")
    print(f"                             R mean {np.nanmean(nR):.4f}")
    print(f"state-specificity criterion (T > 99th pct and R > 0.8): "
          f"{'MET - transient-state marker' if ok else 'NOT met'}")

    zh = pd.read_csv(ROOT / "data/reference/probe_qc/EPIC.hg38.zhou_annotation.tsv.gz",
                     sep="\t", compression="gzip", low_memory=False)
    zh = zh[zh["probe_id"].isin(panel)][["probe_id", "chrm", "beg"]].dropna()
    zh = zh.sort_values(["chrm", "beg"]).reset_index(drop=True)
    cl, cid, prev = [], 0, None
    for _, r in zh.iterrows():
        if not (prev is not None and r["chrm"] == prev["chrm"]
                and abs(r["beg"] - prev["beg"]) < CLUSTER_BP):
            cid += 1
        cl.append(cid); prev = r
    zh["cluster"] = cl
    print(f"\nM3  {len(zh)} panel CpGs with coordinates collapse to "
          f"{zh['cluster'].nunique()} clusters at {CLUSTER_BP} bp")

    pd.DataFrame({"statistic": ["panel_baseline", "panel_intermediate", "panel_ipsc",
                                "T_observed", "R_observed", "T_null_mean", "T_null_sd",
                                "T_null_99pct", "T_empirical_p", "R_null_mean",
                                "n_clusters", "panel_control_beta", "bin_width", "n_draws"],
                  "value": [base, inter, ips, T_obs, R_obs, np.nanmean(nT),
                            np.nanstd(nT, ddof=1), q99, p_T, np.nanmean(nR),
                            zh["cluster"].nunique(), ctrl.reindex(panel).mean(),
                            BIN, N_DRAWS]}).to_csv(
        OUT / "matched_background_state_specificity.csv", index=False)
    zh.to_csv(OUT / "panel_genomic_clusters.csv", index=False)
    print(f"\nwrote {OUT/'matched_background_state_specificity.csv'} and "
          f"{OUT/'panel_genomic_clusters.csv'}")


if __name__ == "__main__":
    main()
