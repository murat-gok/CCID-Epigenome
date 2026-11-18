#!/usr/bin/env python3
"""Regenerate the external validation from GSE142439 and test it against a
matched null.

Two things are wrong with the published version of this result. First, its
provenance file records that it was transcribed from console output rather than
produced by the archived pipeline, so it cannot currently be re-executed.
Second, its exact binomial test assumes a null of p = 0.5, which ignores that
the panel is a set of strongly hypomethylated, high-baseline loci; the honest
null is a background matched on control methylation and on discovery effect
size.

Source: Sarkar, Horvath and Sebastiano (GSE142439), four paired skin-fibroblast
donors, normal versus transiently reprogrammed, EPIC array. The four
vein-endothelial donors are left out: the discovery contrast is fibroblast.
"""

from __future__ import annotations

import gzip
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

UP = Path("/mnt/user-data/uploads/CCID-Epigenome")
OUT = Path("/mnt/user-data/outputs/external")
MATRIX = UP / "data/raw/GSE142439_series_matrix.txt.gz"
RNG = np.random.default_rng(20261118)
N_DRAWS = 10000
BETA_BIN, DELTA_BIN = 0.05, 0.05

FIBRO = {  # person: (normal, treated)
    1: ("GSM4227385", "GSM4227386"), 2: ("GSM4227387", "GSM4227388"),
    3: ("GSM4227389", "GSM4227390"), 4: ("GSM4227391", "GSM4227392"),
}


def external_delta() -> pd.Series:
    """Mean treated-minus-normal difference per locus, over the four donors."""
    cache = OUT / "gse142439_fibroblast_delta.csv.gz"
    if cache.exists():
        return pd.read_csv(cache, index_col=0)["delta"]
    with gzip.open(MATRIX, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
        hdr = [c.strip('"') for c in fh.readline().rstrip("\n").split("\t")]
        pos = {g: hdr.index(g) for p in FIBRO.values() for g in p if g in hdr}
        pairs = [(pos[n], pos[t]) for n, t in FIBRO.values()]
        out = {}
        for line in fh:
            if line.startswith("!"):
                break
            p = line.rstrip("\n").split("\t")
            cpg = p[0].strip('"')
            ds = []
            for ni, ti in pairs:
                a, b = p[ni].strip('"'), p[ti].strip('"')
                if a and b and a not in ("NA", "null") and b not in ("NA", "null"):
                    try:
                        ds.append(float(b) - float(a))
                    except ValueError:
                        pass
            if len(ds) == len(pairs):
                out[cpg] = sum(ds) / len(ds)
    s = pd.Series(out, name="delta")
    s.index.name = "cpg_id"
    OUT.mkdir(parents=True, exist_ok=True)
    s.to_csv(cache, compression="gzip")
    return s


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    panel = pd.read_csv(UP / "results/faz0/testA_panel_floor_anatomy.csv")
    panel_ids = panel["cpg_id"].tolist()
    disc = pd.read_csv(UP / "results/real/safe_osk_target_signature_O1O2_all_loci.csv.gz",
                       usecols=["cpg_id", "target_delta_beta"]).set_index("cpg_id")["target_delta_beta"]
    frozen = set(pd.read_csv(UP / "results/real/safe_osk_target_signature_O1O2.csv.gz")["cpg_id"])
    ctrl = pd.read_csv(UP / "results/faz0/control_beta_all_loci.csv.gz",
                       index_col=0)["control_beta"]

    ext = external_delta()
    print(f"GSE142439 fibroblast donors: {len(FIBRO)} pairs | "
          f"loci with a complete 4-donor delta: {len(ext):,}")

    common = ext.index.intersection(disc.index).intersection(ctrl.index)
    ok = disc.reindex(common).notna() & ctrl.reindex(common).notna() & ext.reindex(common).notna()
    dropped = int((~ok).sum())
    common = common[ok.values]
    if dropped:
        print(f"dropped {dropped:,} loci with a missing discovery delta or control beta")
    print(f"loci shared with the discovery signature: {len(common):,}")

    # ---- observed -------------------------------------------------------
    pin = [c for c in panel_ids if c in common]
    conc = int((np.sign(ext[pin]) == np.sign(disc[pin])).sum())
    rho, rho_p = stats.spearmanr(disc[pin], ext[pin])
    print(f"\npanel loci evaluable in GSE142439: {len(pin)}/16")
    print(f"direction concordant: {conc}/{len(pin)}")
    print(f"effect-size correlation: Spearman rho = {rho:.4f} (p = {rho_p:.3f})")
    p_half = stats.binomtest(conc, len(pin), 0.5, alternative="greater").pvalue
    print(f"exact binomial vs p=0.5 null: p = {p_half:.6f}")

    # global concordance rate across all shared loci - a second naive null
    gl = float((np.sign(ext[common]) == np.sign(disc[common])).mean())
    p_glob = stats.binomtest(conc, len(pin), gl, alternative="greater").pvalue
    print(f"genome-wide concordance rate: {gl:.4f} -> binomial p = {p_glob:.6f}")

    # ---- matched null ---------------------------------------------------
    pool = [c for c in common if c not in frozen and c not in set(panel_ids)]
    pool = pd.Index(pool)
    key = lambda b, d: (int(np.floor(b / BETA_BIN)), int(np.floor(abs(d) / DELTA_BIN)))
    buckets = defaultdict(list)
    for c, b, d in zip(pool, ctrl[pool].values, disc[pool].values):
        buckets[key(b, d)].append(c)

    sign_ok = {c: (np.sign(ext[c]) == np.sign(disc[c])) for c in pool}
    cand_per_slot = []
    for c in pin:
        k = key(ctrl[c], disc[c])
        cand = buckets.get(k, [])
        if len(cand) < 30:                       # widen the effect-size bin only
            cand = sum((buckets.get((k[0], k[1] + o), []) for o in (-1, 0, 1)), [])
        if len(cand) < 30:
            cand = sum((buckets.get((k[0] + o1, k[1] + o2), [])
                        for o1 in (-1, 0, 1) for o2 in (-1, 0, 1)), [])
        cand_per_slot.append(np.array([sign_ok[c2] for c2 in cand], dtype=bool))
        print(f"  {c}: control beta {ctrl[c]:.3f}, discovery delta {disc[c]:+.3f} "
              f"-> {len(cand):,} matched loci, {cand_per_slot[-1].mean():.3f} concordant")

    null = np.zeros(N_DRAWS, dtype=int)
    for arr in cand_per_slot:
        null += RNG.choice(arr, size=N_DRAWS)
    p_matched = (1 + int((null >= conc).sum())) / (1 + N_DRAWS)
    print(f"\nmatched null (n={N_DRAWS}): mean {null.mean():.2f}/{len(pin)} concordant, "
          f"sd {null.std(ddof=1):.2f}")
    print(f"observed {conc}/{len(pin)} -> empirical p = {p_matched:.4g}")

    pd.DataFrame({"statistic": ["n_evaluable", "concordant", "spearman_rho", "spearman_p",
                                "p_null_0.5", "genomewide_concordance", "p_null_genomewide",
                                "matched_null_mean", "matched_null_sd", "p_matched",
                                "n_draws", "beta_bin", "delta_bin"],
                  "value": [len(pin), conc, rho, rho_p, p_half, gl, p_glob,
                            null.mean(), null.std(ddof=1), p_matched,
                            N_DRAWS, BETA_BIN, DELTA_BIN]}).to_csv(
        OUT / "external_validation_matched_null.csv", index=False)
    pd.DataFrame({"cpg_id": pin, "discovery_delta": disc[pin].values,
                  "external_delta": ext[pin].values,
                  "concordant": (np.sign(ext[pin]) == np.sign(disc[pin])).values,
                  "control_beta": ctrl[pin].values}).to_csv(
        OUT / "external_validation_per_locus.csv", index=False)
    print(f"\nwrote {OUT/'external_validation_matched_null.csv'} and per-locus table")


if __name__ == "__main__":
    main()
