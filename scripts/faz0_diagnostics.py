#!/usr/bin/env python3
"""
CCID-Epigenome | FAZ 0 diagnostic gate
======================================

Answers one question before any manuscript revision begins:

    Is the frozen 16-CpG panel a biological signal, or an artefact of the
    discovery ranking statistic?

Background
----------
`src/ccid/signature.py::_summarise_donor_effects` ranks loci with

    rank_score = |mean_delta| * sign_consistency * sqrt(n_donors)
                 / max(standard_error, 0.01)

With n_donors = 2:
    donor_sd (ddof=1) = |d1 - d2| / sqrt(2)
    standard_error    = |d1 - d2| / 2

so the 0.01 floor is reached whenever the two discovery donors agree to
within 0.02 beta. Inside that regime the variance penalty vanishes and the
ranking collapses to |mean_delta| alone, multiplied by a constant ~141.
Every selected locus also has sign_consistency == 1.0 by construction,
because with two donors the only attainable values are 0.5 and 1.0 and the
threshold is 0.67.

This script quantifies how much of the frozen panel lives in that regime,
and whether the panel is simply the extreme tail of a genome-wide
methylation shift.

Required inputs (all produced by scripts/03_build_target_signature.py)
----------------------------------------------------------------------
  --all-loci        <stem>_all_loci.csv.gz
  --donor-effects   <stem>_donor_effects.csv.gz
  --shortlist       data/reference/frozen_primary_shortlist_v090.csv

Optional
--------
  --signature       the frozen top-500 signature CSV (sharpens TEST B)
  --beta            raw beta matrix        (enables TEST D)
  --pairs           <stem>_matched_pairs.csv (enables TEST D)
  --mask            probe mask list(s)    (enables TEST F) - repeatable
  --wcgw            solo-WCGW probe list  (enables TEST F)

Usage
-----
  python faz0_diagnostics.py \
      --all-loci       results/real/.../signature_all_loci.csv.gz \
      --donor-effects  results/real/.../signature_donor_effects.csv.gz \
      --shortlist      data/reference/frozen_primary_shortlist_v090.csv \
      --outdir         results/faz0

Nothing is overwritten; every table is written to --outdir as CSV.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SE_FLOOR_PUBLISHED = 0.01
RULE = "-" * 78


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def banner(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def read_table(path: Path, what: str) -> pd.DataFrame:
    if not path.exists():
        sys.exit(f"[FATAL] {what} not found: {path}")
    frame = pd.read_csv(path)
    if frame.empty:
        sys.exit(f"[FATAL] {what} is empty: {path}")
    return frame


def read_cpg_list(path: Path) -> set[str]:
    """Read a CpG identifier list from a one-per-line or single-column file."""
    frame = pd.read_csv(path, sep=None, engine="python")
    for col in frame.columns:
        values = frame[col].astype(str)
        if values.str.startswith("cg").any():
            return set(values.str.split("_").str[0])
    # fall back: treat the file as a bare list, header included
    raw = [x.strip() for x in path.read_text(encoding="utf-8").splitlines() if x.strip()]
    return {x.split("_")[0] for x in raw if x.startswith("cg")}


def base_id(series: pd.Series) -> pd.Series:
    """Normalise EPICv2 suffixed identifiers (cg00000029_BC21 -> cg00000029)."""
    return series.astype(str).str.split("_").str[0]


def pct_rank(values: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Percentile of each value within the reference distribution (0-100)."""
    reference = np.sort(reference[np.isfinite(reference)])
    if reference.size == 0:
        return np.full(values.shape, np.nan)
    idx = np.searchsorted(reference, values, side="right")
    return 100.0 * idx / reference.size


def recompute_rank(frame: pd.DataFrame, variant: str, floor: float | None) -> pd.Series:
    """Recompute rank_score under an alternative variance treatment."""
    d = frame["target_delta_beta"].abs()
    sc = frame["sign_consistency"]
    n = frame["n_donors"].clip(lower=1)
    se = frame["standard_error"]

    if variant == "floor":
        return d * sc * np.sqrt(n) / np.maximum(se, floor)
    if variant == "nofloor":
        # true t-like statistic; SE==0 -> infinite, keep as +inf so it is visible
        with np.errstate(divide="ignore", invalid="ignore"):
            return d * sc * np.sqrt(n) / se.replace(0.0, np.nan)
    if variant == "delta_only":
        return d
    if variant == "delta_x_sign":
        return d * sc
    raise ValueError(variant)


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return float("nan")
    return len(a & b) / len(a | b)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(
        description="FAZ 0 diagnostic gate for the frozen MPTR CpG panel",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--all-loci", required=True, type=Path)
    ap.add_argument("--donor-effects", required=True, type=Path)
    ap.add_argument("--shortlist", required=True, type=Path)
    ap.add_argument("--signature", type=Path, help="frozen top-500 signature CSV")
    ap.add_argument("--beta", type=Path, help="raw beta matrix (enables TEST D)")
    ap.add_argument("--pairs", type=Path, help="matched pairs CSV (enables TEST D)")
    ap.add_argument("--mask", type=Path, action="append", default=[],
                    help="probe mask list; repeat for several lists")
    ap.add_argument("--wcgw", type=Path, help="solo-WCGW probe list")
    ap.add_argument("--top-n", type=int, default=500, help="frozen list size")
    ap.add_argument("--chunksize", type=int, default=20000, help="beta matrix chunk size")
    ap.add_argument("--outdir", type=Path, default=Path("results/faz0"))
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    verdicts: dict[str, str] = {}
    stats: dict[str, object] = {}

    # ---------------------------------------------------------------- load
    banner("LOADING")
    all_loci = read_table(args.all_loci, "all_loci table")
    needed = {"cpg_id", "target_delta_beta", "standard_error",
              "sign_consistency", "n_donors"}
    missing = needed - set(all_loci.columns)
    if missing:
        sys.exit(f"[FATAL] all_loci is missing columns: {sorted(missing)}")
    all_loci["cpg_base"] = base_id(all_loci["cpg_id"])
    if "rank_score" not in all_loci.columns:
        print("[warn] rank_score absent; recomputing as published")
        all_loci["rank_score"] = recompute_rank(all_loci, "floor", SE_FLOOR_PUBLISHED)

    donor_eff = read_table(args.donor_effects, "donor_effects table")
    id_col = donor_eff.columns[0]
    donor_eff = donor_eff.rename(columns={id_col: "cpg_id"})
    donor_eff["cpg_base"] = base_id(donor_eff["cpg_id"])
    donor_cols = [c for c in donor_eff.columns
                  if c not in {"cpg_id", "cpg_base"}
                  and pd.api.types.is_numeric_dtype(donor_eff[c])]

    shortlist = read_table(args.shortlist, "frozen shortlist")
    if "cpg_id" not in shortlist.columns:
        sys.exit("[FATAL] shortlist lacks a cpg_id column")
    panel = base_id(shortlist["cpg_id"]).tolist()

    print(f"all_loci rows        : {len(all_loci):,}")
    print(f"donor columns        : {donor_cols}")
    print(f"frozen panel size    : {len(panel)}")
    stats["n_all_loci"] = int(len(all_loci))
    stats["donor_columns"] = donor_cols
    stats["panel_size"] = len(panel)

    frozen500: set[str] | None = None
    if args.signature and args.signature.exists():
        sig = pd.read_csv(args.signature)
        if "cpg_id" in sig.columns:
            frozen500 = set(base_id(sig["cpg_id"]))
            print(f"frozen signature size: {len(frozen500)}")

    idx = all_loci.set_index("cpg_base")
    panel_present = [c for c in panel if c in idx.index]
    if len(panel_present) < len(panel):
        print(f"[warn] {len(panel) - len(panel_present)} panel CpGs absent from all_loci")

    # ------------------------------------------------- TEST A: floor anatomy
    banner("TEST A | Donor-level anatomy of the frozen panel")
    print("Is the panel inside the SE-floor regime, where the variance")
    print("penalty vanishes and ranking reduces to |delta| alone?")
    print(RULE)

    de_idx = donor_eff.set_index("cpg_base")
    rows = []
    for cpg in panel_present:
        a = idx.loc[cpg]
        if isinstance(a, pd.DataFrame):
            a = a.iloc[0]
        deltas = {}
        if cpg in de_idx.index:
            row = de_idx.loc[cpg]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            deltas = {c: float(row[c]) for c in donor_cols if np.isfinite(row[c])}
        vals = np.array(list(deltas.values()), dtype=float)
        spread = float(vals.max() - vals.min()) if vals.size >= 2 else np.nan
        se = float(a["standard_error"])
        rows.append({
            "cpg_id": cpg,
            "delta_beta": float(a["target_delta_beta"]),
            "sign_consistency": float(a["sign_consistency"]),
            "n_donors": int(a["n_donors"]),
            **{f"delta_{k}": v for k, v in deltas.items()},
            "donor_spread": spread,
            "standard_error": se,
            "se_at_floor": bool(se < SE_FLOOR_PUBLISHED),
            "rank_score_published": float(a["rank_score"]),
        })

    testA = pd.DataFrame(rows).sort_values("rank_score_published", ascending=False)
    testA.to_csv(args.outdir / "testA_panel_floor_anatomy.csv", index=False)

    n_floor = int(testA["se_at_floor"].sum())
    n_tot = len(testA)
    floor_all = int((all_loci["standard_error"] < SE_FLOOR_PUBLISHED).sum())
    frac_all = floor_all / len(all_loci)

    with pd.option_context("display.width", 200, "display.max_columns", 30):
        cols = ["cpg_id", "delta_beta", "donor_spread", "standard_error",
                "se_at_floor", "rank_score_published"]
        print(testA[cols].to_string(index=False, float_format=lambda v: f"{v: .6f}"))

    print(RULE)
    print(f"Panel CpGs at the SE floor      : {n_floor}/{n_tot}")
    print(f"All loci at the SE floor        : {floor_all:,}/{len(all_loci):,} ({frac_all:.3%})")
    if n_tot:
        enrich = (n_floor / n_tot) / frac_all if frac_all > 0 else float("inf")
        print(f"Enrichment of the panel at floor: {enrich:.1f}x")
        stats["floor_enrichment"] = None if not np.isfinite(enrich) else round(enrich, 2)
    stats["panel_at_floor"] = n_floor
    stats["all_loci_at_floor_fraction"] = round(frac_all, 6)

    if n_tot and n_floor / n_tot >= 0.5:
        verdicts["A"] = ("RED | the floor dominates the panel: ranking inside it "
                         "reduces to |delta| with no variance control")
    elif n_tot and n_floor / n_tot >= 0.25:
        verdicts["A"] = "AMBER | a substantial minority of the panel sits at the floor"
    else:
        verdicts["A"] = "GREEN | the panel is not floor-driven"
    print(f"\nVERDICT A: {verdicts['A']}")

    # ------------------------------------------- TEST B: ranking sensitivity
    banner("TEST B | Sensitivity of the frozen list to the variance treatment")
    print("The top-500 is the input to every downstream stage. If it moves")
    print("when the SE floor moves, everything downstream moves with it.")
    print(RULE)

    variants = [
        ("published (floor 0.01)", "floor", 0.01),
        ("floor 0.02", "floor", 0.02),
        ("floor 0.05", "floor", 0.05),
        ("no floor (true t)", "nofloor", None),
        ("|delta| only", "delta_only", None),
        ("|delta| x sign", "delta_x_sign", None),
    ]

    ref_top: set[str] | None = None
    summary = []
    rank_map: dict[str, pd.Series] = {}
    for label, variant, floor in variants:
        score = recompute_rank(all_loci, variant, floor)
        order = (pd.DataFrame({"cpg_base": all_loci["cpg_base"], "s": score})
                 .sort_values("s", ascending=False, na_position="last"))
        top = set(order["cpg_base"].head(args.top_n))
        rank_map[label] = pd.Series(
            np.arange(1, len(order) + 1), index=order["cpg_base"].values
        )
        if ref_top is None:
            ref_top = top
        kept = sum(1 for c in panel_present if c in top)
        summary.append({
            "variant": label,
            f"jaccard_top{args.top_n}_vs_published": round(jaccard(ref_top, top), 4),
            f"panel_retained_in_top{args.top_n}": f"{kept}/{len(panel_present)}",
            "panel_retained_frac": round(kept / max(len(panel_present), 1), 4),
        })

    testB = pd.DataFrame(summary)
    testB.to_csv(args.outdir / "testB_ranking_sensitivity.csv", index=False)
    print(testB.to_string(index=False))

    # per-candidate rank movement
    move = pd.DataFrame({"cpg_id": panel_present})
    for label in rank_map:
        move[label] = [
            int(rank_map[label].get(c, -1)) for c in panel_present
        ]
    move.to_csv(args.outdir / "testB_panel_rank_movement.csv", index=False)
    print()
    print("Per-candidate rank under each variant (-1 = absent):")
    with pd.option_context("display.width", 220, "display.max_columns", 30):
        print(move.to_string(index=False))

    if frozen500 is not None:
        j_actual = jaccard(frozen500, ref_top)
        print(f"\nSanity check: recomputed published top-{args.top_n} vs the stored "
              f"frozen signature, Jaccard = {j_actual:.4f}")
        stats["sanity_jaccard_vs_stored_signature"] = round(j_actual, 4)
        if j_actual < 0.95:
            print("[warn] recomputation does not reproduce the stored frozen list; "
                  "check that --all-loci matches the run that produced --signature")

    worst = testB.loc[1:, f"jaccard_top{args.top_n}_vs_published"].min()
    stats["min_jaccard_across_variants"] = float(worst)
    if worst < 0.50:
        verdicts["B"] = ("RED | the frozen list is not robust to the variance "
                         "treatment; the 500 is a choice of statistic, not a finding")
    elif worst < 0.75:
        verdicts["B"] = "AMBER | the frozen list shifts materially under reasonable alternatives"
    else:
        verdicts["B"] = "GREEN | the frozen list is stable across variance treatments"
    print(f"\nVERDICT B: {verdicts['B']}")

    # -------------------------------------- TEST C: global direction / drift
    banner("TEST C | Genome-wide delta distribution and the true null direction")
    print("If the methylome shifts globally at day 10, directional concordance")
    print("in any external reprogramming dataset is near-guaranteed and the")
    print("p=0.5 binomial null used for the 14/16 result is invalid.")
    print(RULE)

    d_all = all_loci["target_delta_beta"].to_numpy(dtype=float)
    d_all = d_all[np.isfinite(d_all)]
    hypo_all = float((d_all < 0).mean())

    def describe(name: str, vals: np.ndarray) -> dict:
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            return {"set": name, "n": 0}
        return {
            "set": name,
            "n": int(vals.size),
            "mean_delta": round(float(vals.mean()), 6),
            "median_delta": round(float(np.median(vals)), 6),
            "frac_hypo": round(float((vals < 0).mean()), 4),
            "min": round(float(vals.min()), 6),
            "max": round(float(vals.max()), 6),
        }

    sets = [describe("all loci", d_all)]
    if frozen500 is not None:
        m = all_loci["cpg_base"].isin(frozen500)
        sets.append(describe(f"frozen top-{args.top_n}",
                             all_loci.loc[m, "target_delta_beta"].to_numpy(float)))
    m16 = all_loci["cpg_base"].isin(set(panel_present))
    panel_deltas = all_loci.loc[m16, "target_delta_beta"].to_numpy(float)
    sets.append(describe("frozen panel (16)", panel_deltas))

    testC = pd.DataFrame(sets)
    testC.to_csv(args.outdir / "testC_global_delta_distribution.csv", index=False)
    print(testC.to_string(index=False))

    thr = float(np.nanmin(np.abs(panel_deltas))) if panel_deltas.size else np.nan
    n_extreme = int((np.abs(d_all) >= thr).sum()) if np.isfinite(thr) else -1
    pr = pct_rank(np.abs(panel_deltas), np.abs(d_all)) if panel_deltas.size else np.array([])

    print(RULE)
    print(f"Genome-wide hypomethylated fraction : {hypo_all:.4f}")
    print(f"  -> this, not 0.5, is the null direction probability to use")
    print(f"     for the GSE142439 concordance test (TEST of Faz 2.1)")
    if np.isfinite(thr):
        print(f"\nSmallest |delta| in the panel        : {thr:.4f}")
        print(f"Loci genome-wide with |delta| >= that: {n_extreme:,} "
              f"({n_extreme / len(d_all):.4%})")
        print(f"Panel |delta| percentile range       : "
              f"{pr.min():.3f} - {pr.max():.3f}")

    # what would 14/16 mean under the real null?
    if 0.0 < hypo_all < 1.0:
        from math import comb
        p0 = hypo_all
        p_naive = sum(comb(16, k) * 0.5 ** 16 for k in range(14, 17))
        p_real = sum(comb(16, k) * p0 ** k * (1 - p0) ** (16 - k) for k in range(14, 17))
        print(f"\n14/16 one-sided binomial p, null=0.5      : {p_naive:.6f}")
        print(f"14/16 one-sided binomial p, null={p0:.4f} : {p_real:.6f}")
        print("  (illustrative only - the discovery-set hypo fraction is a proxy;")
        print("   Faz 2.1 must use the external dataset's own global drift and a")
        print("   property-matched permutation null)")
        stats["p_14of16_null_0.5"] = round(p_naive, 8)
        stats["p_14of16_null_global"] = round(p_real, 8)

    stats["genomewide_hypo_fraction"] = round(hypo_all, 6)
    stats["panel_min_abs_delta"] = None if not np.isfinite(thr) else round(thr, 6)
    stats["n_loci_at_or_above_panel_min_abs_delta"] = n_extreme

    if hypo_all >= 0.70 or hypo_all <= 0.30:
        verdicts["C"] = (f"RED | strong genome-wide directional bias "
                         f"(hypo={hypo_all:.3f}); the p=0.5 null is invalid and the "
                         f"14/16 result must be recomputed against a matched null")
    elif hypo_all >= 0.60 or hypo_all <= 0.40:
        verdicts["C"] = f"AMBER | appreciable global drift (hypo={hypo_all:.3f})"
    else:
        verdicts["C"] = f"GREEN | no dominant global direction (hypo={hypo_all:.3f})"
    print(f"\nVERDICT C: {verdicts['C']}")

    # ------------------------------- TEST D: control-sample variance (opt.)
    banner("TEST D | Baseline-sample variance percentile of the panel")
    if not (args.beta and args.pairs and args.beta.exists() and args.pairs.exists()):
        print("SKIPPED - supply --beta and --pairs to run this test.")
        print("It asks whether the panel was selected for technical instability:")
        print("variance of each probe across BASELINE samples only.")
        verdicts["D"] = "SKIPPED"
    else:
        pairs = pd.read_csv(args.pairs)
        base_cols = None
        for cand in ("baseline_sample_id", "control_sample_id", "baseline_sample",
                     "control_sample"):
            if cand in pairs.columns:
                base_cols = cand
                break
        if base_cols is None:
            print(f"[warn] no baseline sample column in {args.pairs.name}; "
                  f"columns = {list(pairs.columns)}")
            verdicts["D"] = "SKIPPED (no baseline column)"
        else:
            samples = sorted(set(pairs[base_cols].astype(str)))
            print(f"baseline samples: {len(samples)}")
            panel_set = set(panel_present)
            var_panel: dict[str, float] = {}
            all_vars: list[np.ndarray] = []
            reader = pd.read_csv(args.beta, chunksize=args.chunksize, sep=None,
                                 engine="python")
            for chunk in reader:
                first = chunk.columns[0]
                chunk["cpg_base"] = base_id(chunk[first])
                use = [s for s in samples if s in chunk.columns]
                if not use:
                    sys.exit("[FATAL] baseline samples absent from the beta matrix")
                vals = chunk[use].apply(pd.to_numeric, errors="coerce")
                v = vals.var(axis=1, ddof=1).to_numpy(dtype=float)
                all_vars.append(v)
                hit = chunk["cpg_base"].isin(panel_set).to_numpy()
                for cpg, vv in zip(chunk.loc[hit, "cpg_base"], v[hit]):
                    var_panel[str(cpg)] = float(vv)
            ref = np.concatenate(all_vars)
            rows = [{"cpg_id": c,
                     "baseline_variance": var_panel.get(c, np.nan),
                     "percentile": float(pct_rank(np.array([var_panel[c]]), ref)[0])
                     if c in var_panel else np.nan}
                    for c in panel_present]
            testD = pd.DataFrame(rows).sort_values("percentile", ascending=False)
            testD.to_csv(args.outdir / "testD_baseline_variance.csv", index=False)
            print(testD.to_string(index=False, float_format=lambda v: f"{v: .6f}"))
            med = float(np.nanmedian(testD["percentile"]))
            print(f"\nMedian baseline-variance percentile of the panel: {med:.2f}")
            stats["panel_median_baseline_variance_percentile"] = round(med, 2)
            if med >= 95:
                verdicts["D"] = ("RED | the panel sits in the top 5% of baseline "
                                 "variance: selected for instability")
            elif med >= 80:
                verdicts["D"] = "AMBER | the panel is variance-enriched at baseline"
            else:
                verdicts["D"] = "GREEN | the panel is not unusually variable at baseline"
            print(f"\nVERDICT D: {verdicts['D']}")

    # ------------------------------------ TEST F: mask / solo-WCGW overlap
    banner("TEST F | Probe-mask and solo-WCGW / PMD overlap")
    if not args.mask and not args.wcgw:
        print("SKIPPED - supply --mask and/or --wcgw.")
        print("Reference resources to fetch:")
        print("  Zhou/Laird/Shen 2017 NAR 45:e22  MASK.general")
        print("    https://zwdzwd.github.io/InfiniumAnnotation")
        print("  Pidsley 2016 Genome Biol 17:208  (cross-reactive + polymorphic)")
        print("  McCartney 2016 Genomics Data 9:22-24")
        print("  Chen 2013 (450K cross-reactive)")
        print("  Zhou 2018 Nat Genet 50:591-602   solo-WCGW probes + common PMDs")
        print("  Endicott 2022 Nat Commun 13:6659 RepliTali CpGs")
        verdicts["F"] = "SKIPPED"
    else:
        rows = []
        panel_set = set(panel_present)
        for mpath in args.mask:
            if not mpath.exists():
                print(f"[warn] mask list not found: {mpath}")
                continue
            masked = read_cpg_list(mpath)
            in_panel = sorted(panel_set & masked)
            in_500 = len(frozen500 & masked) if frozen500 else -1
            rows.append({
                "resource": mpath.name,
                "type": "mask",
                "n_in_resource": len(masked),
                "panel_hits": len(in_panel),
                f"top{args.top_n}_hits": in_500,
                "panel_cpgs": ";".join(in_panel),
            })
            print(f"{mpath.name}: {len(in_panel)}/{len(panel_set)} panel CpGs masked")
            if in_panel:
                print(f"  -> {', '.join(in_panel)}")
        if args.wcgw and args.wcgw.exists():
            wc = read_cpg_list(args.wcgw)
            in_panel = sorted(panel_set & wc)
            in_500 = len(frozen500 & wc) if frozen500 else -1
            base_rate = len(wc & set(all_loci["cpg_base"])) / len(all_loci)
            rows.append({
                "resource": args.wcgw.name,
                "type": "solo-WCGW",
                "n_in_resource": len(wc),
                "panel_hits": len(in_panel),
                f"top{args.top_n}_hits": in_500,
                "panel_cpgs": ";".join(in_panel),
            })
            print(f"\n{args.wcgw.name}: {len(in_panel)}/{len(panel_set)} panel CpGs "
                  f"are solo-WCGW")
            print(f"  genome-wide base rate on this platform: {base_rate:.4%}")
            if len(panel_set):
                obs = len(in_panel) / len(panel_set)
                if base_rate > 0:
                    print(f"  enrichment: {obs / base_rate:.1f}x")
                    stats["wcgw_enrichment"] = round(obs / base_rate, 2)
            if in_panel:
                print(f"  -> {', '.join(in_panel)}")
        if rows:
            pd.DataFrame(rows).to_csv(args.outdir / "testF_mask_wcgw_overlap.csv",
                                      index=False)
        total_masked = sum(r["panel_hits"] for r in rows if r["type"] == "mask")
        wcgw_hits = sum(r["panel_hits"] for r in rows if r["type"] == "solo-WCGW")
        if total_masked >= 0.25 * max(len(panel_set), 1):
            verdicts["F"] = (f"RED | {total_masked} panel CpGs are on exclusion "
                             f"lists; discovery must be rerun after masking")
        elif total_masked > 0:
            verdicts["F"] = f"AMBER | {total_masked} panel CpGs are masked"
        else:
            verdicts["F"] = "GREEN | no panel CpG is on the supplied mask lists"
        if wcgw_hits:
            verdicts["F"] += f" | {wcgw_hits} are solo-WCGW (replication-coupled)"
        print(f"\nVERDICT F: {verdicts['F']}")

    # ------------------------------------------------------------- decision
    banner("FAZ 0 DECISION GATE")
    for key in sorted(verdicts):
        print(f"  TEST {key}: {verdicts[key]}")

    reds = [k for k, v in verdicts.items() if v.startswith("RED")]
    ambers = [k for k, v in verdicts.items() if v.startswith("AMBER")]
    print(RULE)
    if "F" in reds or ("A" in reds and "B" in reds):
        route = ("ROUTE A - artefact-dominated. Rerun discovery with variance "
                 "moderation and probe masking before touching the manuscript.")
    elif reds:
        route = ("ROUTE C - mixed. Separate a defensible sub-panel using the "
                 "independent-dataset reliability filter (GSE179847 FDR and donor "
                 "consistency); move the rest to supplementary.")
    elif ambers:
        route = ("ROUTE C/B - mixed to replication-coupled. Proceed to Faz 1.3 "
                 "baseline comparison and Faz 0.6 PMD annotation before reframing.")
    else:
        route = ("ROUTE B - the panel survives. Reframe around the "
                 "replication-coupled mechanism and proceed to Faz 2.")
    print(route)
    print(RULE)
    print("Regardless of route, these remain mandatory:")
    print("  Faz 2.1  property-matched permutation null for every concordance test")
    print("  Faz 3.1  re-express lifespan slopes per population doubling, not per")
    print("           passage, before any 21% vs 3% oxygen comparison")
    print("  Faz 4.1  rerun scripts/04l_integrate_orthogonal_validation.py; the")
    print("           manuscript currently quotes PRE-POST5 tier counts (12/2/2)")

    stats["verdicts"] = verdicts
    stats["route"] = route
    (args.outdir / "faz0_summary.json").write_text(
        json.dumps(stats, indent=2), encoding="utf-8")
    print(f"\nTables and faz0_summary.json written to: {args.outdir}")


if __name__ == "__main__":
    main()
