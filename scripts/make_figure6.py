#!/usr/bin/env python3
"""Figure 8 for the Functional & Integrative Genomics submission.

Panel a: per-sample mean methylation of the 16-CpG panel and of a background
         sample of loci, by experimental condition.
Panel b: group means by day within each condition, showing that the panel
         recovers over the time course while the background stays flat.

Design constraints:
  - a single y-axis, no dual scales
  - identity carried by marker shape AND fill, never colour alone, so the figure
    survives greyscale printing and colour-vision deficiency
  - no transparency, so the EPS is faithful to the raster versions
  - thin marks, recessive axes, no competing gridlines
  - 174 mm width, well inside the 234 mm height limit
  - vector EPS/PDF plus a 600 dpi LZW TIFF
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

SRC = Path(os.environ.get("FIG8_SRC", "results/faz0/panel_trajectory_per_sample.csv"))
OUT = Path(os.environ.get("FIG8_OUT", "manuscript/figures/submission"))

MM = 1 / 25.4
WIDTH = 174 * MM
HEIGHT = 85 * MM

ORDER = ["baseline", "failed", "reprogrammed_intermediate", "reprogrammed", "ipsc"]
SHORT = {
    "baseline": "Baseline",
    "failed": "Failed",
    "reprogrammed_intermediate": "Intermediate",
    "reprogrammed": "Reprogrammed",
    "ipsc": "iPSC",
}
# panel b has four groups side by side, so it uses abbreviated group labels
TIGHT = {"baseline": "Baseline", "failed": "Failed",
         "reprogrammed_intermediate": "Interm.", "reprogrammed": "Reprog."}
INK = "#1a1a1a"
GREY = "#8c8c8c"

PANEL_KW = dict(marker="o", color=INK, markerfacecolor=INK, markeredgecolor="white",
                markeredgewidth=0.5, markersize=3.6, linestyle="none")
BG_KW = dict(marker="s", color=GREY, markerfacecolor="white", markeredgecolor=GREY,
             markeredgewidth=0.8, markersize=3.2, linestyle="none")
PANEL_LINE_KW = dict(marker="o", color=INK, markerfacecolor=INK, markeredgecolor="white",
                     markeredgewidth=0.5, markersize=4.2, linestyle="none")
BG_LINE_KW = dict(marker="s", color=GREY, markerfacecolor="white", markeredgecolor=GREY,
                  markeredgewidth=0.8, markersize=3.8, linestyle="none")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 7,
    "axes.labelsize": 7.5,
    "axes.titlesize": 9,
    "xtick.labelsize": 6.0,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.2,
    "ytick.major.size": 2.2,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})


def strip(ax, x, vals, kw, jitter=0.075, seed=11):
    """Individual samples as a narrow jittered strip, with a mean bar and 1 SD."""
    rng = np.random.default_rng(seed)
    v = np.asarray(vals, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return
    ax.plot(x + rng.uniform(-jitter, jitter, v.size), v, zorder=3, **kw)
    m = v.mean()
    if v.size > 1:
        sd = v.std(ddof=1)
        ax.plot([x, x], [m - sd, m + sd], color=kw["color"], lw=0.7, zorder=4)
    ax.plot([x - 0.19, x + 0.19], [m, m], color=kw["color"], lw=1.5,
            solid_capstyle="butt", zorder=5)


def main() -> None:
    df = pd.read_csv(SRC)
    df["day"] = pd.to_numeric(df["day"], errors="coerce")
    OUT.mkdir(parents=True, exist_ok=True)

    fig, (axa, axb) = plt.subplots(
        1, 2, figsize=(WIDTH, HEIGHT), gridspec_kw={"width_ratios": [1.0, 1.3]}
    )

    # ---- panel a: by condition -------------------------------------------
    for i, cond in enumerate(ORDER):
        sub = df[df["condition"] == cond]
        strip(axa, i - 0.18, sub["panel_mean"], PANEL_KW)
        strip(axa, i + 0.18, sub["bg_mean"], BG_KW, seed=23)
        axa.annotate(f"n={len(sub)}", (i, 1.025), ha="center", va="bottom",
                     fontsize=6, color="#666666", annotation_clip=False)
    axa.set_xticks(range(len(ORDER)))
    axa.set_xticklabels([SHORT[c] for c in ORDER], rotation=35, ha="right",
                        rotation_mode="anchor")
    axa.set_ylabel("Mean methylation (β)")
    axa.set_xlim(-0.6, len(ORDER) - 0.4)
    axa.set_title("a", loc="left", fontweight="bold", pad=10)

    # ---- panel b: by day --------------------------------------------------
    day_conds = ["baseline", "failed", "reprogrammed_intermediate", "reprogrammed"]
    pos, ticks, labels, seps, spans = 0.0, [], [], [], []
    for ci, cond in enumerate(day_conds):
        sub = df[df["condition"] == cond]
        days = sorted(d for d in sub["day"].dropna().unique())
        xs, pv, bv = [], [], []
        for d in days:
            s = sub[sub["day"] == d]
            xs.append(pos)
            pv.append(s["panel_mean"].mean())
            bv.append(s["bg_mean"].mean())
            ticks.append(pos)
            labels.append(f"{int(d)}")
            pos += 1.0
        axb.plot(xs, pv, lw=1.3, color=INK, zorder=3)
        axb.plot(xs, bv, lw=1.1, color=GREY, zorder=2)
        axb.plot(xs, pv, zorder=5, **PANEL_LINE_KW)
        axb.plot(xs, bv, zorder=4, **BG_LINE_KW)
        spans.append(((xs[0] + xs[-1]) / 2, TIGHT[cond]))
        if ci < len(day_conds) - 1:
            seps.append(pos + 0.05)
        pos += 1.1
    for s in seps:
        axb.axvline(s, color="#dcdcdc", lw=0.5, zorder=1)
    for mid, lab in spans:
        axb.annotate(lab, (mid, 1.025), ha="center", va="bottom",
                     fontsize=6.5, color="#333333", annotation_clip=False)
    axb.set_xticks(ticks)
    axb.set_xticklabels(labels)
    axb.tick_params(axis="x", pad=1.5)
    axb.set_xlabel("Day")
    axb.set_xlim(-0.8, pos - 1.3)
    axb.set_title("b", loc="left", fontweight="bold", pad=10)

    for ax in (axa, axb):
        ax.set_ylim(0, 1.0)
        ax.set_yticks(np.arange(0, 1.01, 0.2))
        ax.spines[["top", "right"]].set_visible(False)
    axa.tick_params(axis="x", length=0)

    handles = [
        Line2D([], [], **PANEL_LINE_KW, label="16-CpG panel"),
        Line2D([], [], **BG_LINE_KW, label="Background loci (n = 8,916)"),
    ]
    axa.legend(handles=handles, loc="lower left", frameon=False,
               handletextpad=0.4, borderaxespad=0.1, labelspacing=0.35)

    fig.subplots_adjust(wspace=0.20)
    for ext, kw in (("eps", {}), ("pdf", {}), ("png", {"dpi": 600}),
                    ("tif", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})):
        fig.savefig(OUT / f"Fig6.{ext}", **kw)
    w, h = fig.get_size_inches()
    plt.close(fig)

    print(f"canvas {w*25.4:.0f} x {h*25.4:.0f} mm")
    g = df.groupby("condition").agg(n=("panel_mean", "size"),
                                    panel=("panel_mean", "mean"),
                                    panel_sd=("panel_mean", "std"),
                                    bg=("bg_mean", "mean"))
    print(g.reindex(ORDER).round(4).to_string())
    for c in day_conds:
        r = df[df["condition"] == c].groupby("day")["panel_mean"].mean().round(3)
        print(f"{c:28s} {dict(r)}")
    for p in sorted(OUT.glob("Fig6.*")):
        print(f"  {p.name}  {p.stat().st_size/1024:.0f} KB")


if __name__ == "__main__":
    main()
