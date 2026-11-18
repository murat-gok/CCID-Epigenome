#!/usr/bin/env python3
"""Figure 7 for the Biogerontology submission.

Panel a: mean trough depth of every evaluable locus by distance to the nearest
         KLF4 or OCT4 binding region.
Panel b: the same gradient computed separately within narrow bands of control
         methylation, which shows that it is not an artefact of the room a
         locus has to lose methylation.

Design constraints follow the other figures of this submission:
  - a single y-axis, no dual scales
  - identity carried by marker shape AND fill, never colour alone, so the
    figure survives greyscale printing and colour-vision deficiency
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

SRC_A = Path(os.environ.get("FIG7_SRC_A", "results/faz0/fig7_bands_overall.csv"))
SRC_B = Path(os.environ.get("FIG7_SRC_B", "results/faz0/fig7_bands_by_control_beta.csv"))
OUT = Path(os.environ.get("FIG7_OUT", "manuscript/figures/submission"))

MM = 1 / 25.4
WIDTH = 174 * MM
HEIGHT = 78 * MM

INK = "#1a1a1a"
GREY = "#8c8c8c"
BAR = "#b4b4b4"

# three strata, separated by marker shape and by fill, so the panel reads in
# greyscale and under colour-vision deficiency alike
STRATA = [
    ("0.6-0.7", dict(marker="o", markerfacecolor=INK, markeredgecolor="white",
                     markeredgewidth=0.5, markersize=4.0, color=INK, lw=1.2)),
    ("0.7-0.8", dict(marker="s", markerfacecolor="white", markeredgecolor=INK,
                     markeredgewidth=0.9, markersize=3.6, color=INK, lw=1.0,
                     dashes=(4, 1.6))),
    ("0.8-0.9", dict(marker="^", markerfacecolor=GREY, markeredgecolor="white",
                     markeredgewidth=0.5, markersize=4.2, color=GREY, lw=1.0,
                     dashes=(1.6, 1.4))),
]

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


def thousands(n: int) -> str:
    return f"{n:,}".replace(",", " ")


def main() -> None:
    a = pd.read_csv(SRC_A)
    b = pd.read_csv(SRC_B)
    order = a["band"].tolist()
    OUT.mkdir(parents=True, exist_ok=True)

    fig, (axa, axb) = plt.subplots(
        1, 2, figsize=(WIDTH, HEIGHT), gridspec_kw={"width_ratios": [1.0, 1.0]}
    )

    # ---- panel a: all evaluable loci --------------------------------------
    x = np.arange(len(a))
    axa.bar(x, a["mean"], width=0.66, facecolor=BAR, edgecolor=INK, linewidth=0.6,
            zorder=3)
    axa.errorbar(x, a["mean"], yerr=a["sem"], fmt="none", ecolor=INK, elinewidth=0.7,
                 capsize=1.6, capthick=0.7, zorder=4)
    axa.axhline(0, color=INK, lw=0.6, zorder=2)
    top = float(a["mean"].max())
    for xi, (m, n) in enumerate(zip(a["mean"], a["n"])):
        off = 0.0022 if m >= 0 else -0.0022
        va = "bottom" if m >= 0 else "top"
        axa.annotate(thousands(int(n)), (xi, m + off), ha="center", va=va,
                     fontsize=5.6, color="#555555", rotation=90,
                     annotation_clip=False)
    axa.set_ylabel("Mean trough depth (Δβ)")
    axa.set_ylim(-0.045, top * 1.55)
    axa.set_xlim(-0.65, len(a) - 0.35)
    axa.set_title("a", loc="left", fontweight="bold", pad=10)

    # ---- panel b: within bands of control methylation ---------------------
    for lab, kw in STRATA:
        s = b[b["beta_band"] == lab].set_index("band").reindex(order).dropna(
            subset=["mean"])
        xs = [order.index(i) for i in s.index]
        axb.plot(xs, s["mean"], zorder=3, **kw)
    axb.axhline(0, color=INK, lw=0.6, zorder=2)
    axb.set_ylabel("Mean trough depth (Δβ)")
    axb.set_ylim(-0.03, 0.205)
    axb.set_xlim(-0.65, len(a) - 0.35)
    axb.set_title("b", loc="left", fontweight="bold", pad=10)

    handles = [Line2D([], [], label=f"Control β {lab}", **kw) for lab, kw in STRATA]
    axb.legend(handles=handles, loc="upper right", frameon=False,
               handletextpad=0.5, handlelength=2.4, borderaxespad=0.2,
               labelspacing=0.35)

    for ax in (axa, axb):
        ax.set_xticks(x)
        ax.set_xticklabels(order, rotation=35, ha="right", rotation_mode="anchor")
        ax.set_xlabel("Distance to nearest KLF4 or OCT4 binding region")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="x", length=0, pad=1.5)

    fig.subplots_adjust(wspace=0.26)
    for ext, kw in (("eps", {}), ("pdf", {}), ("png", {"dpi": 600}),
                    ("tif", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})):
        fig.savefig(OUT / f"Fig7.{ext}", **kw)
    w, h = fig.get_size_inches()
    plt.close(fig)

    print(f"canvas {w*25.4:.0f} x {h*25.4:.0f} mm")
    print(a.round(4).to_string(index=False))
    for p in sorted(OUT.glob("Fig7.*")):
        print(f"  {p.name}  {p.stat().st_size/1024:.0f} KB")


if __name__ == "__main__":
    main()
