#!/usr/bin/env python3
"""Figure 1 for the Functional & Integrative Genomics submission.

The workflow schematic, redrawn as true vector line art. The previous file was a
400 dpi bitmap that came out 233 mm tall at the journal's 174 mm width, one
millimetre inside the height limit, and it carried its own title inside the
image, which the journal does not allow.

Content follows the current manuscript, with one addition: step 8 records the
post hoc challenges (probe quality, genomic domain, state specificity) that the
Results and Discussion now report and that the earlier schematic predates.

Specification: 174 mm wide, Liberation Sans (metrically identical to Arial) at
8-9 pt, all rules 0.7 pt, vector EPS and PDF with Type-42 embedded fonts, and a
1200 dpi TIFF for the bitmap case. No title inside the figure.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUT = Path(os.environ.get("FIG1_OUT", "manuscript/figures/submission"))

MM = 1 / 25.4
W, H = 174, 212                      # mm; the limit is 174 x 234
INK = "#1a1a1a"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
    "text.color": INK,
    "ps.fonttype": 42,
    "pdf.fonttype": 42,
})

# (x, y, w, h) in mm, head, body
FULL = 154.0
LEFT = 8.0
STEPS = [
    ("1. Identity-preserving MPTR state",
     "GSE165177; O1/O2 discovery donors; day 10 selected"),
    ("2. Exploratory CpG ranking",
     "Exact donor/day/batch/stage matching; no locus met the strict false-discovery\n"
     "rule; top 500 frozen"),
    ("3. Sequential three-donor evaluation",
     "O3 frozen-list support; failed-reprogramming controls; time series;\n"
     "leave-one-donor-out rotations"),
    ("4. Coordinate-dependent integration",
     "500 CpGs to 472 harmonized hg38 loci; 28 excluded from coordinate-dependent\n"
     "analyses"),
    (None, None),                     # the 5a / 5b row
    ("6. Integrated robustness assessment",
     "Nine evidence components; score sensitivity, ablation, redundancy and\n"
     "permutation controls; 16 stable CpGs"),
    ("7. Orthogonal validation and biological triage",
     "Independent reprogramming, cellular-aging trajectories, oxygen sensitivity,\n"
     "eQTM, rE2G and guide feasibility; 8 experimental candidates and 2 sentinels"),
    ("8. Post hoc challenges to the frozen panel",
     "Published probe-exclusion resources; partially and highly methylated domain\n"
     "membership; panel methylation across every condition and day of the source series"),
]
SPLIT = ("5a. Editor-network evidence",
         "T-cell and MSC\ndCas9-DNMT3A/3L response",
         "5b. Regulatory context",
         "GENCODE v50, TSS, regions,\nENCODE cCREs and fibroblast rE2G")

BOX_H, SPLIT_H, GAP = 19.0, 21.0, 5.5
NOTE = ("O3 was independent only at the first frozen-list evaluation stage and was not held "
        "out from the final integrated ranking")


def box(ax, x, y, w, h, head, body):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0,rounding_size=2.2",
        linewidth=0.7, edgecolor=INK, facecolor="white", zorder=2))
    ax.text(x + 4, y + h - 5.0, head, fontsize=9, fontweight="bold",
            va="center", ha="left", zorder=3)
    ax.text(x + 4, y + h - 9.6, body, fontsize=8, va="top", ha="left",
            linespacing=1.45, zorder=3)


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7,
        linewidth=0.7, color=INK, shrinkA=0, shrinkB=0, zorder=1))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(W * MM, H * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")

    y = H - 4.0
    centre = LEFT + FULL / 2
    prev_bottom = None
    for head, body in STEPS:
        if head is None:                                   # the 5a / 5b row
            gw = (FULL - 6.0) / 2
            y -= SPLIT_H
            box(ax, LEFT, y, gw, SPLIT_H, SPLIT[0], SPLIT[1])
            box(ax, LEFT + gw + 6.0, y, gw, SPLIT_H, SPLIT[2], SPLIT[3])
            arrow(ax, centre, prev_bottom, centre, y + SPLIT_H)
            left_c, right_c = LEFT + gw / 2, LEFT + gw + 6.0 + gw / 2
            nxt = y - GAP
            arrow(ax, left_c, y, centre - 14, nxt)
            arrow(ax, right_c, y, centre + 14, nxt)
            prev_bottom = None
            y -= GAP
            continue
        y -= BOX_H
        box(ax, LEFT, y, FULL, BOX_H, head, body)
        if prev_bottom is not None:
            arrow(ax, centre, prev_bottom, centre, y + BOX_H)
        prev_bottom = y
        y -= GAP

    ax.text(LEFT, y + GAP - 3.0, NOTE, fontsize=8, style="italic",
            va="top", ha="left", linespacing=1.45)

    for ext, kw in (("eps", {}), ("pdf", {}), ("png", {"dpi": 600}),
                    ("tif", {"dpi": 1200,
                             "pil_kwargs": {"compression": "tiff_lzw"}})):
        fig.savefig(OUT / f"Fig1.{ext}", **kw)
    plt.close(fig)
    print(f"canvas {W} x {H} mm (limit 174 x 234); bottom of content y={y:.1f} mm")
    for p in sorted(OUT.glob("Fig1.*")):
        print(f"  {p.name}  {p.stat().st_size/1024:.0f} KB")


if __name__ == "__main__":
    main()
