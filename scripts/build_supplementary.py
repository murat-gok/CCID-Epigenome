#!/usr/bin/env python3
"""Build the two supplementary workbooks the manuscript cites as Online
Resource 1 and Online Resource 2.

Online Resource 1  machine-readable outputs of the workflow stages
Online Resource 2  candidate-level probe quality and genomic-domain annotation

Each workbook opens with a metadata sheet carrying the article title, the
journal, the author and the corresponding address, as the journal asks, and a
caption describing the item. Run from the repository root.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd

ROOT = Path(".").resolve()
OUT = ROOT / "manuscript" / "supplementary"
QC = ROOT / "data" / "reference" / "probe_qc"

TITLE = ("CCID-Epigenome: an evidence-convergence workflow for nominating DNA methylation "
         "editing targets, and the transient nature of partial-reprogramming candidates")
JOURNAL = "Functional & Integrative Genomics"
AUTHOR = "Murat Gok"
AFFIL = ("Department of Computer Engineering, Yalova University, Yalova, Turkiye; "
         "murat.gok@yalova.edu.tr; ORCID 0000-0003-2261-9288")
REPO = ("https://github.com/murat-gok/CCID-Epigenome; archived at Zenodo, "
        "doi:10.5281/zenodo.23765257 (v0.9.1)")


def meta_sheet(caption, contents):
    rows = [("Article title", TITLE), ("Journal", JOURNAL), ("Author", AUTHOR),
            ("Corresponding author", AFFIL), ("Software", REPO), ("", ""),
            ("Caption", caption), ("", "")]
    rows += [("Sheet", "Contents")] + list(contents)
    return pd.DataFrame(rows, columns=["Field", "Value"])


def write(path, sheets):
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        for name, df in sheets:
            df.to_excel(xl, sheet_name=name[:31], index=False)
            ws = xl.sheets[name[:31]]
            for col in ws.columns:
                width = max((len(str(c.value)) for c in col if c.value is not None),
                            default=10)
                ws.column_dimensions[col[0].column_letter].width = min(max(width + 2, 10), 70)
    print(f"  wrote {path.name} ({path.stat().st_size/1024:.0f} KB, {len(sheets)} sheets)")


# --------------------------------------------------------------- resource 1
def resource_1():
    faz = ROOT / "results" / "faz0"
    summary = json.loads((faz / "faz0_summary.json").read_text())
    flat = []

    def walk(prefix, obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                walk(f"{prefix}.{k}" if prefix else k, v)
        elif isinstance(obj, list):
            flat.append((prefix, "; ".join(map(str, obj))))
        else:
            flat.append((prefix, obj))

    walk("", summary)
    sheets = [
        ("Diagnostic summary", pd.DataFrame(flat, columns=["Key", "Value"])),
        ("Rank statistic anatomy", pd.read_csv(faz / "testA_panel_floor_anatomy.csv")),
        ("Variance floor sensitivity", pd.read_csv(faz / "testB_ranking_sensitivity.csv")),
        ("Panel rank movement", pd.read_csv(faz / "testB_panel_rank_movement.csv")),
        ("Global effect distribution", pd.read_csv(faz / "testC_global_delta_distribution.csv")),
        ("Panel trajectory", pd.read_csv(faz / "panel_trajectory_per_sample.csv")),
        ("Orthogonal evidence overview",
         pd.read_csv(ROOT / "data/derived/orthogonal_results_overview.csv")),
        ("Natural aging audit",
         pd.read_csv(ROOT / "data/derived/natural_aging_candidate_audit_v090.csv")),
    ]
    caption = ("Online Resource 1. Machine-readable outputs of the CCID-Epigenome workflow "
               "stages: the rank-statistic anatomy and variance-floor sensitivity analyses, "
               "the genome-wide effect distribution, per-sample panel and background "
               "methylation across every condition of GSE165179, and the orthogonal "
               "validation summaries.")
    contents = [(n, f"{len(df):,} rows x {len(df.columns)} columns") for n, df in sheets]
    write(OUT / "ESM_1.xlsx", [("Metadata", meta_sheet(caption, contents))] + sheets)


# --------------------------------------------------------------- resource 2
def resource_2():
    faz = pd.read_csv(ROOT / "results/faz0/testA_panel_floor_anatomy.csv")
    panel = faz["cpg_id"].tolist()
    moves = pd.read_csv(ROOT / "results/faz0/testB_panel_rank_movement.csv")

    # Zhou annotation: coordinates, mapping quality, masks, domains
    zh = pd.read_csv(QC / "EPIC.hg38.zhou_annotation.tsv.gz", sep="\t",
                     compression="gzip", low_memory=False)
    zh = zh[zh["probe_id"].isin(panel)]

    # Pidsley: every sheet of the EPIC off-target workbook, as a flat probe list
    pid = set()
    with open(QC / "13059_2016_1066_MOESM1_ESM.csv", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            tok = line.split(",")[0].strip().strip('"')
            if tok.startswith("cg") or tok.startswith("ch."):
                pid.add(tok)

    # McCartney: polymorphic / off-target lists
    mcc = set()
    for name in ("1-s2.0-S221359601630071X-mmc2.txt", "1-s2.0-S221359601630071X-mmc3.txt"):
        with open(QC / name, encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                tok = line.split("\t")[0].split(",")[0].strip().strip('"')
                if tok.startswith("cg") or tok.startswith("ch."):
                    mcc.add(tok)

    # Control methylation at the 16 loci. This must be the negative controls matched
    # to the day-10 discovery contrast, which is the quantity the manuscript reports;
    # averaging every baseline sample of the series gives a different number.
    pairs = pd.read_csv(ROOT / "results/real/safe_osk_target_signature_O1O2_matched_pairs.csv")
    ctrl = sorted({s for ids in pairs["baseline_sample_ids"]
                   for s in str(ids).split(";") if s.strip()})
    keep = set(panel)
    rows = {}
    with gzip.open(ROOT / "data/processed/fibroblast_beta.csv.gz", "rt") as fh:
        header = fh.readline().rstrip("\n").split(",")
        idx = [i for i, c in enumerate(header) if c in set(ctrl)]
        for line in fh:
            parts = line.rstrip("\n").split(",")
            if parts[0] in keep:
                vals = [float(parts[i]) for i in idx
                        if parts[i] not in ("", "NA", "nan")]
                rows[parts[0]] = sum(vals) / len(vals) if vals else float("nan")
            if len(rows) == len(keep):
                break
    ctrl_mean = pd.Series(rows, name="control_mean_beta")

    audit = pd.read_csv(ROOT / "data/derived/natural_aging_candidate_audit_v090.csv")

    t = faz.merge(zh, left_on="cpg_id", right_on="probe_id", how="left")
    t = t.merge(moves, on="cpg_id", how="left")
    t = t.merge(ctrl_mean, left_on="cpg_id", right_index=True, how="left")
    t = t.merge(audit[["cpg_id", "candidate_identity", "robust_rank", "candidate_rank",
                       "natural_aging_evidence_class"]], on="cpg_id", how="left")
    t["flagged_pidsley_2016"] = t["cpg_id"].isin(pid)
    t["flagged_mccartney_2016"] = t["cpg_id"].isin(mcc)
    t["masked_zhou_2017_M_general"] = t["M_general"].astype("boolean")
    t["flagged_by_any_resource"] = (t["flagged_pidsley_2016"] | t["flagged_mccartney_2016"]
                                    | t["masked_zhou_2017_M_general"].fillna(False))

    cols = {
        "cpg_id": "CpG identifier",
        "candidate_identity": "Locus identity",
        "robust_rank": "Computational robustness rank",
        "candidate_rank": "Biological priority rank",
        "chrm": "Chromosome (hg38)",
        "beg": "Position (hg38)",
        "mapQ": "Mapping quality",
        "flagged_pidsley_2016": "Flagged by Pidsley et al. 2016",
        "flagged_mccartney_2016": "Flagged by McCartney et al. 2016",
        "masked_zhou_2017_M_general": "Masked by Zhou et al. 2017 (M_general)",
        "flagged_by_any_resource": "Flagged by any resource",
        "commonPMD": "In a common partially methylated domain",
        "commonHMD": "In a common highly methylated domain",
        "control_mean_beta": "Mean methylation in matched controls",
        "delta_beta": "Discovery mean delta beta",
        "delta_O1": "Donor O1 delta beta",
        "delta_O2": "Donor O2 delta beta",
        "sign_consistency": "Donor sign consistency",
        "donor_spread": "Donor spread",
        "standard_error": "Donor-level standard error",
        "se_at_floor": "Standard error at the 0.01 floor",
        "rank_score_published": "Published rank score",
        "published (floor 0.01)": "Rank, published floor 0.01",
        "floor 0.02": "Rank, floor 0.02",
        "floor 0.05": "Rank, floor 0.05",
        "no floor (true t)": "Rank, no floor",
        "natural_aging_evidence_class": "Natural-aging evidence class",
    }
    out = t[[c for c in cols if c in t.columns]].rename(columns=cols)
    out = out.sort_values("Computational robustness rank")

    notes = pd.DataFrame([
        ("Pidsley et al. 2016", "Genome Biol 17:208; EPIC off-target and polymorphic probe lists"),
        ("McCartney et al. 2016", "Genomics Data 9:22-24; polymorphic and off-target probe lists"),
        ("Zhou et al. 2017", "Nucleic Acids Res 45:e22; M_general mapping, repeat and SNP mask"),
        ("Zhou et al. 2018", "Nat Genet 50:591-602; common partially and highly methylated domains"),
        ("Platform-wide rates",
         "union of exclusion resources 5.13% of 864,719 evaluable loci; Zhou mask alone 5.08%; "
         "common PMD 24.3%; common HMD 39.7%"),
        ("Variance floor",
         "the discovery rank divides the mean donor effect by the donor-level standard error, "
         "floored at 0.01; ranks under alternative floors are given for comparison"),
        ("Control methylation",
         "mean beta across the negative-control samples matched to the day-10 discovery "
         "contrast for donors O1 and O2"),
    ], columns=["Item", "Description"])

    caption = ("Online Resource 2. Candidate-level probe quality and genomic-domain "
               "annotation for the 16 shortlisted CpGs: coordinates and mapping quality, "
               "flagging by three published exclusion resources, partially and highly "
               "methylated domain membership, control methylation, discovery effect "
               "components, and rank under alternative variance floors.")
    sheets = [("Candidate annotation", out), ("Sources and notes", notes)]
    contents = [(n, f"{len(df):,} rows x {len(df.columns)} columns") for n, df in sheets]
    write(OUT / "ESM_2.xlsx", [("Metadata", meta_sheet(caption, contents))] + sheets)

    # the checks the manuscript's claims rest on
    print("  verification:")
    print(f"    candidates: {len(out)} (expected 16)")
    print(f"    flagged by any resource: {int(out['Flagged by any resource'].sum())} (expected 0)")
    print(f"    in a common PMD: {int(out['In a common partially methylated domain'].sum())} "
          "(expected 0)")
    print(f"    in a common HMD: {int(out['In a common highly methylated domain'].sum())} "
          "(expected 13)")
    print(f"    mapping quality 60: {int((out['Mapping quality'] == 60).sum())} (expected 15)")
    print(f"    mean control methylation: {out['Mean methylation in matched controls'].mean():.3f} "
          "(expected 0.766)")
    print(f"    standard error at floor: {int(out['Standard error at the 0.01 floor'].sum())} "
          "(expected 14)")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    print("Online Resource 1")
    resource_1()
    print("Online Resource 2")
    resource_2()
