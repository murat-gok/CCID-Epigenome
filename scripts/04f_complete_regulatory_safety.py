from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from ccid.regulatory_safety import (
    DEFAULT_V081_CONFIG,
    load_cancer_gene_census,
    load_enhancer_gene_links,
    load_literature_evidence,
    load_regulatory_intervals,
    load_tss_reference,
    score_regulatory_safety_completion,
    write_regulatory_safety_outputs,
)


def _optional_frame(path: str, loader):
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        print(f"Optional input not found; marking evidence as not supplied: {p}")
        return None
    return loader(p)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Complete v0.8.0 candidate scoring with traceable regulatory, literature and safety evidence."
    )
    parser.add_argument("--v080-dir", default="results/real/genomic_causality_v080")
    parser.add_argument("--gencode-gtf", default="data/reference/v081/gencode.v50.annotation.gtf.gz")
    parser.add_argument("--ccre", default="")
    parser.add_argument("--ccre-biosamples", default="")
    parser.add_argument("--enhancer-gene-links", default="")
    parser.add_argument(
        "--literature-evidence",
        default="data/reference/literature_evidence_seed_v081.csv",
    )
    parser.add_argument("--cancer-gene-census", default="")
    parser.add_argument("--reference-manifest", default="data/reference/v081/reference_manifest_v081.json")
    parser.add_argument("--config-json", default="configs/regulatory_safety_v081.json")
    parser.add_argument("--shortlist-size", type=int, default=100)
    parser.add_argument("--output-dir", default="results/real/regulatory_safety_v081")
    args = parser.parse_args()

    v080 = Path(args.v080_dir)
    cpg_path = v080 / "candidate_cpg_annotation_causality_scores.csv.gz"
    region_path = v080 / "candidate_region_annotation_causality_scores.csv"
    if not cpg_path.exists() or not region_path.exists():
        raise FileNotFoundError("Run v0.8.0 first; required CpG/region score files were not found")
    cpg = pd.read_csv(cpg_path, low_memory=False)
    regions = pd.read_csv(region_path, low_memory=False)

    tss = _optional_frame(args.gencode_gtf, load_tss_reference)

    ccre_paths = [x for x in args.ccre.split(",") if x]
    if not ccre_paths:
        reference_dir = Path(args.gencode_gtf).parent
        ccre_paths = [str(x) for x in sorted(reference_dir.glob("encode_ccre_*.bed.gz"))]
    biosamples = [x for x in args.ccre_biosamples.split(",") if x]
    regulatory = []
    for i, path in enumerate(ccre_paths):
        p = Path(path)
        if not p.exists():
            print(f"Optional cCRE input not found; skipping: {p}")
            continue
        biosample = biosamples[i] if i < len(biosamples) else p.stem
        regulatory.append(load_regulatory_intervals(p, biosample=biosample))

    enhancer = _optional_frame(args.enhancer_gene_links, load_enhancer_gene_links)
    literature_paths = [x for x in args.literature_evidence.split(",") if x and Path(x).exists()]
    literature = load_literature_evidence(literature_paths) if literature_paths else pd.DataFrame()
    if args.cancer_gene_census and Path(args.cancer_gene_census).exists():
        cgc = load_cancer_gene_census(args.cancer_gene_census)
        literature = pd.concat([literature, cgc], ignore_index=True) if not literature.empty else cgc

    config = dict(DEFAULT_V081_CONFIG)
    config_path = Path(args.config_json)
    if config_path.exists():
        config.update(json.loads(config_path.read_text(encoding="utf-8")))
    config["shortlist_size"] = args.shortlist_size

    result = score_regulatory_safety_completion(
        cpg,
        regions,
        tss_reference=tss,
        regulatory_intervals=regulatory,
        enhancer_gene_links=enhancer,
        literature_evidence=literature,
        config=config,
    )
    reference_manifest = Path(args.reference_manifest)
    if reference_manifest.exists():
        result.manifest["reference_download_manifest"] = json.loads(reference_manifest.read_text(encoding="utf-8"))
    write_regulatory_safety_outputs(result, args.output_dir)
    print("Candidate CpGs:", result.summary["n_candidate_cpgs"])
    print("TSS annotated:", result.summary["n_tss_annotated"])
    print("Regulatory overlaps:", result.summary["n_regulatory_overlap"])
    print("Literature reviewed:", result.summary["n_literature_reviewed"])
    print("Safety reviewed:", result.summary["n_safety_reviewed"])
    print("Manual review queue:", result.summary["n_manual_review_queue"])
    print("P1 candidates:", result.summary["n_priority_P1"])
    print("Results:", args.output_dir)
    print(json.dumps(result.summary["top_cpgs"][:10], indent=2))


if __name__ == "__main__":
    main()
