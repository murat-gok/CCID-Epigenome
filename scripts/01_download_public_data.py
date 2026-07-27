#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

DATASETS = {
    "gse165177_rna": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165177/suppl/GSE165177_Log2_RPM_Transient_reprogramming.txt.gz",
        "filename": "GSE165177_Log2_RPM_Transient_reprogramming.txt.gz",
        "role": "Fibroblast expression, MPTR part 1 (required)",
    },
    "gse165177_rna_part2": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165177/suppl/GSE165177_Log2_RPM_Transient_reprogramming_part2_170621.txt.gz",
        "filename": "GSE165177_Log2_RPM_Transient_reprogramming_part2_170621.txt.gz",
        "role": "Fibroblast expression, MPTR part 2 (required for the complete trajectory)",
    },
    "gse165178_sendai_methylation": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165178/suppl/GSE165178_Matrix_processed_sendai.txt.gz",
        "filename": "GSE165178_Matrix_processed_sendai.txt.gz",
        "role": "Full Sendai reprogramming methylation trajectory",
    },
    "gse165179_mptr_methylation": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165179/suppl/GSE165179_Matrix_processed_transient.txt.gz",
        "filename": "GSE165179_Matrix_processed_transient.txt.gz",
        "role": "Primary MPTR methylation target signature (large, ~530 MB compressed)",
    },
    "gse269758_tcell_bystander": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269758/suppl/GSE269758_T_cell_hyper_beta_detp.csv.gz",
        "filename": "GSE269758_T_cell_hyper_beta_detp.csv.gz",
        "role": "External T-cell bystander validation, replicates 1-2",
    },
    "gse269758_tcell_bystander_r3": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269758/suppl/GSE269758_T_cell_R3_Hyperlibrary_beta_detp.csv.gz",
        "filename": "GSE269758_T_cell_R3_Hyperlibrary_beta_detp.csv.gz",
        "role": "External T-cell bystander validation, replicate 3 (required)",
    },
    "gse285251_msc_bystander": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285251/suppl/GSE285251_MSC_Hyperlibrary_beta_detp.csv.gz",
        "filename": "GSE285251_MSC_Hyperlibrary_beta_detp.csv.gz",
        "role": "External MSC epigenome-editing/bystander validation",
    },
}


def download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0:
        print(f"SKIP {target.name}: already exists")
        return
    partial = target.with_suffix(target.suffix + ".part")
    print(f"DOWNLOAD {url}\n      -> {target}")
    try:
        urllib.request.urlretrieve(url, partial)
        partial.replace(target)
    except Exception:
        partial.unlink(missing_ok=True)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description="Download public datasets for the CCID project")
    parser.add_argument(
        "--datasets",
        default="gse165177_rna,gse165177_rna_part2,gse269758_tcell_bystander,gse269758_tcell_bystander_r3,gse285251_msc_bystander",
        help="Comma-separated dataset keys, or 'all'",
    )
    parser.add_argument("--output", default="data/raw")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()

    if args.list:
        for key, item in DATASETS.items():
            print(f"{key}: {item['role']}")
        return

    keys = list(DATASETS) if args.datasets == "all" else args.datasets.split(",")
    unknown = sorted(set(keys).difference(DATASETS))
    if unknown:
        sys.exit(f"Unknown dataset keys: {unknown}")
    output = Path(args.output)
    for key in keys:
        item = DATASETS[key]
        download(item["url"], output / item["filename"])


if __name__ == "__main__":
    main()
