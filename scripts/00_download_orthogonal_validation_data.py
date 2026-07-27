#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

DATASETS = {
    "gse142439_external_reprogramming": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE142nnn/GSE142439/matrix/GSE142439_series_matrix.txt.gz",
        "filename": "GSE142439_series_matrix.txt.gz",
        "role": "Independent paired OSKMLN transient-reprogramming EPIC matrix (4 skin-fibroblast donors)",
        "set": "core",
    },
    "gse131280_replicative_aging": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131280/matrix/GSE131280_series_matrix.txt.gz",
        "filename": "GSE131280_series_matrix.txt.gz",
        "role": "EPIC methylation across replicative lifespan in primary human fibroblasts",
        "set": "core",
    },
    "gse22595_chronological_aging": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22595/matrix/GSE22595_series_matrix.txt.gz",
        "filename": "GSE22595_series_matrix.txt.gz",
        "role": "Cross-sectional young-versus-old dermal fibroblast HM27 methylation",
        "set": "core",
    },
    "gse52025_fibroblast_methylation": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52025/matrix/GSE52025_series_matrix.txt.gz",
        "filename": "GSE52025_series_matrix.txt.gz",
        "role": "Fibroblast HM450 methylation in 62 unrelated individuals for eQTM analysis",
        "set": "core",
    },
    "gse53243_fibroblast_expression": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/GSE53243/matrix/GSE53243_series_matrix.txt.gz",
        "filename": "GSE53243_series_matrix.txt.gz",
        "role": "Matched fibroblast expression matrix for eQTM analysis",
        "set": "core",
    },
    "gpl6883_expression_annotation": {
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6883&targ=self&view=data&form=text",
        "filename": "GPL6883_data.txt",
        "role": "Illumina HumanRef-8 v3 probe-to-gene platform table (NCBI GEO Accession Display)",
        "set": "core",
    },
    "gse179847_multidonor_metadata": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179847/matrix/GSE179847_series_matrix.txt.gz",
        "filename": "GSE179847_series_matrix.txt.gz",
        "role": "GSE179847 sample metadata Series Matrix (small; no deposited beta matrix)",
        "set": "extended",
    },
    "gse179847_multidonor_lifespan": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179847/suppl/GSE179847_Cell_lifespan_DNAm_processed_matrix.csv.gz",
        "filename": "GSE179847_Cell_lifespan_DNAm_processed_matrix.csv.gz",
        "role": "Optional 2-GB normalized beta matrix for multidonor longitudinal fibroblast aging",
        "set": "extended",
    },
}


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, target: Path, retries: int = 5) -> dict[str, object]:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 0:
        print(f"SKIP {target.name}: already exists ({target.stat().st_size:,} bytes)")
        return {
            "url": url,
            "path": str(target),
            "size_bytes": target.stat().st_size,
            "sha256": sha256_file(target),
            "status": "existing",
        }
    partial = target.with_suffix(target.suffix + ".part")
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        existing = partial.stat().st_size if partial.exists() else 0
        headers = {"User-Agent": "CCID-Epigenome/0.9.0.post4"}
        if existing:
            headers["Range"] = f"bytes={existing}-"
        request = urllib.request.Request(url, headers=headers)
        mode = "ab" if existing else "wb"
        print(
            f"DOWNLOAD {url}\n      -> {target} (attempt {attempt}/{retries}, "
            f"resume={existing:,} bytes)"
        )
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                status = getattr(response, "status", None)
                # If the server ignored a Range request, restart rather than append a full file.
                if existing and status != 206:
                    existing = 0
                    mode = "wb"
                with partial.open(mode) as handle:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        handle.write(chunk)
            partial.replace(target)
            return {
                "url": url,
                "path": str(target),
                "size_bytes": target.stat().st_size,
                "sha256": sha256_file(target),
                "status": "downloaded",
            }
        except Exception as exc:
            last_error = exc
            # Keep the .part file so the next attempt or a later invocation can resume.
            if attempt < retries:
                time.sleep(2**attempt)
    raise RuntimeError(
        f"Download failed after {retries} attempts: {url}. Partial file retained at {partial}"
    ) from last_error


def main() -> None:
    parser = argparse.ArgumentParser(description="Download orthogonal-validation datasets for CCID-Epigenome")
    parser.add_argument("--datasets", default="core", help="core, extended, all, or comma-separated keys")
    parser.add_argument("--output", default="data/raw/orthogonal_validation")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()

    if args.list:
        for key, item in DATASETS.items():
            print(f"{key} [{item['set']}]: {item['role']}")
        return

    if args.datasets in {"core", "extended"}:
        keys = [k for k, v in DATASETS.items() if v["set"] == args.datasets]
    elif args.datasets == "all":
        keys = list(DATASETS)
    else:
        keys = [x.strip() for x in args.datasets.split(",") if x.strip()]
    unknown = sorted(set(keys).difference(DATASETS))
    if unknown:
        sys.exit(f"Unknown dataset keys: {unknown}")
    output = Path(args.output)
    manifest: list[dict[str, object]] = []
    for key in keys:
        item = DATASETS[key]
        record = download(item["url"], output / item["filename"])
        record.update({"dataset_key": key, "role": item["role"], "set": item["set"]})
        manifest.append(record)
    manifest_path = output / "orthogonal_validation_download_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("Manifest:", manifest_path)


if __name__ == "__main__":
    main()
