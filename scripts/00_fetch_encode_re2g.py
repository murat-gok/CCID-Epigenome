#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import urllib.parse
import urllib.request
from pathlib import Path

DEFAULTS = {
    "ENCSR360BQF": "fibroblast_skin_scalp",
    "ENCSR886NXW": "fibroblast_skin_back",
}


def fetch_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "CCID-Epigenome/0.9.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.load(response)


def download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "CCID-Epigenome/0.9.0"})
    with urllib.request.urlopen(request, timeout=300) as response, partial.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    partial.replace(target)


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch ENCODE-rE2G fibroblast association-graph files.")
    parser.add_argument("--accessions", default=",".join(DEFAULTS))
    parser.add_argument("--output", default="data/raw/orthogonal_validation/encode_re2g")
    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    manifest = []
    for accession in [x.strip() for x in args.accessions.split(",") if x.strip()]:
        query = urllib.parse.urlencode(
            {
                "type": "File",
                "dataset": f"/annotations/{accession}/",
                "format": "json",
                "limit": "all",
            }
        )
        payload = fetch_json(f"https://www.encodeproject.org/search/?{query}")
        files = payload.get("@graph", [])
        released = [f for f in files if f.get("status") == "released"]
        candidates = sorted(
            released,
            key=lambda f: (
                str(f.get("output_type", "")).lower() != "association graph",
                str(f.get("file_format", "")).lower() not in {"bed", "tsv"},
                -int(f.get("file_size") or 0),
            ),
        )
        selected = candidates[0] if candidates else None
        if selected is None:
            manifest.append({"accession": accession, "status": "no_released_file_found"})
            print(f"WARN {accession}: no released file was returned by the ENCODE API")
            continue
        href = selected.get("href")
        if not href:
            manifest.append({"accession": accession, "status": "file_missing_href"})
            continue
        source_url = urllib.parse.urljoin("https://www.encodeproject.org", href)
        extension = Path(urllib.parse.urlparse(source_url).path).suffixes
        suffix = "".join(extension[-2:]) if extension[-1:] == [".gz"] else "".join(extension[-1:])
        target = output / f"{accession}_{DEFAULTS.get(accession, accession)}{suffix}"
        print(f"DOWNLOAD {source_url}\n      -> {target}")
        download(source_url, target)
        converted = None
        if target.suffix.lower() in {".bigbed", ".bb"}:
            converter = shutil.which("bigBedToBed")
            if converter:
                converted = target.with_suffix(".bed")
                subprocess.run([converter, str(target), str(converted)], check=True)
        manifest.append(
            {
                "accession": accession,
                "status": "downloaded",
                "file_accession": selected.get("accession"),
                "file_format": selected.get("file_format"),
                "output_type": selected.get("output_type"),
                "source_url": source_url,
                "path": str(converted or target),
            }
        )
    (output / "encode_re2g_download_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print("Manifest:", output / "encode_re2g_download_manifest.json")


if __name__ == "__main__":
    main()
