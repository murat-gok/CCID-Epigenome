from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path


GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz"
ENCODE_ANNOTATIONS = {
    "agnostic": "ENCSR800VNX",
    "activated_t_cell": "ENCSR441PXE",
    "wi38_fibroblast": "ENCSR084ZCE",
}


def _download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "ccid-epigenome/0.8.1"})
    with urllib.request.urlopen(request, timeout=120) as response, target.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def _get_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "ccid-epigenome/0.8.1"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def _encode_preferred_bed(annotation: str) -> tuple[str, str]:
    payload = _get_json(f"https://www.encodeproject.org/annotations/{annotation}/?format=json")
    files = payload.get("files", [])
    candidates = []
    for item in files:
        if isinstance(item, str):
            item = _get_json(f"https://www.encodeproject.org{item}?format=json")
        accession = item.get("accession") or str(item.get("@id", "")).strip("/").split("/")[-1]
        if not accession:
            continue
        if item.get("status") not in {None, "released"}:
            continue
        if item.get("file_format") != "bed":
            continue
        if item.get("output_type") != "candidate Cis-Regulatory Elements":
            continue
        preferred = bool(item.get("preferred_default"))
        candidates.append((preferred, accession))
    if not candidates:
        search = _get_json(
            "https://www.encodeproject.org/search/?type=File&format=json&limit=all"
            f"&dataset=/annotations/{annotation}/&file_format=bed&output_type=candidate%20Cis-Regulatory%20Elements"
        )
        for item in search.get("@graph", []):
            accession = item.get("accession")
            if accession and item.get("status") == "released":
                candidates.append((bool(item.get("preferred_default")), accession))
    if not candidates:
        raise RuntimeError(f"No released cCRE BED found for {annotation}")
    candidates.sort(reverse=True)
    accession = candidates[0][1]
    return accession, f"https://www.encodeproject.org/files/{accession}/@@download/{accession}.bed.gz"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch current official GENCODE and selected ENCODE cCRE references for v0.8.1")
    parser.add_argument("--output-dir", default="data/reference/v081")
    parser.add_argument("--skip-gencode", action="store_true")
    parser.add_argument("--skip-encode", action="store_true")
    args = parser.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, object] = {"package_version": "0.8.1", "files": []}

    if not args.skip_gencode:
        target = out / "gencode.v50.annotation.gtf.gz"
        _download(GENCODE_URL, target)
        with gzip.open(target, "rt", encoding="utf-8", errors="replace") as handle:
            first_data = next(line for line in handle if not line.startswith("#"))
        if len(first_data.split("\t")) < 9:
            raise RuntimeError("Downloaded GENCODE file failed GTF validation")
        manifest["files"].append({"name": target.name, "source": GENCODE_URL, "type": "GENCODE_GTF", "size_bytes": target.stat().st_size, "sha256": _sha256(target)})
        print("Downloaded", target)

    if not args.skip_encode:
        for label, annotation in ENCODE_ANNOTATIONS.items():
            accession, url = _encode_preferred_bed(annotation)
            target = out / f"encode_ccre_{label}_{accession}.bed.gz"
            _download(url, target)
            manifest["files"].append(
                {"name": target.name, "source": url, "annotation_accession": annotation, "file_accession": accession, "type": "ENCODE_cCRE", "size_bytes": target.stat().st_size, "sha256": _sha256(target)}
            )
            print("Downloaded", target)

    (out / "reference_manifest_v081.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
