#!/usr/bin/env python3
"""Apply CCID-Epigenome v0.9.0.post5 metadata/longitudinal hotfix.

This patch is intentionally applied in place to an existing v0.9.0.post4 tree so
that unrelated local project files and post4 changes are preserved.
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

POST5 = "0.9.0.post5"


def replace_top_level_function(source: str, name: str, replacement: str) -> str:
    pattern = re.compile(
        rf"(?ms)^def {re.escape(name)}\(.*?(?=^def [A-Za-z_]\w*\(|\Z)"
    )
    match = pattern.search(source)
    if not match:
        raise RuntimeError(f"Top-level function not found: {name}")
    return source[: match.start()] + replacement.rstrip() + "\n\n\n" + source[match.end() :]


METADATA_FRAME = r'''def _metadata_frame(metadata: Mapping[str, list[str]], sample_columns: Sequence[str]) -> pd.DataFrame:
    """Build sample metadata while preserving GEO array identifiers.

    GSE179847 uses ``!Sample_description`` for the Sentrix array/position ID
    present in the deposited processed-matrix column names.  Earlier versions
    discarded that row, preventing the processed beta matrix from being linked
    back to its biological sample title.
    """
    n = len(sample_columns)
    titles = metadata.get("!Sample_title", list(sample_columns))
    accessions = metadata.get("!Sample_geo_accession", [""] * n)
    descriptions = metadata.get("!Sample_description", [""] * n)
    if len(titles) != n:
        titles = list(sample_columns)
    if len(accessions) != n:
        accessions = [""] * n
    if len(descriptions) != n:
        descriptions = [""] * n
    frame = pd.DataFrame(
        {
            "matrix_column": list(sample_columns),
            "sample_title": [_strip_geo_value(x) for x in titles],
            "sample_geo_accession": [_strip_geo_value(x) for x in accessions],
            "sample_description": [_strip_geo_value(x) for x in descriptions],
        }
    )
    # Preserve characteristic rows without assuming a fixed schema.
    char_rows = [values for key, values in metadata.items() if key.startswith("!Sample_characteristics")]
    for idx, values in enumerate(char_rows, start=1):
        if len(values) == n:
            frame[f"characteristic_{idx}"] = [_strip_geo_value(x) for x in values]
    return frame'''


GSE_READER = r'''def read_gse179847_processed_matrix(
    path: str | Path,
    *,
    target_ids: Iterable[str],
    metadata_matrix_file: str | Path | None = None,
) -> GeoSeriesMatrix:
    """Stream the GSE179847 processed matrix and retain beta columns only.

    The deposited file alternates one beta-value column and one detection-p
    column for each of 479 samples.  Only columns ending in `` beta`` are
    retained.  Their suffix is removed so the resulting matrix columns are the
    Sentrix array/position identifiers stored in GEO ``!Sample_description``.
    """
    wanted = {str(x).strip() for x in target_ids}
    rows: list[list[float]] = []
    row_ids: list[str] = []
    header: list[str] | None = None
    beta_positions: list[int] = []
    matrix_columns: list[str] = []

    with _open_text(path) as handle:
        reader = csv.reader(handle, delimiter=",", quotechar='"')
        for fields in reader:
            if not fields:
                continue
            if header is None:
                header = [_strip_geo_value(v) for v in fields]
                for position, column in enumerate(header[1:]):
                    if re.search(r"\s+beta$", column, flags=re.I):
                        beta_positions.append(position)
                        matrix_columns.append(re.sub(r"\s+beta$", "", column, flags=re.I).strip())
                if not beta_positions:
                    # Backward-compatible test/legacy path: earlier releases also
                    # accepted a conventional GEO series matrix here. Re-open it
                    # with the tab-delimited streaming reader.
                    return read_geo_series_matrix(path, target_ids=wanted)
                continue

            row_id = _strip_geo_value(fields[0])
            match = _CPG_RE.search(row_id)
            normalized_id = match.group(0) if match else row_id
            if normalized_id not in wanted:
                continue

            raw_values = fields[1:]
            numeric: list[float] = []
            for position in beta_positions:
                value = raw_values[position] if position < len(raw_values) else ""
                numeric.append(pd.to_numeric(_strip_geo_value(value), errors="coerce"))
            rows.append(numeric)
            row_ids.append(normalized_id)
            if len(set(row_ids)) >= len(wanted):
                break

    if header is None or len(header) < 2:
        raise ValueError(f"No CSV matrix header was found in {path}")

    values = pd.DataFrame(rows, index=row_ids, columns=matrix_columns, dtype=float)
    if not values.empty:
        values = values[~values.index.duplicated(keep="first")]
    values.index.name = "feature_id"

    sample_metadata = pd.DataFrame(
        {
            "matrix_column": matrix_columns,
            "sample_title": matrix_columns,
            "sample_geo_accession": [""] * len(matrix_columns),
            "sample_description": matrix_columns,
        }
    )

    if metadata_matrix_file is not None and Path(metadata_matrix_file).exists():
        metadata_matrix = read_geo_series_matrix(metadata_matrix_file, target_ids=[])
        deposited = metadata_matrix.metadata.copy()
        by_identifier: dict[str, tuple[str, str, str]] = {}

        for row in deposited.itertuples(index=False):
            title = str(getattr(row, "sample_title", ""))
            accession = str(getattr(row, "sample_geo_accession", ""))
            description = str(getattr(row, "sample_description", ""))
            matrix_col = str(getattr(row, "matrix_column", ""))
            payload = (title, accession, description)
            for key in (description, accession, matrix_col):
                normalized_key = _strip_geo_value(key)
                if normalized_key:
                    by_identifier[normalized_key] = payload

        mapped_titles: list[str] = []
        mapped_accessions: list[str] = []
        mapped_descriptions: list[str] = []
        mapped_count = 0
        for column in matrix_columns:
            payload = by_identifier.get(column)
            if payload is None:
                mapped_titles.append(column)
                mapped_accessions.append("")
                mapped_descriptions.append(column)
                continue
            title, accession, description = payload
            mapped_titles.append(title or column)
            mapped_accessions.append(accession)
            mapped_descriptions.append(description or column)
            mapped_count += 1

        if mapped_count == 0:
            raise ValueError(
                "GSE179847 metadata was loaded, but no processed-matrix beta column "
                "matched GEO !Sample_description values"
            )
        sample_metadata["sample_title"] = mapped_titles
        sample_metadata["sample_geo_accession"] = mapped_accessions
        sample_metadata["sample_description"] = mapped_descriptions

    return GeoSeriesMatrix(values=values, metadata=sample_metadata, series_metadata={})'''


TITLE_PARSER = r'''def _parse_multidonor_lifespan_title(
    title: str,
) -> tuple[str | None, str | None, float | None, float | None, bool]:
    """Parse healthy Normal-Control GSE179847 longitudinal sample titles.

    Returns donor, experimental series (sp1/sp2/...), oxygen percentage,
    passage, and eligibility.  Oxygen is retained explicitly so 21% O2 can be
    used as the primary analysis without mixing it with low-oxygen series.
    """
    normalized = str(title).strip().replace(".", "_")
    match = re.search(
        r"(?:^|_)(?:Methylation_)?(HC\d+)_(sp\d+)_Normal_Control_ox(\d+)_p(\d+)(?:$|_)",
        normalized,
        flags=re.I,
    )
    if not match:
        return None, None, None, None, False
    return (
        match.group(1).upper(),
        match.group(2).lower(),
        float(match.group(3)),
        float(match.group(4)),
        True,
    )'''


MULTIDONOR_STATS = r'''def _multidonor_lifespan_statistics(
    matrix: GeoSeriesMatrix,
    candidate_ids: Sequence[str],
) -> pd.DataFrame:
    """Estimate robust passage slopes without pooling distinct culture series.

    Primary estimates use healthy Normal-Control samples cultured at 21% O2.
    A Theil-Sen slope is first estimated for each donor x study-part series.
    Series slopes are then collapsed within donor by the median, followed by a
    median across donors.  This prevents separate experimental series from
    being treated as one continuous passage trajectory.
    """
    meta = matrix.metadata.copy()
    parsed = meta["sample_title"].map(_parse_multidonor_lifespan_title)
    meta[
        [
            "lifespan_donor",
            "lifespan_series",
            "oxygen_percent",
            "passage",
            "eligible_control",
        ]
    ] = pd.DataFrame(parsed.tolist(), index=meta.index)
    meta = meta[
        meta["eligible_control"].fillna(False)
        & meta["passage"].notna()
        & meta["oxygen_percent"].notna()
    ].copy()
    if meta.empty:
        raise ValueError("No healthy Normal-Control lifespan samples were identified in GSE179847")

    primary_meta = meta[np.isclose(pd.to_numeric(meta["oxygen_percent"], errors="coerce"), 21.0)].copy()
    if primary_meta.empty:
        raise ValueError("No 21% O2 healthy Normal-Control lifespan samples were identified in GSE179847")

    def summarize_for_oxygen(cpg: str, selected_meta: pd.DataFrame, prefix: str) -> dict[str, object]:
        series_slopes_by_donor: dict[str, list[float]] = {}
        centered_x: list[float] = []
        centered_y: list[float] = []
        series_details: dict[str, object] = {}
        n_series = 0

        for (donor, series), group in selected_meta.groupby(
            ["lifespan_donor", "lifespan_series"], dropna=True
        ):
            cols = group["matrix_column"].tolist()
            passages = pd.to_numeric(group["passage"], errors="coerce").to_numpy(dtype=float)
            beta = pd.to_numeric(matrix.values.loc[cpg, cols].to_numpy(), errors="coerce")
            valid = np.isfinite(passages) & np.isfinite(beta)
            passages, beta = passages[valid], beta[valid]
            if len(passages) < 3 or len(np.unique(passages)) < 3:
                continue

            series_frame = pd.DataFrame({"passage": passages, "beta": beta}).groupby(
                "passage", as_index=False
            )["beta"].mean()
            if len(series_frame) < 3:
                continue

            slope = float(stats.theilslopes(series_frame["beta"], series_frame["passage"])[0])
            donor_key = str(donor)
            series_key = str(series)
            series_slopes_by_donor.setdefault(donor_key, []).append(slope)
            series_details[f"{prefix}_{donor_key.lower()}_{series_key.lower()}_slope"] = slope
            n_series += 1

            px = series_frame["passage"].to_numpy(dtype=float)
            by = series_frame["beta"].to_numpy(dtype=float)
            centered_x.extend((px - np.mean(px)).tolist())
            centered_y.extend((by - np.mean(by)).tolist())

        donor_slopes: list[float] = []
        donor_details: dict[str, object] = {}
        for donor, slopes in sorted(series_slopes_by_donor.items()):
            donor_slope = float(np.median(slopes))
            donor_slopes.append(donor_slope)
            donor_details[f"{prefix}_{donor.lower()}_slope"] = donor_slope
            donor_details[f"{prefix}_{donor.lower()}_n_series"] = int(len(slopes))

        if donor_slopes:
            median_slope = float(np.median(donor_slopes))
            dominant_sign = int(np.sign(median_slope))
            nonzero = [value for value in donor_slopes if value != 0]
            sign_consistency = (
                float(np.mean([np.sign(value) == dominant_sign for value in nonzero]))
                if nonzero and dominant_sign != 0
                else float("nan")
            )
            sign_successes = int(
                sum(np.sign(value) == dominant_sign for value in nonzero)
            ) if dominant_sign != 0 else 0
            sign_p = _binomial_p(sign_successes, len(nonzero))
        else:
            median_slope = sign_consistency = sign_p = float("nan")

        if len(centered_x) >= 6:
            corr = stats.spearmanr(centered_x, centered_y, nan_policy="omit")
            within_rho, within_p = float(corr.statistic), float(corr.pvalue)
        else:
            within_rho = within_p = float("nan")

        return {
            f"{prefix}_available": bool(donor_slopes),
            f"{prefix}_n_series": int(n_series),
            f"{prefix}_n_donors": int(len(donor_slopes)),
            f"{prefix}_median_theil_sen_slope": median_slope,
            f"{prefix}_donor_sign_consistency": sign_consistency,
            f"{prefix}_donor_sign_binomial_p": sign_p,
            f"{prefix}_within_donor_spearman_rho": within_rho,
            f"{prefix}_within_donor_spearman_p": within_p,
            **donor_details,
            **series_details,
        }

    records: list[dict[str, object]] = []
    oxygen_values = pd.to_numeric(meta["oxygen_percent"], errors="coerce")
    low_oxygen_meta = meta[np.isclose(oxygen_values, 3.0)].copy()

    for cpg in candidate_ids:
        if cpg not in matrix.values.index:
            records.append({"cpg_id": cpg, "multidonor_available": False})
            continue

        primary = summarize_for_oxygen(cpg, primary_meta, "multidonor")
        sensitivity: dict[str, object] = {}
        if not low_oxygen_meta.empty:
            sensitivity = summarize_for_oxygen(cpg, low_oxygen_meta, "multidonor_ox3")
        records.append(
            {
                "cpg_id": cpg,
                "multidonor_primary_oxygen_percent": 21,
                **primary,
                **sensitivity,
            }
        )

    return pd.DataFrame(records)'''


def update_version_file(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    updated, count = re.subn(
        r'__version__\s*=\s*"[^"]+"', f'__version__ = "{POST5}"', text, count=1
    )
    if count != 1:
        raise RuntimeError(f"Could not update version in {path}")
    path.write_text(updated, encoding="utf-8")


def update_pyproject(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    updated, count = re.subn(
        r'(?m)^version\s*=\s*"[^"]+"', f'version = "{POST5}"', text, count=1
    )
    if count != 1:
        raise RuntimeError(f"Could not update project version in {path}")
    path.write_text(updated, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    root = Path(args.project_root).resolve()
    source_path = root / "src" / "ccid" / "orthogonal_validation.py"
    init_path = root / "src" / "ccid" / "__init__.py"
    pyproject_path = root / "pyproject.toml"
    required = [source_path, init_path, pyproject_path]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Project files not found: {missing}")

    init_text = init_path.read_text(encoding="utf-8")
    if POST5 in init_text:
        print(f"CCID-Epigenome {POST5} is already applied.")
        return
    if "0.9.0.post4" not in init_text and not args.force:
        raise RuntimeError(
            "This patch expects v0.9.0.post4. Re-run with --force only after verifying the tree."
        )

    for path in required:
        backup = path.with_suffix(path.suffix + ".post4.bak")
        if not backup.exists():
            shutil.copy2(path, backup)

    source = source_path.read_text(encoding="utf-8")
    for name, replacement in [
        ("_metadata_frame", METADATA_FRAME),
        ("read_gse179847_processed_matrix", GSE_READER),
        ("_parse_multidonor_lifespan_title", TITLE_PARSER),
        ("_multidonor_lifespan_statistics", MULTIDONOR_STATS),
    ]:
        source = replace_top_level_function(source, name, replacement)

    compile(source, str(source_path), "exec")
    source_path.write_text(source, encoding="utf-8")
    update_version_file(init_path)
    update_pyproject(pyproject_path)

    changelog = root / "CHANGELOG.md"
    entry = (
        "\n## 0.9.0.post5\n"
        "- Maps GSE179847 processed beta columns through GEO !Sample_description.\n"
        "- Excludes detection-p columns from the methylation matrix.\n"
        "- Estimates 21% O2 slopes per donor x study-part series before donor aggregation.\n"
        "- Adds separate 3% O2 sensitivity fields when available.\n"
    )
    if changelog.exists() and "## 0.9.0.post5" not in changelog.read_text(encoding="utf-8"):
        with changelog.open("a", encoding="utf-8") as handle:
            handle.write(entry)

    print(f"Applied CCID-Epigenome {POST5} to {root}")
    print("Backups use the suffix .post4.bak")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
