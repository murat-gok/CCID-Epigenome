from __future__ import annotations

import csv
import gzip
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd


ANNOTATION_END_CANDIDATES = ("Distance", "distance")
_CPG_RE = re.compile(r"(cg\d+)", re.IGNORECASE)
_DAY_RE = re.compile(r"(?:^|[_\s])d?(\d+)(?:days?)?(?=$|[_\s])", re.IGNORECASE)
_EXP_RE = re.compile(r"(?:^|_)exp(\d+)(?:_|$)", re.IGNORECASE)
_DONOR_RE = re.compile(r"^(O\d+|Y\d+|N\d+|MSC\d+)", re.IGNORECASE)
_REPLICATE_RE = re.compile(r"(?:^|_)R(\d+)(?:_|$)", re.IGNORECASE)


@dataclass(frozen=True)
class MatrixPreprocessResult:
    matrix_path: Path
    metadata_path: Path
    sample_qc_path: Path
    summary_path: Path
    n_input_loci: int
    n_output_loci: int
    n_samples: int


def _normalise_sample_name(name: str) -> str:
    return re.sub(r"\s+", "_", str(name).strip())


def parse_reprogramming_sample_name(sample_name: str, *, source_series: str) -> dict[str, object]:
    """Parse MPTR/Sendai sample names into analysis-ready metadata.

    The parser intentionally maps failed reprogramming samples to condition="failed"
    so the safe-window selector does not accidentally treat them as successful
    reprogramming samples merely because their raw name contains "reprogram".
    """

    raw = str(sample_name).strip()
    norm = _normalise_sample_name(raw)
    lower = norm.lower()

    donor_match = _DONOR_RE.search(norm)
    donor = donor_match.group(1).upper() if donor_match else "unknown"

    day_match = _DAY_RE.search(norm)
    if day_match:
        day: float | int | None = int(day_match.group(1))
    elif re.search(r"(?:^|_)fib(?:_|$)", lower):
        day = 0
    else:
        day = None

    exp_match = _EXP_RE.search(norm)
    batch = f"exp{exp_match.group(1)}" if exp_match else "unspecified"

    if "negative_control" in lower or re.search(r"(?:^|_)fib(?:_|$)", lower):
        condition = "baseline"
    elif "failed_to_transiently_reprogram" in lower or "failing_to_transiently_reprogram" in lower:
        condition = "failed"
    elif "transiently_reprogrammed" in lower:
        condition = "reprogrammed"
    elif "transient_reprogramming_intermediate" in lower:
        condition = "reprogrammed_intermediate"
    elif "ipsc" in lower:
        condition = "ipsc"
    elif "ssea4" in lower:
        condition = "sendai_reprogramming_ssea4"
    elif "cd13" in lower:
        condition = "sendai_reprogramming_cd13"
    else:
        condition = "unknown"

    stage = "intermediate" if "intermediate" in lower else "endpoint"
    if condition == "baseline" and day == 0:
        stage = "baseline"
    if condition == "ipsc":
        stage = "pluripotent"

    if donor.startswith("O"):
        age_group = "older"
    elif donor.startswith("Y") or donor.startswith("N"):
        age_group = "younger"
    else:
        age_group = "unknown"

    return {
        "sample_id": raw,
        "donor_id": donor,
        "day": day,
        "condition": condition,
        "cell_type": "dermal_fibroblast",
        "batch": batch,
        "stage": stage,
        "age_group": age_group,
        "source_series": source_series,
        "raw_sample_name": raw,
    }


def parse_bystander_sample_name(sample_name: str, *, cell_type: str, source_series: str) -> dict[str, object]:
    raw = str(sample_name).strip()
    stem = re.sub(r"\.(beta|detp)$", "", raw, flags=re.IGNORECASE)
    lower = stem.lower()

    if "scramble" in lower or "scrdnmt" in lower:
        treatment = "control"
    elif "hyper" in lower or "target" in lower:
        treatment = "targeted"
    else:
        treatment = "unknown"

    day_match = _DAY_RE.search(stem)
    day = int(day_match.group(1)) if day_match else None

    replicate_match = _REPLICATE_RE.search(stem)
    donor_match = _DONOR_RE.search(stem)
    donor = donor_match.group(1).upper() if donor_match else "unknown"
    if donor.startswith("MSC"):
        replicate = donor
    elif replicate_match:
        replicate = f"R{replicate_match.group(1)}"
    else:
        replicate = stem

    editor = "dCas9-DNMT3A" if ("dnmt" in lower or "hyper" in lower) else "unknown"

    return {
        "sample_id": stem,
        "donor_id": donor,
        "day": day,
        "condition": treatment,
        "treatment": treatment,
        "replicate": replicate,
        "cell_type": cell_type,
        "editor": editor,
        "batch": source_series,
        "source_series": source_series,
        "raw_sample_name": raw,
    }


def _read_expression_file(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    annotation_end = None
    for candidate in ANNOTATION_END_CANDIDATES:
        if candidate in frame.columns:
            annotation_end = frame.columns.get_loc(candidate)
            break
    if annotation_end is None:
        raise ValueError(f"Could not locate the final annotation column ('Distance') in {path}")

    sample_columns = list(frame.columns[annotation_end + 1 :])
    if not sample_columns:
        raise ValueError(f"No expression sample columns found in {path}")

    gene_column = "Probe" if "Probe" in frame.columns else frame.columns[0]
    genes = frame[gene_column].astype(str).str.strip()
    expression = frame[sample_columns].apply(pd.to_numeric, errors="coerce")
    expression.index = genes
    expression = expression.loc[expression.index.notna() & (expression.index != "")]
    expression = expression.groupby(level=0, sort=False).mean()
    return expression.astype(np.float32)


def preprocess_fibroblast_expression(
    input_paths: Sequence[str | Path],
    output_expression: str | Path,
    output_metadata: str | Path,
    output_qc: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    paths = [Path(p) for p in input_paths]
    if not paths:
        raise ValueError("At least one expression input file is required")
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing expression input files: {missing}")

    matrices: list[pd.DataFrame] = []
    for path in paths:
        matrices.append(_read_expression_file(path))

    all_genes = matrices[0].index
    for matrix in matrices[1:]:
        all_genes = all_genes.union(matrix.index, sort=False)

    aligned = [matrix.reindex(all_genes) for matrix in matrices]
    expression = pd.concat(aligned, axis=1)

    duplicate_samples = expression.columns[expression.columns.duplicated()].unique().tolist()
    if duplicate_samples:
        collapsed: dict[str, pd.Series] = {}
        for sample in pd.unique(expression.columns):
            sub = expression.loc[:, expression.columns == sample]
            collapsed[str(sample)] = sub.mean(axis=1)
        expression = pd.DataFrame(collapsed, index=expression.index)

    metadata = pd.DataFrame(
        [parse_reprogramming_sample_name(c, source_series="GSE165177") for c in expression.columns]
    )

    expression_path = Path(output_expression)
    metadata_path = Path(output_metadata)
    qc_path = Path(output_qc)
    for path in (expression_path, metadata_path, qc_path):
        path.parent.mkdir(parents=True, exist_ok=True)

    expression.to_csv(expression_path, compression="infer")
    metadata.to_csv(metadata_path, index=False)

    available_gene_counts = expression.notna().sum(axis=0)
    qc = {
        "source_files": [str(p) for p in paths],
        "n_genes": int(expression.shape[0]),
        "n_samples": int(expression.shape[1]),
        "duplicate_samples_collapsed": duplicate_samples,
        "n_samples_by_condition": metadata["condition"].value_counts(dropna=False).to_dict(),
        "n_samples_by_donor": metadata["donor_id"].value_counts(dropna=False).to_dict(),
        "min_nonmissing_genes_per_sample": int(available_gene_counts.min()),
        "median_nonmissing_genes_per_sample": float(available_gene_counts.median()),
        "max_nonmissing_genes_per_sample": int(available_gene_counts.max()),
    }
    qc_path.write_text(json.dumps(qc, indent=2, default=str))
    return expression, metadata, qc


def read_paired_beta_header(path: str | Path) -> tuple[list[str], list[int], list[int], str]:
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", errors="replace") as handle:
        header = next(csv.reader(handle))
    if len(header) < 3:
        raise ValueError(f"Unexpected methylation matrix header in {path}")
    id_column = header[0].strip().strip('"')

    beta_positions: list[int] = []
    det_positions: list[int] = []
    sample_names: list[str] = []
    i = 1
    while i < len(header):
        beta_name = header[i].strip().strip('"')
        if i + 1 >= len(header):
            raise ValueError(f"Unpaired beta column at position {i} in {path}")
        det_name = header[i + 1].strip().strip('"').lower()
        if "detection" not in det_name and "detp" not in det_name:
            raise ValueError(
                f"Expected a detection-p column after {beta_name!r}, found {header[i + 1]!r}"
            )
        beta_positions.append(i)
        det_positions.append(i + 1)
        sample_names.append(beta_name)
        i += 2
    return sample_names, beta_positions, det_positions, id_column


def _append_csv_gzip(frame: pd.DataFrame, path: Path, *, first: bool) -> None:
    # Appending creates concatenated gzip members, which gzip/pandas read transparently.
    frame.to_csv(
        path,
        mode="wt" if first else "at",
        header=first,
        index=True,
        compression="gzip",
        float_format="%.7g",
    )


def preprocess_paired_beta_matrix(
    input_path: str | Path,
    output_beta: str | Path,
    output_metadata: str | Path,
    output_sample_qc: str | Path,
    output_summary: str | Path,
    *,
    source_series: str,
    cell_type: str = "dermal_fibroblast",
    detection_p_threshold: float = 0.01,
    max_locus_missing: float = 0.20,
    chunksize: int = 5000,
) -> MatrixPreprocessResult:
    """Stream a beta/detection-p paired GEO matrix into a clean beta matrix."""

    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(input_path)
    output_beta = Path(output_beta)
    output_metadata = Path(output_metadata)
    output_sample_qc = Path(output_sample_qc)
    output_summary = Path(output_summary)
    for path in (output_beta, output_metadata, output_sample_qc, output_summary):
        path.parent.mkdir(parents=True, exist_ok=True)

    sample_names, beta_positions, det_positions, _ = read_paired_beta_header(input_path)
    n_columns = 1 + 2 * len(sample_names)

    metadata_rows = []
    for sample in sample_names:
        if cell_type == "dermal_fibroblast":
            metadata_rows.append(parse_reprogramming_sample_name(sample, source_series=source_series))
        else:
            metadata_rows.append(
                parse_bystander_sample_name(sample, cell_type=cell_type, source_series=source_series)
            )
    metadata = pd.DataFrame(metadata_rows)
    metadata.to_csv(output_metadata, index=False)

    total_loci = 0
    output_loci = 0
    detection_fail = np.zeros(len(sample_names), dtype=np.int64)
    out_of_range = np.zeros(len(sample_names), dtype=np.int64)
    nonmissing = np.zeros(len(sample_names), dtype=np.int64)
    duplicate_loci = 0
    seen_loci: set[str] = set()

    output_beta.unlink(missing_ok=True)
    reader = pd.read_csv(
        input_path,
        compression="infer",
        header=None,
        skiprows=1,
        names=list(range(n_columns)),
        chunksize=chunksize,
        low_memory=False,
    )
    first_write = True
    for chunk in reader:
        ids = chunk.iloc[:, 0].astype(str).str.strip().str.strip('"')
        beta = chunk.iloc[:, beta_positions].apply(pd.to_numeric, errors="coerce")
        det = chunk.iloc[:, det_positions].apply(pd.to_numeric, errors="coerce")
        beta.columns = sample_names
        det.columns = sample_names

        # pandas 3.x may expose a read-only NumPy view here.  We mask invalid
        # beta values in-place below, so request an explicit writable copy.
        beta_values = beta.to_numpy(dtype=np.float32, copy=True)
        det_values = det.to_numpy(dtype=np.float32, copy=False)
        detection_bad = ~np.isfinite(det_values) | (det_values > detection_p_threshold)
        range_bad = ~np.isfinite(beta_values) | (beta_values < 0.0) | (beta_values > 1.0)
        bad = detection_bad | range_bad
        beta_values[bad] = np.nan

        detection_fail += detection_bad.sum(axis=0)
        out_of_range += range_bad.sum(axis=0)
        nonmissing += np.isfinite(beta_values).sum(axis=0)
        total_loci += len(chunk)

        missing_fraction = np.mean(~np.isfinite(beta_values), axis=1)
        keep = missing_fraction <= max_locus_missing
        if not np.any(keep):
            continue

        kept_ids = ids.loc[keep].to_numpy(dtype=str)
        duplicate_loci += sum(locus in seen_loci for locus in kept_ids)
        seen_loci.update(kept_ids)

        clean = pd.DataFrame(beta_values[keep], index=kept_ids, columns=sample_names)
        clean.index.name = "locus_id"
        _append_csv_gzip(clean, output_beta, first=first_write)
        first_write = False
        output_loci += int(keep.sum())

    if first_write:
        raise ValueError("No loci passed preprocessing filters")

    denom = max(total_loci, 1)
    sample_qc = metadata[["sample_id", "donor_id", "condition", "cell_type", "source_series"]].copy()
    sample_qc["total_loci"] = total_loci
    sample_qc["valid_loci"] = nonmissing
    sample_qc["valid_fraction"] = nonmissing / denom
    sample_qc["detection_failure_fraction"] = detection_fail / denom
    sample_qc["out_of_range_fraction"] = out_of_range / denom
    sample_qc.to_csv(output_sample_qc, index=False)

    summary = {
        "source_file": str(input_path),
        "source_series": source_series,
        "n_input_loci": total_loci,
        "n_output_loci": output_loci,
        "n_samples": len(sample_names),
        "detection_p_threshold": detection_p_threshold,
        "max_locus_missing": max_locus_missing,
        "duplicate_loci_seen": duplicate_loci,
        "minimum_sample_valid_fraction": float(sample_qc["valid_fraction"].min()),
        "median_sample_valid_fraction": float(sample_qc["valid_fraction"].median()),
        "conditions": metadata["condition"].value_counts(dropna=False).to_dict(),
        "donors": metadata["donor_id"].value_counts(dropna=False).to_dict(),
    }
    output_summary.write_text(json.dumps(summary, indent=2, default=str))
    return MatrixPreprocessResult(
        output_beta,
        output_metadata,
        output_sample_qc,
        output_summary,
        total_loci,
        output_loci,
        len(sample_names),
    )


def _extract_base_cpg(value: object) -> str:
    text = str(value).strip().strip('"')
    match = _CPG_RE.search(text)
    return match.group(1).lower() if match else text


def _read_bystander_wide(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    raw = pd.read_csv(path, compression="infer", low_memory=False)
    id_candidates = [c for c in raw.columns if str(c).lower() in {"probeid", "id_ref", "locus_id"}]
    if not id_candidates:
        raise ValueError(f"No probe ID column found in {path}")
    id_col = id_candidates[0]

    beta_cols = [c for c in raw.columns if str(c).lower().endswith(".beta")]
    det_cols = [c for c in raw.columns if str(c).lower().endswith((".detp", ".detp".lower()))]
    if not beta_cols:
        raise ValueError(f"No .beta columns found in {path}")

    beta = raw[beta_cols].apply(pd.to_numeric, errors="coerce")
    beta.columns = [re.sub(r"\.beta$", "", str(c), flags=re.IGNORECASE) for c in beta_cols]

    det_map: dict[str, str] = {}
    for col in det_cols:
        stem = re.sub(r"\.detp$", "", str(col), flags=re.IGNORECASE)
        det_map[stem] = col

    for sample in beta.columns:
        det_col = det_map.get(sample)
        if det_col is None:
            continue
        det = pd.to_numeric(raw[det_col], errors="coerce")
        beta.loc[det > 0.01, sample] = np.nan

    beta[(beta < 0) | (beta > 1)] = np.nan
    beta.index = raw[id_col].map(_extract_base_cpg)
    duplicate_count = int(beta.index.duplicated().sum())
    if duplicate_count:
        beta = beta.groupby(level=0, sort=False).mean()

    cell_type = "mesenchymal_stromal_cell" if "MSC" in path.name.upper() else "t_cell"
    source_series_match = re.search(r"GSE\d+", path.name, re.IGNORECASE)
    source_series = source_series_match.group(0).upper() if source_series_match else path.stem
    metadata = pd.DataFrame(
        [parse_bystander_sample_name(c, cell_type=cell_type, source_series=source_series) for c in beta.columns]
    )
    qc = {
        "source_file": str(path),
        "n_loci": int(beta.shape[0]),
        "n_samples": int(beta.shape[1]),
        "duplicate_base_cpgs_collapsed": duplicate_count,
        "minimum_sample_valid_fraction": float(beta.notna().mean(axis=0).min()),
        "median_sample_valid_fraction": float(beta.notna().mean(axis=0).median()),
    }
    return beta.astype(np.float32), metadata, qc


def preprocess_bystander_matrices(
    tcell_inputs: Sequence[str | Path],
    msc_input: str | Path,
    output_dir: str | Path,
    *,
    max_locus_missing: float = 0.20,
) -> dict[str, dict[str, object]]:
    output_dir = Path(output_dir)
    qc_dir = output_dir / "qc"
    output_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    tcell_matrices: list[pd.DataFrame] = []
    tcell_meta: list[pd.DataFrame] = []
    tcell_qc: list[dict[str, object]] = []
    for raw_path in [Path(p) for p in tcell_inputs]:
        if not raw_path.exists():
            raise FileNotFoundError(raw_path)
        beta, meta, qc = _read_bystander_wide(raw_path)
        tcell_matrices.append(beta)
        tcell_meta.append(meta)
        tcell_qc.append(qc)

    tcell = pd.concat(tcell_matrices, axis=1, join="outer")
    if tcell.columns.duplicated().any():
        tcell = tcell.T.groupby(level=0, sort=False).mean().T
    tcell = tcell.loc[tcell.isna().mean(axis=1) <= max_locus_missing]
    tcell_metadata = pd.concat(tcell_meta, ignore_index=True).drop_duplicates("sample_id")
    tcell_metadata = tcell_metadata[tcell_metadata["sample_id"].isin(tcell.columns)]
    tcell.to_csv(output_dir / "tcell_beta.csv.gz", compression="gzip", float_format="%.7g")
    tcell_metadata.to_csv(output_dir / "tcell_metadata.csv", index=False)

    msc_path = Path(msc_input)
    if not msc_path.exists():
        raise FileNotFoundError(msc_path)
    msc, msc_metadata, msc_qc = _read_bystander_wide(msc_path)
    msc = msc.loc[msc.isna().mean(axis=1) <= max_locus_missing]
    msc.to_csv(output_dir / "msc_beta.csv.gz", compression="gzip", float_format="%.7g")
    msc_metadata.to_csv(output_dir / "msc_metadata.csv", index=False)

    summary = {
        "tcell": {
            "n_loci": int(tcell.shape[0]),
            "n_samples": int(tcell.shape[1]),
            "source_qc": tcell_qc,
            "treatment_counts": tcell_metadata["treatment"].value_counts().to_dict(),
            "replicate_counts": tcell_metadata["replicate"].value_counts().to_dict(),
        },
        "msc": {
            "n_loci": int(msc.shape[0]),
            "n_samples": int(msc.shape[1]),
            "source_qc": msc_qc,
            "treatment_counts": msc_metadata["treatment"].value_counts().to_dict(),
            "replicate_counts": msc_metadata["replicate"].value_counts().to_dict(),
        },
        "max_locus_missing": max_locus_missing,
    }
    (qc_dir / "bystander_preprocessing_summary.json").write_text(
        json.dumps(summary, indent=2, default=str)
    )
    return summary


def build_overlap_report(
    expression_metadata_path: str | Path,
    methylation_metadata_path: str | Path,
    output_path: str | Path,
) -> dict[str, object]:
    expression_meta = pd.read_csv(expression_metadata_path)
    methylation_meta = pd.read_csv(methylation_metadata_path)
    expression_ids = set(expression_meta["sample_id"].astype(str))
    methylation_ids = set(methylation_meta["sample_id"].astype(str))
    overlap = sorted(expression_ids.intersection(methylation_ids))

    overlap_meta = methylation_meta[methylation_meta["sample_id"].isin(overlap)].copy()
    donor_condition = (
        overlap_meta.groupby(["donor_id", "condition"]).size().unstack(fill_value=0).to_dict("index")
        if not overlap_meta.empty
        else {}
    )
    eligible_donors = []
    for donor, counts in donor_condition.items():
        has_baseline = counts.get("baseline", 0) > 0
        has_reprogrammed = counts.get("reprogrammed", 0) + counts.get("reprogrammed_intermediate", 0) > 0
        if has_baseline and has_reprogrammed:
            eligible_donors.append(str(donor))

    report = {
        "n_expression_samples": len(expression_ids),
        "n_methylation_samples": len(methylation_ids),
        "n_overlapping_samples": len(overlap),
        "overlap_fraction_expression": len(overlap) / max(len(expression_ids), 1),
        "overlap_fraction_methylation": len(overlap) / max(len(methylation_ids), 1),
        "eligible_paired_donors": eligible_donors,
        "n_eligible_paired_donors": len(eligible_donors),
        "donor_condition_counts": donor_condition,
        "overlapping_sample_ids": overlap,
    }
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, default=str))
    return report
