from __future__ import annotations

import csv
import gzip
import io
import json
import math
import re
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


_CPG_RE = re.compile(r"cg\d+", re.I)


@dataclass(frozen=True)
class GeoSeriesMatrix:
    values: pd.DataFrame
    metadata: pd.DataFrame
    series_metadata: dict[str, list[str]]


def _open_text(path: str | Path):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace", newline="")
    return p.open("rt", encoding="utf-8", errors="replace", newline="")


def _strip_geo_value(value: str) -> str:
    return str(value).strip().strip('"')


def _metadata_frame(metadata: Mapping[str, list[str]], sample_columns: Sequence[str]) -> pd.DataFrame:
    n = len(sample_columns)
    titles = metadata.get("!Sample_title", list(sample_columns))
    accessions = metadata.get("!Sample_geo_accession", [""] * n)
    if len(titles) != n:
        titles = list(sample_columns)
    if len(accessions) != n:
        accessions = [""] * n
    frame = pd.DataFrame(
        {
            "matrix_column": list(sample_columns),
            "sample_title": [_strip_geo_value(x) for x in titles],
            "sample_geo_accession": [_strip_geo_value(x) for x in accessions],
        }
    )
    # Preserve characteristic rows without assuming a fixed schema.
    char_rows = [values for key, values in metadata.items() if key.startswith("!Sample_characteristics")]
    for idx, values in enumerate(char_rows, start=1):
        if len(values) == n:
            frame[f"characteristic_{idx}"] = [_strip_geo_value(x) for x in values]
    return frame


def read_geo_series_matrix(
    path: str | Path,
    *,
    target_ids: Iterable[str] | None = None,
) -> GeoSeriesMatrix:
    """Read a GEO Series Matrix, optionally retaining only selected row IDs.

    The parser is streaming and therefore suitable for EPIC-scale matrices when
    target_ids contains a small frozen candidate set.
    """

    wanted = None if target_ids is None else {str(x).strip() for x in target_ids}
    metadata: dict[str, list[str]] = {}
    rows: list[list[object]] = []
    row_ids: list[str] = []
    header: list[str] | None = None
    in_table = False

    with _open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t", quotechar='"')
        for fields in reader:
            if not fields:
                continue
            first = fields[0].strip()
            if not in_table and first.startswith("!"):
                metadata[first] = [_strip_geo_value(v) for v in fields[1:]]
                if first == "!series_matrix_table_begin":
                    in_table = True
                continue
            if not in_table:
                continue
            if first == "!series_matrix_table_end":
                break
            if header is None:
                header = [_strip_geo_value(v) for v in fields]
                continue
            row_id = _strip_geo_value(fields[0])
            if wanted is not None and row_id not in wanted:
                continue
            numeric = [pd.to_numeric(_strip_geo_value(v), errors="coerce") for v in fields[1:]]
            if header is not None and len(numeric) < len(header) - 1:
                numeric.extend([np.nan] * (len(header) - 1 - len(numeric)))
            rows.append(numeric[: len(header) - 1])
            row_ids.append(row_id)

    if header is None:
        raise ValueError(f"No series matrix table was found in {path}")
    sample_columns = header[1:]
    values = pd.DataFrame(rows, index=row_ids, columns=sample_columns, dtype=float)
    values.index.name = "feature_id"
    sample_metadata = _metadata_frame(metadata, sample_columns)
    return GeoSeriesMatrix(values=values, metadata=sample_metadata, series_metadata=metadata)


def read_gse179847_processed_matrix(
    path: str | Path,
    *,
    target_ids: Iterable[str],
    metadata_matrix_file: str | Path | None = None,
) -> GeoSeriesMatrix:
    """Stream the large GSE179847 processed beta CSV and retain target CpGs only.

    The GEO Series Matrix for GSE179847 contains sample metadata but not the
    deposited 866k-by-479 processed beta matrix.  The latter is distributed as
    ``GSE179847_Cell_lifespan_DNAm_processed_matrix.csv.gz``.  This reader keeps
    only the prespecified CpG rows, so memory use is proportional to the frozen
    panel rather than the full two-gigabyte matrix.
    """

    wanted = {str(x).strip() for x in target_ids}
    rows: list[list[float]] = []
    row_ids: list[str] = []
    header: list[str] | None = None

    with _open_text(path) as handle:
        reader = csv.reader(handle, delimiter=",", quotechar='"')
        for fields in reader:
            if not fields:
                continue
            if header is None:
                header = [_strip_geo_value(v) for v in fields]
                continue
            row_id = _strip_geo_value(fields[0])
            match = _CPG_RE.search(row_id)
            normalized_id = match.group(0) if match else row_id
            if normalized_id not in wanted:
                continue
            numeric = [pd.to_numeric(_strip_geo_value(v), errors="coerce") for v in fields[1:]]
            if len(numeric) < len(header) - 1:
                numeric.extend([np.nan] * (len(header) - 1 - len(numeric)))
            rows.append(numeric[: len(header) - 1])
            row_ids.append(normalized_id)
            if len(set(row_ids)) >= len(wanted):
                # All frozen candidates have been found; no need to scan further.
                break

    if header is None or len(header) < 2:
        raise ValueError(f"No CSV matrix header was found in {path}")

    sample_columns = [_strip_geo_value(v) for v in header[1:]]
    values = pd.DataFrame(rows, index=row_ids, columns=sample_columns, dtype=float)
    if not values.empty:
        values = values[~values.index.duplicated(keep="first")]
    values.index.name = "feature_id"

    sample_metadata = pd.DataFrame(
        {
            "matrix_column": sample_columns,
            "sample_title": sample_columns,
            "sample_geo_accession": [""] * len(sample_columns),
        }
    )

    if metadata_matrix_file is not None and Path(metadata_matrix_file).exists():
        metadata_matrix = read_geo_series_matrix(metadata_matrix_file, target_ids=[])
        deposited = metadata_matrix.metadata.copy()
        title_map: dict[str, str] = {}
        accession_map: dict[str, str] = {}
        for row in deposited.itertuples(index=False):
            matrix_col = str(row.matrix_column)
            accession = str(row.sample_geo_accession)
            title = str(row.sample_title)
            title_map[matrix_col] = title
            if accession:
                title_map[accession] = title
                accession_map[accession] = accession
        sample_metadata["sample_title"] = [title_map.get(c, c) for c in sample_columns]
        sample_metadata["sample_geo_accession"] = [accession_map.get(c, "") for c in sample_columns]

    return GeoSeriesMatrix(values=values, metadata=sample_metadata, series_metadata={})


def _find_column(frame: pd.DataFrame, options: Sequence[str], *, required: bool = True) -> str | None:
    lower = {str(c).lower(): str(c) for c in frame.columns}
    for option in options:
        if option.lower() in lower:
            return lower[option.lower()]
    if required:
        raise ValueError(f"Required column was not found. Tried: {options}")
    return None


def load_candidate_directions(
    candidate_file: str | Path,
    *,
    frozen_shortlist_file: str | Path | None = None,
) -> pd.DataFrame:
    candidates = pd.read_csv(candidate_file, compression="infer", low_memory=False)
    cpg_col = _find_column(candidates, ["cpg_id", "probe_id", "ID_REF"])
    delta_col = _find_column(
        candidates,
        [
            "target_delta_beta",
            "discovery_delta_beta",
            "analysis_discovery",
            "mean_delta",
            "delta_beta",
        ],
    )
    rank_col = _find_column(
        candidates,
        ["manuscript_shortlist_rank", "robust_rank", "candidate_rank", "rank"],
        required=False,
    )
    keep = [cpg_col, delta_col] + ([rank_col] if rank_col else [])
    out = candidates[keep].copy()
    out = out.rename(columns={cpg_col: "cpg_id", delta_col: "discovery_delta_beta"})
    if rank_col:
        out = out.rename(columns={rank_col: "candidate_rank"})
    else:
        out["candidate_rank"] = np.arange(1, len(out) + 1)
    out["cpg_id"] = out["cpg_id"].astype(str)
    out["discovery_delta_beta"] = pd.to_numeric(out["discovery_delta_beta"], errors="coerce")
    out = out.dropna(subset=["discovery_delta_beta"]).drop_duplicates("cpg_id")

    if frozen_shortlist_file is not None:
        frozen = pd.read_csv(frozen_shortlist_file)
        frozen_cpg = _find_column(frozen, ["cpg_id", "probe_id"])
        frozen = frozen.rename(columns={frozen_cpg: "cpg_id"})
        frozen["cpg_id"] = frozen["cpg_id"].astype(str)
        out = frozen.merge(out, on="cpg_id", how="left", suffixes=("", "_candidate"))
        if "candidate_rank" not in out.columns:
            if "robust_rank" in out.columns:
                out["candidate_rank"] = pd.to_numeric(out["robust_rank"], errors="coerce")
            elif "biological_shortlist_rank" in out.columns:
                out["candidate_rank"] = pd.to_numeric(out["biological_shortlist_rank"], errors="coerce")
            elif "candidate_rank_candidate" in out.columns:
                out["candidate_rank"] = pd.to_numeric(out["candidate_rank_candidate"], errors="coerce")
            else:
                out["candidate_rank"] = np.arange(1, len(out) + 1)
        elif "candidate_rank_candidate" in out.columns:
            out["candidate_rank"] = pd.to_numeric(out["candidate_rank"], errors="coerce").fillna(
                pd.to_numeric(out["candidate_rank_candidate"], errors="coerce")
            )
        missing = out.loc[out["discovery_delta_beta"].isna(), "cpg_id"].tolist()
        if missing:
            raise ValueError(f"Frozen candidates missing discovery direction: {missing}")
    out["discovery_direction"] = np.sign(out["discovery_delta_beta"]).astype(int)
    return out.sort_values("candidate_rank").reset_index(drop=True)


def _bootstrap_mean_ci(values: Sequence[float], *, iterations: int, seed: int) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) < 2:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    boot = np.mean(rng.choice(arr, size=(iterations, len(arr)), replace=True), axis=1)
    return tuple(np.quantile(boot, [0.025, 0.975]).astype(float))


def _safe_wilcoxon(values: Sequence[float]) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) < 2 or np.allclose(arr, 0):
        return float("nan")
    try:
        return float(stats.wilcoxon(arr, alternative="two-sided", zero_method="wilcox").pvalue)
    except ValueError:
        return float("nan")


def _binomial_p(successes: int, trials: int) -> float:
    if trials <= 0:
        return float("nan")
    return float(stats.binomtest(successes, trials, 0.5, alternative="greater").pvalue)


def _paired_reprogramming_metadata(meta: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in meta.itertuples(index=False):
        title = str(row.sample_title)
        lower = title.lower()
        if "skinfibroblast" not in lower.replace(" ", ""):
            continue
        person = re.search(r"person\s*(\d+)", title, flags=re.I)
        condition = "treated" if "treated" in lower else "normal" if "normal" in lower else "unknown"
        rows.append(
            {
                "matrix_column": row.matrix_column,
                "sample_title": title,
                "sample_geo_accession": row.sample_geo_accession,
                "person_id": person.group(1) if person else "unknown",
                "condition": condition,
            }
        )
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError("No SkinFibroblasts samples were found in the external reprogramming matrix")
    return frame


def validate_external_reprogramming(
    candidate_file: str | Path,
    frozen_shortlist_file: str | Path,
    series_matrix_file: str | Path,
    *,
    bootstrap_iterations: int = 5000,
    seed: int = 20260724,
) -> tuple[pd.DataFrame, dict[str, object]]:
    candidates = load_candidate_directions(candidate_file, frozen_shortlist_file=frozen_shortlist_file)
    matrix = read_geo_series_matrix(series_matrix_file, target_ids=candidates["cpg_id"])
    meta = _paired_reprogramming_metadata(matrix.metadata)
    paired_people = sorted(
        person
        for person, group in meta.groupby("person_id")
        if {"normal", "treated"}.issubset(set(group["condition"]))
    )
    if not paired_people:
        raise ValueError("No normal/treated fibroblast pairs were identified")

    records: list[dict[str, object]] = []
    for idx, candidate in candidates.iterrows():
        cpg = candidate["cpg_id"]
        if cpg not in matrix.values.index:
            records.append(
                {
                    **candidate.to_dict(),
                    "external_available": False,
                    "n_pairs": 0,
                }
            )
            continue
        deltas: list[float] = []
        pair_values: dict[str, float] = {}
        for person in paired_people:
            group = meta[meta["person_id"] == person]
            normal_cols = group.loc[group["condition"] == "normal", "matrix_column"].tolist()
            treated_cols = group.loc[group["condition"] == "treated", "matrix_column"].tolist()
            normal = pd.to_numeric(matrix.values.loc[cpg, normal_cols], errors="coerce").mean()
            treated = pd.to_numeric(matrix.values.loc[cpg, treated_cols], errors="coerce").mean()
            delta = float(treated - normal) if np.isfinite(normal) and np.isfinite(treated) else float("nan")
            pair_values[f"person_{person}_delta"] = delta
            if np.isfinite(delta):
                deltas.append(delta)
        mean_delta = float(np.mean(deltas)) if deltas else float("nan")
        median_delta = float(np.median(deltas)) if deltas else float("nan")
        desired_sign = int(candidate["discovery_direction"])
        same = [int(np.sign(v) == desired_sign) for v in deltas if v != 0]
        ci_low, ci_high = _bootstrap_mean_ci(
            deltas,
            iterations=bootstrap_iterations,
            seed=seed + idx,
        )
        records.append(
            {
                **candidate.to_dict(),
                "external_available": True,
                "n_pairs": len(deltas),
                "external_mean_delta_beta": mean_delta,
                "external_median_delta_beta": median_delta,
                "external_bootstrap_ci_low": ci_low,
                "external_bootstrap_ci_high": ci_high,
                "external_pair_sign_consistency": float(np.mean(same)) if same else float("nan"),
                "external_mean_direction_concordant": bool(np.sign(mean_delta) == desired_sign)
                if np.isfinite(mean_delta) and mean_delta != 0
                else False,
                "external_paired_wilcoxon_p": _safe_wilcoxon(deltas),
                **pair_values,
            }
        )

    result = pd.DataFrame(records)
    available = result[result["external_available"] & result["external_mean_delta_beta"].notna()].copy()
    successes = int(available["external_mean_direction_concordant"].sum()) if len(available) else 0
    rho, rho_p = (float("nan"), float("nan"))
    if len(available) >= 3:
        rho_res = stats.spearmanr(
            available["discovery_delta_beta"],
            available["external_mean_delta_beta"],
            nan_policy="omit",
        )
        rho, rho_p = float(rho_res.statistic), float(rho_res.pvalue)
    summary = {
        "dataset": "GSE142439",
        "design": "four paired human skin fibroblast donors, normal versus four-day OSKMLN treatment",
        "n_frozen_candidates": int(len(candidates)),
        "n_candidates_available": int(len(available)),
        "n_paired_donors": int(len(paired_people)),
        "direction_concordant_candidates": successes,
        "direction_concordance_fraction": float(successes / len(available)) if len(available) else None,
        "direction_concordance_binomial_p": _binomial_p(successes, len(available)),
        "effect_spearman": rho,
        "effect_spearman_p": rho_p,
        "interpretation": "Independent transient-reprogramming direction validation; no candidate weights or ranks are refit.",
    }
    return result, summary


def _parse_day_title(title: str) -> tuple[str, float | None]:
    lower = title.lower()
    condition = "normal_glucose" if "_ng_" in lower else "high_glucose" if "_hg_" in lower else "unknown"
    match = re.search(r"_(\d+)days", lower)
    return condition, float(match.group(1)) if match else None


def _parse_age_title(title: str) -> tuple[float | None, str | None, int | None]:
    age_match = re.search(r"_(\d+)yrs", title, flags=re.I)
    passage_match = re.search(r"_P(\d+)", title, flags=re.I)
    # Site is the second token after HDF- in GSE22595 titles.
    site_match = re.search(r"HDF-([^_]+)", title, flags=re.I)
    return (
        float(age_match.group(1)) if age_match else None,
        site_match.group(1).lower() if site_match else None,
        int(passage_match.group(1)) if passage_match else None,
    )


def _slope_statistics(x: Sequence[float], y: Sequence[float]) -> dict[str, float]:
    xv = np.asarray(x, dtype=float)
    yv = np.asarray(y, dtype=float)
    valid = np.isfinite(xv) & np.isfinite(yv)
    xv, yv = xv[valid], yv[valid]
    if len(xv) < 3 or len(np.unique(xv)) < 2:
        return {
            "n": int(len(xv)),
            "spearman_rho": float("nan"),
            "spearman_p": float("nan"),
            "theil_sen_slope": float("nan"),
            "slope_ci_low": float("nan"),
            "slope_ci_high": float("nan"),
        }
    corr = stats.spearmanr(xv, yv, nan_policy="omit")
    slope, intercept, low, high = stats.theilslopes(yv, xv, alpha=0.95)
    return {
        "n": int(len(xv)),
        "spearman_rho": float(corr.statistic),
        "spearman_p": float(corr.pvalue),
        "theil_sen_slope": float(slope),
        "slope_ci_low": float(low),
        "slope_ci_high": float(high),
    }


def _parse_multidonor_lifespan_title(title: str) -> tuple[str | None, float | None, bool]:
    """Parse untreated healthy-control samples from GSE179847 titles."""
    normalized = str(title).strip().replace(".", "_")
    match = re.search(
        r"(?:^|_)(?:Methylation_)?(HC\d+)_sp\d+_Normal_Control_ox21_p(\d+)(?:$|_)",
        normalized,
        flags=re.I,
    )
    if not match:
        return None, None, False
    return match.group(1).upper(), float(match.group(2)), True


def _multidonor_lifespan_statistics(
    matrix: GeoSeriesMatrix,
    candidate_ids: Sequence[str],
) -> pd.DataFrame:
    meta = matrix.metadata.copy()
    parsed = meta["sample_title"].map(_parse_multidonor_lifespan_title)
    meta[["lifespan_donor", "passage", "eligible_control"]] = pd.DataFrame(
        parsed.tolist(), index=meta.index
    )
    meta = meta[meta["eligible_control"].fillna(False) & meta["passage"].notna()].copy()
    if meta.empty:
        raise ValueError("No untreated healthy-control lifespan samples were identified in GSE179847")

    records: list[dict[str, object]] = []
    for cpg in candidate_ids:
        if cpg not in matrix.values.index:
            records.append({"cpg_id": cpg, "multidonor_available": False})
            continue
        donor_slopes: list[float] = []
        centered_x: list[float] = []
        centered_y: list[float] = []
        donor_details: dict[str, object] = {}
        for donor, group in meta.groupby("lifespan_donor"):
            cols = group["matrix_column"].tolist()
            passages = pd.to_numeric(group["passage"], errors="coerce").to_numpy(dtype=float)
            beta = pd.to_numeric(matrix.values.loc[cpg, cols].to_numpy(), errors="coerce")
            valid = np.isfinite(passages) & np.isfinite(beta)
            passages, beta = passages[valid], beta[valid]
            if len(passages) < 3 or len(np.unique(passages)) < 2:
                continue
            # Collapse duplicate passages before estimating a donor-specific robust slope.
            donor_frame = pd.DataFrame({"passage": passages, "beta": beta}).groupby(
                "passage", as_index=False
            )["beta"].mean()
            if len(donor_frame) < 3:
                continue
            slope = float(stats.theilslopes(donor_frame["beta"], donor_frame["passage"])[0])
            donor_slopes.append(slope)
            donor_details[f"{donor.lower()}_slope"] = slope
            px = donor_frame["passage"].to_numpy(dtype=float)
            by = donor_frame["beta"].to_numpy(dtype=float)
            centered_x.extend((px - np.mean(px)).tolist())
            centered_y.extend((by - np.mean(by)).tolist())
        if donor_slopes:
            median_slope = float(np.median(donor_slopes))
            dominant_sign = int(np.sign(median_slope))
            nonzero = [v for v in donor_slopes if v != 0]
            sign_consistency = (
                float(np.mean([np.sign(v) == dominant_sign for v in nonzero])) if nonzero else np.nan
            )
            sign_successes = int(sum(np.sign(v) == dominant_sign for v in nonzero))
            sign_p = _binomial_p(sign_successes, len(nonzero))
        else:
            median_slope = sign_consistency = sign_p = float("nan")
        if len(centered_x) >= 6:
            corr = stats.spearmanr(centered_x, centered_y, nan_policy="omit")
            within_rho, within_p = float(corr.statistic), float(corr.pvalue)
        else:
            within_rho = within_p = float("nan")
        records.append(
            {
                "cpg_id": cpg,
                "multidonor_available": bool(donor_slopes),
                "multidonor_n_donors": int(len(donor_slopes)),
                "multidonor_median_theil_sen_slope": median_slope,
                "multidonor_donor_sign_consistency": sign_consistency,
                "multidonor_donor_sign_binomial_p": sign_p,
                "multidonor_within_donor_spearman_rho": within_rho,
                "multidonor_within_donor_spearman_p": within_p,
                **donor_details,
            }
        )
    return pd.DataFrame(records)


def validate_natural_aging(
    candidate_file: str | Path,
    frozen_shortlist_file: str | Path,
    *,
    replicative_matrix_file: str | Path | None = None,
    chronological_matrix_file: str | Path | None = None,
    multidonor_matrix_file: str | Path | None = None,
    multidonor_metadata_file: str | Path | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    candidates = load_candidate_directions(candidate_file, frozen_shortlist_file=frozen_shortlist_file)
    result = candidates.copy()

    if replicative_matrix_file is not None:
        rep = read_geo_series_matrix(replicative_matrix_file, target_ids=candidates["cpg_id"])
        rep_meta = rep.metadata.copy()
        parsed = rep_meta["sample_title"].map(_parse_day_title)
        rep_meta[["aging_condition", "culture_day"]] = pd.DataFrame(parsed.tolist(), index=rep_meta.index)
        rep_meta = rep_meta[(rep_meta["aging_condition"] == "normal_glucose") & rep_meta["culture_day"].notna()]
        rep_rows: list[dict[str, object]] = []
        for cpg in candidates["cpg_id"]:
            if cpg not in rep.values.index:
                rep_rows.append({"cpg_id": cpg, "replicative_available": False})
                continue
            data = pd.DataFrame(
                {
                    "day": rep_meta["culture_day"].astype(float),
                    "beta": pd.to_numeric(
                        rep.values.loc[cpg, rep_meta["matrix_column"].tolist()].to_numpy(), errors="coerce"
                    ),
                }
            )
            # Average technical replicates at the same day before trend estimation.
            data = data.groupby("day", as_index=False)["beta"].mean()
            stat = _slope_statistics(data["day"], data["beta"])
            rep_rows.append(
                {
                    "cpg_id": cpg,
                    "replicative_available": True,
                    **{f"replicative_{k}": v for k, v in stat.items()},
                }
            )
        result = result.merge(pd.DataFrame(rep_rows), on="cpg_id", how="left")

    if chronological_matrix_file is not None:
        age = read_geo_series_matrix(chronological_matrix_file, target_ids=candidates["cpg_id"])
        age_meta = age.metadata.copy()
        parsed = age_meta["sample_title"].map(_parse_age_title)
        age_meta[["chronological_age", "anatomical_site", "passage"]] = pd.DataFrame(
            parsed.tolist(), index=age_meta.index
        )
        # Cross-sectional comparison is restricted to passage 3; the P21 sample is a culture-aging sample.
        age_meta = age_meta[(age_meta["passage"] == 3) & age_meta["chronological_age"].notna()]
        age_rows: list[dict[str, object]] = []
        for cpg in candidates["cpg_id"]:
            if cpg not in age.values.index:
                age_rows.append({"cpg_id": cpg, "chronological_available": False})
                continue
            y = pd.to_numeric(
                age.values.loc[cpg, age_meta["matrix_column"].tolist()].to_numpy(), errors="coerce"
            )
            stat = _slope_statistics(age_meta["chronological_age"].astype(float), y)
            site_data = age_meta[["anatomical_site", "chronological_age"]].copy()
            site_data["beta"] = y
            site_data["age_group"] = np.where(site_data["chronological_age"] < 30, "young", "old")
            site_deltas: list[float] = []
            site_detail: dict[str, object] = {}
            for site, group in site_data.groupby("anatomical_site", dropna=True):
                if not {"young", "old"}.issubset(set(group["age_group"])):
                    continue
                young_mean = pd.to_numeric(
                    group.loc[group["age_group"] == "young", "beta"], errors="coerce"
                ).mean()
                old_mean = pd.to_numeric(
                    group.loc[group["age_group"] == "old", "beta"], errors="coerce"
                ).mean()
                delta = float(old_mean - young_mean)
                if np.isfinite(delta):
                    site_deltas.append(delta)
                    site_detail[f"chronological_site_{site}_old_minus_young"] = delta
            site_median = float(np.median(site_deltas)) if site_deltas else float("nan")
            site_sign = int(np.sign(site_median)) if np.isfinite(site_median) and site_median != 0 else 0
            nonzero = [v for v in site_deltas if v != 0]
            site_consistency = (
                float(np.mean([np.sign(v) == site_sign for v in nonzero]))
                if nonzero and site_sign != 0
                else float("nan")
            )
            age_rows.append(
                {
                    "cpg_id": cpg,
                    "chronological_available": True,
                    **{f"chronological_{k}": v for k, v in stat.items()},
                    "chronological_site_matched_n_sites": int(len(site_deltas)),
                    "chronological_site_median_old_minus_young": site_median,
                    "chronological_site_sign_consistency": site_consistency,
                    "chronological_site_wilcoxon_p": _safe_wilcoxon(site_deltas),
                    **site_detail,
                }
            )
        result = result.merge(pd.DataFrame(age_rows), on="cpg_id", how="left")

    if multidonor_matrix_file is not None:
        multi_path = Path(multidonor_matrix_file)
        if "processed_matrix.csv" in multi_path.name.lower() or multi_path.name.lower().endswith(".csv.gz"):
            multi = read_gse179847_processed_matrix(
                multi_path,
                target_ids=candidates["cpg_id"],
                metadata_matrix_file=multidonor_metadata_file,
            )
        else:
            multi = read_geo_series_matrix(multi_path, target_ids=candidates["cpg_id"])
        if multi.values.empty:
            raise ValueError(
                "The GSE179847 input contains no frozen CpG rows. Use the deposited "
                "GSE179847_Cell_lifespan_DNAm_processed_matrix.csv.gz file, not the small "
                "Series Matrix metadata file."
            )
        multi_rows = _multidonor_lifespan_statistics(multi, candidates["cpg_id"].tolist())
        result = result.merge(multi_rows, on="cpg_id", how="left")
        slope_col = "multidonor_median_theil_sen_slope"
        consistency_col = "multidonor_donor_sign_consistency"
        p_col = "multidonor_within_donor_spearman_p"
        if slope_col not in result.columns:
            result[slope_col] = np.nan
        if consistency_col not in result.columns:
            result[consistency_col] = np.nan
        if p_col not in result.columns:
            result[p_col] = np.nan
        result["multidonor_reversal_concordant"] = (
            np.sign(pd.to_numeric(result[slope_col], errors="coerce"))
            == -np.sign(result["discovery_delta_beta"])
        ) & (pd.to_numeric(result[consistency_col], errors="coerce") >= 0.67)
        result["multidonor_within_donor_spearman_fdr"] = np.nan
        valid_multi = result[p_col].notna()
        if valid_multi.any():
            result.loc[valid_multi, "multidonor_within_donor_spearman_fdr"] = multipletests(
                result.loc[valid_multi, p_col], method="fdr_bh"
            )[1]

    for prefix in ["replicative", "chronological"]:
        slope_col = f"{prefix}_theil_sen_slope"
        if slope_col in result.columns:
            direction_value = pd.to_numeric(result[slope_col], errors="coerce")
            if prefix == "chronological" and "chronological_site_median_old_minus_young" in result.columns:
                site_value = pd.to_numeric(
                    result["chronological_site_median_old_minus_young"], errors="coerce"
                )
                direction_value = site_value.where(site_value.notna(), direction_value)
            result[f"{prefix}_reversal_concordant"] = (
                np.sign(direction_value) == -np.sign(result["discovery_delta_beta"])
            )
            if prefix == "chronological" and "chronological_site_sign_consistency" in result.columns:
                consistency = pd.to_numeric(
                    result["chronological_site_sign_consistency"], errors="coerce"
                )
                enough_sites = pd.to_numeric(
                    result["chronological_site_matched_n_sites"], errors="coerce"
                ).ge(2)
                result[f"{prefix}_reversal_concordant"] &= (~enough_sites) | consistency.ge(0.50)
            p_col = f"{prefix}_spearman_p"
            fdr_col = f"{prefix}_spearman_fdr"
            result[fdr_col] = np.nan
            valid = result[p_col].notna()
            if valid.any():
                result.loc[valid, fdr_col] = multipletests(result.loc[valid, p_col], method="fdr_bh")[1]

    reversal_columns = [c for c in result.columns if c.endswith("_reversal_concordant")]
    if reversal_columns:
        result["natural_aging_reversal_support_count"] = result[reversal_columns].fillna(False).sum(axis=1)
        available_cols = [c.replace("_reversal_concordant", "_available") for c in reversal_columns]
        result["natural_aging_datasets_available"] = sum(
            result[c].fillna(False).astype(int) for c in available_cols if c in result.columns
        )

    summary: dict[str, object] = {
        "n_frozen_candidates": int(len(candidates)),
        "replicative_dataset": "GSE131280" if replicative_matrix_file is not None else None,
        "chronological_dataset": "GSE22595" if chronological_matrix_file is not None else None,
        "multidonor_dataset": "GSE179847" if multidonor_matrix_file is not None else None,
        "interpretation": "A candidate is directionally supported when the partial-reprogramming delta opposes the aging slope.",
    }
    for prefix in ["replicative", "chronological", "multidonor"]:
        col = f"{prefix}_reversal_concordant"
        avail = f"{prefix}_available"
        if col in result.columns:
            subset = result[result[avail].fillna(False)]
            successes = int(subset[col].fillna(False).sum())
            summary[f"{prefix}_available_candidates"] = int(len(subset))
            summary[f"{prefix}_reversal_concordant_candidates"] = successes
            summary[f"{prefix}_binomial_p"] = _binomial_p(successes, len(subset))
    return result, summary


def read_platform_annotation(path: str | Path) -> pd.DataFrame:
    """Read a GEO platform annotation table from plain TSV or SOFT text.

    GEO's ``view=data&form=text`` endpoint returns a SOFT document containing
    metadata lines before the tabular section.  Passing that whole file to
    :func:`pandas.read_csv` produces variable-field parser errors.  This reader
    extracts the text between ``!platform_table_begin`` and
    ``!platform_table_end`` when those markers are present, while retaining
    compatibility with ordinary ``.txt``/``.gz`` annotation tables.
    """
    from io import StringIO

    annotation_path = Path(path)
    opener = gzip.open if annotation_path.suffix.lower() == ".gz" else open
    with opener(annotation_path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    begin = next(
        (i for i, line in enumerate(lines) if line.strip().lower() == "!platform_table_begin"),
        None,
    )
    if begin is not None:
        end = next(
            (
                i
                for i in range(begin + 1, len(lines))
                if lines[i].strip().lower() == "!platform_table_end"
            ),
            len(lines),
        )
        table_lines = lines[begin + 1 : end]
    else:
        # Plain exported tables may have leading comments but no SOFT markers.
        header_index = next(
            (
                i
                for i, line in enumerate(lines)
                if "\t" in line
                and line.lstrip('\ufeff"').upper().startswith(("ID\t", "ID_REF\t", "PROBE_ID\t", "ILMNID\t"))
            ),
            None,
        )
        if header_index is None:
            raise ValueError(f"Could not locate a platform annotation table in {annotation_path}")
        table_lines = lines[header_index:]

    table_text = "".join(line for line in table_lines if line.strip())
    if not table_text.strip():
        raise ValueError(f"Platform annotation table is empty in {annotation_path}")

    frame = pd.read_csv(StringIO(table_text), sep="\t", low_memory=False)
    frame.columns = [str(column).strip().strip('\ufeff"') for column in frame.columns]
    id_col = _find_column(frame, ["ID", "ID_REF", "Probe_Id", "IlmnID"])
    gene_col = _find_column(frame, ["Gene Symbol", "Gene symbol", "Symbol", "GENE_SYMBOL"])
    out = frame[[id_col, gene_col]].copy()
    out.columns = ["probe_id", "gene_symbol"]
    out["probe_id"] = out["probe_id"].astype(str).str.strip().str.strip('\"')
    out["gene_symbol"] = out["gene_symbol"].fillna("").astype(str)
    out["gene_symbol"] = out["gene_symbol"].str.split(r"[;,|]")
    out = out.explode("gene_symbol")
    out["gene_symbol"] = out["gene_symbol"].astype(str).str.strip().str.strip('\"').str.upper()
    return out[out["gene_symbol"].ne("")].drop_duplicates()


def _normalise_fibroblast_sample(title: str) -> str:
    return re.sub(r"\s*\[gene expression\]\s*$", "", str(title), flags=re.I).strip()


def validate_fibroblast_eqtm(
    candidate_file: str | Path,
    frozen_shortlist_file: str | Path,
    methylation_matrix_file: str | Path,
    expression_matrix_file: str | Path,
    expression_annotation_file: str | Path,
    gene_map_file: str | Path,
) -> tuple[pd.DataFrame, dict[str, object]]:
    candidates = load_candidate_directions(candidate_file, frozen_shortlist_file=frozen_shortlist_file)
    gene_map = pd.read_csv(gene_map_file)
    cpg_col = _find_column(gene_map, ["cpg_id", "probe_id"])
    direct_col = _find_column(gene_map, ["direct_gene", "gene"], required=False)
    nearest_col = _find_column(gene_map, ["nearest_tss_gene", "nearest_gene"], required=False)
    pairs: list[dict[str, str]] = []
    for _, row in gene_map.iterrows():
        cpg = str(row[cpg_col])
        for link_type, col in [("direct_gene", direct_col), ("nearest_tss_gene", nearest_col)]:
            if col is None:
                continue
            value = row.get(col)
            if value is None or pd.isna(value):
                continue
            for gene in re.split(r"[;,|]", str(value)):
                gene = gene.strip().upper()
                if gene and not gene.startswith("ENSG"):
                    pairs.append({"cpg_id": cpg, "gene": gene, "gene_link_type": link_type})
    pair_frame = pd.DataFrame(pairs).drop_duplicates()
    if not pair_frame.empty:
        pair_frame["_link_priority"] = pair_frame["gene_link_type"].map(
            {"direct_gene": 0, "nearest_tss_gene": 1}
        ).fillna(2)
        pair_frame = (
            pair_frame.sort_values("_link_priority")
            .drop_duplicates(["cpg_id", "gene"], keep="first")
            .drop(columns="_link_priority")
        )
    pair_frame = pair_frame[pair_frame["cpg_id"].isin(candidates["cpg_id"])]
    if pair_frame.empty:
        return pair_frame, {"n_tested_pairs": 0, "interpretation": "No array-compatible gene symbols were available."}

    methyl = read_geo_series_matrix(methylation_matrix_file, target_ids=pair_frame["cpg_id"].unique())
    annotation = read_platform_annotation(expression_annotation_file)
    wanted_genes = set(pair_frame["gene"])
    wanted_probes = annotation.loc[annotation["gene_symbol"].isin(wanted_genes), "probe_id"].unique()
    expr = read_geo_series_matrix(expression_matrix_file, target_ids=wanted_probes)

    meth_meta = methyl.metadata.copy()
    expr_meta = expr.metadata.copy()
    meth_meta["sample_key"] = meth_meta["sample_title"].map(_normalise_fibroblast_sample)
    expr_meta["sample_key"] = expr_meta["sample_title"].map(_normalise_fibroblast_sample)
    common = sorted(set(meth_meta["sample_key"]) & set(expr_meta["sample_key"]))
    if len(common) < 10:
        raise ValueError(f"Only {len(common)} matched fibroblast samples were found between methylation and expression")
    meth_col = meth_meta.set_index("sample_key")["matrix_column"].to_dict()
    expr_col = expr_meta.set_index("sample_key")["matrix_column"].to_dict()

    probe_gene = annotation[annotation["probe_id"].isin(expr.values.index)].copy()
    rows: list[dict[str, object]] = []
    for pair in pair_frame.itertuples(index=False):
        if pair.cpg_id not in methyl.values.index:
            rows.append({**pair._asdict(), "methylation_available": False, "expression_available": False})
            continue
        probes = probe_gene.loc[probe_gene["gene_symbol"] == pair.gene, "probe_id"].tolist()
        probes = [p for p in probes if p in expr.values.index]
        if not probes:
            rows.append({**pair._asdict(), "methylation_available": True, "expression_available": False})
            continue
        meth_values = np.array([methyl.values.loc[pair.cpg_id, meth_col[k]] for k in common], dtype=float)
        expr_values = np.vstack(
            [np.array([expr.values.loc[p, expr_col[k]] for k in common], dtype=float) for p in probes]
        )
        gene_expression = np.nanmean(expr_values, axis=0)
        valid = np.isfinite(meth_values) & np.isfinite(gene_expression)
        if valid.sum() < 10:
            rho = pvalue = float("nan")
        else:
            corr = stats.spearmanr(meth_values[valid], gene_expression[valid], nan_policy="omit")
            rho, pvalue = float(corr.statistic), float(corr.pvalue)
        rows.append(
            {
                **pair._asdict(),
                "methylation_available": True,
                "expression_available": True,
                "n_matched_samples": int(valid.sum()),
                "n_expression_probes": int(len(probes)),
                "eqtm_spearman_rho": rho,
                "eqtm_p": pvalue,
            }
        )
    result = pd.DataFrame(rows)
    result["eqtm_fdr"] = np.nan
    valid = result["eqtm_p"].notna() if "eqtm_p" in result.columns else pd.Series(False, index=result.index)
    if valid.any():
        result.loc[valid, "eqtm_fdr"] = multipletests(result.loc[valid, "eqtm_p"], method="fdr_bh")[1]
    summary = {
        "datasets": ["GSE52025", "GSE53243"],
        "n_matched_fibroblast_samples": int(len(common)),
        "n_tested_pairs": int(valid.sum()),
        "n_nominal_p_lt_0_05": int((result.get("eqtm_p", pd.Series(dtype=float)) < 0.05).sum()),
        "n_fdr_lt_0_10": int((result.get("eqtm_fdr", pd.Series(dtype=float)) < 0.10).sum()),
        "interpretation": "Cross-sectional fibroblast eQTM support is associative and does not establish mediation or causality.",
    }
    return result, summary


def _open_text_auto(path: str | Path):
    """Open plain or gzip-compressed text without loading the file into memory."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def _normalise_re2g_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def _resolve_re2g_columns(header: Sequence[str] | None) -> dict[str, int | None]:
    """Resolve ENCODE/ABC element-gene-link columns from a header when available."""
    if header is None:
        return {
            "chromosome": 0,
            "start": 1,
            "end": 2,
            "name": 3,
            "target_gene": 3,
            "score": 4,
        }
    normalised = {_normalise_re2g_header(value): idx for idx, value in enumerate(header)}

    def first(*aliases: str) -> int | None:
        for alias in aliases:
            idx = normalised.get(_normalise_re2g_header(alias))
            if idx is not None:
                return idx
        return None

    columns = {
        "chromosome": first("chr", "chrom", "chromosome"),
        "start": first("start", "chromStart"),
        "end": first("end", "chromEnd"),
        "name": first("name", "element", "element_name"),
        "target_gene": first("TargetGene", "target_gene", "gene", "gene_symbol"),
        # ENCODE rE2G files expose a final combined Score. Fall back to ABC/power-law.
        "score": first("Score", "ABC.Score", "powerlaw.Score", "ABCScore", "powerlawScore"),
    }
    required = ["chromosome", "start", "end", "target_gene"]
    missing = [key for key in required if columns[key] is None]
    if missing:
        raise ValueError(f"rE2G header is missing required columns {missing}: {header}")
    return columns


def _stream_re2g_overlaps(
    path: str | Path,
    *,
    biosample: str,
    candidates_by_chromosome: dict[str, list[tuple[str, int]]],
) -> list[dict[str, object]]:
    """Stream a large BED-like rE2G file and retain only overlaps with frozen CpGs.

    ENCODE rE2G/ABC files can contain millions of rows. Reading both files into
    pandas at once can exhaust WSL memory. This parser keeps only the header and
    overlap records for the small frozen candidate panel.
    """
    records: list[dict[str, object]] = []
    header: list[str] | None = None
    columns: dict[str, int | None] | None = None
    with _open_text_auto(path) as handle:
        for raw_line in handle:
            if not raw_line.strip():
                continue
            if raw_line.startswith("#"):
                possible = raw_line[1:].rstrip("\n\r").split("\t")
                normalised = {_normalise_re2g_header(value) for value in possible}
                if {"chr", "start", "end"}.issubset(normalised) or {
                    "chromosome", "start", "end"
                }.issubset(normalised):
                    header = possible
                    columns = _resolve_re2g_columns(header)
                continue
            fields = raw_line.rstrip("\n\r").split("\t")
            if columns is None:
                columns = _resolve_re2g_columns(header)
            max_index = max(idx for idx in columns.values() if idx is not None)
            if len(fields) <= max_index:
                continue
            chrom = str(fields[int(columns["chromosome"])]).strip()
            chrom = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)
            candidate_points = candidates_by_chromosome.get(chrom)
            if not candidate_points:
                continue
            try:
                start0 = int(float(fields[int(columns["start"])]))
                end0 = int(float(fields[int(columns["end"])]))
            except (TypeError, ValueError):
                continue
            if end0 <= start0:
                continue
            # BED intervals are zero-based half-open. Candidate positions are one-based.
            matched_cpgs = [
                cpg_id
                for cpg_id, position_1based in candidate_points
                if start0 <= int(position_1based) - 1 < end0
            ]
            if not matched_cpgs:
                continue
            target_gene = str(fields[int(columns["target_gene"])]).strip().upper()
            if not target_gene or target_gene in {".", "NA", "NAN"}:
                target_gene = ""
            score = np.nan
            score_idx = columns.get("score")
            if score_idx is not None and int(score_idx) < len(fields):
                try:
                    score = float(fields[int(score_idx)])
                except (TypeError, ValueError):
                    score = np.nan
            name = ""
            name_idx = columns.get("name")
            if name_idx is not None and int(name_idx) < len(fields):
                name = str(fields[int(name_idx)]).strip()
            for cpg_id in matched_cpgs:
                records.append(
                    {
                        "cpg_id": cpg_id,
                        "re2g_linked": True,
                        "re2g_target_gene": target_gene or np.nan,
                        "re2g_biosample": biosample,
                        "re2g_score": score,
                        "re2g_element": name or f"chr{chrom}:{start0 + 1}-{end0}",
                    }
                )
    return records


def parse_bedlike_re2g(path: str | Path, *, biosample: str) -> pd.DataFrame:
    """Parse a BED-like rE2G file.

    This compatibility helper remains suitable for small files. Large ENCODE
    files should be processed through :func:`overlap_candidates_with_re2g`,
    which streams them without materialising the complete table.
    """
    header: list[str] | None = None
    rows: list[dict[str, object]] = []
    with _open_text_auto(path) as handle:
        for raw_line in handle:
            if not raw_line.strip():
                continue
            if raw_line.startswith("#"):
                possible = raw_line[1:].rstrip("\n\r").split("\t")
                normalised = {_normalise_re2g_header(value) for value in possible}
                if "start" in normalised and "end" in normalised:
                    header = possible
                continue
            fields = raw_line.rstrip("\n\r").split("\t")
            columns = _resolve_re2g_columns(header)
            max_index = max(idx for idx in columns.values() if idx is not None)
            if len(fields) <= max_index:
                continue
            try:
                start0 = int(float(fields[int(columns["start"])]))
                end0 = int(float(fields[int(columns["end"])]))
            except (TypeError, ValueError):
                continue
            score = np.nan
            score_idx = columns.get("score")
            if score_idx is not None:
                try:
                    score = float(fields[int(score_idx)])
                except (TypeError, ValueError):
                    score = np.nan
            rows.append(
                {
                    "chromosome": re.sub(
                        r"^chr", "", str(fields[int(columns["chromosome"])]), flags=re.IGNORECASE
                    ),
                    "start": start0 + 1,
                    "end": end0,
                    "target_gene": str(fields[int(columns["target_gene"])]).strip().upper(),
                    "biosample": biosample,
                    "re2g_score": score,
                }
            )
    return pd.DataFrame(rows)


def overlap_candidates_with_re2g(
    candidate_file: str | Path,
    frozen_shortlist_file: str | Path,
    re2g_files: Sequence[str | Path],
    *,
    biosample_names: Sequence[str] | None = None,
) -> pd.DataFrame:
    candidates = pd.read_csv(candidate_file, compression="infer", low_memory=False)
    frozen = load_candidate_directions(candidate_file, frozen_shortlist_file=frozen_shortlist_file)
    candidates = candidates[candidates[_find_column(candidates, ["cpg_id"])].astype(str).isin(frozen["cpg_id"])].copy()
    cpg_col = _find_column(candidates, ["cpg_id"])
    chr_col = _find_column(candidates, ["chromosome", "chrom", "chr"])
    pos_col = _find_column(candidates, ["position", "coordinate", "pos", "mapinfo"])
    candidates = candidates.rename(columns={cpg_col: "cpg_id", chr_col: "chromosome", pos_col: "position"})
    candidates["chromosome"] = candidates["chromosome"].astype(str).str.replace("^chr", "", regex=True)
    candidates["position"] = pd.to_numeric(candidates["position"], errors="coerce")
    candidates = candidates.dropna(subset=["position"]).copy()
    candidates["position"] = candidates["position"].astype(int)

    names = list(biosample_names or [Path(p).stem for p in re2g_files])
    if len(names) != len(re2g_files):
        raise ValueError("biosample_names must match re2g_files")
    candidates_by_chromosome: dict[str, list[tuple[str, int]]] = {}
    for row in candidates[["cpg_id", "chromosome", "position"]].itertuples(index=False):
        candidates_by_chromosome.setdefault(str(row.chromosome), []).append((str(row.cpg_id), int(row.position)))

    records: list[dict[str, object]] = []
    for path, name in zip(re2g_files, names):
        records.extend(
            _stream_re2g_overlaps(
                path,
                biosample=name,
                candidates_by_chromosome=candidates_by_chromosome,
            )
        )
    linked_ids = {str(record["cpg_id"]) for record in records}
    for cpg_id in candidates["cpg_id"].astype(str):
        if cpg_id not in linked_ids:
            records.append({"cpg_id": cpg_id, "re2g_linked": False})
    result = pd.DataFrame(records)
    if result.empty:
        return pd.DataFrame(columns=["cpg_id", "re2g_linked"])
    dedupe_columns = [
        column
        for column in [
            "cpg_id",
            "re2g_linked",
            "re2g_target_gene",
            "re2g_biosample",
            "re2g_score",
            "re2g_element",
        ]
        if column in result.columns
    ]
    return result.drop_duplicates(subset=dedupe_columns).reset_index(drop=True)


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1]


def enumerate_spcas9_guides(
    sequence: str,
    *,
    sequence_start_1based: int,
    cpg_position_1based: int,
    max_distance: int = 150,
) -> pd.DataFrame:
    seq = sequence.upper()
    records: list[dict[str, object]] = []
    # Plus-strand guide: 20 nt immediately upstream of NGG.
    for pam_start in range(20, len(seq) - 2):
        pam = seq[pam_start : pam_start + 3]
        if len(pam) == 3 and pam[1:] == "GG":
            protospacer = seq[pam_start - 20 : pam_start]
            if re.fullmatch(r"[ACGT]{20}", protospacer):
                guide_center = sequence_start_1based + pam_start - 10
                records.append(
                    {
                        "strand": "+",
                        "protospacer": protospacer,
                        "pam": pam,
                        "guide_center": guide_center,
                        "distance_to_cpg": abs(guide_center - cpg_position_1based),
                    }
                )
    # Minus-strand PAM is CCN on the plus-strand sequence.
    for pam_start in range(0, len(seq) - 23):
        pam_plus = seq[pam_start : pam_start + 3]
        if len(pam_plus) == 3 and pam_plus[:2] == "CC":
            protospacer_plus = seq[pam_start + 3 : pam_start + 23]
            if re.fullmatch(r"[ACGT]{20}", protospacer_plus):
                protospacer = reverse_complement(protospacer_plus)
                pam = reverse_complement(pam_plus)
                guide_center = sequence_start_1based + pam_start + 13
                records.append(
                    {
                        "strand": "-",
                        "protospacer": protospacer,
                        "pam": pam,
                        "guide_center": guide_center,
                        "distance_to_cpg": abs(guide_center - cpg_position_1based),
                    }
                )
    frame = pd.DataFrame(records)
    if frame.empty:
        return frame
    frame["gc_fraction"] = frame["protospacer"].map(lambda s: (s.count("G") + s.count("C")) / 20)
    frame["homopolymer_4"] = frame["protospacer"].str.contains(r"A{4}|C{4}|G{4}|T{4}", regex=True)
    frame["poly_t_4"] = frame["protospacer"].str.contains("TTTT", regex=False)
    frame["basic_sequence_pass"] = (
        frame["gc_fraction"].between(0.35, 0.70)
        & ~frame["homopolymer_4"]
        & ~frame["poly_t_4"]
        & (frame["distance_to_cpg"] <= max_distance)
    )
    return frame.sort_values(["basic_sequence_pass", "distance_to_cpg"], ascending=[False, True]).reset_index(drop=True)


def fetch_ucsc_sequence(chromosome: str, start0: int, end0: int, *, genome: str = "hg38") -> str:
    query = urllib.parse.urlencode(
        {"genome": genome, "chrom": chromosome, "start": int(start0), "end": int(end0)}
    )
    url = f"https://api.genome.ucsc.edu/getData/sequence?{query}"
    with urllib.request.urlopen(url, timeout=60) as response:
        payload = json.load(response)
    dna = str(payload.get("dna", "")).upper()
    if not dna:
        raise ValueError(f"No sequence returned for {chromosome}:{start0}-{end0}")
    return dna


def design_candidate_guides(
    candidate_file: str | Path,
    frozen_shortlist_file: str | Path,
    *,
    flank: int = 250,
    max_distance: int = 150,
    genome: str = "hg38",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = pd.read_csv(candidate_file, compression="infer", low_memory=False)
    frozen = load_candidate_directions(candidate_file, frozen_shortlist_file=frozen_shortlist_file)
    cpg_col = _find_column(candidates, ["cpg_id"])
    chr_col = _find_column(candidates, ["chromosome", "chrom", "chr"])
    pos_col = _find_column(candidates, ["position", "coordinate", "pos", "mapinfo"])
    candidates = candidates.rename(columns={cpg_col: "cpg_id", chr_col: "chromosome", pos_col: "position"})
    candidates = candidates[candidates["cpg_id"].astype(str).isin(frozen["cpg_id"])].copy()
    candidates["position"] = pd.to_numeric(candidates["position"], errors="coerce")

    all_guides: list[pd.DataFrame] = []
    summary_rows: list[dict[str, object]] = []
    for row in candidates.itertuples(index=False):
        chrom = str(row.chromosome)
        if not chrom.lower().startswith("chr"):
            chrom = f"chr{chrom}"
        position = int(row.position)
        start0 = max(0, position - 1 - flank)
        end0 = position - 1 + flank + 1
        error = ""
        try:
            sequence = fetch_ucsc_sequence(chrom, start0, end0, genome=genome)
            guides = enumerate_spcas9_guides(
                sequence,
                sequence_start_1based=start0 + 1,
                cpg_position_1based=position,
                max_distance=max_distance,
            )
        except Exception as exc:  # Network and coordinate failures are recorded per candidate.
            guides = pd.DataFrame()
            error = f"{type(exc).__name__}: {exc}"
        if not guides.empty:
            guides.insert(0, "cpg_id", row.cpg_id)
            guides.insert(1, "chromosome", chrom)
            guides.insert(2, "cpg_position", position)
            all_guides.append(guides)
        summary_rows.append(
            {
                "cpg_id": row.cpg_id,
                "chromosome": chrom,
                "position": position,
                "n_spcas9_ngg_guides": int(len(guides)),
                "n_basic_sequence_pass": int(guides["basic_sequence_pass"].sum()) if len(guides) else 0,
                "nearest_guide_distance": float(guides["distance_to_cpg"].min()) if len(guides) else np.nan,
                "guide_feasible": bool(guides["basic_sequence_pass"].any()) if len(guides) else False,
                "sequence_fetch_error": error,
            }
        )
    guide_frame = pd.concat(all_guides, ignore_index=True) if all_guides else pd.DataFrame()
    return guide_frame, pd.DataFrame(summary_rows)


def integrate_orthogonal_evidence(
    frozen_shortlist_file: str | Path,
    external_reprogramming_file: str | Path,
    natural_aging_file: str | Path,
    *,
    eqtm_file: str | Path | None = None,
    re2g_file: str | Path | None = None,
    guide_summary_file: str | Path | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    frozen = pd.read_csv(frozen_shortlist_file)
    external = pd.read_csv(external_reprogramming_file)
    aging = pd.read_csv(natural_aging_file)
    out = frozen.merge(external, on="cpg_id", how="left", suffixes=("", "_external"))
    out = out.merge(aging, on="cpg_id", how="left", suffixes=("", "_aging"))

    if eqtm_file is not None and Path(eqtm_file).exists():
        eqtm = pd.read_csv(eqtm_file)
        eqtm_summary = (
            eqtm.groupby("cpg_id", as_index=False)
            .agg(
                eqtm_min_fdr=("eqtm_fdr", "min"),
                eqtm_max_abs_rho=("eqtm_spearman_rho", lambda s: pd.to_numeric(s, errors="coerce").abs().max()),
                eqtm_supported_genes=("gene", lambda s: ";".join(sorted(set(map(str, s))))),
            )
        )
        out = out.merge(eqtm_summary, on="cpg_id", how="left")
    if re2g_file is not None and Path(re2g_file).exists():
        re2g = pd.read_csv(re2g_file)
        re2g_summary = (
            re2g.groupby("cpg_id", as_index=False)
            .agg(
                re2g_linked=("re2g_linked", "max"),
                re2g_target_genes=("re2g_target_gene", lambda s: ";".join(sorted(set(str(x) for x in s.dropna())))),
                re2g_biosamples=("re2g_biosample", lambda s: ";".join(sorted(set(str(x) for x in s.dropna())))),
            )
        )
        out = out.merge(re2g_summary, on="cpg_id", how="left")
    if guide_summary_file is not None and Path(guide_summary_file).exists():
        guides = pd.read_csv(guide_summary_file)
        out = out.merge(guides, on="cpg_id", how="left")

    bool_cols = [
        "external_mean_direction_concordant",
        "replicative_reversal_concordant",
        "chronological_reversal_concordant",
        "multidonor_reversal_concordant",
        "re2g_linked",
        "guide_feasible",
    ]
    for col in bool_cols:
        if col not in out.columns:
            out[col] = False
        out[col] = out[col].fillna(False).astype(bool)
    out["eqtm_supported"] = out.get("eqtm_min_fdr", pd.Series(np.nan, index=out.index)).lt(0.10)
    out["orthogonal_support_count"] = out[
        [
            "external_mean_direction_concordant",
            "replicative_reversal_concordant",
            "chronological_reversal_concordant",
            "multidonor_reversal_concordant",
            "eqtm_supported",
            "re2g_linked",
            "guide_feasible",
        ]
    ].sum(axis=1)

    def tier(row: pd.Series) -> str:
        panel_role = str(row.get("panel_role", "")).lower()
        if "sentinel" in panel_role:
            return "S - safety sentinel"
        if row["external_mean_direction_concordant"] and row["orthogonal_support_count"] >= 4:
            return "V1 - multiply supported experimental candidate"
        if row["external_mean_direction_concordant"] and row["orthogonal_support_count"] >= 2:
            return "V2 - directionally replicated candidate"
        if row["orthogonal_support_count"] >= 2:
            return "V3 - mechanistically supported but direction-incomplete"
        return "V4 - insufficient orthogonal support"

    out["orthogonal_validation_tier"] = out.apply(tier, axis=1)
    rank_sort = "robust_rank" if "robust_rank" in out.columns else "biological_shortlist_rank"
    out = out.sort_values(["orthogonal_validation_tier", rank_sort], kind="stable")
    summary = {
        "n_candidates": int(len(out)),
        "tier_counts": out["orthogonal_validation_tier"].value_counts().to_dict(),
        "principle": "The frozen panel is evaluated without refitting the original evidence weights or ranks.",
        "scope": "Orthogonal in-silico validation does not demonstrate locus-level causal rejuvenation.",
    }
    return out.reset_index(drop=True), summary
