from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import json
import math
import re

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


DEFAULT_WEIGHTS: dict[str, float] = {
    "discovery_strength_score": 0.10,
    "donor_stability_score": 0.18,
    "o3_support_score": 0.14,
    "successful_failed_specificity_score": 0.12,
    "temporal_consistency_score": 0.12,
    "editor_responsiveness_score": 0.10,
    "regional_support_score": 0.08,
    "genomic_context_score": 0.06,
    "biological_relevance_score": 0.10,
}

RISK_CATEGORIES = {
    "oncogenic_risk",
    "pluripotency_reprogramming_risk",
    "tumor_suppressor_risk",
    "developmental_identity_risk",
}


@dataclass(frozen=True)
class GenomicCausalityResult:
    cpg_scores: pd.DataFrame
    region_scores: pd.DataFrame
    gene_scores: pd.DataFrame
    category_enrichment: pd.DataFrame
    safety_flags: pd.DataFrame
    evidence_component_summary: pd.DataFrame
    shortlist: pd.DataFrame
    summary: dict[str, object]
    manifest: dict[str, object]


def _jsonable(value: object) -> object:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def _split_tokens(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    return sorted({x.strip().upper() for x in re.split(r"[;,|]", str(value)) if x.strip()})


def _first_existing(frame: pd.DataFrame, names: Sequence[str]) -> str | None:
    lower = {str(c).lower(): c for c in frame.columns}
    for name in names:
        if name.lower() in lower:
            return str(lower[name.lower()])
    return None


def _normalise_cpg(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    col = _first_existing(out, ["cpg_id", "probeid", "probe_id", "id_ref", "ilmnid"])
    if col is None:
        raise ValueError("Input table does not contain a CpG/probe identifier column")
    if col != "cpg_id":
        out = out.rename(columns={col: "cpg_id"})
    out["cpg_id"] = out["cpg_id"].astype(str).str.replace(r"_(BC|TC)\\d+$", "", regex=True)
    return out.drop_duplicates("cpg_id", keep="first")


def _percentile(values: pd.Series, *, higher_is_better: bool = True) -> pd.Series:
    x = pd.to_numeric(values, errors="coerce")
    if x.notna().sum() <= 1:
        return pd.Series(np.where(x.notna(), 0.5, np.nan), index=x.index, dtype=float)
    ranks = x.rank(method="average", pct=True, ascending=higher_is_better)
    # pandas ascending=True gives largest rank to the largest value; this is desired.
    return ranks.astype(float)


def _weighted_available_score(frame: pd.DataFrame, weights: Mapping[str, float]) -> tuple[pd.Series, pd.Series]:
    score = np.zeros(len(frame), dtype=float)
    denom = np.zeros(len(frame), dtype=float)
    available_count = np.zeros(len(frame), dtype=float)
    for col, weight in weights.items():
        if col not in frame.columns:
            continue
        values = pd.to_numeric(frame[col], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(values)
        score[valid] += weight * np.clip(values[valid], 0.0, 1.0)
        denom[valid] += weight
        available_count[valid] += 1.0
    out = np.divide(score, denom, out=np.full(len(frame), np.nan), where=denom > 0)
    completeness = available_count / max(len(weights), 1)
    return pd.Series(out, index=frame.index), pd.Series(completeness, index=frame.index)


def load_gene_evidence(
    seed_path: str | Path,
    external_paths: Sequence[str | Path] = (),
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for source, path_like in [("seed", seed_path), *[("external", x) for x in external_paths]]:
        path = Path(path_like)
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path)
        required = {"gene", "category"}
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"Gene evidence file {path} is missing {sorted(missing)}")
        frame = frame.copy()
        frame["gene"] = frame["gene"].astype(str).str.upper().str.strip()
        frame["category"] = frame["category"].astype(str).str.strip()
        if "relevance_weight" not in frame.columns:
            frame["relevance_weight"] = 0.0
        if "risk_weight" not in frame.columns:
            frame["risk_weight"] = 0.0
        if "evidence_level" not in frame.columns:
            frame["evidence_level"] = "unspecified"
        if "source_id" not in frame.columns:
            frame["source_id"] = path.name
        if "source_note" not in frame.columns:
            frame["source_note"] = ""
        frame["evidence_source_class"] = source
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)
    combined["relevance_weight"] = pd.to_numeric(combined["relevance_weight"], errors="coerce").fillna(0).clip(0, 1)
    combined["risk_weight"] = pd.to_numeric(combined["risk_weight"], errors="coerce").fillna(0).clip(0, 1)
    combined = combined[combined["gene"].ne("")].copy()
    # Preserve source provenance while using the strongest evidence for scoring.
    return combined.drop_duplicates(["gene", "category", "source_id"], keep="first")


def _merge_optional_cpg_evidence(
    cpg: pd.DataFrame,
    successful_vs_failed: pd.DataFrame | None,
    time_series: pd.DataFrame | None,
) -> pd.DataFrame:
    out = cpg.copy()
    if successful_vs_failed is not None and not successful_vs_failed.empty:
        svf = _normalise_cpg(successful_vs_failed)
        discovery_col = _first_existing(svf, ["discovery_delta_beta", "target_delta_beta"])
        validation_col = _first_existing(svf, ["validation_delta_beta", "successful_vs_failed_delta", "delta_beta"])
        direction_col = _first_existing(svf, ["same_direction", "centered_same_direction"])
        retention_col = _first_existing(svf, ["effect_retention"])
        if validation_col is not None:
            d = pd.to_numeric(svf[discovery_col], errors="coerce") if discovery_col else pd.Series(np.nan, index=svf.index)
            v = pd.to_numeric(svf[validation_col], errors="coerce")
            if direction_col:
                same = svf[direction_col].fillna(False).astype(bool).astype(float)
            else:
                same = (np.sign(d) == np.sign(v)).astype(float)
            if retention_col:
                retention = pd.to_numeric(svf[retention_col], errors="coerce")
            else:
                retention = v.abs() / d.abs().clip(lower=1e-12)
            svf_small = pd.DataFrame(
                {
                    "cpg_id": svf["cpg_id"],
                    "successful_failed_delta_beta": v,
                    "successful_failed_same_direction": same,
                    "successful_failed_effect_retention": retention,
                    "successful_failed_specificity_score": 0.60 * same + 0.40 * (retention.clip(0, 1.5) / 1.5),
                }
            ).groupby("cpg_id", as_index=False).mean(numeric_only=True)
            out = out.merge(svf_small, on="cpg_id", how="left")
    if "successful_failed_specificity_score" not in out.columns:
        out["successful_failed_specificity_score"] = np.nan

    if time_series is not None and not time_series.empty:
        ts = _normalise_cpg(time_series)
        success_col = _first_existing(ts, ["success_sign_consistency"])
        failed_col = _first_existing(ts, ["failed_sign_consistency"])
        o3_col = _first_existing(ts, ["o3_success_sign_consistency"])
        earliest_col = _first_existing(ts, ["earliest_o3_supported_day"])
        if success_col:
            success = pd.to_numeric(ts[success_col], errors="coerce")
            failed = pd.to_numeric(ts[failed_col], errors="coerce") if failed_col else pd.Series(np.nan, index=ts.index)
            o3 = pd.to_numeric(ts[o3_col], errors="coerce") if o3_col else success
            contrast = (success - failed).clip(lower=0) if failed_col else pd.Series(np.nan, index=ts.index)
            temporal = 0.55 * success + 0.30 * o3
            temporal += np.where(contrast.notna(), 0.15 * contrast, 0.0)
            ts_small = pd.DataFrame(
                {
                    "cpg_id": ts["cpg_id"],
                    "success_sign_consistency": success,
                    "failed_sign_consistency": failed,
                    "o3_success_sign_consistency": o3,
                    "temporal_consistency_score": temporal.clip(0, 1),
                }
            )
            if earliest_col:
                ts_small["earliest_o3_supported_day"] = pd.to_numeric(ts[earliest_col], errors="coerce")
            out = out.merge(ts_small.drop_duplicates("cpg_id"), on="cpg_id", how="left")
    if "temporal_consistency_score" not in out.columns:
        out["temporal_consistency_score"] = np.nan
    return out


def _annotate_gene_evidence(cpg: pd.DataFrame, evidence: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    for row in cpg[["cpg_id", "gene"]].itertuples(index=False):
        genes = _split_tokens(row.gene)
        for gene in genes:
            rows.append({"cpg_id": row.cpg_id, "gene_symbol": gene})
    map_frame = pd.DataFrame(rows, columns=["cpg_id", "gene_symbol"])
    if map_frame.empty:
        out = cpg.copy()
        out["annotated_gene_count"] = 0
        out["biological_relevance_score"] = 0.0
        out["safety_risk_score"] = 0.0
        out["biological_categories"] = ""
        out["risk_categories"] = ""
        return out, pd.DataFrame()

    joined = map_frame.merge(evidence, left_on="gene_symbol", right_on="gene", how="left")
    joined["relevance_weight"] = pd.to_numeric(joined["relevance_weight"], errors="coerce").fillna(0.0)
    joined["risk_weight"] = pd.to_numeric(joined["risk_weight"], errors="coerce").fillna(0.0)
    joined["is_risk_category"] = joined["category"].isin(RISK_CATEGORIES) | (joined["risk_weight"] > 0)

    summaries: list[dict[str, object]] = []
    for cpg_id, group in joined.groupby("cpg_id", sort=False):
        positive = group[group["relevance_weight"] > 0]
        risk = group[group["risk_weight"] > 0]
        rel_values = positive["relevance_weight"].to_numpy(dtype=float)
        relevance = float(0.60 * rel_values.max() + 0.40 * rel_values.mean()) if len(rel_values) else 0.0
        risk_score = float(risk["risk_weight"].max()) if len(risk) else 0.0
        summaries.append(
            {
                "cpg_id": cpg_id,
                "annotated_gene_count": int(group["gene_symbol"].nunique()),
                "biological_relevance_score": relevance,
                "safety_risk_score": risk_score,
                "biological_categories": ";".join(sorted(positive["category"].dropna().astype(str).unique())),
                "risk_categories": ";".join(sorted(risk["category"].dropna().astype(str).unique())),
                "evidence_genes": ";".join(sorted(positive["gene_symbol"].dropna().astype(str).unique())),
                "risk_genes": ";".join(sorted(risk["gene_symbol"].dropna().astype(str).unique())),
            }
        )
    out = cpg.merge(pd.DataFrame(summaries), on="cpg_id", how="left")
    for col, default in [
        ("annotated_gene_count", 0),
        ("biological_relevance_score", 0.0),
        ("safety_risk_score", 0.0),
        ("biological_categories", ""),
        ("risk_categories", ""),
        ("evidence_genes", ""),
        ("risk_genes", ""),
    ]:
        out[col] = out[col].fillna(default)
    return out, joined


def _genomic_context_scores(frame: pd.DataFrame, region_scores: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    group_text = out.get("gene_group", pd.Series("", index=out.index)).fillna("").astype(str).str.upper()
    island_text = out.get("island_relation", pd.Series("", index=out.index)).fillna("").astype(str).str.upper()
    promoter = group_text.str.contains(r"TSS200|TSS1500|5'UTR|5UTR|1STEXON", regex=True).astype(float)
    island = np.select(
        [
            island_text.str.contains("ISLAND", regex=False),
            island_text.str.contains("SHORE", regex=False),
            island_text.str.contains("SHELF", regex=False),
            island_text.str.contains("OPEN", regex=False),
        ],
        [1.0, 0.8, 0.5, 0.2],
        default=0.0,
    )
    region = region_scores.copy()
    if not region.empty:
        region["region_responsiveness_percentile"] = _percentile(region["regional_responsiveness_score"])
        reg_cols = [
            "unique_region_id",
            "n_candidate_cpgs",
            "region_class",
            "regional_responsiveness_score",
            "region_responsiveness_percentile",
        ]
        out = out.merge(region[reg_cols], on="unique_region_id", how="left", suffixes=("", "_region"))
    else:
        out["n_candidate_cpgs"] = 1
        out["region_class"] = "singleton_locus"
        out["regional_responsiveness_score"] = np.nan
        out["region_responsiveness_percentile"] = np.nan
    multi = (pd.to_numeric(out["n_candidate_cpgs"], errors="coerce").fillna(1).clip(upper=3) - 1) / 2
    out["promoter_proximal"] = promoter.astype(bool)
    out["cpg_context_score"] = island.astype(float)
    out["genomic_context_score"] = 0.40 * promoter + 0.25 * island + 0.35 * multi
    out["regional_support_score"] = (
        0.55 * pd.to_numeric(out["region_responsiveness_percentile"], errors="coerce").fillna(0.0)
        + 0.45 * multi
    )
    return out


def score_candidate_cpgs(
    membership: pd.DataFrame,
    region_scores: pd.DataFrame,
    gene_evidence: pd.DataFrame,
    *,
    successful_vs_failed: pd.DataFrame | None = None,
    time_series: pd.DataFrame | None = None,
    weights: Mapping[str, float] = DEFAULT_WEIGHTS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cpg = _normalise_cpg(membership)
    if "gene" not in cpg.columns:
        cpg["gene"] = ""
    if "target_delta_beta" not in cpg.columns:
        cpg["target_delta_beta"] = np.nan
    if "n_lodo_rotations" not in cpg.columns:
        cpg["n_lodo_rotations"] = np.nan

    abs_target = pd.to_numeric(cpg["target_delta_beta"], errors="coerce").abs()
    effect_rank = _percentile(abs_target)
    if "rank_score" in cpg.columns:
        rank_score = _percentile(pd.to_numeric(cpg["rank_score"], errors="coerce"))
        cpg["discovery_strength_score"] = pd.concat([effect_rank, rank_score], axis=1).mean(axis=1)
    else:
        cpg["discovery_strength_score"] = effect_rank

    inferred_lodo = pd.Series(np.nan, index=cpg.index, dtype=float)
    if "in_stable_core_3of3" in cpg.columns:
        inferred_lodo = inferred_lodo.mask(cpg["in_stable_core_3of3"].fillna(False), 3.0)
    if "in_consensus_at_least_2of3" in cpg.columns:
        inferred_lodo = inferred_lodo.mask(
            inferred_lodo.isna() & cpg["in_consensus_at_least_2of3"].fillna(False), 2.0
        )
    lodo = pd.to_numeric(cpg["n_lodo_rotations"], errors="coerce").fillna(inferred_lodo)
    cpg["donor_stability_score"] = (lodo / 3.0).clip(0, 1)
    cpg["n_lodo_rotations"] = lodo

    tier1 = cpg.get("in_o3_tier1", pd.Series(False, index=cpg.index)).fillna(False).astype(float)
    frozen = cpg.get("in_frozen_500", pd.Series(True, index=cpg.index)).fillna(False).astype(float)
    cpg = _merge_optional_cpg_evidence(cpg, successful_vs_failed, time_series)
    specificity = pd.to_numeric(cpg["successful_failed_specificity_score"], errors="coerce")
    cpg["o3_support_score"] = pd.concat(
        [tier1.rename("tier1"), (0.45 * frozen + 0.55 * specificity).rename("specificity")], axis=1
    ).max(axis=1, skipna=True)

    for cell in ["tcell", "msc"]:
        delta_col = f"{cell}_delta"
        consistency_col = f"{cell}_consistency"
        if delta_col not in cpg.columns:
            cpg[delta_col] = np.nan
        if consistency_col not in cpg.columns:
            cpg[consistency_col] = np.nan
        magnitude = _percentile(pd.to_numeric(cpg[delta_col], errors="coerce").abs())
        consistency = pd.to_numeric(cpg[consistency_col], errors="coerce").clip(0, 1)
        cpg[f"{cell}_responsiveness_score"] = 0.40 * magnitude + 0.60 * consistency
    cpg["editor_responsiveness_score"] = cpg[
        ["tcell_responsiveness_score", "msc_responsiveness_score"]
    ].max(axis=1, skipna=True)
    cpg["cell_context_divergence"] = (
        cpg["tcell_responsiveness_score"] - cpg["msc_responsiveness_score"]
    ).abs()

    cpg = _genomic_context_scores(cpg, region_scores)
    cpg, gene_map = _annotate_gene_evidence(cpg, gene_evidence)

    raw, completeness = _weighted_available_score(cpg, weights)
    cpg["causal_evidence_score"] = 100.0 * raw
    cpg["evidence_completeness"] = completeness
    cpg["safety_adjusted_priority_score"] = (
        cpg["causal_evidence_score"]
        * (1.0 - 0.45 * cpg["safety_risk_score"].clip(0, 1))
        * (0.75 + 0.25 * cpg["evidence_completeness"])
    )
    cpg["causal_evidence_grade"] = np.select(
        [
            (cpg["causal_evidence_score"] >= 75) & (cpg["evidence_completeness"] >= 0.70),
            (cpg["causal_evidence_score"] >= 60) & (cpg["evidence_completeness"] >= 0.55),
            cpg["causal_evidence_score"] >= 45,
        ],
        ["E1_high_priority", "E2_supported", "E3_exploratory"],
        default="E4_insufficient",
    )
    cpg["experimental_priority_tier"] = np.select(
        [
            (cpg["safety_adjusted_priority_score"] >= 65) & (cpg["safety_risk_score"] < 0.50),
            (cpg["safety_adjusted_priority_score"] >= 52) & (cpg["safety_risk_score"] < 0.70),
            cpg["safety_adjusted_priority_score"] >= 38,
        ],
        ["P1", "P2", "P3"],
        default="P4",
    )
    target = pd.to_numeric(cpg["target_delta_beta"], errors="coerce")
    cpg["directional_editor_hypothesis"] = np.select(
        [target < 0, target > 0],
        ["TET_like_hypomethylating_editor", "DNMT3A_like_hypermethylating_editor"],
        default="undetermined",
    )
    cpg["interpretation"] = "triangulated_prioritisation_not_causal_proof"
    cpg = cpg.sort_values(
        ["safety_adjusted_priority_score", "causal_evidence_score"], ascending=False
    ).reset_index(drop=True)
    cpg["global_priority_rank"] = np.arange(1, len(cpg) + 1)
    return cpg, gene_map


def score_regions(cpg_scores: pd.DataFrame, region_scores: pd.DataFrame) -> pd.DataFrame:
    if region_scores.empty:
        return pd.DataFrame()
    rows: list[dict[str, object]] = []
    for region_id, group in cpg_scores.groupby("unique_region_id", dropna=True, sort=True):
        base = region_scores[region_scores["unique_region_id"] == region_id]
        region_response = float(base["regional_responsiveness_score"].iloc[0]) if len(base) else np.nan
        row: dict[str, object] = {
            "unique_region_id": region_id,
            "chromosome": str(group["chromosome"].iloc[0]) if "chromosome" in group.columns else "",
            "start": int(pd.to_numeric(group["position"], errors="coerce").min()),
            "end": int(pd.to_numeric(group["position"], errors="coerce").max()),
            "n_candidate_cpgs": int(group["cpg_id"].nunique()),
            "candidate_cpg_ids": ";".join(sorted(group["cpg_id"].astype(str).unique())),
            "genes": ";".join(sorted({g for value in group["gene"] for g in _split_tokens(value)})),
            "median_causal_evidence_score": float(group["causal_evidence_score"].median()),
            "max_causal_evidence_score": float(group["causal_evidence_score"].max()),
            "median_safety_adjusted_priority": float(group["safety_adjusted_priority_score"].median()),
            "max_safety_risk_score": float(group["safety_risk_score"].max()),
            "mean_biological_relevance_score": float(group["biological_relevance_score"].mean()),
            "mean_evidence_completeness": float(group["evidence_completeness"].mean()),
            "n_p1_cpgs": int((group["experimental_priority_tier"] == "P1").sum()),
            "regional_responsiveness_score": region_response,
            "region_class": "multi_cpg_region" if group["cpg_id"].nunique() >= 2 else "singleton_locus",
        }
        for col in [
            "in_stable_core_3of3",
            "in_consensus_at_least_2of3",
            "in_o3_tier1",
            "in_frozen_500",
        ]:
            if col in group.columns:
                row[col] = bool(group[col].fillna(False).any())
        rows.append(row)
    out = pd.DataFrame(rows)
    response_pct = _percentile(out["regional_responsiveness_score"])
    multi_bonus = (out["n_candidate_cpgs"].clip(upper=4) - 1) / 3
    out["region_causal_priority_score"] = (
        0.55 * out["median_safety_adjusted_priority"]
        + 0.25 * out["max_causal_evidence_score"]
        + 20.0 * (0.12 * response_pct.fillna(0) + 0.08 * multi_bonus)
    ) * (1.0 - 0.30 * out["max_safety_risk_score"].clip(0, 1))
    out["region_priority_tier"] = np.select(
        [
            (out["region_causal_priority_score"] >= 65) & (out["max_safety_risk_score"] < 0.5),
            (out["region_causal_priority_score"] >= 50) & (out["max_safety_risk_score"] < 0.7),
            out["region_causal_priority_score"] >= 35,
        ],
        ["R1", "R2", "R3"],
        default="R4",
    )
    out["interpretation"] = "region_priority_not_direct_editability_or_causal_proof"
    return out.sort_values("region_causal_priority_score", ascending=False).reset_index(drop=True)


def score_genes(cpg_scores: pd.DataFrame, gene_evidence: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for row in cpg_scores.itertuples(index=False):
        for gene in _split_tokens(getattr(row, "gene", "")):
            rows.append(
                {
                    "gene": gene,
                    "cpg_id": row.cpg_id,
                    "causal_evidence_score": row.causal_evidence_score,
                    "safety_adjusted_priority_score": row.safety_adjusted_priority_score,
                    "safety_risk_score": row.safety_risk_score,
                    "experimental_priority_tier": row.experimental_priority_tier,
                    "biological_relevance_score": row.biological_relevance_score,
                }
            )
    if not rows:
        return pd.DataFrame()
    long = pd.DataFrame(rows)
    evidence_categories = (
        gene_evidence.groupby("gene")["category"]
        .agg(lambda x: ";".join(sorted(set(x.astype(str)))))
        .rename("evidence_categories")
    )
    out = (
        long.groupby("gene", as_index=False)
        .agg(
            n_candidate_cpgs=("cpg_id", "nunique"),
            median_causal_evidence_score=("causal_evidence_score", "median"),
            max_causal_evidence_score=("causal_evidence_score", "max"),
            max_safety_adjusted_priority=("safety_adjusted_priority_score", "max"),
            max_safety_risk_score=("safety_risk_score", "max"),
            mean_biological_relevance_score=("biological_relevance_score", "mean"),
            n_p1_cpgs=("experimental_priority_tier", lambda x: int((x == "P1").sum())),
        )
        .merge(evidence_categories, left_on="gene", right_index=True, how="left")
    )
    return out.sort_values(["max_safety_adjusted_priority", "max_causal_evidence_score"], ascending=False)


def read_gmt(paths: Sequence[str | Path]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for path_like in paths:
        path = Path(path_like)
        if not path.exists():
            raise FileNotFoundError(path)
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                set_name, description, *genes = parts
                for gene in genes:
                    if gene.strip():
                        rows.append(
                            {
                                "category": set_name,
                                "gene": gene.strip().upper(),
                                "source": path.name,
                                "description": description,
                            }
                        )
    return pd.DataFrame(rows)


def category_enrichment(
    cpg_scores: pd.DataFrame,
    annotation: pd.DataFrame,
    gene_evidence: pd.DataFrame,
    *,
    gmt: pd.DataFrame | None = None,
) -> pd.DataFrame:
    universe = {g for value in annotation.get("gene", pd.Series(dtype=str)).fillna("") for g in _split_tokens(value)}
    category_rows = gene_evidence[["category", "gene"]].drop_duplicates().copy()
    category_rows["source"] = "gene_evidence_panel"
    if gmt is not None and not gmt.empty:
        category_rows = pd.concat(
            [category_rows, gmt[["category", "gene", "source"]].drop_duplicates()], ignore_index=True
        )
    candidate_all = {g for value in cpg_scores["gene"].fillna("") for g in _split_tokens(value)} & universe
    candidate_high = {
        g
        for value in cpg_scores.loc[
            cpg_scores["experimental_priority_tier"].isin(["P1", "P2"]), "gene"
        ].fillna("")
        for g in _split_tokens(value)
    } & universe
    rows: list[dict[str, object]] = []
    M = len(universe)
    if M == 0:
        return pd.DataFrame()
    for scope, genes in [("all_candidates", candidate_all), ("priority_P1_P2", candidate_high)]:
        N = len(genes)
        if N == 0:
            continue
        for (category, source), group in category_rows.groupby(["category", "source"], sort=True):
            set_genes = set(group["gene"].astype(str).str.upper()) & universe
            n = len(set_genes)
            if n < 2:
                continue
            overlap = genes & set_genes
            k = len(overlap)
            p = float(stats.hypergeom.sf(k - 1, M, n, N)) if k else 1.0
            rows.append(
                {
                    "scope": scope,
                    "category": category,
                    "source": source,
                    "universe_gene_count": M,
                    "candidate_gene_count": N,
                    "category_gene_count": n,
                    "overlap_count": k,
                    "overlap_genes": ";".join(sorted(overlap)),
                    "fold_enrichment": (k / N) / (n / M) if N and n else np.nan,
                    "p_value": p,
                }
            )
    enrichment_columns = [
        "scope", "category", "source", "universe_gene_count", "candidate_gene_count",
        "category_gene_count", "overlap_count", "overlap_genes", "fold_enrichment",
        "p_value", "bh_q_global", "bh_q_within_scope",
    ]
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=enrichment_columns)
    out["bh_q_global"] = multipletests(out["p_value"].fillna(1), method="fdr_bh")[1]
    out["bh_q_within_scope"] = np.nan
    for scope, idx in out.groupby("scope").groups.items():
        out.loc[idx, "bh_q_within_scope"] = multipletests(
            out.loc[idx, "p_value"].fillna(1), method="fdr_bh"
        )[1]
    return out.sort_values(["bh_q_global", "fold_enrichment"], ascending=[True, False])


def build_safety_flags(cpg_scores: pd.DataFrame, gene_map: pd.DataFrame) -> pd.DataFrame:
    if gene_map.empty:
        return pd.DataFrame()
    risk = gene_map[(gene_map["risk_weight"] > 0) | gene_map["category"].isin(RISK_CATEGORIES)].copy()
    if risk.empty:
        return pd.DataFrame()
    cols = [
        "cpg_id",
        "gene_symbol",
        "category",
        "risk_weight",
        "evidence_level",
        "source_id",
        "source_note",
    ]
    risk = risk[[c for c in cols if c in risk.columns]].drop_duplicates()
    attach = cpg_scores[
        [
            "cpg_id",
            "unique_region_id",
            "target_delta_beta",
            "causal_evidence_score",
            "safety_adjusted_priority_score",
            "experimental_priority_tier",
        ]
    ]
    risk = risk.merge(attach, on="cpg_id", how="left")
    risk["flag_interpretation"] = "proximity_flag_requires_manual_experimental_safety_review"
    return risk.sort_values(["risk_weight", "safety_adjusted_priority_score"], ascending=False)


def build_component_summary(cpg_scores: pd.DataFrame, weights: Mapping[str, float]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for col, weight in weights.items():
        if col not in cpg_scores.columns:
            continue
        x = pd.to_numeric(cpg_scores[col], errors="coerce")
        rows.append(
            {
                "component": col,
                "weight": weight,
                "n_available": int(x.notna().sum()),
                "availability_fraction": float(x.notna().mean()),
                "median_score": float(x.median()) if x.notna().any() else np.nan,
                "mean_score": float(x.mean()) if x.notna().any() else np.nan,
                "q25": float(x.quantile(0.25)) if x.notna().any() else np.nan,
                "q75": float(x.quantile(0.75)) if x.notna().any() else np.nan,
            }
        )
    return pd.DataFrame(rows)


def run_genomic_causality_scoring(
    *,
    membership: pd.DataFrame,
    region_scores: pd.DataFrame,
    annotation: pd.DataFrame,
    gene_evidence: pd.DataFrame,
    successful_vs_failed: pd.DataFrame | None = None,
    time_series: pd.DataFrame | None = None,
    gmt: pd.DataFrame | None = None,
    weights: Mapping[str, float] = DEFAULT_WEIGHTS,
    shortlist_size: int = 100,
) -> GenomicCausalityResult:
    membership = _normalise_cpg(membership)
    annotation = _normalise_cpg(annotation)
    annotation_keep = [
        c for c in ["cpg_id", "chromosome", "position", "strand", "gene", "gene_group", "island_relation", "genome_build", "annotation_source"]
        if c in annotation.columns
    ]
    merge = membership.merge(annotation[annotation_keep], on="cpg_id", how="left", suffixes=("", "_annotation"))
    for col in ["chromosome", "position", "gene", "gene_group", "island_relation", "strand"]:
        alt = f"{col}_annotation"
        if alt in merge.columns:
            if col not in merge.columns:
                merge[col] = merge[alt]
            else:
                merge[col] = merge[col].where(merge[col].notna() & merge[col].astype(str).ne(""), merge[alt])
            merge = merge.drop(columns=alt)

    cpg_scores, gene_map = score_candidate_cpgs(
        merge,
        region_scores,
        gene_evidence,
        successful_vs_failed=successful_vs_failed,
        time_series=time_series,
        weights=weights,
    )
    regions = score_regions(cpg_scores, region_scores)
    genes = score_genes(cpg_scores, gene_evidence)
    enrichment = category_enrichment(cpg_scores, annotation, gene_evidence, gmt=gmt)
    safety = build_safety_flags(cpg_scores, gene_map)
    components = build_component_summary(cpg_scores, weights)
    shortlist = cpg_scores.head(max(int(shortlist_size), 1)).copy()

    summary = {
        "analysis_status": "genomic_annotation_biological_relevance_causality_evidence_scoring",
        "interpretation_scope": (
            "Scores triangulate discovery, donor stability, O3/temporal specificity, editor responsiveness, "
            "genomic context and heuristic biological evidence. They prioritise experiments and do not prove "
            "that editing a CpG causes rejuvenation or safety."
        ),
        "n_candidate_cpgs": int(len(cpg_scores)),
        "n_annotated_to_gene": int(cpg_scores["annotated_gene_count"].gt(0).sum()),
        "n_unique_regions": int(len(regions)),
        "n_multicpg_regions": int((regions["n_candidate_cpgs"] >= 2).sum()) if not regions.empty else 0,
        "n_priority_P1": int((cpg_scores["experimental_priority_tier"] == "P1").sum()),
        "n_priority_P2": int((cpg_scores["experimental_priority_tier"] == "P2").sum()),
        "n_safety_flagged_cpgs": int(safety["cpg_id"].nunique()) if not safety.empty else 0,
        "top_cpgs": cpg_scores.head(20)[
            [
                "cpg_id",
                "gene",
                "unique_region_id",
                "causal_evidence_score",
                "safety_adjusted_priority_score",
                "experimental_priority_tier",
                "directional_editor_hypothesis",
            ]
        ].to_dict("records"),
        "top_regions": regions.head(15)[
            [
                "unique_region_id",
                "genes",
                "n_candidate_cpgs",
                "region_causal_priority_score",
                "region_priority_tier",
            ]
        ].to_dict("records") if not regions.empty else [],
        "weights": dict(weights),
    }
    manifest = {
        "package_version": "0.8.0",
        "candidate_lists_frozen_before_v080": True,
        "causal_language_guardrail": "prioritisation_not_causal_proof",
        "gene_panel_guardrail": "bundled_seed_panel_is_heuristic_and_must_be_supplemented_by_traceable_literature_for_manuscript_claims",
        "safety_guardrail": "risk_flags_are_proximity_flags_not_clinical_safety_assessments",
        "score_range": "0_to_100",
        "missing_component_policy": "renormalise_over_available_prespecified_components_and_report_completeness",
        "weights": dict(weights),
        "editor_direction_policy": "derived_only_from_safe_reprogramming_target_delta_sign",
        "recommended_next_stage": "sequence_level_sgRNA_editor_design_and_off_target_safety_screening",
    }
    return GenomicCausalityResult(
        cpg_scores=cpg_scores,
        region_scores=regions,
        gene_scores=genes,
        category_enrichment=enrichment,
        safety_flags=safety,
        evidence_component_summary=components,
        shortlist=shortlist,
        summary=_jsonable(summary),
        manifest=_jsonable(manifest),
    )


def write_genomic_causality_outputs(result: GenomicCausalityResult, output_dir: str | Path) -> None:
    out = Path(output_dir)
    figures = out / "figures"
    out.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    result.cpg_scores.to_csv(out / "candidate_cpg_annotation_causality_scores.csv.gz", index=False, compression="gzip")
    result.region_scores.to_csv(out / "candidate_region_annotation_causality_scores.csv", index=False)
    result.gene_scores.to_csv(out / "gene_evidence_scores.csv", index=False)
    result.category_enrichment.to_csv(out / "biological_category_enrichment.csv", index=False)
    result.safety_flags.to_csv(out / "safety_flags.csv", index=False)
    result.evidence_component_summary.to_csv(out / "evidence_component_summary.csv", index=False)
    result.shortlist.to_csv(out / "top_candidate_shortlist.csv", index=False)
    (out / "genomic_causality_v080_summary.json").write_text(
        json.dumps(result.summary, indent=2), encoding="utf-8"
    )
    (out / "genomic_causality_v080_manifest.json").write_text(
        json.dumps(result.manifest, indent=2), encoding="utf-8"
    )
    _write_figures(result, figures)


def _write_figures(result: GenomicCausalityResult, figures: Path) -> None:
    import matplotlib
    # Figures are written to disk and never displayed, so force the
    # non-interactive backend. Without this, matplotlib may select an
    # interactive backend (e.g. TkAgg) and fail on machines with a
    # broken or absent Tcl/Tk, and in any headless environment.
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cpg = result.cpg_scores
    if not cpg.empty:
        plt.figure(figsize=(8, 5))
        plt.hist(cpg["causal_evidence_score"].dropna(), bins=25, alpha=0.75, label="Evidence")
        plt.hist(cpg["safety_adjusted_priority_score"].dropna(), bins=25, alpha=0.55, label="Safety-adjusted")
        plt.xlabel("Score")
        plt.ylabel("CpG count")
        plt.title("Candidate evidence and safety-adjusted priority")
        plt.legend()
        plt.tight_layout()
        plt.savefig(figures / "cpg_score_distribution.png", dpi=220)
        plt.close()

        components = list(DEFAULT_WEIGHTS)
        top = cpg.head(min(30, len(cpg))).set_index("cpg_id")
        matrix = top[[c for c in components if c in top.columns]].to_numpy(dtype=float)
        plt.figure(figsize=(11, max(5, len(top) * 0.24)))
        plt.imshow(matrix, aspect="auto", vmin=0, vmax=1)
        plt.colorbar(label="Component score")
        plt.yticks(np.arange(len(top)), top.index)
        plt.xticks(np.arange(matrix.shape[1]), [c.replace("_score", "") for c in components if c in top.columns], rotation=50, ha="right")
        plt.title("Top CpG evidence-component matrix")
        plt.tight_layout()
        plt.savefig(figures / "top_cpg_evidence_heatmap.png", dpi=220)
        plt.close()

        plt.figure(figsize=(7, 5))
        plt.scatter(cpg["safety_risk_score"], cpg["causal_evidence_score"], s=18, alpha=0.65)
        plt.xlabel("Safety proximity-risk score")
        plt.ylabel("Causal evidence score")
        plt.title("Evidence versus safety proximity risk")
        plt.tight_layout()
        plt.savefig(figures / "evidence_vs_safety_risk.png", dpi=220)
        plt.close()

    if not result.region_scores.empty:
        top_regions = result.region_scores.head(20)
        plt.figure(figsize=(9, 6))
        plt.barh(top_regions["unique_region_id"][::-1], top_regions["region_causal_priority_score"][::-1])
        plt.xlabel("Region causal-priority score")
        plt.title("Top candidate regions")
        plt.tight_layout()
        plt.savefig(figures / "top_candidate_regions.png", dpi=220)
        plt.close()

    if not result.category_enrichment.empty:
        enrich = result.category_enrichment[result.category_enrichment["scope"] == "priority_P1_P2"].head(20).copy()
        if len(enrich):
            values = -np.log10(enrich["bh_q_within_scope"].clip(lower=1e-12))
            plt.figure(figsize=(9, 6))
            plt.barh(enrich["category"][::-1], values[::-1])
            plt.xlabel("-log10(BH q within scope)")
            plt.title("Biological-category enrichment of priority candidates")
            plt.tight_layout()
            plt.savefig(figures / "biological_category_enrichment.png", dpi=220)
            plt.close()
