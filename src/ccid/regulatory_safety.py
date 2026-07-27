from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import gzip
import json
import math
import re
from urllib.parse import quote_plus

import numpy as np
import pandas as pd


DEFAULT_V081_CONFIG: dict[str, object] = {
    "experimental_weight": 0.72,
    "regulatory_weight": 0.18,
    "literature_relevance_weight": 0.10,
    "risk_penalty": 0.50,
    "uncertainty_penalty": 0.25,
    "nearest_tss_max_bp": 1_000_000,
    "linked_gene_max_bp": 100_000,
    "promoter_bp": 2_000,
    "proximal_bp": 10_000,
    "shortlist_size": 100,
}

REGULATORY_CLASS_SCORE: dict[str, float] = {
    "PLS": 1.00,
    "PROMOTER": 1.00,
    "PROMOTER-LIKE": 1.00,
    "PELS": 0.82,
    "PROXIMAL ENHANCER": 0.82,
    "DELS": 0.76,
    "DISTAL ENHANCER": 0.76,
    "CA-H3K4ME3": 0.70,
    "DNASE-H3K4ME3": 0.70,
    "CA-TF": 0.62,
    "DNASE-ONLY": 0.55,
    "CA": 0.55,
    "CTCF-ONLY": 0.42,
    "CA-CTCF": 0.42,
    "CTCF": 0.42,
    "TF-ONLY": 0.35,
}

RELATIONSHIP_WEIGHT: dict[str, float] = {
    "direct_epic_gene": 1.00,
    "promoter_nearest_tss": 0.95,
    "proximal_nearest_tss": 0.82,
    "distal_nearest_tss": 0.55,
    "enhancer_link": 0.90,
}

RISK_DOMAINS = {
    "oncogenic_risk",
    "tumor_suppressor_risk",
    "pluripotency_reprogramming_risk",
    "developmental_identity_risk",
    "genome_stability_risk",
    "essential_gene_risk",
}


@dataclass(frozen=True)
class RegulatorySafetyResult:
    cpg_scores: pd.DataFrame
    region_scores: pd.DataFrame
    gene_links: pd.DataFrame
    literature_audit: pd.DataFrame
    safety_audit: pd.DataFrame
    evidence_status_summary: pd.DataFrame
    manual_review_queue: pd.DataFrame
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


def _normalise_chromosome(value: object) -> str:
    text = str(value).strip()
    text = re.sub(r"^chr", "", text, flags=re.I)
    if text == "23":
        return "X"
    if text == "24":
        return "Y"
    if text.upper() in {"M", "MT", "MITO"}:
        return "MT"
    return text.upper()


def _split_genes(value: object) -> list[str]:
    if value is None:
        return []
    try:
        missing = pd.isna(value)
        if isinstance(missing, (bool, np.bool_)) and bool(missing):
            return []
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "<na>", "na", "null"}:
        return []
    return sorted({x.strip().upper() for x in re.split(r"[;,|]", text) if x.strip()})


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "reviewed"}


def _open_text(path: str | Path):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return p.open("rt", encoding="utf-8", errors="replace")


def _read_delimited(path: str | Path, *, header: bool = True) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    suffixes = "".join(p.suffixes).lower()
    if suffixes.endswith((".csv", ".csv.gz")):
        return pd.read_csv(p, compression="infer", low_memory=False, header=0 if header else None)
    return pd.read_csv(p, sep="\t", compression="infer", comment="#", low_memory=False, header=0 if header else None)


def parse_gtf_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in text.strip().strip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        match = re.match(r"([^\s]+)\s+[\"]?(.+?)[\"]?$", part)
        if match:
            attrs[match.group(1)] = match.group(2).strip().strip('"')
    return attrs


def load_tss_reference(path: str | Path) -> pd.DataFrame:
    """Load a GENCODE-style GTF or a precomputed TSS TSV/CSV.

    The returned coordinate is 1-based and compatible with EPIC annotation positions.
    """
    p = Path(path)
    name = p.name.lower()
    if name.endswith((".gtf", ".gtf.gz", ".gff", ".gff3", ".gff.gz", ".gff3.gz")):
        rows: list[dict[str, object]] = []
        with _open_text(p) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue
                attrs = parse_gtf_attributes(fields[8])
                strand = fields[6]
                start = int(fields[3])
                end = int(fields[4])
                tss = start if strand == "+" else end
                symbol = attrs.get("gene_name") or attrs.get("Name") or attrs.get("gene_id", "")
                rows.append(
                    {
                        "chromosome": _normalise_chromosome(fields[0]),
                        "tss": tss,
                        "gene": str(symbol).upper(),
                        "gene_id": attrs.get("gene_id", ""),
                        "gene_type": attrs.get("gene_type") or attrs.get("gene_biotype") or "",
                        "strand": strand,
                        "tss_source": p.name,
                    }
                )
        frame = pd.DataFrame(rows)
    else:
        frame = _read_delimited(p)
        rename: dict[str, str] = {}
        lower = {str(c).lower(): c for c in frame.columns}
        for canonical, options in {
            "chromosome": ["chromosome", "chrom", "chr"],
            "tss": ["tss", "position", "start"],
            "gene": ["gene", "gene_symbol", "symbol"],
            "gene_id": ["gene_id", "ensembl_gene_id"],
            "gene_type": ["gene_type", "gene_biotype", "biotype"],
            "strand": ["strand"],
        }.items():
            for option in options:
                if option in lower:
                    rename[lower[option]] = canonical
                    break
        frame = frame.rename(columns=rename)
        missing = {"chromosome", "tss", "gene"}.difference(frame.columns)
        if missing:
            raise ValueError(f"TSS reference is missing {sorted(missing)}")
        if "gene_id" not in frame.columns:
            frame["gene_id"] = ""
        if "gene_type" not in frame.columns:
            frame["gene_type"] = ""
        if "strand" not in frame.columns:
            frame["strand"] = ""
        frame["tss_source"] = p.name
    if frame.empty:
        raise ValueError(f"No gene/TSS records were parsed from {p}")
    frame["chromosome"] = frame["chromosome"].map(_normalise_chromosome)
    frame["tss"] = pd.to_numeric(frame["tss"], errors="coerce")
    frame["gene"] = frame["gene"].fillna("").astype(str).str.upper().str.strip()
    frame = frame.dropna(subset=["tss"])
    frame = frame[frame["gene"].ne("")].copy()
    frame["tss"] = frame["tss"].astype(int)
    return frame.drop_duplicates(["chromosome", "tss", "gene", "strand"]).reset_index(drop=True)


def _guess_ccre_class(row: Sequence[object]) -> str:
    tokens = [str(x).upper() for x in row if x is not None and str(x) != "nan"]
    joined = "|".join(tokens)
    for key in ["PLS", "PELS", "DELS", "CA-H3K4ME3", "CA-CTCF", "CA-TF", "CTCF-ONLY", "DNASE-ONLY", "TF-ONLY"]:
        if key in joined:
            return key
    if "PROMOTER" in joined:
        return "PLS"
    if "PROXIMAL" in joined and "ENHANCER" in joined:
        return "PELS"
    if "DISTAL" in joined and "ENHANCER" in joined:
        return "DELS"
    return "UNCLASSIFIED"


def load_regulatory_intervals(path: str | Path, *, biosample: str = "agnostic") -> pd.DataFrame:
    """Load BED/TSV regulatory intervals.

    BED coordinates are converted from 0-based half-open to 1-based inclusive.
    """
    p = Path(path)
    frame = _read_delimited(p, header=False)
    if frame.shape[1] < 3:
        raise ValueError(f"Regulatory interval file requires at least three columns: {p}")
    rows: list[dict[str, object]] = []
    for values in frame.itertuples(index=False, name=None):
        try:
            start0 = int(values[1])
            end0 = int(values[2])
        except (TypeError, ValueError):
            continue
        name = str(values[3]) if len(values) > 3 else ""
        reg_class = _guess_ccre_class(values)
        rows.append(
            {
                "chromosome": _normalise_chromosome(values[0]),
                "start": start0 + 1,
                "end": end0,
                "regulatory_id": name,
                "regulatory_class": reg_class,
                "biosample": biosample,
                "regulatory_source": p.name,
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError(f"No regulatory intervals were parsed from {p}")
    return out


def load_enhancer_gene_links(path: str | Path) -> pd.DataFrame:
    frame = _read_delimited(path)
    lower = {str(c).lower(): c for c in frame.columns}
    rename: dict[str, str] = {}
    choices = {
        "chromosome": ["chromosome", "chrom", "chr"],
        "start": ["start", "enhancer_start"],
        "end": ["end", "enhancer_end"],
        "target_gene": ["target_gene", "gene", "gene_symbol"],
        "confidence": ["confidence", "score", "abc_score", "link_score"],
        "source_id": ["source_id", "source", "dataset"],
    }
    for canonical, options in choices.items():
        for option in options:
            if option in lower:
                rename[lower[option]] = canonical
                break
    frame = frame.rename(columns=rename)
    missing = {"chromosome", "start", "end", "target_gene"}.difference(frame.columns)
    if missing:
        raise ValueError(f"Enhancer-gene link table is missing {sorted(missing)}")
    if "confidence" not in frame.columns:
        frame["confidence"] = 1.0
    if "source_id" not in frame.columns:
        frame["source_id"] = Path(path).name
    frame["chromosome"] = frame["chromosome"].map(_normalise_chromosome)
    frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
    frame["end"] = pd.to_numeric(frame["end"], errors="coerce")
    frame["target_gene"] = frame["target_gene"].fillna("").astype(str).str.upper().str.strip()
    frame["confidence"] = pd.to_numeric(frame["confidence"], errors="coerce").clip(0, 1)
    return frame.dropna(subset=["start", "end"])[frame["target_gene"].ne("")].copy()


def load_literature_evidence(paths: Sequence[str | Path]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    required = {"gene", "evidence_role", "evidence_domain", "evidence_status", "evidence_strength", "source_id"}
    for path_like in paths:
        path = Path(path_like)
        frame = _read_delimited(path)
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"Literature evidence file {path} is missing {sorted(missing)}")
        frame = frame.copy()
        frame["evidence_file"] = path.name
        frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=sorted(required))
    out = pd.concat(frames, ignore_index=True)
    out["gene"] = out["gene"].fillna("").astype(str).str.upper().str.strip()
    out["evidence_role"] = out["evidence_role"].fillna("").astype(str).str.lower().str.strip()
    out["evidence_domain"] = out["evidence_domain"].fillna("").astype(str).str.lower().str.strip()
    out["evidence_status"] = out["evidence_status"].fillna("").astype(str).str.lower().str.strip()
    out["evidence_strength"] = pd.to_numeric(out["evidence_strength"], errors="coerce").clip(0, 1)
    if "reviewed" not in out.columns:
        out["reviewed"] = True
    out["reviewed"] = out["reviewed"].map(_truthy)
    for col in ["evidence_level", "source_type", "source_year", "source_url_or_doi", "source_note"]:
        if col not in out.columns:
            out[col] = ""
    return out[out["gene"].ne("") & out["reviewed"]].drop_duplicates(
        ["gene", "evidence_role", "evidence_domain", "source_id"], keep="first"
    )



def load_cancer_gene_census(path: str | Path) -> pd.DataFrame:
    """Convert a user-supplied COSMIC Cancer Gene Census table into risk evidence.

    COSMIC files are not redistributed. Users must obtain and use them under the
    applicable COSMIC licence. Common CGC column names are detected.
    """
    frame = _read_delimited(path)
    lower = {str(c).lower(): c for c in frame.columns}
    gene_col = next((lower[x] for x in ["gene symbol", "gene_symbol", "gene", "symbol"] if x in lower), None)
    role_col = next((lower[x] for x in ["role in cancer", "role_in_cancer", "role"] if x in lower), None)
    tier_col = next((lower[x] for x in ["tier", "cgc tier", "cgc_tier"] if x in lower), None)
    if gene_col is None:
        raise ValueError("Cancer Gene Census table requires a gene-symbol column")
    rows: list[dict[str, object]] = []
    for _, original in frame.iterrows():
        gene = str(original[gene_col]).upper().strip()
        if not gene or gene == "NAN":
            continue
        role_value = str(original[role_col]) if role_col is not None and pd.notna(original[role_col]) else ""
        tier_value = pd.to_numeric(original[tier_col], errors="coerce") if tier_col is not None else np.nan
        role_upper = role_value.upper()
        if "TSG" in role_upper or "TUMOR SUPPRESSOR" in role_upper:
            domain = "tumor_suppressor_risk"
        else:
            domain = "oncogenic_risk"
        strength = 1.0 if tier_value == 1 else 0.85 if tier_value == 2 else 0.75
        rows.append(
            {
                "gene": gene,
                "evidence_role": "risk",
                "evidence_domain": domain,
                "evidence_status": "positive",
                "evidence_strength": strength,
                "evidence_level": "curated_cancer_gene_catalog",
                "source_type": "curated_database",
                "source_id": f"COSMIC_CGC:{Path(path).name}",
                "source_year": "",
                "source_url_or_doi": "https://cancer.sanger.ac.uk/census",
                "reviewed": True,
                "source_note": (
                    f"User-supplied Cancer Gene Census role={role_value}; proximity is a safety flag, "
                    "not proof that CpG editing alters the gene."
                ),
                "evidence_file": Path(path).name,
            }
        )
    return pd.DataFrame(rows)


def annotate_nearest_tss(cpg: pd.DataFrame, tss: pd.DataFrame, *, max_bp: int) -> pd.DataFrame:
    out = cpg.copy()
    defaults = {
        "nearest_tss_gene": "",
        "nearest_tss_gene_id": "",
        "nearest_tss_gene_type": "",
        "nearest_tss_position": np.nan,
        "nearest_tss_distance_bp": np.nan,
        "nearest_tss_abs_distance_bp": np.nan,
        "nearest_tss_strand": "",
    }
    for col, default in defaults.items():
        out[col] = default
    for chrom, idx in out.groupby("chromosome").groups.items():
        ref = tss[tss["chromosome"] == chrom].sort_values("tss")
        if ref.empty:
            continue
        ref_pos = ref["tss"].to_numpy(dtype=int)
        for i in idx:
            position = pd.to_numeric(out.at[i, "position"], errors="coerce")
            if not np.isfinite(position):
                continue
            j = int(np.searchsorted(ref_pos, int(position)))
            candidates = [x for x in [j - 1, j] if 0 <= x < len(ref)]
            best = min(candidates, key=lambda k: abs(int(position) - ref_pos[k]))
            record = ref.iloc[best]
            signed = int(position) - int(record["tss"])
            if abs(signed) > max_bp:
                continue
            out.at[i, "nearest_tss_gene"] = record["gene"]
            out.at[i, "nearest_tss_gene_id"] = record.get("gene_id", "")
            out.at[i, "nearest_tss_gene_type"] = record.get("gene_type", "")
            out.at[i, "nearest_tss_position"] = int(record["tss"])
            out.at[i, "nearest_tss_distance_bp"] = signed
            out.at[i, "nearest_tss_abs_distance_bp"] = abs(signed)
            out.at[i, "nearest_tss_strand"] = record.get("strand", "")
    return out


def annotate_interval_overlaps(cpg: pd.DataFrame, intervals: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = cpg.copy()
    links: list[dict[str, object]] = []
    for chrom, idx in out.groupby("chromosome").groups.items():
        ref = intervals[intervals["chromosome"] == chrom]
        if ref.empty:
            continue
        starts = ref["start"].to_numpy(dtype=int)
        ends = ref["end"].to_numpy(dtype=int)
        for i in idx:
            pos = pd.to_numeric(out.at[i, "position"], errors="coerce")
            if not np.isfinite(pos):
                continue
            hit = ref[(starts <= int(pos)) & (ends >= int(pos))]
            for record in hit.itertuples(index=False):
                links.append(
                    {
                        "cpg_id": out.at[i, "cpg_id"],
                        "regulatory_id": record.regulatory_id,
                        "regulatory_class": record.regulatory_class,
                        "biosample": record.biosample,
                        "regulatory_source": record.regulatory_source,
                    }
                )
    overlap = pd.DataFrame(links)
    if overlap.empty:
        out["regulatory_overlap_count"] = 0
        out["regulatory_classes"] = ""
        out["regulatory_biosamples"] = ""
        out["regulatory_class_score"] = 0.0
        return out, overlap
    overlap["regulatory_class_score"] = overlap["regulatory_class"].map(
        lambda x: REGULATORY_CLASS_SCORE.get(str(x).upper(), 0.20)
    )
    summary = overlap.groupby("cpg_id", as_index=False).agg(
        regulatory_overlap_count=("regulatory_id", "size"),
        regulatory_classes=("regulatory_class", lambda x: ";".join(sorted(set(map(str, x))))),
        regulatory_biosamples=("biosample", lambda x: ";".join(sorted(set(map(str, x))))),
        regulatory_class_score=("regulatory_class_score", "max"),
    )
    out = out.merge(summary, on="cpg_id", how="left")
    out["regulatory_overlap_count"] = out["regulatory_overlap_count"].fillna(0).astype(int)
    out["regulatory_classes"] = out["regulatory_classes"].fillna("")
    out["regulatory_biosamples"] = out["regulatory_biosamples"].fillna("")
    out["regulatory_class_score"] = out["regulatory_class_score"].fillna(0.0)
    return out, overlap


def annotate_enhancer_links(cpg: pd.DataFrame, enhancer_links: pd.DataFrame | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = cpg.copy()
    if enhancer_links is None or enhancer_links.empty:
        out["enhancer_linked_genes"] = ""
        out["enhancer_link_max_confidence"] = np.nan
        out["enhancer_link_count"] = 0
        return out, pd.DataFrame()
    rows: list[dict[str, object]] = []
    for chrom, idx in out.groupby("chromosome").groups.items():
        ref = enhancer_links[enhancer_links["chromosome"] == chrom]
        if ref.empty:
            continue
        for i in idx:
            pos = pd.to_numeric(out.at[i, "position"], errors="coerce")
            if not np.isfinite(pos):
                continue
            hit = ref[(ref["start"] <= int(pos)) & (ref["end"] >= int(pos))]
            for record in hit.itertuples(index=False):
                rows.append(
                    {
                        "cpg_id": out.at[i, "cpg_id"],
                        "gene": record.target_gene,
                        "confidence": record.confidence,
                        "source_id": record.source_id,
                    }
                )
    links = pd.DataFrame(rows)
    if links.empty:
        out["enhancer_linked_genes"] = ""
        out["enhancer_link_max_confidence"] = 0.0
        out["enhancer_link_count"] = 0
        return out, links
    summary = links.groupby("cpg_id", as_index=False).agg(
        enhancer_linked_genes=("gene", lambda x: ";".join(sorted(set(map(str, x))))),
        enhancer_link_max_confidence=("confidence", "max"),
        enhancer_link_count=("gene", "size"),
    )
    out = out.merge(summary, on="cpg_id", how="left")
    out["enhancer_linked_genes"] = out["enhancer_linked_genes"].fillna("")
    out["enhancer_link_max_confidence"] = pd.to_numeric(
        out["enhancer_link_max_confidence"], errors="coerce"
    ).fillna(0.0)
    out["enhancer_link_count"] = out["enhancer_link_count"].fillna(0).astype(int)
    return out, links


def build_gene_links(cpg: pd.DataFrame, enhancer_links: pd.DataFrame | None, *, config: Mapping[str, object]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    promoter_bp = int(config["promoter_bp"])
    proximal_bp = int(config["proximal_bp"])
    max_bp = int(config["linked_gene_max_bp"])
    for row in cpg.itertuples(index=False):
        cpg_id = str(row.cpg_id)
        for gene in _split_genes(getattr(row, "gene", "")):
            rows.append(
                {
                    "cpg_id": cpg_id,
                    "gene": gene,
                    "relationship": "direct_epic_gene",
                    "relationship_weight": RELATIONSHIP_WEIGHT["direct_epic_gene"],
                    "distance_bp": np.nan,
                    "source_id": "EPICv2_annotation",
                }
            )
        nearest = str(getattr(row, "nearest_tss_gene", "") or "").upper()
        distance = pd.to_numeric(getattr(row, "nearest_tss_abs_distance_bp", np.nan), errors="coerce")
        if nearest and np.isfinite(distance) and distance <= max_bp:
            if distance <= promoter_bp:
                relation = "promoter_nearest_tss"
                base = RELATIONSHIP_WEIGHT[relation]
            elif distance <= proximal_bp:
                relation = "proximal_nearest_tss"
                base = RELATIONSHIP_WEIGHT[relation]
            else:
                relation = "distal_nearest_tss"
                base = RELATIONSHIP_WEIGHT[relation] * math.exp(-float(distance) / max(max_bp, 1))
            rows.append(
                {
                    "cpg_id": cpg_id,
                    "gene": nearest,
                    "relationship": relation,
                    "relationship_weight": float(base),
                    "distance_bp": float(distance),
                    "source_id": "GENCODE_nearest_TSS",
                }
            )
    if enhancer_links is not None and not enhancer_links.empty:
        for row in enhancer_links.itertuples(index=False):
            rows.append(
                {
                    "cpg_id": str(row.cpg_id),
                    "gene": str(row.gene).upper(),
                    "relationship": "enhancer_link",
                    "relationship_weight": float(RELATIONSHIP_WEIGHT["enhancer_link"] * row.confidence),
                    "distance_bp": np.nan,
                    "source_id": row.source_id,
                }
            )
    links = pd.DataFrame(rows)
    if links.empty:
        return pd.DataFrame(columns=["cpg_id", "gene", "relationship", "relationship_weight", "distance_bp", "source_id"])
    links["relationship_weight"] = pd.to_numeric(links["relationship_weight"], errors="coerce").clip(0, 1)
    links = links.sort_values("relationship_weight", ascending=False).drop_duplicates(
        ["cpg_id", "gene", "relationship"], keep="first"
    )
    return links.reset_index(drop=True)


def _evidence_status(group: pd.DataFrame, role: str) -> str:
    subset = group[group["evidence_role"] == role]
    if subset.empty:
        return "not_reviewed"
    statuses = set(subset["evidence_status"].astype(str))
    if role == "risk":
        if statuses & {"positive", "risk", "mixed"}:
            return "risk_signal" if "mixed" not in statuses else "mixed"
        if statuses & {"negative", "no_known_signal", "clear"}:
            return "no_known_signal"
    else:
        if "mixed" in statuses or ({"positive", "negative"}.issubset(statuses)):
            return "mixed"
        if "positive" in statuses:
            return "positive"
        if "negative" in statuses:
            return "negative"
    return "reviewed_unclear"


def annotate_literature_and_safety(
    cpg: pd.DataFrame,
    gene_links: pd.DataFrame,
    evidence: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    out = cpg.copy()
    if gene_links.empty:
        for col in ["literature_relevance_score", "safety_risk_score_v081"]:
            out[col] = np.nan
        out["literature_evidence_status"] = "not_reviewed"
        out["safety_evidence_status"] = "not_reviewed"
        out["literature_source_ids"] = ""
        out["safety_source_ids"] = ""
        out["linked_gene_count"] = 0
        return out, pd.DataFrame(), pd.DataFrame()

    joined = gene_links.merge(evidence, on="gene", how="left")
    joined["weighted_strength"] = joined["relationship_weight"] * pd.to_numeric(
        joined["evidence_strength"], errors="coerce"
    )
    summaries: list[dict[str, object]] = []
    for cpg_id, group in joined.groupby("cpg_id", sort=False):
        reviewed = group[group["source_id_y"].notna()].copy() if "source_id_y" in group.columns else group.iloc[0:0]
        relevance = reviewed[reviewed["evidence_role"] == "relevance"]
        risk = reviewed[reviewed["evidence_role"] == "risk"]
        protective = reviewed[reviewed["evidence_role"] == "protective"]

        rel_status = _evidence_status(reviewed, "relevance")
        risk_status = _evidence_status(reviewed, "risk")
        if relevance.empty:
            rel_score = np.nan
        else:
            positive = relevance[relevance["evidence_status"].isin(["positive", "mixed"])]
            negative = relevance[relevance["evidence_status"] == "negative"]
            pos = positive["weighted_strength"].max() if not positive.empty else 0.0
            neg = negative["weighted_strength"].max() if not negative.empty else 0.0
            rel_score = float(np.clip(pos - 0.5 * neg, 0, 1))
        if risk.empty:
            risk_score = np.nan
        else:
            positive_risk = risk[risk["evidence_status"].isin(["positive", "risk", "mixed"])]
            clear = risk[risk["evidence_status"].isin(["negative", "no_known_signal", "clear"])]
            if not positive_risk.empty:
                risk_score = float(positive_risk["weighted_strength"].max())
            elif not clear.empty:
                risk_score = 0.0
            else:
                risk_score = np.nan
        if not protective.empty and np.isfinite(risk_score):
            risk_score = float(np.clip(risk_score - 0.20 * protective["weighted_strength"].max(), 0, 1))

        summaries.append(
            {
                "cpg_id": cpg_id,
                "linked_gene_count": int(group["gene"].nunique()),
                "reviewed_linked_gene_count": int(reviewed["gene"].nunique()),
                "literature_relevance_score": rel_score,
                "literature_evidence_status": rel_status,
                "safety_risk_score_v081": risk_score,
                "safety_evidence_status": risk_status,
                "literature_source_ids": ";".join(sorted(set(relevance.get("source_id_y", pd.Series(dtype=str)).dropna().astype(str)))),
                "safety_source_ids": ";".join(sorted(set(risk.get("source_id_y", pd.Series(dtype=str)).dropna().astype(str)))),
                "linked_genes": ";".join(
                    sorted(
                        {
                            gene
                            for value in group["gene"].tolist()
                            for gene in _split_genes(value)
                        }
                    )
                ),
            }
        )
    out = out.merge(pd.DataFrame(summaries), on="cpg_id", how="left")
    out["linked_gene_count"] = out["linked_gene_count"].fillna(0).astype(int)
    out["reviewed_linked_gene_count"] = out["reviewed_linked_gene_count"].fillna(0).astype(int)
    out["literature_evidence_status"] = out["literature_evidence_status"].fillna("not_reviewed")
    out["safety_evidence_status"] = out["safety_evidence_status"].fillna("not_reviewed")
    out["literature_source_ids"] = out["literature_source_ids"].fillna("")
    out["safety_source_ids"] = out["safety_source_ids"].fillna("")

    literature_audit = joined[joined.get("evidence_role", pd.Series(index=joined.index)).eq("relevance")].copy()
    safety_audit = joined[joined.get("evidence_role", pd.Series(index=joined.index)).isin(["risk", "protective"])].copy()
    return out, literature_audit, safety_audit


def score_regulatory_safety_completion(
    v080_cpg_scores: pd.DataFrame,
    v080_region_scores: pd.DataFrame,
    *,
    tss_reference: pd.DataFrame | None = None,
    regulatory_intervals: Sequence[pd.DataFrame] = (),
    enhancer_gene_links: pd.DataFrame | None = None,
    literature_evidence: pd.DataFrame | None = None,
    config: Mapping[str, object] = DEFAULT_V081_CONFIG,
) -> RegulatorySafetyResult:
    cfg = {**DEFAULT_V081_CONFIG, **dict(config)}
    cpg = v080_cpg_scores.copy()
    required = {"cpg_id", "chromosome", "position"}
    missing = required.difference(cpg.columns)
    if missing:
        raise ValueError(f"v0.8.0 CpG score table is missing {sorted(missing)}")
    cpg["chromosome"] = cpg["chromosome"].map(_normalise_chromosome)
    cpg["position"] = pd.to_numeric(cpg["position"], errors="coerce")
    if "causal_evidence_score" in cpg.columns:
        experimental = pd.to_numeric(cpg["causal_evidence_score"], errors="coerce") / 100.0
    elif "safety_adjusted_priority_score" in cpg.columns:
        experimental = pd.to_numeric(cpg["safety_adjusted_priority_score"], errors="coerce") / 100.0
    else:
        raise ValueError("v0.8.0 CpG table requires causal_evidence_score or safety_adjusted_priority_score")
    cpg["integrated_experimental_evidence_score"] = 100.0 * experimental.clip(0, 1)

    tss_supplied = tss_reference is not None and not tss_reference.empty
    if tss_supplied:
        cpg = annotate_nearest_tss(cpg, tss_reference, max_bp=int(cfg["nearest_tss_max_bp"]))
    else:
        for col, default in {
            "nearest_tss_gene": "",
            "nearest_tss_gene_id": "",
            "nearest_tss_gene_type": "",
            "nearest_tss_position": np.nan,
            "nearest_tss_distance_bp": np.nan,
            "nearest_tss_abs_distance_bp": np.nan,
            "nearest_tss_strand": "",
        }.items():
            cpg[col] = default

    reg_supplied = bool(regulatory_intervals)
    overlap_frames: list[pd.DataFrame] = []
    if reg_supplied:
        all_intervals = pd.concat(list(regulatory_intervals), ignore_index=True)
        cpg, overlaps = annotate_interval_overlaps(cpg, all_intervals)
        overlap_frames.append(overlaps)
    else:
        cpg["regulatory_overlap_count"] = 0
        cpg["regulatory_classes"] = ""
        cpg["regulatory_biosamples"] = ""
        cpg["regulatory_class_score"] = np.nan

    enhancer_supplied = enhancer_gene_links is not None and not enhancer_gene_links.empty
    cpg, enhancer_overlap = annotate_enhancer_links(cpg, enhancer_gene_links)
    gene_links = build_gene_links(cpg, enhancer_overlap, config=cfg)

    evidence = literature_evidence if literature_evidence is not None else pd.DataFrame()
    cpg, literature_audit, safety_audit = annotate_literature_and_safety(cpg, gene_links, evidence)

    distance = pd.to_numeric(cpg["nearest_tss_abs_distance_bp"], errors="coerce")
    tss_score = np.select(
        [distance <= 200, distance <= int(cfg["promoter_bp"]), distance <= int(cfg["proximal_bp"]), distance <= 100_000],
        [1.0, 0.90, 0.70, 0.40],
        default=np.where(distance.notna(), 0.15, np.nan),
    )
    cpg["nearest_tss_regulatory_score"] = tss_score
    reg_parts = pd.DataFrame(
        {
            "tss": cpg["nearest_tss_regulatory_score"] if tss_supplied else pd.Series(np.nan, index=cpg.index),
            "ccre": cpg["regulatory_class_score"] if reg_supplied else pd.Series(np.nan, index=cpg.index),
            "enhancer": cpg["enhancer_link_max_confidence"] if enhancer_supplied else pd.Series(np.nan, index=cpg.index),
        }
    )
    cpg["regulatory_evidence_score"] = reg_parts.mean(axis=1, skipna=True)
    cpg["regulatory_evidence_completeness"] = reg_parts.notna().sum(axis=1) / 3.0
    cpg["regulatory_annotation_status"] = np.select(
        [cpg["regulatory_evidence_completeness"] >= 0.99, cpg["regulatory_evidence_completeness"] >= 0.34],
        ["complete", "partial"],
        default="not_supplied",
    )

    literature_available = cpg["literature_relevance_score"].notna()
    safety_available = cpg["safety_risk_score_v081"].notna()
    cpg["literature_review_completeness"] = literature_available.astype(float)
    cpg["safety_review_completeness"] = safety_available.astype(float)
    cpg["safety_uncertainty_score"] = 1.0 - (
        0.55 * cpg["safety_review_completeness"]
        + 0.30 * cpg["regulatory_evidence_completeness"]
        + 0.15 * gene_links.groupby("cpg_id")["gene"].nunique().reindex(cpg["cpg_id"]).fillna(0).gt(0).to_numpy(dtype=float)
    )
    cpg["safety_uncertainty_score"] = cpg["safety_uncertainty_score"].clip(0, 1)

    weights = {
        "experimental": float(cfg["experimental_weight"]),
        "regulatory": float(cfg["regulatory_weight"]),
        "literature": float(cfg["literature_relevance_weight"]),
    }
    numer = weights["experimental"] * experimental.to_numpy(dtype=float)
    denom = np.full(len(cpg), weights["experimental"], dtype=float)
    reg = pd.to_numeric(cpg["regulatory_evidence_score"], errors="coerce").to_numpy(dtype=float)
    valid_reg = np.isfinite(reg)
    numer[valid_reg] += weights["regulatory"] * reg[valid_reg]
    denom[valid_reg] += weights["regulatory"]
    lit = pd.to_numeric(cpg["literature_relevance_score"], errors="coerce").to_numpy(dtype=float)
    valid_lit = np.isfinite(lit)
    numer[valid_lit] += weights["literature"] * lit[valid_lit]
    denom[valid_lit] += weights["literature"]
    base = np.divide(numer, denom, out=np.full(len(cpg), np.nan), where=denom > 0)
    cpg["literature_grounded_priority_score"] = 100.0 * base
    cpg["overall_evidence_completeness_v081"] = (
        weights["experimental"]
        + weights["regulatory"] * cpg["regulatory_evidence_completeness"]
        + weights["literature"] * cpg["literature_review_completeness"]
    ) / sum(weights.values())

    risk = pd.to_numeric(cpg["safety_risk_score_v081"], errors="coerce")
    risk_factor = np.where(risk.notna(), 1.0 - float(cfg["risk_penalty"]) * risk.fillna(0), 1.0)
    uncertainty_factor = 1.0 - float(cfg["uncertainty_penalty"]) * cpg["safety_uncertainty_score"]
    cpg["risk_uncertainty_adjusted_priority_score"] = (
        cpg["literature_grounded_priority_score"] * risk_factor * uncertainty_factor
    )

    cpg["evidence_readiness_grade"] = np.select(
        [
            (cpg["risk_uncertainty_adjusted_priority_score"] >= 65)
            & (cpg["regulatory_evidence_completeness"] >= 2 / 3)
            & safety_available
            & (risk < 0.50),
            cpg["risk_uncertainty_adjusted_priority_score"] >= 58,
            cpg["risk_uncertainty_adjusted_priority_score"] >= 45,
        ],
        ["E1_design_ready", "E2_high_priority_review_required", "E3_exploratory"],
        default="E4_insufficient",
    )
    cpg["experimental_priority_tier_v081"] = np.select(
        [
            cpg["evidence_readiness_grade"] == "E1_design_ready",
            cpg["evidence_readiness_grade"] == "E2_high_priority_review_required",
            cpg["evidence_readiness_grade"] == "E3_exploratory",
        ],
        ["P1", "P2", "P3"],
        default="P4",
    )
    cpg["manual_review_required"] = (
        (cpg["integrated_experimental_evidence_score"] >= 60)
        & ((cpg["literature_evidence_status"] == "not_reviewed") | (cpg["safety_evidence_status"] == "not_reviewed"))
    )
    cpg["interpretation_v081"] = "literature_grounded_regulatory_priority_not_causal_or_clinical_safety_proof"
    cpg = cpg.sort_values(
        ["risk_uncertainty_adjusted_priority_score", "integrated_experimental_evidence_score"],
        ascending=False,
    ).reset_index(drop=True)
    cpg["global_priority_rank_v081"] = np.arange(1, len(cpg) + 1)

    region = _score_regions_v081(cpg, v080_region_scores)
    status = _build_status_summary(cpg)
    manual = cpg[cpg["manual_review_required"]].copy()
    manual["recommended_review_action"] = np.select(
        [
            (manual["literature_evidence_status"] == "not_reviewed") & (manual["safety_evidence_status"] == "not_reviewed"),
            manual["safety_evidence_status"] == "not_reviewed",
        ],
        ["review_biological_relevance_and_safety", "review_safety_evidence"],
        default="review_biological_relevance",
    )
    def _review_gene(row: pd.Series) -> str:
        genes = _split_genes(row.get("linked_genes", ""))
        if genes:
            return genes[0]
        for col in ["nearest_tss_gene", "gene"]:
            value = str(row.get(col, "") or "").strip().upper()
            if value and value != "NAN":
                return value
        return ""
    manual["primary_review_gene"] = manual.apply(_review_gene, axis=1)
    manual["pubmed_relevance_query"] = manual["primary_review_gene"].map(
        lambda gene: "https://pubmed.ncbi.nlm.nih.gov/?term=" + quote_plus(
            f'{gene} AND (aging OR senescence OR rejuvenation OR reprogramming OR DNA methylation)'
        ) if gene else ""
    )
    manual["pubmed_safety_query"] = manual["primary_review_gene"].map(
        lambda gene: "https://pubmed.ncbi.nlm.nih.gov/?term=" + quote_plus(
            f'{gene} AND (cancer OR oncogene OR tumor suppressor OR pluripotency OR genome stability)'
        ) if gene else ""
    )
    shortlist = cpg.head(max(int(cfg["shortlist_size"]), 1)).copy()

    summary = {
        "analysis_status": "literature_grounded_regulatory_annotation_and_safety_evidence_completion",
        "interpretation_scope": (
            "v0.8.1 separates not-reviewed evidence from reviewed negative evidence, adds nearest-TSS and cCRE/enhancer context, "
            "and penalises safety uncertainty. Scores prioritise experiments and do not prove causal rejuvenation or clinical safety."
        ),
        "n_candidate_cpgs": len(cpg),
        "n_tss_annotated": int(cpg["nearest_tss_gene"].ne("").sum()),
        "n_regulatory_overlap": int(cpg["regulatory_overlap_count"].gt(0).sum()),
        "n_enhancer_linked": int(cpg["enhancer_link_count"].gt(0).sum()),
        "n_literature_reviewed": int(cpg["literature_evidence_status"].ne("not_reviewed").sum()),
        "n_safety_reviewed": int(cpg["safety_evidence_status"].ne("not_reviewed").sum()),
        "n_safety_risk_signal": int(cpg["safety_evidence_status"].isin(["risk_signal", "mixed"]).sum()),
        "n_manual_review_queue": len(manual),
        "n_priority_P1": int((cpg["experimental_priority_tier_v081"] == "P1").sum()),
        "n_priority_P2": int((cpg["experimental_priority_tier_v081"] == "P2").sum()),
        "top_cpgs": cpg.head(20)[
            [
                "cpg_id",
                "gene",
                "nearest_tss_gene",
                "unique_region_id",
                "integrated_experimental_evidence_score",
                "regulatory_evidence_score",
                "literature_relevance_score",
                "safety_risk_score_v081",
                "safety_uncertainty_score",
                "risk_uncertainty_adjusted_priority_score",
                "evidence_readiness_grade",
            ]
        ].to_dict("records"),
        "top_regions": region.head(15)[
            ["unique_region_id", "genes", "n_candidate_cpgs", "region_priority_score_v081", "region_priority_tier_v081"]
        ].to_dict("records") if not region.empty else [],
    }
    manifest = {
        "package_version": "0.8.1",
        "candidate_lists_frozen_before_v081": True,
        "unknown_vs_zero_policy": "not_reviewed_is_NaN; zero_is_used_only_for_explicit_reviewed_negative_or_no_known_signal_evidence",
        "coordinate_policy": "GRCh38; EPIC positions and GENCODE TSS are 1-based; BED intervals converted from 0-based half-open",
        "regulatory_sources_supplied": {
            "tss": tss_supplied,
            "ccre_or_interval": reg_supplied,
            "enhancer_gene_links": enhancer_supplied,
        },
        "literature_guardrail": "only reviewed rows with traceable source_id contribute to literature or safety scores",
        "safety_guardrail": "absence_of_a_risk_signal_is_not_safety_proof; unreviewed evidence receives uncertainty penalty",
        "causal_language_guardrail": "prioritisation_not_causal_proof",
        "recommended_next_stage": "candidate_specific_manual_literature_review_then_sequence_level_sgRNA_design",
        "config": dict(cfg),
    }
    return RegulatorySafetyResult(
        cpg_scores=cpg,
        region_scores=region,
        gene_links=gene_links,
        literature_audit=literature_audit,
        safety_audit=safety_audit,
        evidence_status_summary=status,
        manual_review_queue=manual,
        shortlist=shortlist,
        summary=_jsonable(summary),
        manifest=_jsonable(manifest),
    )


def _score_regions_v081(cpg: pd.DataFrame, base_regions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for region_id, group in cpg.groupby("unique_region_id", dropna=False):
        if pd.isna(region_id):
            continue
        rows.append(
            {
                "unique_region_id": region_id,
                "genes": ";".join(
                    sorted(
                        {
                            gene
                            for value in group.get(
                                "linked_genes",
                                pd.Series("", index=group.index, dtype="object"),
                            ).tolist()
                            for gene in _split_genes(value)
                        }
                    )
                ),
                "n_candidate_cpgs": int(len(group)),
                "median_priority_score_v081": float(group["risk_uncertainty_adjusted_priority_score"].median()),
                "max_priority_score_v081": float(group["risk_uncertainty_adjusted_priority_score"].max()),
                "mean_regulatory_evidence_score": float(group["regulatory_evidence_score"].mean()),
                "mean_safety_uncertainty_score": float(group["safety_uncertainty_score"].mean()),
                "n_p1_cpgs_v081": int((group["experimental_priority_tier_v081"] == "P1").sum()),
                "n_manual_review_required": int(group["manual_review_required"].sum()),
            }
        )
    region = pd.DataFrame(rows)
    if region.empty:
        return region
    if not base_regions.empty and "unique_region_id" in base_regions.columns:
        keep = [c for c in ["unique_region_id", "chromosome", "start", "end", "n_candidate_cpgs", "region_class", "regional_responsiveness_score"] if c in base_regions.columns]
        base = base_regions[keep].drop_duplicates("unique_region_id")
        region = region.merge(base, on="unique_region_id", how="left", suffixes=("", "_base"))
        if "n_candidate_cpgs_base" in region.columns:
            region["n_candidate_cpgs"] = region["n_candidate_cpgs_base"].fillna(region["n_candidate_cpgs"])
            region = region.drop(columns="n_candidate_cpgs_base")
    multi_bonus = np.minimum(pd.to_numeric(region["n_candidate_cpgs"], errors="coerce").fillna(1), 3) - 1
    multi_bonus = multi_bonus / 2.0
    response = pd.to_numeric(region.get("regional_responsiveness_score", np.nan), errors="coerce")
    response_rank = response.rank(pct=True).fillna(0.0)
    region["region_priority_score_v081"] = (
        0.68 * region["median_priority_score_v081"]
        + 12.0 * multi_bonus
        + 8.0 * response_rank
        - 8.0 * region["mean_safety_uncertainty_score"]
    )
    region["region_priority_tier_v081"] = np.select(
        [region["region_priority_score_v081"] >= 65, region["region_priority_score_v081"] >= 52, region["region_priority_score_v081"] >= 40],
        ["R1", "R2", "R3"],
        default="R4",
    )
    region["interpretation_v081"] = "regulatory_literature_safety_completed_region_priority_not_causal_proof"
    return region.sort_values("region_priority_score_v081", ascending=False).reset_index(drop=True)


def _build_status_summary(cpg: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for column in [
        "regulatory_annotation_status",
        "literature_evidence_status",
        "safety_evidence_status",
        "evidence_readiness_grade",
        "experimental_priority_tier_v081",
    ]:
        for status, count in cpg[column].value_counts(dropna=False).items():
            rows.append(
                {
                    "dimension": column,
                    "status": str(status),
                    "count": int(count),
                    "fraction": float(count / max(len(cpg), 1)),
                }
            )
    return pd.DataFrame(rows)


def write_regulatory_safety_outputs(result: RegulatorySafetyResult, output_dir: str | Path) -> None:
    out = Path(output_dir)
    figures = out / "figures"
    out.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    result.cpg_scores.to_csv(out / "candidate_cpg_regulatory_safety_v081.csv.gz", index=False, compression="gzip")
    result.region_scores.to_csv(out / "candidate_region_regulatory_safety_v081.csv", index=False)
    result.gene_links.to_csv(out / "candidate_gene_links_v081.csv.gz", index=False, compression="gzip")
    result.literature_audit.to_csv(out / "literature_evidence_audit.csv", index=False)
    result.safety_audit.to_csv(out / "safety_evidence_audit.csv", index=False)
    result.evidence_status_summary.to_csv(out / "evidence_status_summary.csv", index=False)
    result.manual_review_queue.to_csv(out / "manual_literature_safety_review_queue.csv", index=False)
    result.shortlist.to_csv(out / "top_candidate_shortlist_v081.csv", index=False)
    (out / "regulatory_safety_v081_summary.json").write_text(json.dumps(result.summary, indent=2), encoding="utf-8")
    (out / "regulatory_safety_v081_manifest.json").write_text(json.dumps(result.manifest, indent=2), encoding="utf-8")
    _write_figures(result, figures)


def _write_figures(result: RegulatorySafetyResult, figures: Path) -> None:
    import matplotlib.pyplot as plt

    cpg = result.cpg_scores
    if not cpg.empty:
        distances = pd.to_numeric(cpg["nearest_tss_abs_distance_bp"], errors="coerce").dropna()
        if len(distances):
            plt.figure(figsize=(8, 5))
            plt.hist(np.log10(distances + 1), bins=30)
            plt.xlabel("log10(nearest TSS distance + 1 bp)")
            plt.ylabel("CpG count")
            plt.title("Nearest-TSS distance distribution")
            plt.tight_layout()
            plt.savefig(figures / "nearest_tss_distance.png", dpi=220)
            plt.close()

        plt.figure(figsize=(7, 5))
        plt.scatter(cpg["safety_uncertainty_score"], cpg["integrated_experimental_evidence_score"], s=18, alpha=0.65)
        plt.xlabel("Safety uncertainty score")
        plt.ylabel("Integrated experimental evidence score")
        plt.title("Experimental evidence versus safety uncertainty")
        plt.tight_layout()
        plt.savefig(figures / "experimental_evidence_vs_safety_uncertainty.png", dpi=220)
        plt.close()

        status = cpg["regulatory_annotation_status"].value_counts()
        plt.figure(figsize=(7, 4))
        plt.bar(status.index.astype(str), status.values)
        plt.ylabel("CpG count")
        plt.title("Regulatory annotation completeness")
        plt.tight_layout()
        plt.savefig(figures / "regulatory_annotation_status.png", dpi=220)
        plt.close()

        top = cpg.head(min(25, len(cpg)))
        plt.figure(figsize=(9, 7))
        plt.barh(top["cpg_id"][::-1], top["risk_uncertainty_adjusted_priority_score"][::-1])
        plt.xlabel("Risk/uncertainty-adjusted priority")
        plt.title("Top v0.8.1 candidate CpGs")
        plt.tight_layout()
        plt.savefig(figures / "top_candidate_cpgs_v081.png", dpi=220)
        plt.close()

    if not result.region_scores.empty:
        top_regions = result.region_scores.head(20)
        plt.figure(figsize=(9, 6))
        plt.barh(top_regions["unique_region_id"][::-1], top_regions["region_priority_score_v081"][::-1])
        plt.xlabel("Region priority score v0.8.1")
        plt.title("Top literature/regulatory-completed regions")
        plt.tight_layout()
        plt.savefig(figures / "top_candidate_regions_v081.png", dpi=220)
        plt.close()
