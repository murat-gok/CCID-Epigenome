from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ccid.regulatory_safety import (
    _split_genes,
    annotate_nearest_tss,
    load_cancer_gene_census,
    load_literature_evidence,
    load_regulatory_intervals,
    load_tss_reference,
    parse_gtf_attributes,
    score_regulatory_safety_completion,
    write_regulatory_safety_outputs,
)


def _cpg() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cpg_id": ["cg1", "cg2", "cg3", "cg4"],
            "chromosome": ["1", "1", "2", "3"],
            "position": [1000, 5200, 1000, 9000],
            "gene": ["SIRT1", "", "TP53", ""],
            "unique_region_id": ["chr1:1000-1000", "chr1:5200-5200", "chr2:1000-1000", "chr3:9000-9000"],
            "causal_evidence_score": [75.0, 70.0, 76.0, 65.0],
            "safety_adjusted_priority_score": [75.0, 70.0, 76.0, 65.0],
            "target_delta_beta": [-0.6, -0.5, -0.7, -0.4],
        }
    )


def _regions() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "unique_region_id": ["chr1:1000-1000", "chr1:5200-5200", "chr2:1000-1000", "chr3:9000-9000"],
            "chromosome": ["1", "1", "2", "3"],
            "start": [1000, 5200, 1000, 9000],
            "end": [1000, 5200, 1000, 9000],
            "n_candidate_cpgs": [1, 1, 1, 1],
            "region_class": ["singleton_locus"] * 4,
            "regional_responsiveness_score": [0.04, 0.03, 0.02, 0.01],
        }
    )


def _tss() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chromosome": ["1", "1", "2"],
            "tss": [900, 5000, 1100],
            "gene": ["SIRT1", "GENE2", "TP53"],
            "gene_id": ["ENSG1", "ENSG2", "ENSG3"],
            "gene_type": ["protein_coding"] * 3,
            "strand": ["+", "+", "-"],
            "tss_source": ["test"] * 3,
        }
    )


def _regulatory(tmp_path: Path) -> pd.DataFrame:
    bed = tmp_path / "ccre.bed"
    bed.write_text(
        "chr1\t950\t1050\tEH1\t0\t.\t0\t0\t0\tPLS\n"
        "chr1\t5100\t5300\tEH2\t0\t.\t0\t0\t0\tdELS\n"
        "chr2\t950\t1050\tEH3\t0\t.\t0\t0\t0\tpELS\n",
        encoding="utf-8",
    )
    return load_regulatory_intervals(bed, biosample="test")


def _literature(tmp_path: Path) -> pd.DataFrame:
    path = tmp_path / "literature.csv"
    pd.DataFrame(
        {
            "gene": ["SIRT1", "TP53", "GENE2"],
            "evidence_role": ["relevance", "risk", "risk"],
            "evidence_domain": ["aging_longevity", "tumor_suppressor_risk", "oncogenic_risk"],
            "evidence_status": ["positive", "positive", "no_known_signal"],
            "evidence_strength": [0.8, 0.95, 0.0],
            "source_id": ["PMID:1", "PMID:2", "REVIEW:1"],
            "reviewed": [True, True, True],
        }
    ).to_csv(path, index=False)
    return load_literature_evidence([path])


def test_split_genes_ignores_pandas3_missing_scalars() -> None:
    assert _split_genes(np.nan) == []
    assert _split_genes(pd.NA) == []
    assert _split_genes("nan") == []
    assert _split_genes("SSH1;tmcc1|SSH1") == ["SSH1", "TMCC1"]


def test_gtf_attribute_parser_and_tss_loading(tmp_path: Path) -> None:
    attrs = parse_gtf_attributes('gene_id "ENSG1"; gene_name "ABC"; gene_type "protein_coding";')
    assert attrs["gene_name"] == "ABC"
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\ttest\tgene\t100\t200\t.\t+\t.\tgene_id "ENSG1"; gene_name "ABC"; gene_type "protein_coding";\n'
        'chr1\ttest\tgene\t400\t500\t.\t-\t.\tgene_id "ENSG2"; gene_name "XYZ"; gene_type "lncRNA";\n',
        encoding="utf-8",
    )
    tss = load_tss_reference(gtf)
    assert tss.set_index("gene").loc["ABC", "tss"] == 100
    assert tss.set_index("gene").loc["XYZ", "tss"] == 500


def test_unknown_is_not_zero_and_reviewed_no_signal_is_zero(tmp_path: Path) -> None:
    result = score_regulatory_safety_completion(
        _cpg(),
        _regions(),
        tss_reference=_tss(),
        regulatory_intervals=[_regulatory(tmp_path)],
        literature_evidence=_literature(tmp_path),
        config={"shortlist_size": 4},
    )
    scored = result.cpg_scores.set_index("cpg_id")
    assert scored.loc["cg1", "literature_relevance_score"] > 0
    assert np.isnan(scored.loc["cg4", "literature_relevance_score"])
    assert np.isnan(scored.loc["cg4", "safety_risk_score_v081"])
    assert scored.loc["cg2", "safety_risk_score_v081"] == 0.0
    assert scored.loc["cg3", "safety_risk_score_v081"] > 0.8
    assert scored.loc["cg3", "risk_uncertainty_adjusted_priority_score"] < scored.loc["cg3", "literature_grounded_priority_score"]
    assert scored.loc["cg4", "manual_review_required"]


def test_regulatory_annotation_and_outputs(tmp_path: Path) -> None:
    result = score_regulatory_safety_completion(
        _cpg(),
        _regions(),
        tss_reference=_tss(),
        regulatory_intervals=[_regulatory(tmp_path)],
        literature_evidence=_literature(tmp_path),
        config={"shortlist_size": 3},
    )
    scored = result.cpg_scores.set_index("cpg_id")
    assert scored.loc["cg1", "nearest_tss_gene"] == "SIRT1"
    assert scored.loc["cg1", "regulatory_classes"] == "PLS"
    assert scored.loc["cg2", "regulatory_classes"] == "DELS"
    assert result.summary["n_candidate_cpgs"] == 4
    assert len(result.shortlist) == 3
    assert not result.manual_review_queue.empty

    out = tmp_path / "out"
    write_regulatory_safety_outputs(result, out)
    expected = [
        "candidate_cpg_regulatory_safety_v081.csv.gz",
        "candidate_region_regulatory_safety_v081.csv",
        "candidate_gene_links_v081.csv.gz",
        "literature_evidence_audit.csv",
        "safety_evidence_audit.csv",
        "evidence_status_summary.csv",
        "manual_literature_safety_review_queue.csv",
        "top_candidate_shortlist_v081.csv",
        "regulatory_safety_v081_summary.json",
        "regulatory_safety_v081_manifest.json",
        "figures/nearest_tss_distance.png",
        "figures/experimental_evidence_vs_safety_uncertainty.png",
        "figures/top_candidate_regions_v081.png",
    ]
    for name in expected:
        assert (out / name).exists(), name


def test_missing_regulatory_inputs_are_partial_not_false_negative(tmp_path: Path) -> None:
    result = score_regulatory_safety_completion(
        _cpg(),
        _regions(),
        literature_evidence=_literature(tmp_path),
    )
    assert result.cpg_scores["regulatory_evidence_score"].isna().all()
    assert (result.cpg_scores["regulatory_annotation_status"] == "not_supplied").all()
    assert (result.cpg_scores["regulatory_evidence_completeness"] == 0).all()


def test_cancer_gene_census_adapter(tmp_path: Path) -> None:
    path = tmp_path / "cgc.tsv"
    pd.DataFrame(
        {
            "Gene Symbol": ["TP53", "MYC"],
            "Role in Cancer": ["TSG", "oncogene"],
            "Tier": [1, 2],
        }
    ).to_csv(path, sep="\t", index=False)
    evidence = load_cancer_gene_census(path)
    assert set(evidence["gene"]) == {"TP53", "MYC"}
    assert evidence.set_index("gene").loc["TP53", "evidence_domain"] == "tumor_suppressor_risk"
    assert evidence.set_index("gene").loc["MYC", "evidence_strength"] == 0.85
