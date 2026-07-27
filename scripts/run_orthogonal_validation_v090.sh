#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CANDIDATE_FILE="${CANDIDATE_FILE:-results/real/regulatory_safety_v081/candidate_cpg_regulatory_safety_v081.csv.gz}"
FROZEN="data/reference/frozen_primary_shortlist_v090.csv"
RAW="data/raw/orthogonal_validation"
OUT="results/real/orthogonal_validation_v090"

python scripts/00_download_orthogonal_validation_data.py --datasets core --output "$RAW"

MULTI_ARGS=()
if [[ "${RUN_EXTENDED:-0}" == "1" ]]; then
  python scripts/00_download_orthogonal_validation_data.py --datasets extended --output "$RAW"
  MULTI_ARGS=(
    --multidonor-matrix "$RAW/GSE179847_Cell_lifespan_DNAm_processed_matrix.csv.gz"
    --multidonor-metadata "$RAW/GSE179847_series_matrix.txt.gz"
  )
fi

python scripts/04h_validate_external_reprogramming.py \
  --candidate-file "$CANDIDATE_FILE" \
  --frozen-shortlist "$FROZEN" \
  --series-matrix "$RAW/GSE142439_series_matrix.txt.gz" \
  --output-dir "$OUT/external_reprogramming"

python scripts/04i_validate_natural_aging_direction.py \
  --candidate-file "$CANDIDATE_FILE" \
  --frozen-shortlist "$FROZEN" \
  --replicative-matrix "$RAW/GSE131280_series_matrix.txt.gz" \
  --chronological-matrix "$RAW/GSE22595_series_matrix.txt.gz" \
  "${MULTI_ARGS[@]}" \
  --output-dir "$OUT/natural_aging"

python scripts/04j_validate_fibroblast_eqtm.py \
  --candidate-file "$CANDIDATE_FILE" \
  --frozen-shortlist "$FROZEN" \
  --methylation-matrix "$RAW/GSE52025_series_matrix.txt.gz" \
  --expression-matrix "$RAW/GSE53243_series_matrix.txt.gz" \
  --expression-annotation "$RAW/GPL6883_data.txt" \
  --gene-map "$FROZEN" \
  --output-dir "$OUT/eqtm"

python scripts/04k_assess_regulatory_links_and_editability.py \
  --candidate-file "$CANDIDATE_FILE" \
  --frozen-shortlist "$FROZEN" \
  --output-dir "$OUT/regulatory_editability"

RE2G_ARGS=()
if [[ -f "$OUT/regulatory_editability/fibroblast_re2g_candidate_links_v090.csv" ]]; then
  RE2G_ARGS=(
    --re2g-file "$OUT/regulatory_editability/fibroblast_re2g_candidate_links_v090.csv"
  )
fi

python scripts/04l_integrate_orthogonal_validation.py \
  --frozen-shortlist "$FROZEN" \
  --external-file "$OUT/external_reprogramming/external_reprogramming_validation_v090.csv.gz" \
  --aging-file "$OUT/natural_aging/natural_aging_direction_validation_v090.csv.gz" \
  --eqtm-file "$OUT/eqtm/fibroblast_eqtm_validation_v090.csv" \
  "${RE2G_ARGS[@]}" \
  --guide-summary-file "$OUT/regulatory_editability/candidate_guide_feasibility_summary_v090.csv" \
  --output-dir "$OUT/integrated"

echo "Orthogonal validation completed: $OUT"
