#!/usr/bin/env bash
set -euo pipefail

WORK_DIR="/public1/home/stu_gaoruiting/project/Ssin/HiFiasm/merge/2024-8-28/pop/Ai-input"
MARKER_SOURCE="${WORK_DIR}/results/task1_gwas/task1_GZZTF_gwas_full.all_markers.tsv.gz"
TAG="${TAG:-task1_gp_full}"
MAX_MARKERS="${MAX_MARKERS:-0}"
MAX_LINES_PER_FILE="${MAX_LINES_PER_FILE:-0}"

mkdir -p "${WORK_DIR}/results/task1_genomic_prediction" \
         "${WORK_DIR}/tmp/task1_genomic_prediction" \
         "${WORK_DIR}/logs" \
         "${WORK_DIR}/metadata"

echo "[build GP marker matrix] $(date)"
python "${WORK_DIR}/programs/task1_build_gp_matrix.py" \
  --work-dir "${WORK_DIR}" \
  --marker-source "${MARKER_SOURCE}" \
  --tag "${TAG}" \
  --max-markers "${MAX_MARKERS}" \
  --max-lines-per-file "${MAX_LINES_PER_FILE}" \
  --progress-every 10

echo "[run genomic prediction CV] $(date)"
module load apps/R/4.4.3
Rscript "${WORK_DIR}/programs/task1_genomic_prediction.R" \
  --work-dir "${WORK_DIR}" \
  --tag "${TAG}" \
  --folds 5 \
  --repeats 5 \
  --seed 20260509

echo "[done] $(date)"
