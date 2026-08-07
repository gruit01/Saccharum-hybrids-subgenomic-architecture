#!/usr/bin/env bash
set -euo pipefail

WORK_DIR="/public1/home/stu_gaoruiting/project/Ssin/HiFiasm/merge/2024-8-28/pop/Ai-input"
TAG="${TAG:-task1_gp_full}"
OUT_TAG="${OUT_TAG:-task1_ridge_elasticnet_full}"
FOLDS="${FOLDS:-5}"
REPEATS="${REPEATS:-5}"
RIDGE_TOP_K="${RIDGE_TOP_K:-5000}"
ENET_TOP_K="${ENET_TOP_K:-1000}"
ENET_ALPHA="${ENET_ALPHA:-0.5}"

mkdir -p "${WORK_DIR}/results/task1_ridge_elasticnet" "${WORK_DIR}/logs"

echo "[run Ridge / ElasticNet CV] $(date)"
module load apps/R/4.4.3
Rscript "${WORK_DIR}/programs/task1_ridge_elasticnet.R" \
  --work-dir "${WORK_DIR}" \
  --tag "${TAG}" \
  --out-tag "${OUT_TAG}" \
  --folds "${FOLDS}" \
  --repeats "${REPEATS}" \
  --ridge-top-k "${RIDGE_TOP_K}" \
  --enet-top-k "${ENET_TOP_K}" \
  --enet-alpha "${ENET_ALPHA}" \
  --seed 20260509

echo "[done] $(date)"
