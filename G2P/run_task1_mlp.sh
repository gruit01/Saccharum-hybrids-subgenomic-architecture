#!/usr/bin/env bash
set -euo pipefail

WORK_DIR="/public1/home/stu_gaoruiting/project/Ssin/HiFiasm/merge/2024-8-28/pop/Ai-input"
TAG="${TAG:-task1_gp_full}"
OUT_TAG="${OUT_TAG:-task1_mlp_full}"
FOLDS="${FOLDS:-5}"
REPEATS="${REPEATS:-5}"
PCA_TOP_K="${PCA_TOP_K:-5000}"
PCA_N="${PCA_N:-30}"
PCA_HIDDEN="${PCA_HIDDEN:-16}"
DIRECT_TOP_K="${DIRECT_TOP_K:-200}"
DIRECT_HIDDEN="${DIRECT_HIDDEN:-16}"
EPOCHS="${EPOCHS:-400}"
LR="${LR:-0.01}"
WD="${WD:-0.001}"

mkdir -p "${WORK_DIR}/results/task1_mlp" "${WORK_DIR}/logs"

echo "[run PCA/screened-feature MLP baselines] $(date)"
module load apps/R/4.4.3
Rscript "${WORK_DIR}/programs/task1_mlp_baseline.R" \
  --work-dir "${WORK_DIR}" \
  --tag "${TAG}" \
  --out-tag "${OUT_TAG}" \
  --folds "${FOLDS}" \
  --repeats "${REPEATS}" \
  --pca-top-k "${PCA_TOP_K}" \
  --pca-n "${PCA_N}" \
  --pca-hidden "${PCA_HIDDEN}" \
  --direct-top-k "${DIRECT_TOP_K}" \
  --direct-hidden "${DIRECT_HIDDEN}" \
  --epochs "${EPOCHS}" \
  --learning-rate "${LR}" \
  --weight-decay "${WD}" \
  --seed 20260509

echo "[done] $(date)"
