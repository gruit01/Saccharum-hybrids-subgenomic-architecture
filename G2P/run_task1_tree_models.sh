#!/usr/bin/env bash
set -euo pipefail

WORK_DIR="/public1/home/stu_gaoruiting/project/Ssin/HiFiasm/merge/2024-8-28/pop/Ai-input"
TAG="${TAG:-task1_gp_full}"
OUT_TAG="${OUT_TAG:-task1_tree_models_full}"
FOLDS="${FOLDS:-5}"
REPEATS="${REPEATS:-5}"
RF_TOP_K="${RF_TOP_K:-1000}"
RF_TREES="${RF_TREES:-80}"
GBDT_TOP_K="${GBDT_TOP_K:-500}"
GBDT_ROUNDS="${GBDT_ROUNDS:-80}"
GBDT_LR="${GBDT_LR:-0.05}"

mkdir -p "${WORK_DIR}/results/task1_tree_models" "${WORK_DIR}/logs"

echo "[run XGBoost-like / RandomForest-like tree baselines] $(date)"
module load apps/R/4.4.3
Rscript "${WORK_DIR}/programs/task1_tree_models.R" \
  --work-dir "${WORK_DIR}" \
  --tag "${TAG}" \
  --out-tag "${OUT_TAG}" \
  --folds "${FOLDS}" \
  --repeats "${REPEATS}" \
  --rf-top-k "${RF_TOP_K}" \
  --rf-trees "${RF_TREES}" \
  --gbdt-top-k "${GBDT_TOP_K}" \
  --gbdt-rounds "${GBDT_ROUNDS}" \
  --gbdt-lr "${GBDT_LR}" \
  --seed 20260509

echo "[done] $(date)"
