#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop(paste("Missing value for", flag))
  args[[idx + 1]]
}

work_dir <- get_arg("--work-dir", "/public1/home/stu_gaoruiting/project/Ssin/HiFiasm/merge/2024-8-28/pop/Ai-input")
tag <- get_arg("--tag", "task1_gp_full")
out_tag <- get_arg("--out-tag", "task1_tree_models_full")
folds <- as.integer(get_arg("--folds", "5"))
repeats <- as.integer(get_arg("--repeats", "5"))
seed <- as.integer(get_arg("--seed", "20260509"))
rf_top_k <- as.integer(get_arg("--rf-top-k", "1000"))
rf_trees <- as.integer(get_arg("--rf-trees", "80"))
rf_maxdepth <- as.integer(get_arg("--rf-maxdepth", "5"))
gbdt_top_k <- as.integer(get_arg("--gbdt-top-k", "500"))
gbdt_rounds <- as.integer(get_arg("--gbdt-rounds", "80"))
gbdt_lr <- as.numeric(get_arg("--gbdt-lr", "0.05"))
gbdt_maxdepth <- as.integer(get_arg("--gbdt-maxdepth", "2"))
selected_report_n <- as.integer(get_arg("--selected-report-n", "50"))

suppressPackageStartupMessages(library(rpart))

metadata_dir <- file.path(work_dir, "metadata")
result_dir <- file.path(work_dir, "results", "task1_tree_models")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

read_manifest <- function(path) {
  tab <- read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  out <- as.list(tab$value)
  names(out) <- tab$key
  out
}

pearson <- function(obs, pred) {
  if (length(obs) < 3 || sd(obs) == 0 || sd(pred) == 0) return(NA_real_)
  as.numeric(cor(obs, pred))
}
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))
mae <- function(obs, pred) mean(abs(obs - pred))

make_folds <- function(n, k, seed_value) {
  set.seed(seed_value)
  order <- sample(seq_len(n))
  fold_id <- integer(n)
  fold_id[order] <- rep(seq_len(k), length.out = n)
  fold_id
}

screen_markers <- function(M, y, train_idx, top_k) {
  y_train <- y[train_idx]
  yc <- y_train - mean(y_train)
  Xtr <- M[train_idx, , drop = FALSE]
  xy <- as.numeric(crossprod(Xtr, yc))
  x2 <- colSums(Xtr^2)
  y2 <- sum(yc^2)
  score <- abs(xy) / sqrt(pmax(x2 * y2, .Machine$double.eps))
  score[!is.finite(score)] <- 0
  k <- min(top_k, length(score))
  idx <- order(score, decreasing = TRUE)[seq_len(k)]
  list(idx = idx, score = score[idx])
}

standardize_train_test <- function(X_train, X_test) {
  mu <- colMeans(X_train)
  Xc <- sweep(X_train, 2, mu, "-")
  sdv <- sqrt(colSums(Xc^2) / pmax(nrow(Xc) - 1, 1))
  keep <- is.finite(sdv) & sdv > 0
  Xtr <- sweep(Xc[, keep, drop = FALSE], 2, sdv[keep], "/")
  Xte <- sweep(sweep(X_test[, keep, drop = FALSE], 2, mu[keep], "-"), 2, sdv[keep], "/")
  Xtr[!is.finite(Xtr)] <- 0
  Xte[!is.finite(Xte)] <- 0
  list(train = Xtr, test = Xte, keep = keep)
}

fit_rf_like <- function(X_train, y_train, X_test, n_trees, maxdepth, seed_value) {
  set.seed(seed_value)
  n <- nrow(X_train)
  p <- ncol(X_train)
  mtry <- max(1, floor(sqrt(p)))
  pred <- rep(0, nrow(X_test))
  ctrl <- rpart.control(
    maxdepth = maxdepth,
    minsplit = max(8, floor(n * 0.08)),
    minbucket = max(3, floor(n * 0.03)),
    cp = 0.0005,
    xval = 0
  )
  for (tree_i in seq_len(n_trees)) {
    boot <- sample(seq_len(n), size = n, replace = TRUE)
    feat <- sample(seq_len(p), size = mtry, replace = FALSE)
    train_df <- data.frame(y = y_train[boot], X_train[boot, feat, drop = FALSE], check.names = FALSE)
    names(train_df) <- c("y", paste0("V", seq_along(feat)))
    fit <- rpart(y ~ ., data = train_df, method = "anova", control = ctrl)
    test_df <- data.frame(X_test[, feat, drop = FALSE], check.names = FALSE)
    names(test_df) <- paste0("V", seq_along(feat))
    pred <- pred + as.numeric(predict(fit, newdata = test_df))
  }
  pred / n_trees
}

fit_gbdt_like <- function(X_train, y_train, X_test, n_rounds, learning_rate, maxdepth, seed_value) {
  set.seed(seed_value)
  n <- nrow(X_train)
  p <- ncol(X_train)
  colsample <- max(1, min(p, floor(sqrt(p) * 3)))
  pred_train <- rep(mean(y_train), n)
  pred_test <- rep(mean(y_train), nrow(X_test))
  ctrl <- rpart.control(
    maxdepth = maxdepth,
    minsplit = max(8, floor(n * 0.08)),
    minbucket = max(3, floor(n * 0.03)),
    cp = 0,
    xval = 0
  )
  for (round_i in seq_len(n_rounds)) {
    residual <- y_train - pred_train
    feat <- sample(seq_len(p), size = colsample, replace = FALSE)
    train_df <- data.frame(y = residual, X_train[, feat, drop = FALSE], check.names = FALSE)
    names(train_df) <- c("y", paste0("V", seq_along(feat)))
    fit <- rpart(y ~ ., data = train_df, method = "anova", control = ctrl)
    train_pred_df <- data.frame(X_train[, feat, drop = FALSE], check.names = FALSE)
    test_pred_df <- data.frame(X_test[, feat, drop = FALSE], check.names = FALSE)
    names(train_pred_df) <- paste0("V", seq_along(feat))
    names(test_pred_df) <- paste0("V", seq_along(feat))
    pred_train <- pred_train + learning_rate * as.numeric(predict(fit, newdata = train_pred_df))
    pred_test <- pred_test + learning_rate * as.numeric(predict(fit, newdata = test_pred_df))
  }
  pred_test
}

add_pred_rows <- function(rows, trait, model, rep_i, fold, samples, idx, obs, pred, top_k) {
  start <- length(rows) + 1
  for (j in seq_along(idx)) {
    i <- idx[[j]]
    rows[[start + j - 1]] <- data.frame(
      trait = trait,
      model = model,
      repeat_id = rep_i,
      fold = fold,
      sample_id = samples$sample_id[[i]],
      observed = obs[[j]],
      predicted = pred[[j]],
      top_k = top_k,
      stringsAsFactors = FALSE
    )
  }
  rows
}

manifest <- read_manifest(file.path(metadata_dir, paste0(tag, ".manifest.tsv")))
n_samples <- as.integer(manifest$n_samples)
n_markers <- as.integer(manifest$n_markers)
matrix_file <- manifest$matrix_file
samples_file <- manifest$samples_file
markers_file <- manifest$markers_file

cat("Loading matrix:", matrix_file, "\n")
cat("n_samples=", n_samples, " n_markers=", n_markers, "\n", sep = "")
x <- readBin(matrix_file, what = "numeric", n = n_samples * n_markers, size = 4)
M <- matrix(x, nrow = n_samples, ncol = n_markers, byrow = TRUE)
rm(x)
gc()

samples <- read.table(samples_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
markers <- read.table(markers_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")

cat("Global unsupervised imputation, centering, scaling\n")
obs <- !is.na(M)
col_n <- colSums(obs)
mu <- colSums(M, na.rm = TRUE) / pmax(col_n, 1)
M <- sweep(M, 2, mu, "-")
M[!obs] <- 0
ss <- colSums(M^2)
sdv <- sqrt(ss / pmax(col_n - 1, 1))
keep <- is.finite(sdv) & sdv > 0 & col_n >= 3
M <- sweep(M[, keep, drop = FALSE], 2, sdv[keep], "/")
M[!is.finite(M)] <- 0
markers <- markers[keep, , drop = FALSE]
cat("markers_after_sd_filter=", ncol(M), "\n", sep = "")
rm(obs, col_n, mu, ss, sdv, keep)
gc()

traits <- c("GZZTF", "TCD", "ZZZTF","HYTF","ZGD","SCD","ZLCD","ZXWF")
pred_rows <- list()
fold_rows <- list()
selected_rows <- list()
fold_row_i <- 1
selected_i <- 1

for (trait_i in seq_along(traits)) {
  trait <- traits[[trait_i]]
  y <- as.numeric(samples[[trait]])
  n <- length(y)
  cat("Trait=", trait, "\n", sep = "")
  for (rep_i in seq_len(repeats)) {
    fold_id <- make_folds(n, folds, seed + 1000 * trait_i + rep_i)
    for (fold in seq_len(folds)) {
      test_idx <- which(fold_id == fold)
      train_idx <- which(fold_id != fold)
      y_train <- y[train_idx]
      y_test <- y[test_idx]

      screen_rf <- screen_markers(M, y, train_idx, rf_top_k)
      if (selected_report_n > 0) {
        report_n <- min(selected_report_n, length(screen_rf$idx))
        for (rank_i in seq_len(report_n)) {
          m_idx <- screen_rf$idx[[rank_i]]
          selected_rows[[selected_i]] <- data.frame(
            trait = trait,
            model = "RandomForest_like_bagged_rpart",
            repeat_id = rep_i,
            fold = fold,
            rank = rank_i,
            marker = markers$marker[[m_idx]],
            score = screen_rf$score[[rank_i]],
            stringsAsFactors = FALSE
          )
          selected_i <- selected_i + 1
        }
      }
      st <- standardize_train_test(M[train_idx, screen_rf$idx, drop = FALSE], M[test_idx, screen_rf$idx, drop = FALSE])
      pred <- fit_rf_like(st$train, y_train, st$test, rf_trees, rf_maxdepth, seed + trait_i * 10000 + rep_i * 100 + fold)
      model_name <- paste0("RandomForest_like_bagged_rpart_top", rf_top_k, "_trees", rf_trees)
      fold_rows[[fold_row_i]] <- data.frame(
        trait = trait, model = model_name, repeat_id = rep_i, fold = fold,
        n_train = length(train_idx), n_test = length(test_idx), top_k = rf_top_k,
        pearson = pearson(y_test, pred), rmse = rmse(y_test, pred), mae = mae(y_test, pred),
        stringsAsFactors = FALSE
      )
      fold_row_i <- fold_row_i + 1
      pred_rows <- add_pred_rows(pred_rows, trait, model_name, rep_i, fold, samples, test_idx, y_test, pred, rf_top_k)

      screen_gbdt <- if (gbdt_top_k == rf_top_k) screen_rf else screen_markers(M, y, train_idx, gbdt_top_k)
      if (selected_report_n > 0) {
        report_n <- min(selected_report_n, length(screen_gbdt$idx))
        for (rank_i in seq_len(report_n)) {
          m_idx <- screen_gbdt$idx[[rank_i]]
          selected_rows[[selected_i]] <- data.frame(
            trait = trait,
            model = "GBDT_XGBoost_like_rpart",
            repeat_id = rep_i,
            fold = fold,
            rank = rank_i,
            marker = markers$marker[[m_idx]],
            score = screen_gbdt$score[[rank_i]],
            stringsAsFactors = FALSE
          )
          selected_i <- selected_i + 1
        }
      }
      st <- standardize_train_test(M[train_idx, screen_gbdt$idx, drop = FALSE], M[test_idx, screen_gbdt$idx, drop = FALSE])
      pred <- fit_gbdt_like(st$train, y_train, st$test, gbdt_rounds, gbdt_lr, gbdt_maxdepth, seed + trait_i * 20000 + rep_i * 100 + fold)
      model_name <- paste0("GBDT_XGBoost_like_rpart_top", gbdt_top_k, "_rounds", gbdt_rounds)
      fold_rows[[fold_row_i]] <- data.frame(
        trait = trait, model = model_name, repeat_id = rep_i, fold = fold,
        n_train = length(train_idx), n_test = length(test_idx), top_k = gbdt_top_k,
        pearson = pearson(y_test, pred), rmse = rmse(y_test, pred), mae = mae(y_test, pred),
        stringsAsFactors = FALSE
      )
      fold_row_i <- fold_row_i + 1
      pred_rows <- add_pred_rows(pred_rows, trait, model_name, rep_i, fold, samples, test_idx, y_test, pred, gbdt_top_k)
      cat("trait=", trait, " repeat=", rep_i, " fold=", fold, " done\n", sep = "")
    }
  }
}

pred_tab <- do.call(rbind, pred_rows)
fold_tab <- do.call(rbind, fold_rows)
selected_tab <- if (length(selected_rows)) do.call(rbind, selected_rows) else data.frame()

summary_rows <- list()
si <- 1
for (trait in unique(pred_tab$trait)) {
  for (model in unique(pred_tab$model)) {
    idx <- pred_tab$trait == trait & pred_tab$model == model
    summary_rows[[si]] <- data.frame(
      tag = out_tag,
      source_matrix_tag = tag,
      trait = trait,
      model = model,
      n_predictions = sum(idx),
      repeats = repeats,
      folds = folds,
      pearson = pearson(pred_tab$observed[idx], pred_tab$predicted[idx]),
      rmse = rmse(pred_tab$observed[idx], pred_tab$predicted[idx]),
      mae = mae(pred_tab$observed[idx], pred_tab$predicted[idx]),
      stringsAsFactors = FALSE
    )
    si <- si + 1
  }
}
summary_tab <- do.call(rbind, summary_rows)

write.table(pred_tab, file.path(result_dir, paste0(out_tag, ".cv_predictions.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fold_tab, file.path(result_dir, paste0(out_tag, ".cv_fold_metrics.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_tab, file.path(result_dir, paste0(out_tag, ".cv_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(selected_tab, file.path(result_dir, paste0(out_tag, ".selected_markers_top", selected_report_n, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote tree model outputs for tag=", out_tag, "\n", sep = "")
