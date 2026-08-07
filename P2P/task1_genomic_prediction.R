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
folds <- as.integer(get_arg("--folds", "5"))
repeats <- as.integer(get_arg("--repeats", "5"))
seed <- as.integer(get_arg("--seed", "20260509"))

metadata_dir <- file.path(work_dir, "metadata")
result_dir <- file.path(work_dir, "results", "task1_genomic_prediction")
tmp_dir <- file.path(work_dir, "tmp", "task1_genomic_prediction")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

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

predict_kernel <- function(K, y, train_idx, test_idx, lambda) {
  y_mean <- mean(y[train_idx])
  y_center <- y[train_idx] - y_mean
  A <- K[train_idx, train_idx, drop = FALSE] + diag(lambda, length(train_idx))
  alpha <- solve(A, y_center)
  pred <- K[test_idx, train_idx, drop = FALSE] %*% alpha + y_mean
  as.numeric(pred)
}

inner_select <- function(K, y, train_idx, lambdas, k_inner = 4, seed_value = 1) {
  if (length(train_idx) < k_inner) k_inner <- max(2, length(train_idx))
  fold_id <- make_folds(length(train_idx), k_inner, seed_value)
  scores <- numeric(length(lambdas))
  for (li in seq_along(lambdas)) {
    lambda <- lambdas[[li]]
    pred_all <- rep(NA_real_, length(train_idx))
    for (fold in seq_len(k_inner)) {
      inner_test_local <- which(fold_id == fold)
      inner_train_local <- which(fold_id != fold)
      inner_test <- train_idx[inner_test_local]
      inner_train <- train_idx[inner_train_local]
      pred_all[inner_test_local] <- predict_kernel(K, y, inner_train, inner_test, lambda)
    }
    scores[[li]] <- rmse(y[train_idx], pred_all)
  }
  lambdas[[which.min(scores)]]
}

run_cv <- function(K, y, samples, trait, model, lambdas, repeats, folds, seed_offset) {
  rows_pred <- list()
  rows_fold <- list()
  row_i <- 1
  fold_i <- 1
  n <- length(y)
  for (rep_i in seq_len(repeats)) {
    fold_id <- make_folds(n, folds, seed_offset + rep_i)
    for (fold in seq_len(folds)) {
      test_idx <- which(fold_id == fold)
      train_idx <- which(fold_id != fold)
      best_lambda <- inner_select(K, y, train_idx, lambdas, k_inner = min(4, folds - 1), seed_value = seed_offset + 1000 * rep_i + fold)
      pred <- predict_kernel(K, y, train_idx, test_idx, best_lambda)
      rows_fold[[fold_i]] <- data.frame(
        trait = trait,
        model = model,
        repeat_id = rep_i,
        fold = fold,
        n_train = length(train_idx),
        n_test = length(test_idx),
        lambda = best_lambda,
        pearson = pearson(y[test_idx], pred),
        rmse = rmse(y[test_idx], pred),
        mae = mae(y[test_idx], pred),
        stringsAsFactors = FALSE
      )
      fold_i <- fold_i + 1
      for (j in seq_along(test_idx)) {
        idx <- test_idx[[j]]
        rows_pred[[row_i]] <- data.frame(
          trait = trait,
          model = model,
          repeat_id = rep_i,
          fold = fold,
          sample_id = samples$sample_id[[idx]],
          observed = y[[idx]],
          predicted = pred[[j]],
          lambda = best_lambda,
          stringsAsFactors = FALSE
        )
        row_i <- row_i + 1
      }
    }
  }
  list(pred = do.call(rbind, rows_pred), fold = do.call(rbind, rows_fold))
}

manifest <- read_manifest(file.path(metadata_dir, paste0(tag, ".manifest.tsv")))
n_samples <- as.integer(manifest$n_samples)
n_markers <- as.integer(manifest$n_markers)
matrix_file <- manifest$matrix_file
samples_file <- manifest$samples_file

cat("Loading matrix:", matrix_file, "\n")
cat("n_samples=", n_samples, " n_markers=", n_markers, "\n", sep = "")
x <- readBin(matrix_file, what = "numeric", n = n_samples * n_markers, size = 4)
M <- matrix(x, nrow = n_samples, ncol = n_markers, byrow = TRUE)
rm(x)
gc()

samples <- read.table(samples_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
traits <- c("GZZTF", "TCD", "ZZZTF","HYTF","ZGD","SCD","ZLCD","ZXWF")

cat("Centering and scaling markers\n")
obs <- !is.na(M)
col_n <- colSums(obs)
mu <- colSums(M, na.rm = TRUE) / pmax(col_n, 1)
M <- sweep(M, 2, mu, "-")
M[!obs] <- 0
ss <- colSums(M^2)
sdv <- sqrt(ss / pmax(col_n - 1, 1))
keep <- is.finite(sdv) & sdv > 0 & col_n >= 3
cat("markers_before=", n_markers, " markers_after_sd_filter=", sum(keep), "\n", sep = "")
M <- sweep(M[, keep, drop = FALSE], 2, sdv[keep], "/")
M[!is.finite(M)] <- 0
rm(obs, col_n, mu, ss, sdv)
gc()

cat("Computing linear genomic relationship matrix\n")
K_linear <- tcrossprod(M) / ncol(M)
diag_mean <- mean(diag(K_linear))
if (is.finite(diag_mean) && diag_mean > 0) K_linear <- K_linear / diag_mean

d <- diag(K_linear)
D2 <- outer(d, d, "+") - 2 * K_linear
D2[D2 < 0] <- 0
median_d2 <- median(D2[upper.tri(D2)], na.rm = TRUE)
if (!is.finite(median_d2) || median_d2 <= 0) median_d2 <- 1
K_rbf <- exp(-D2 / median_d2)

write.table(K_linear, file.path(tmp_dir, paste0(tag, ".linear_G.tsv")), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(K_rbf, file.path(tmp_dir, paste0(tag, ".rkhs_rbf_K.tsv")), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

kernel_summary <- data.frame(
  tag = tag,
  n_samples = n_samples,
  markers_before = n_markers,
  markers_after_sd_filter = ncol(M),
  linear_diag_min = min(diag(K_linear)),
  linear_diag_mean = mean(diag(K_linear)),
  linear_diag_max = max(diag(K_linear)),
  rbf_median_d2 = median_d2,
  stringsAsFactors = FALSE
)
write.table(kernel_summary, file.path(result_dir, paste0(tag, ".kernel_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

lambda_grid <- c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100)
all_pred <- list()
all_fold <- list()
block <- 1

for (trait in traits) {
  y <- as.numeric(samples[[trait]])
  cat("Running CV for trait=", trait, "\n", sep = "")
  out_linear <- run_cv(K_linear, y, samples, trait, "RRBLUP_GBLUP_BayesianRidge_linear", lambda_grid, repeats, folds, seed + block * 100)
  all_pred[[length(all_pred) + 1]] <- out_linear$pred
  all_fold[[length(all_fold) + 1]] <- out_linear$fold
  block <- block + 1

  out_rbf <- run_cv(K_rbf, y, samples, trait, "RKHS_GBLUP_RBF", lambda_grid, repeats, folds, seed + block * 100)
  all_pred[[length(all_pred) + 1]] <- out_rbf$pred
  all_fold[[length(all_fold) + 1]] <- out_rbf$fold
  block <- block + 1
}

pred_tab <- do.call(rbind, all_pred)
fold_tab <- do.call(rbind, all_fold)

summary_rows <- list()
si <- 1
for (trait in unique(pred_tab$trait)) {
  for (model in unique(pred_tab$model)) {
    idx <- pred_tab$trait == trait & pred_tab$model == model
    summary_rows[[si]] <- data.frame(
      tag = tag,
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

write.table(pred_tab, file.path(result_dir, paste0(tag, ".cv_predictions.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fold_tab, file.path(result_dir, paste0(tag, ".cv_fold_metrics.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_tab, file.path(result_dir, paste0(tag, ".cv_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote genomic prediction outputs for tag=", tag, "\n", sep = "")
