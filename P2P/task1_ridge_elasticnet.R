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
out_tag <- get_arg("--out-tag", "task1_ridge_elasticnet_full")
folds <- as.integer(get_arg("--folds", "5"))
repeats <- as.integer(get_arg("--repeats", "5"))
seed <- as.integer(get_arg("--seed", "20260509"))
ridge_top_k <- as.integer(get_arg("--ridge-top-k", "5000"))
enet_top_k <- as.integer(get_arg("--enet-top-k", "1000"))
enet_alpha <- as.numeric(get_arg("--enet-alpha", "0.5"))
selected_report_n <- as.integer(get_arg("--selected-report-n", "50"))

metadata_dir <- file.path(work_dir, "metadata")
result_dir <- file.path(work_dir, "results", "task1_ridge_elasticnet")
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

soft_threshold <- function(z, gamma) {
  sign(z) * pmax(abs(z) - gamma, 0)
}

standardize_train_test <- function(X_train, X_test) {
  mu <- colMeans(X_train)
  Xc <- sweep(X_train, 2, mu, "-")
  sdv <- sqrt(colSums(Xc^2) / pmax(nrow(Xc) - 1, 1))
  keep <- is.finite(sdv) & sdv > 0
  Xtr <- sweep(Xc[, keep, drop = FALSE], 2, sdv[keep], "/")
  Xte <- sweep(sweep(X_test[, keep, drop = FALSE], 2, mu[keep], "-"), 2, sdv[keep], "/")
  list(train = Xtr, test = Xte, keep = keep)
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

predict_ridge_kernel <- function(K_train, K_test_train, y_train, train_local, test_local, lambda) {
  y_mean <- mean(y_train[train_local])
  yc <- y_train[train_local] - y_mean
  A <- K_train[train_local, train_local, drop = FALSE] + diag(lambda, length(train_local))
  alpha <- solve(A, yc)
  as.numeric(K_test_train[test_local, train_local, drop = FALSE] %*% alpha + y_mean)
}

select_ridge_lambda <- function(K, y, lambdas, seed_value) {
  n <- length(y)
  inner_k <- min(4, max(2, n - 1))
  fold_id <- make_folds(n, inner_k, seed_value)
  scores <- rep(NA_real_, length(lambdas))
  for (li in seq_along(lambdas)) {
    pred <- rep(NA_real_, n)
    for (fold in seq_len(inner_k)) {
      te <- which(fold_id == fold)
      tr <- which(fold_id != fold)
      pred[te] <- predict_ridge_kernel(K, K, y, tr, te, lambdas[[li]])
    }
    scores[[li]] <- rmse(y, pred)
  }
  lambdas[[which.min(scores)]]
}

fit_enet_fista <- function(X, y, lambda, alpha = 0.5, max_iter = 350, tol = 1e-5) {
  n <- nrow(X)
  p <- ncol(X)
  beta <- rep(0, p)
  z <- beta
  tval <- 1
  # Conservative Lipschitz bound; cheap and stable for small selected matrices.
  L0 <- sum(X^2) / n
  L <- L0 + lambda * (1 - alpha)
  if (!is.finite(L) || L <= 0) L <- 1
  for (iter in seq_len(max_iter)) {
    beta_old <- beta
    grad <- as.numeric(crossprod(X, (X %*% z - y))) / n + lambda * (1 - alpha) * z
    beta <- soft_threshold(z - grad / L, lambda * alpha / L)
    t_new <- (1 + sqrt(1 + 4 * tval^2)) / 2
    z <- beta + ((tval - 1) / t_new) * (beta - beta_old)
    if (max(abs(beta - beta_old)) < tol) break
    tval <- t_new
  }
  beta
}

select_enet_lambda <- function(X, y, alpha, seed_value) {
  n <- length(y)
  inner_k <- min(4, max(2, n - 1))
  fold_id <- make_folds(n, inner_k, seed_value)
  y_center <- y - mean(y)
  lambda_max <- max(abs(as.numeric(crossprod(X, y_center)))) / (n * max(alpha, 1e-6))
  if (!is.finite(lambda_max) || lambda_max <= 0) lambda_max <- 1
  lambdas <- lambda_max * c(0.3, 0.1, 0.03, 0.01, 0.003)
  scores <- rep(NA_real_, length(lambdas))
  for (li in seq_along(lambdas)) {
    pred <- rep(NA_real_, n)
    for (fold in seq_len(inner_k)) {
      te <- which(fold_id == fold)
      tr <- which(fold_id != fold)
      y_mean <- mean(y[tr])
      beta <- fit_enet_fista(X[tr, , drop = FALSE], y[tr] - y_mean, lambdas[[li]], alpha = alpha)
      pred[te] <- as.numeric(X[te, , drop = FALSE] %*% beta + y_mean)
    }
    scores[[li]] <- rmse(y, pred)
  }
  lambdas[[which.min(scores)]]
}

add_pred_rows <- function(rows, trait, model, rep_i, fold, samples, idx, obs, pred, lambda, top_k) {
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
      lambda = lambda,
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
ridge_lambdas <- c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100)

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

      # Ridge with screened top markers.
      screen_ridge <- screen_markers(M, y, train_idx, ridge_top_k)
      if (selected_report_n > 0) {
        report_n <- min(selected_report_n, length(screen_ridge$idx))
        for (rank_i in seq_len(report_n)) {
          m_idx <- screen_ridge$idx[[rank_i]]
          selected_rows[[selected_i]] <- data.frame(
            trait = trait,
            model = "Ridge_screened",
            repeat_id = rep_i,
            fold = fold,
            rank = rank_i,
            marker = markers$marker[[m_idx]],
            score = screen_ridge$score[[rank_i]],
            stringsAsFactors = FALSE
          )
          selected_i <- selected_i + 1
        }
      }
      Xtr0 <- M[train_idx, screen_ridge$idx, drop = FALSE]
      Xte0 <- M[test_idx, screen_ridge$idx, drop = FALSE]
      st <- standardize_train_test(Xtr0, Xte0)
      Xtr <- st$train
      Xte <- st$test
      Ktr <- tcrossprod(Xtr) / ncol(Xtr)
      Kte <- Xte %*% t(Xtr) / ncol(Xtr)
      best_lambda <- select_ridge_lambda(Ktr, y_train, ridge_lambdas, seed + rep_i * 100 + fold)
      pred <- predict_ridge_kernel(Ktr, Kte, y_train, seq_along(train_idx), seq_along(test_idx), best_lambda)
      model_name <- paste0("Ridge_screened_top", ridge_top_k)
      fold_rows[[fold_row_i]] <- data.frame(
        trait = trait, model = model_name, repeat_id = rep_i, fold = fold,
        n_train = length(train_idx), n_test = length(test_idx), lambda = best_lambda,
        top_k = ridge_top_k, pearson = pearson(y_test, pred), rmse = rmse(y_test, pred), mae = mae(y_test, pred),
        stringsAsFactors = FALSE
      )
      fold_row_i <- fold_row_i + 1
      pred_rows <- add_pred_rows(pred_rows, trait, model_name, rep_i, fold, samples, test_idx, y_test, pred, best_lambda, ridge_top_k)

      # ElasticNet with a smaller screened set.
      screen_enet <- if (enet_top_k == ridge_top_k) screen_ridge else screen_markers(M, y, train_idx, enet_top_k)
      if (selected_report_n > 0) {
        report_n <- min(selected_report_n, length(screen_enet$idx))
        for (rank_i in seq_len(report_n)) {
          m_idx <- screen_enet$idx[[rank_i]]
          selected_rows[[selected_i]] <- data.frame(
            trait = trait,
            model = paste0("ElasticNet_alpha", enet_alpha),
            repeat_id = rep_i,
            fold = fold,
            rank = rank_i,
            marker = markers$marker[[m_idx]],
            score = screen_enet$score[[rank_i]],
            stringsAsFactors = FALSE
          )
          selected_i <- selected_i + 1
        }
      }
      Xtr0 <- M[train_idx, screen_enet$idx, drop = FALSE]
      Xte0 <- M[test_idx, screen_enet$idx, drop = FALSE]
      st <- standardize_train_test(Xtr0, Xte0)
      Xtr <- st$train
      Xte <- st$test
      best_lambda <- select_enet_lambda(Xtr, y_train, enet_alpha, seed + rep_i * 1000 + fold)
      y_mean <- mean(y_train)
      beta <- fit_enet_fista(Xtr, y_train - y_mean, best_lambda, alpha = enet_alpha)
      pred <- as.numeric(Xte %*% beta + y_mean)
      model_name <- paste0("ElasticNet_alpha", enet_alpha, "_screened_top", enet_top_k)
      fold_rows[[fold_row_i]] <- data.frame(
        trait = trait, model = model_name, repeat_id = rep_i, fold = fold,
        n_train = length(train_idx), n_test = length(test_idx), lambda = best_lambda,
        top_k = enet_top_k, pearson = pearson(y_test, pred), rmse = rmse(y_test, pred), mae = mae(y_test, pred),
        stringsAsFactors = FALSE
      )
      fold_row_i <- fold_row_i + 1
      pred_rows <- add_pred_rows(pred_rows, trait, model_name, rep_i, fold, samples, test_idx, y_test, pred, best_lambda, enet_top_k)
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

cat("Wrote Ridge/ElasticNet outputs for tag=", out_tag, "\n", sep = "")
