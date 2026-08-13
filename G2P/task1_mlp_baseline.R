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
out_tag <- get_arg("--out-tag", "task1_mlp_full")
folds <- as.integer(get_arg("--folds", "5"))
repeats <- as.integer(get_arg("--repeats", "5"))
seed <- as.integer(get_arg("--seed", "20260509"))
pca_top_k <- as.integer(get_arg("--pca-top-k", "5000"))
pca_n <- as.integer(get_arg("--pca-n", "30"))
pca_hidden <- as.integer(get_arg("--pca-hidden", "16"))
direct_top_k <- as.integer(get_arg("--direct-top-k", "200"))
direct_hidden <- as.integer(get_arg("--direct-hidden", "16"))
epochs <- as.integer(get_arg("--epochs", "400"))
learning_rate <- as.numeric(get_arg("--learning-rate", "0.01"))
weight_decay <- as.numeric(get_arg("--weight-decay", "0.001"))
patience <- as.integer(get_arg("--patience", "12"))
selected_report_n <- as.integer(get_arg("--selected-report-n", "50"))

metadata_dir <- file.path(work_dir, "metadata")
result_dir <- file.path(work_dir, "results", "task1_mlp")
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

make_pca_features <- function(X_train, X_test, npc) {
  npc <- min(npc, nrow(X_train) - 1, ncol(X_train))
  if (npc < 1) stop("No PCA components available")
  sv <- svd(X_train, nu = npc, nv = npc)
  train_scores <- sv$u[, seq_len(npc), drop = FALSE] %*% diag(sv$d[seq_len(npc)], nrow = npc)
  test_scores <- X_test %*% sv$v[, seq_len(npc), drop = FALSE]
  st <- standardize_train_test(train_scores, test_scores)
  list(train = st$train, test = st$test, npc = ncol(st$train))
}

mlp_predict <- function(X, W1, b1, W2, b2) {
  H <- tanh(sweep(X %*% W1, 2, b1, "+"))
  as.numeric(H %*% W2 + b2)
}

fit_mlp_predict <- function(X_train, y_train, X_test, hidden, seed_value,
                            epochs = 400, learning_rate = 0.01,
                            weight_decay = 0.001, patience = 12) {
  set.seed(seed_value)
  n <- nrow(X_train)
  p <- ncol(X_train)
  if (p < 1) stop("No MLP input features")
  val_n <- max(12, round(0.2 * n))
  val_idx <- sort(sample(seq_len(n), val_n))
  fit_idx <- setdiff(seq_len(n), val_idx)

  y_mean <- mean(y_train)
  y_sd <- sd(y_train)
  if (!is.finite(y_sd) || y_sd == 0) y_sd <- 1
  y_scaled <- (y_train - y_mean) / y_sd

  X_fit <- X_train[fit_idx, , drop = FALSE]
  y_fit <- y_scaled[fit_idx]
  X_val <- X_train[val_idx, , drop = FALSE]
  y_val <- y_scaled[val_idx]

  W1 <- matrix(rnorm(p * hidden, sd = sqrt(1 / max(1, p))), nrow = p, ncol = hidden)
  b1 <- rep(0, hidden)
  W2 <- matrix(rnorm(hidden, sd = sqrt(1 / max(1, hidden))), nrow = hidden, ncol = 1)
  b2 <- 0

  mW1 <- W1 * 0; vW1 <- W1 * 0
  mb1 <- b1 * 0; vb1 <- b1 * 0
  mW2 <- W2 * 0; vW2 <- W2 * 0
  mb2 <- 0; vb2 <- 0
  beta1 <- 0.9; beta2 <- 0.999; eps <- 1e-8

  best <- list(W1 = W1, b1 = b1, W2 = W2, b2 = b2)
  best_loss <- Inf
  stale <- 0
  eval_every <- 10
  step <- 0

  for (epoch in seq_len(epochs)) {
    step <- step + 1
    H <- tanh(sweep(X_fit %*% W1, 2, b1, "+"))
    pred <- as.numeric(H %*% W2 + b2)
    dY <- 2 * (pred - y_fit) / length(y_fit)
    dW2 <- t(H) %*% matrix(dY, ncol = 1) + 2 * weight_decay * W2
    db2 <- sum(dY)
    dH <- matrix(dY, ncol = 1) %*% t(W2)
    dZ <- dH * (1 - H^2)
    dW1 <- t(X_fit) %*% dZ + 2 * weight_decay * W1
    db1 <- colSums(dZ)

    mW1 <- beta1 * mW1 + (1 - beta1) * dW1
    vW1 <- beta2 * vW1 + (1 - beta2) * (dW1^2)
    mb1 <- beta1 * mb1 + (1 - beta1) * db1
    vb1 <- beta2 * vb1 + (1 - beta2) * (db1^2)
    mW2 <- beta1 * mW2 + (1 - beta1) * dW2
    vW2 <- beta2 * vW2 + (1 - beta2) * (dW2^2)
    mb2 <- beta1 * mb2 + (1 - beta1) * db2
    vb2 <- beta2 * vb2 + (1 - beta2) * (db2^2)

    lr_t <- learning_rate * sqrt(1 - beta2^step) / (1 - beta1^step)
    W1 <- W1 - lr_t * mW1 / (sqrt(vW1) + eps)
    b1 <- b1 - lr_t * mb1 / (sqrt(vb1) + eps)
    W2 <- W2 - lr_t * mW2 / (sqrt(vW2) + eps)
    b2 <- b2 - lr_t * mb2 / (sqrt(vb2) + eps)

    if (epoch %% eval_every == 0) {
      val_pred <- mlp_predict(X_val, W1, b1, W2, b2)
      val_loss <- mean((val_pred - y_val)^2) + weight_decay * (sum(W1^2) + sum(W2^2))
      if (is.finite(val_loss) && val_loss < best_loss - 1e-6) {
        best_loss <- val_loss
        best <- list(W1 = W1, b1 = b1, W2 = W2, b2 = b2)
        stale <- 0
      } else {
        stale <- stale + 1
      }
      if (stale >= patience) break
    }
  }

  pred_scaled <- mlp_predict(X_test, best$W1, best$b1, best$W2, best$b2)
  list(pred = pred_scaled * y_sd + y_mean, best_val_loss = best_loss, epochs_used = epoch)
}

add_pred_rows <- function(rows, trait, model, rep_i, fold, samples, idx, obs, pred, top_k, extra) {
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
      extra = extra,
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

      screen <- screen_markers(M, y, train_idx, pca_top_k)
      if (selected_report_n > 0) {
        report_n <- min(selected_report_n, length(screen$idx))
        for (rank_i in seq_len(report_n)) {
          m_idx <- screen$idx[[rank_i]]
          selected_rows[[selected_i]] <- data.frame(
            trait = trait,
            model = "PCA_MLP_and_Screened_MLP_screening",
            repeat_id = rep_i,
            fold = fold,
            rank = rank_i,
            marker = markers$marker[[m_idx]],
            score = screen$score[[rank_i]],
            stringsAsFactors = FALSE
          )
          selected_i <- selected_i + 1
        }
      }

      # PCA + small MLP.
      st_pca <- standardize_train_test(M[train_idx, screen$idx, drop = FALSE], M[test_idx, screen$idx, drop = FALSE])
      pc <- make_pca_features(st_pca$train, st_pca$test, pca_n)
      fit <- fit_mlp_predict(
        pc$train, y_train, pc$test, hidden = pca_hidden,
        seed_value = seed + trait_i * 10000 + rep_i * 100 + fold,
        epochs = epochs, learning_rate = learning_rate,
        weight_decay = weight_decay, patience = patience
      )
      model_name <- paste0("PCA_MLP_screened_top", pca_top_k, "_pc", pc$npc, "_hidden", pca_hidden)
      pred <- fit$pred
      fold_rows[[fold_row_i]] <- data.frame(
        trait = trait, model = model_name, repeat_id = rep_i, fold = fold,
        n_train = length(train_idx), n_test = length(test_idx), top_k = pca_top_k,
        extra = paste0("pc=", pc$npc, ";epochs=", fit$epochs_used, ";val_loss=", signif(fit$best_val_loss, 5)),
        pearson = pearson(y_test, pred), rmse = rmse(y_test, pred), mae = mae(y_test, pred),
        stringsAsFactors = FALSE
      )
      fold_row_i <- fold_row_i + 1
      pred_rows <- add_pred_rows(pred_rows, trait, model_name, rep_i, fold, samples, test_idx, y_test, pred, pca_top_k, paste0("pc=", pc$npc))

      # Direct screened-feature MLP.
      direct_idx <- screen$idx[seq_len(min(direct_top_k, length(screen$idx)))]
      st_direct <- standardize_train_test(M[train_idx, direct_idx, drop = FALSE], M[test_idx, direct_idx, drop = FALSE])
      fit <- fit_mlp_predict(
        st_direct$train, y_train, st_direct$test, hidden = direct_hidden,
        seed_value = seed + trait_i * 20000 + rep_i * 100 + fold,
        epochs = epochs, learning_rate = learning_rate,
        weight_decay = weight_decay, patience = patience
      )
      model_name <- paste0("Screened_MLP_top", length(direct_idx), "_hidden", direct_hidden)
      pred <- fit$pred
      fold_rows[[fold_row_i]] <- data.frame(
        trait = trait, model = model_name, repeat_id = rep_i, fold = fold,
        n_train = length(train_idx), n_test = length(test_idx), top_k = length(direct_idx),
        extra = paste0("epochs=", fit$epochs_used, ";val_loss=", signif(fit$best_val_loss, 5)),
        pearson = pearson(y_test, pred), rmse = rmse(y_test, pred), mae = mae(y_test, pred),
        stringsAsFactors = FALSE
      )
      fold_row_i <- fold_row_i + 1
      pred_rows <- add_pred_rows(pred_rows, trait, model_name, rep_i, fold, samples, test_idx, y_test, pred, length(direct_idx), "")
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

cat("Wrote MLP outputs for tag=", out_tag, "\n", sep = "")
