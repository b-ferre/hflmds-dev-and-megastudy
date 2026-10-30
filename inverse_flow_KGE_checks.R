## investigate how inverse KGE looks for my models based on thw knowledge
## that this does better to emphasize low flows

library(stringr)
library(Matrix)
library(hydroGOF)

## helper functions directly from predictive run
get_design_matrix <- function(x, y, first_day_index,
                              m, max_lag, ars = 0,
                              remove_na_rows = TRUE) {

  if (length(x) != length(y) || length(y) < max_lag + ars) {
    stop("x and y must have equal lengths and enough observations for the lags.")              ## nolint
  }

  ## lags vary from 0 to ars + (max_lag - 1), so the first y we can predict is
  ## at the index ars + max_lag
  target_vec <- y[(max_lag + ars) : length(y)]
  nrow <- length(target_vec)

  ## first, build the first (max_lag + ars) * (m + 1) columns using x values
  ncol_x <- (max_lag + ars) * (m + 1)
  i_x <- sapply(seq_len(nrow), function(q) rep(q, (max_lag + ars)))
  j_x <- sapply(seq_len(nrow), function(q) {
    qstar <- (((q - 1) + (max_lag + ars) + (first_day_index - 1) - 1) %% (m + 1)) + 1           ## nolint
    seq(from = ((((max_lag + ars) * (qstar - 1))) %% ncol_x) + 1,
        to = ((((max_lag + ars) * (qstar)) - 1) %% ncol_x) + 1,
        by = 1)
  })
  x_x <- sapply(seq_len(nrow), function(q) {
    # qstar <- (((q - 1) + (max_lag + ars) + (first_day_index - 1) - 1) %% (m + 1)) + 1           ## nolint
    # max_lag_i <- dvec[qstar] + ars    ## nolint
    c(rev(x[q:(q + (max_lag + ars) - 1)])) # [seq_len(max_lag_i)], rep(0, (max_lag + ars) - max_lag_i))    ## nolint
  })

  ## then, if ars > 0, build additional y columns and append them
  i_y <- c()
  j_y <- c()
  x_y <- c()
  ncol_y <- ars
  if (ars > 0) {
    i_y <- sapply(seq_len(nrow), function(q) rep(q, ncol_y))
    j_y <- rep(seq(from = ncol_x + 1, to = ncol_x + ars, by = 1), nrow)
    x_y <- sapply(seq_len(nrow), function(q) {
      rev(y[(max_lag + q - 1):(max_lag + q + ars - 2)])
    })
  }

  mat <- sparseMatrix(
    i = c(i_x, i_y),
    j = c(j_x, j_y),
    x = c(x_x, x_y),
    dims = c(nrow, ncol_x + ncol_y)
  )
  if (remove_na_rows) {
    bad_x <- is.na(x_x) | !is.finite(x_x)
    bad_y <- if (ars > 0) is.na(x_y) | !is.finite(x_y) else logical(0)

    na_mat_rows <- unique(c(i_x[bad_x], i_y[bad_y]))
    na_target_rows <- which(!is.finite(target_vec))

    comb_na_rows <- union(na_mat_rows, na_target_rows)
    rows_to_keep <- setdiff(seq_len(nrow), comb_na_rows)

    mat <- mat[rows_to_keep, , drop = FALSE]
    target_vec <- target_vec[rows_to_keep]
  }

  list(design_matrix = mat, y_trunc = target_vec)
}
calculate_predictive_kge <- function(bvec, x, y,
                                     first_day_index,
                                     m, max_lag, ars = 0) {
  ## build design matrix and y_trunc
  fit_mats <- get_design_matrix(x = x, y = y,
                                first_day_index = first_day_index,
                                m = m, max_lag = max_lag, ars = ars)

  pred_y <- as.numeric(fit_mats$design_matrix %*% bvec)
  y_trunc <- fit_mats$y_trunc

  KGE(sim = as.numeric(pred_y), obs = as.numeric(y_trunc))
}

beta_dir <- "./results"
scores <- readRDS("./clustering/precomputed_score_df.rds")
get_gc <- function(fn) str_remove_all(fn, "[^\\d]")
all_beta_files <- list.files(beta_dir)

scores_ext <- data.frame(
  gridcode = to_check,
  true_pred_kge = rep(NA, n),
  check_pred_kge = rep(NA, n),
  check_inf_kge = rep(NA, n)
)

p <- txtProgressBar(min = 0, max = n, style = 3)
for (gc_i in to_check) {
  ## get predictive B(s, t) estimate
  load(paste0(beta_dir, "pred_model_", dszn_key, "_grid", gc_i, ".rda"))
  pred_beta <- model$pred_beta
  rm(model)

  ## get inferential B(s, t) estimate
  load(paste0(beta_dir, "inf_model_", dszn_key, "_grid", gc_i, ".rda"))
  inf_beta <- model$pred_beta
  rm(model)

  ## get data vectors
  x <- as.numeric(t(readRDS(paste0("./data/processed/grid", gc_i, "_rain.rds"))))
  y <- as.numeric(t(readRDS(paste0("./data/processed/grid", gc_i, "_streamflow.rds"))))

  ## code directly from predictive run
  m <- 364
  max_lag <- 150
  trainfrac <- 0.8
  n_train <- floor(trainfrac * length(x) / (m + 1)) * (m + 1)
  x_train <- x[1 : n_train]
  y_train <- y[1 : n_train]
  x_test <- x[(n_train + 1) : length(x)]
  y_test <- y[(n_train + 1) : length(y)]
  first_day_test <- (n_train %% (m + 1)) + 1

  pred_kge_known <- known_scores$pred_kge[known_scores$model_identifier == paste0(gc_i, "_", dszn_key)]
  pred_kge_check <- calculate_predictive_kge(
    bvec = pred_beta,
    x = x_test, y = y_test,
    first_day_index = first_day_test,
    m = 364, max_lag = 150, ars = 0
  )

  inf_kge_check <- calculate_predictive_kge(
    bvec = inf_beta,
    x = x_test, y = y_test,
    first_day_index = first_day_test,
    m = 364, max_lag = 150, ars = 0
  )

  #cat("GRID:", gc_i, "KNOWN PRED KGE:", pred_kge_known, "CHECK PRED KGE:", pred_kge_check, "INF KGE:", inf_kge_check, "\n")
  check_results[check_results$grid == gc_i, ] <- list(
    grid = gc_i,
    true_pred_kge = pred_kge_known,
    check_pred_kge = pred_kge_check,
    check_inf_kge = inf_kge_check
  )
  setTxtProgressBar(p, which(to_check == gc_i))
}

check_results$inf_advantage <- check_results$check_inf_kge - check_results$check_pred_kge
median(check_results$inf_advantage, na.rm = TRUE)
sum(check_results$inf_advantage > 0, na.rm = TRUE) / sum(!is.na(check_results$inf_advantage))
