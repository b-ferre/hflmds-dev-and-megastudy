library(Matrix)
library(GPfit)
library(hydroGOF)

## Fork-based parallelism is available on Fir (Linux). Keep one core available
## when n_cores is omitted, and fall back to serial evaluation off-cluster.
parallel_lapply <- function(x, fun, n_cores = 1L) {
  n_cores <- as.integer(n_cores)
  if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 1L) n_cores <- 1L
  n_cores <- max(1L, min(n_cores, length(x)))
  if (n_cores == 1L || .Platform$OS.type == "windows") {
    lapply(x, fun)
  } else {
    parallel::mclapply(x, fun, mc.cores = n_cores, mc.preschedule = TRUE)
  }
}

## helper for checking if given set of arcfs are stationary
## TODO: document properly
is_stationary <- function(phi, tolerance = 1e-6) {
  if (length(phi) == 0L) {
    return(TRUE)
  }

  roots <- polyroot(c(1, -phi))
  all(is.finite(Mod(roots))) &&
    all(Mod(roots) > 1 + tolerance)
}

### get design matrix #############
##
## build a design matrix X from x (and y, if autocorrelated error
## whitening is being used) which is masked by dvec and constructed
## such that Xb ~ y
##
###################################
#'
#' @param x is the flattened time series of covariate values
#' @param y is the flattened time series of target values
#' @param first_day_index is the day of the year index (in 1:365) of the first
#' entry in the provided x and y time series, used to order the columns of the
#' design matrix accordingly
#'
#' @param m is the length of the cycle which t is indexed to minus 1
#' @param max_lag is the number of effective lags s to consider (starting at 0)
#' @param ars is the degree of the ar process used to model the residuals;
#' ars = 0 indicates that no autocorrelated residual correction is being used
#'
#' @param remove_na_rows boolean which dictates whether the returned design
#' matrix should be trimmed to disinclude any rows where either the target
#' or covariates are NA
#'
#' @param dvec if provided, corresponds to a 365-length vector that lists the
#' number of effective lags s to consider (starting at 0) on day i of the cycle
#'
#' @return a list containing design_matrix, a matrix built from lagged values of
#' x (and y if ars !=0) and masked according to dvec, and a truncated y vector
#' such that design matrix %*% bvec ~ y
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

### get regularization matrices ###
##
## bvec :
##   - [(M + 1) * (max_lags + ars)] irf parameters indexed primarily by
##   increasing s and secondarily by increasing t
##   - ars ar coefficients tacked on to the end
##
## D_h : applies first difference penalty on the t axis of the irf parameters
## and does so cyclically as the final day of the year should be similar to
## the first day of the year
##
## D_v : applies second derivative Tikhonov-Phillips regularization along the
## s axis of the irf-related parameters
##
###################################
#'
#' @param m is the length of the cycle which t is indexed to minus 1
#' @param max_lag is the number of effective lags s to consider (starting at 0)
#' @param ars is the degree of the ar process used to model the residuals;
#' ars = 0 indicates that no autocorrelated residual correction is being used
#' @param d_v_diff is a number (either 1 or 2) representing which
#' differences should be penalized along the lag axis. 1, 2 correspond to the
#' first and second differences, respectively
#'
#' @returns a list of d_h and d_v - the day-wise and lag-wise regularization matrices     ## nolint
get_regularization_mats <- function(m, max_lag, ars = 0, d_v_diff = 1) {

  ## build d_h
  nrow_h <- (m + 1) * (max_lag + ars)
  ncol_h <- ((m + 1) * (max_lag + ars)) + ars
  i_h <- as.vector(sapply(seq_len(nrow_h), function(q) rep(q, 2)))
  j_h <- as.vector(sapply(seq_len(nrow_h), function(q) {
    c(
      as.integer(q),
      as.integer((max_lag + ars + q - 1) %% (ncol_h - ars)) + 1
    )
  }))
  x_h <- as.vector(sapply(seq_len(nrow_h), function(q) c(-1, 1)))
  d_h <- sparseMatrix(i = i_h, j = j_h, x = x_h, dims = c(nrow_h, ncol_h))

  ## build d_v (2nd diff)
  if (d_v_diff == 2) {
    num_interior_lags <- (max_lag + ars - 2)
    nrow_v <- num_interior_lags * (m + 1)
    ncol_v <- ((m + 1) * (max_lag + ars)) + ars ## 0-pad for AR terms
    i_v <- as.vector(sapply(seq_len(nrow_v), function(q) rep(q, 3)))
    j_v <- as.vector(sapply(seq_len(nrow_v), function(q) {
      qstar <- ((q - 1) %% num_interior_lags) + 1
      block_idx <- as.integer(ceiling(q / num_interior_lags))
      start_idx <- ((block_idx - 1) * (max_lag + ars)) + qstar
      c(start_idx,
        start_idx + 1,
        start_idx + 2
      )
    }))
    x_v <- as.vector(sapply(seq_len(nrow_v), function(q) c(1, -2, 1)))
    d_v <- sparseMatrix(i = i_v, j = j_v, x = x_v, dims = c(nrow_v, ncol_v))
  } else {
    stopifnot(d_v_diff == 1)
    num_interior_lags <- (max_lag + ars - 1)
    nrow_v <- (num_interior_lags + 1) * (m + 1)
    ncol_v <- ((m + 1) * (max_lag + ars)) + ars ## 0-pad for AR terms
    i_v <- as.vector(unlist(sapply(seq_len(nrow_v), function(q) {
      if ((q %% (num_interior_lags + 1)) == 0) rep(q, 1) else rep(q, 2)
    })))
    j_v <- as.vector(unlist(sapply(seq_len(nrow_v), function(q) {
      if ((q %% (num_interior_lags + 1)) == 0) c(q) else c(q, q + 1)
    })))
    x_v <- as.vector(unlist(sapply(seq_len(nrow_v), function(q) {
      if ((q %% (num_interior_lags + 1)) == 0) c(1) else c(-1, 1)
    })))
    d_v <- sparseMatrix(i = i_v, j = j_v, x = x_v, dims = c(nrow_v, ncol_v))
  }

  ## return list of regularization matrices
  list(d_h = d_h, d_v = d_v)
}

### perform a model fit ###########
##
## find the optimal set of b coefficients that satisfy:
## min_b (||Xb - y|| + w_h||D_h*b|| + w_v||D_v*b||)
## where X is a design matrix built from x and dvec (and y, if
## autocorrelated error whitening is being used)
##
###################################
#'
#' @param m is the length of the cycle which t is indexed to minus 1
#' @param max_lag is the number of effective lags s to consider (starting at 0)
#' @param ars is the degree of the ar process used to model the residuals;
#' ars = 0 indicates that no autocorrelated residual correction is being used
#'
#' @param w_h is the logarithm of the hyperparameter that tunes day-wise
#' smoothness regularization
#' @param w_v is the logarithm of the hyperparameter that tunes lag-wise
#' smoothness regularization
#'
#' @param x is the flattened time series of covariate values to train on
#' @param y is the flattened time series of target values to train on
#' @param first_day_index is the day of the year index (in 1:365) of the first
#' entry in the provided x and y time series, used to order the columns of the
#' design matrix accordingly
#'
#' @param dvec if provided, corresponds to the sequence of lags after which the
#' covariate is assumed to have no effect on the target variable and zeroed out
perform_fit <- function(m, max_lag, ars = 0,
                        w_h, w_v,
                        x, y, first_day_index, dvec = NULL) {

  ## get regularization matrices
  regmats <- get_regularization_mats(m, max_lag, ars)
  d_h <- regmats$d_h
  d_v <- regmats$d_v

  fit_mats <- get_design_matrix(x = x, y = y,
                                first_day_index = first_day_index,
                                m = m, max_lag = max_lag,
                                ars = ars)
  design_mat <- fit_mats$design_matrix
  y_trunc <- fit_mats$y_trunc

  # [Z.T @ Z + exp(w_h) * D_h.T @ D_h + exp(w_v) * D_v.T @ D_v] @ b = Z.T @ y
  lhs <- crossprod(design_mat) +
    exp(w_h) * crossprod(d_h) +
    exp(w_v) * crossprod(d_v)
  rhs <- crossprod(design_mat, y_trunc)

  ## fit all non-zeroed B(s, t) indices
  keep <- seq_len(ncol(design_mat))
  if (!is.null(dvec)) {
    keep <- c(unlist(Map(function(i, d) (i - 1) * (max_lag + ars) + seq_len(d + ars), ## nolint
                         seq_along(dvec), dvec)),
              (m + 1) * (max_lag + ars) + seq_len(ars))
  }
  fit <- numeric(ncol(design_mat))
  fit[keep] <- Matrix::solve(lhs[keep, keep, drop = FALSE], rhs[keep])
  fit
}

### predict y_hat given x and b/B coeffs ##
##
## generate a prediction of y_hat given x
## DO NOT prune out NA predictions as this
## destroys the continuity of yhat which
## breaks residual calculation
##
###################################
#'
#' @param bvec is either a flat vec of Beta(s, t) parameters indexed primarily
#' by increasing lag and secondarily by increasing day OR a sequence of b
#' parameters that consist of (max_lag + ars) beta-esque parameters for each
#' day of the cycle and ars autoregressive coefficients at the end
#'
#' @param m is the length of the cycle which t is indexed to minus 1
#' @param max_lag is the number of effective lags s to consider (starting at 0)
#' @param ars is the degree of the ar process used to model the residuals;
#' ars = 0 indicates that no autocorrelated residual correction is being used
#'
#' @param x is the flattened time series of covariate values to test on
#' @param y is the flattened time series of target values to test on
#' @param first_day_index is the day of the year index (in 1:365) of the first
#' entry in the provided x and y time series, used to order the columns of the
#' design matrix accordingly
#'
#' @param dvec if provided, corresponds to the sequence of lags after which the
#' covariate is assumed to have no effect on the target variable and zeroed out
#'
#' @return a list of y_hat and y_true
calculate_yhat <- function(bvec, x, y,
                           first_day_index,
                           m, max_lag, ars = 0) {
  design <- get_design_matrix(x = x, y = y,
                              first_day_index = first_day_index,
                              m = m, max_lag = max_lag,
                              ars = ars,
                              remove_na_rows = FALSE)
  design_mat <- design$design_matrix
  true_y <- design$y_trunc

  pred_y <- as.numeric(design_mat %*% bvec)

  list(y_hat = pred_y, y_true = true_y)
}

### calculate KGE GOF of given model ##
##
## calculate the predictive r2 of a set of b/beta coefficients
## when using x (and maybe y) to predict y
##
###################################
#'
#' @param bvec is either a flat vec of Beta(s, t) parameters indexed primarily
#' by increasing lag and secondarily by increasing day OR a sequence of b
#' parameters that consist of (max_lag + ars) beta-esque parameters for each
#' day of the cycle and ars autoregressive coefficients at the end
#'
#' @param m is the length of the cycle which t is indexed to minus 1
#' @param max_lag is the number of effective lags s to consider (starting at 0)
#' @param ars is the degree of the ar process used to model the residuals;
#' ars = 0 indicates that no autocorrelated residual correction is being used
#'
#' @param x is the flattened time series of covariate values to test on
#' @param y is the flattened time series of target values to test on of the same
#' length as the eventual yhat estimate
#' @param first_day_index is the day of the year index (in 1:365) of the first
#' entry in the provided x and y time series, used to order the columns of the
#' design matrix accordingly
#'
#' @param dvec if provided, corresponds to the sequence of lags after which the
#' covariate is assumed to have no effect on the target variable and zeroed out
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

### calculate r2 GOF of given model ##
##
## calculate the predictive r2 of a set of b/beta coefficients
## when using x (and maybe y) to predict y
##
###################################
#'
#' @param bvec is either a flat vec of Beta(s, t) parameters indexed primarily
#' by increasing lag and secondarily by increasing day OR a sequence of b
#' parameters that consist of (max_lag + ars) beta-esque parameters for each
#' day of the cycle and ars autoregressive coefficients at the end
#'
#' @param m is the length of the cycle which t is indexed to minus 1
#' @param max_lag is the number of effective lags s to consider (starting at 0)
#' @param ars is the degree of the ar process used to model the residuals;
#' ars = 0 indicates that no autocorrelated residual correction is being used
#'
#' @param x is the flattened time series of covariate values to test on
#' @param y is the flattened time series of target values to test on of the same
#' length as the eventual yhat estimate
#' @param first_day_index is the day of the year index (in 1:365) of the first
#' entry in the provided x and y time series, used to order the columns of the
#' design matrix accordingly
#'
#' @param dvec if provided, corresponds to the sequence of lags after which the
#' covariate is assumed to have no effect on the target variable and zeroed out
calculate_predictive_r2 <- function(bvec, x, y,
                                    first_day_index,
                                    m, max_lag, ars = 0) {
  ## build design matrix and y_trunc
  fit_mats <- get_design_matrix(x = x, y = y,
                                first_day_index = first_day_index,
                                m = m, max_lag = max_lag, ars = ars)

  pred_y <- as.numeric(fit_mats$design_matrix %*% bvec)
  y_trunc <- fit_mats$y_trunc

  1 - (sum((pred_y - y_trunc) ** 2) / sum((y_trunc - mean(y_trunc)) ** 2))
}

### get beta convolution matrix #
##
## build a matrix M such that Mbeta = b
##
###################################
#'
#' @param m is the length of the cycle which t is indexed to minus 1
#' @param max_lag is the number of effective lags s to consider (starting at 0)
#' @param arcfs are the actual estimated ar coefficients used for deconvolution
get_conv_mtx <- function(arcfs, m, max_lag) {

  ## FIXME: actually check this code instead of blindly trusting chatslop

  ars <- length(arcfs)

  max_elements <- (max_lag * (m + 1)) + (max_lag * (m + 1) * ars)
  row_idx <- integer(max_elements)
  col_idx <- integer(max_elements)
  vals <- numeric(max_elements)
  counter <- 1

  for (k in 1:max_lag) {
    for (d in 1:(m + 1)) {
      row_idx[counter] <- (k - 1) * (m + 1) + d
      col_idx[counter] <- (k - 1) * (m + 1) + d
      vals[counter] <- 1
      counter <- counter + 1
    }
  }

  if (max_lag > 1) {
    for (k in 2:max_lag) {
      for (i in 1:min(k - 1, ars)) {
        row_base <- (k - 1) * (m + 1)
        col_base <- (k - 1 - i) * (m + 1)
        for (d in 1:(m + 1)) {
          row_idx[counter] <- row_base + d
          c_day <- (d - 1 - i) %% (m + 1) + 1
          col_idx[counter] <- col_base + c_day
          vals[counter] <- -arcfs[i]
          counter <- counter + 1
        }
      }
    }
  }

  actual_len <- counter - 1
  m_matrix <- sparseMatrix(i = row_idx[1:actual_len],
                           j = col_idx[1:actual_len],
                           x = vals[1:actual_len],
                           dims = c(max_lag * (m + 1), max_lag * (m + 1)))
  m_matrix
}

### deconvolve b to beta ##########
##
## deconvolve the provided vector of estimated b coefficients into the
## implicit B(s, t) coefficients they represent by solving M %*% beta = b
##
###################################
#'
#' @param bvec is a flattened vector of b coefficients representing
#' convolutions of underlying B(s, t) coefficients with AR coeffs
#' @param m is the length of the cycle which t is indexed to minus 1
#' @param max_lag is the number of effective lags s to consider (starting at 0)
#' @param ars is a positive integer representing the degree of the ar
#' process used to model the residuals
#' @return a list of pred_beta (the IRF coeffs) and pred_phi (the AR coeffs)
deconvolve_b_to_beta <- function(bvec, m, max_lag, ars) {

  a <- tail(bvec, ars)
  b_trunc <- head(bvec, -ars)
  b_blocks <- matrix(b_trunc, nrow = max_lag + ars, ncol = m + 1)
  b_core <- b_blocks[1:max_lag, , drop = FALSE]
  b_core_vec <- as.vector(t(b_core))
  m_matrix <- get_conv_mtx(arcfs = a, m = m, max_lag = max_lag)
  beta_vec <- as.vector(Matrix::solve(m_matrix, b_core_vec))
  beta_matrix <- matrix(beta_vec, nrow = max_lag, ncol = m + 1, byrow = TRUE)
  beta_vec <- as.vector(beta_matrix)

  list(pred_beta = beta_vec, pred_phi = a)
}

### SAFELY decode beta from b #####
##
## deconvolve the provided vector of estimated b coefficients into the
## implicit B(s, t) coefficients they represent by solving M %*% beta = b
## after checking stationarity of estimated phi coefficients. Then, if
## provided dvec is not null, re-enforce mask (for each day t, set all
## coefficients B(s, t) such that s > dvec[t] to 0)
## Throws an error if phi is nonstationary.
##
###################################
#'
#' @param bvec is a flattened vector of b coefficients representing
#' convolutions of underlying B(s, t) coefficients with AR coeffs
#' @param m is the length of the cycle which t is indexed to minus 1
#' @param max_lag is the number of effective lags s to consider (starting at 0)
#' @param ars is a positive integer representing the degree of the ar
#' process used to model the residuals
#' @param dvec is the sequence of maximum effective lags for each day t
#'
#' @return a lag-minor indexed vector of B(s, t) coefficients
get_beta_from_b_safely <- function(bvec, m, max_lag, ars, dvec = NULL) {
  if (ars > 0) {
    if (is.null(dvec)) {
      dvec <- rep(max_lag, (m + 1))
    }
    dec <- deconvolve_b_to_beta(bvec = bvec,
                                m = m, max_lag = max_lag,
                                ars = ars)
    if (!is_stationary(dec$pred_phi)) stop("Nonstationary phi found by get_beta_from_b_safely().")                                    ## nolint
    betavec <- dec$pred_beta
    betavec[rep(seq_len(max_lag), m + 1) >
              rep(dvec, each = max_lag)] <- 0
    betavec
  }
}

### hyperparameter optimization ###
##
## tests different values of w_h and w_v and optimize for performance
## on a validation set. begins by performing a grid search over
## w_h_range x w_v_range, and then uses a gaussian process to refine
## the estimates of the optimal hyperparameter values
##
###################################
#'
#' @param m is the length of the cycle which t is indexed to minus 1
#' @param max_lag is the number of effective lags s to consider (starting at 0)
#' @param ars is the degree of the ar process used to model the residuals
#'
#' @param n_grid the number of points to use in the grid search stage IN
#' EACH DIMENSION (total points checked = n_grid ** 2)
#' @param n_gp the number of points to use in the gaussian process stage
#'
#' @param w_h_range the logarithm of the range of values to consider for w_h
#' i.e. a parameter value of c(-1, 3) checks w_h values in exp(-1) to exp(3)
#' @param w_v_range the logarithm of the range of values to consider for w_v
#' i.e. a parameter value of c(-1, 3) checks w_v values in exp(-1) to exp(3)
#'
#' @param x is the full useable time series for the covariate (flattened) and
#' the first entry in x MUST align with the first day of the cycle in t
#' @param y is the full useable series for the target variable (flattened) and
#' the first entry in y MUST align with the first day of the cycle in t
#' @param dvec if provided, corresponds to the sequence of lags after which the
#' covariate is assumed to have no effect on the target variable and zeroed out
#'
#' @param datasplit is a vector of the form (training frac, validation frac)
#' which dictates how x and y should be partitioned
#'
#' @param verbose is a boolean dictating whether a progress bar for the
#' two rounds of optimization should be displayed
#'
#' TODO: finish documenting
#' TODO: finish making completely observable post-mortem
#'
#' @returns the optimal values of log(w_h) and log(w_v)
optimize_hyperparameters <- function(m, max_lag, ars,
                                     n_grid, n_gp, w_h_range, w_v_range,
                                     x, y, dvec = NULL,
                                     datasplit = c(0.75, 0.25),
                                     verbose = TRUE, n_cores = 1L,
                                     k_folds_for_cv = 1,
                                     d_v_diff = 1,
                                     large_hp_pseudoprior = FALSE) {

  ## assign a default dvec if none was provided
  if (is.null(dvec)) {
    dvec <- rep(max_lag, m + 1)
  }

  print("running optimize_hyperparameters...")

  ## build entire design matrix using all data which will later
  ## be split into train/val splits
  design_mat_full <- get_design_matrix(
    x, y, first_day_index = 1, m = m,
    max_lag = max_lag, ars = ars, remove_na_rows = TRUE
  )

  ## precompute all training/validation data pairs
  default_val_starts <- floor(seq(from = 1 - datasplit[2], to = 0, by = -datasplit[2]) * nrow(design_mat_full)) ## nolint
  val_length_d <- floor(datasplit[2] * nrow(design_mat_full))
  train_length_d <- nrow(design_mat_full) - val_length_d
  train_val_splits <- lapply(
    seq_len(k_folds_for_cv),
    function(fold_index) {
      if (fold_index <= length(default_val_starts)) {
        val_start <- default_val_starts[fold_index]
      } else {
        val_start <- floor(runif(1, min = 0, max = 1) * nrow(design_mat_full))
      }
      val_end <- val_start + val_length_d - 1
      val_indexes <- ((seq(from = val_start, to = val_end, by = 1) - 1) %% nrow(design_mat_full)) + 1         ## nolint

      train_start <- val_end + 1
      train_end <- train_start + train_length_d - 1
      train_indexes <- ((seq(from = train_start, to = train_end, by = 1) - 1) %% nrow(design_mat_full)) + 1   ## nolint

      list(
        train_mat = design_mat_full$design_matrix[train_indexes, , drop = FALSE],                             ## nolint
        train_y = design_mat_full$y_trunc[train_indexes],
        val_mat = design_mat_full$design_matrix[val_indexes, , drop = FALSE],
        val_y_trunc = design_mat_full$y_trunc[val_indexes]
      )
    }
  )

  ## precalculate our dvec column mask
  keep_coef <- seq_len(ncol(design_mat_full$design_matrix))
  if (!is.null(dvec)) {
    keep_coef <- c(
      unlist(Map(
        function(i, d) (i - 1L) * (max_lag + ars) + seq_len(d + ars),
        seq_along(dvec), dvec
      )),
      (m + 1L) * (max_lag + ars) + seq_len(ars)
    )
  }

  ## build our regularization matrices & precompute penalty terms
  regmats <- get_regularization_mats(m, max_lag, ars, d_v_diff = d_v_diff)
  h_penalty <- crossprod(regmats$d_h[, keep_coef, drop = FALSE])
  v_penalty <- crossprod(regmats$d_v[, keep_coef, drop = FALSE])

  ## for each fold, mask training matrix and precompute fold-specific
  ## terms in each linear regression to reduce solving overhead
  train_val_splits <- lapply(train_val_splits, function(fold) {
    train_mat_kept <- fold$train_mat[, keep_coef, drop = FALSE]

    fold$xtx <- crossprod(train_mat_kept)
    fold$xty <- crossprod(train_mat_kept, fold$train_y)
    fold$val_tss <- sum((fold$val_y_trunc - mean(fold$val_y_trunc))^2)

    fold
  })

  ## helper function for testing a candidate pair of w_h/w_v values across
  ## all available train/val splits. report mean validation error and
  ## standard error
  score_hp <- function(w_h, w_v) {
    tryCatch({
      ## evaluate this hyperparameter pair across all folds
      fold_r2s <- vapply(train_val_splits, function(fold) {

        ## Combine pre-computed matrices for this specific fold
        lhs <- fold$xtx + exp(w_h) * h_penalty + exp(w_v) * v_penalty

        fit <- numeric(ncol(design_mat_full$design_matrix))
        fit[keep_coef] <- Matrix::solve(lhs, fold$xty)

        if (ars > 0) {
          fit <- get_beta_from_b_safely(fit,
                                        m = m,
                                        max_lag = max_lag,
                                        ars = ars,
                                        dvec = dvec)
        }

        pred_y <- as.numeric(fold$val_mat %*% fit)
        r2 <- 1 - sum((pred_y - fold$val_y_trunc)^2) / fold$val_tss

        r2
      }, FUN.VALUE = numeric(1))

      ## filter all fits to only include the ones that didn't fail
      valid_r2s <- fold_r2s[is.finite(fold_r2s)]
      k_valid <- length(valid_r2s)

      ## calculate the mean validation r2 across all successful fits
      mean_r2 <- if (k_valid > 0) mean(valid_r2s) else NA_real_

      ## calculate the standard error of the r2 across all successful fits
      se_r2 <- if (k_valid > 1) sd(valid_r2s) / sqrt(k_valid) else NA_real_

      ## return values
      list(
        mean_val_r2 = mean_r2,
        se_val_r2 = se_r2
      )
    }, error = function(e) {
      warning(paste0("score_hp failed with w_h/w_v = ", w_h, "/", w_v, "\n"))
      NA_real_
    })
  }

  ############################## GRID SEARCH ##################################

  ## define a grid to search over
  grid_w_h <- seq(from = w_h_range[1], to = w_h_range[2], length.out = n_grid)              ## nolint
  grid_w_v <- seq(from = w_v_range[1], to = w_v_range[2], length.out = n_grid)              ## nolint
  grid_cells <- expand.grid(grid_w_h, grid_w_v)
  names(grid_cells) <- c("w_h", "w_v")

  ## run score_hp over the grid
  grid_scores_list <- parallel_lapply(seq_len(nrow(grid_cells)), function(i) {
    res <- score_hp(grid_cells$w_h[i], grid_cells$w_v[i])
    if (verbose) {
      cat("grid search iteration", i, "out of", nrow(grid_cells),
          "complete; R2 is", res["r2"], "(SE:", res["se"], ")\n",
          "w_h is ", grid_cells$w_h[i], ", and w_v is ", grid_cells$w_v[i],
          "\n")
    }
    res
  }, n_cores)

  ## extract means and SEs from grid fits
  grid_r2s <- vapply(grid_scores_list, function(x) x["r2"], numeric(1))
  grid_ses <- vapply(grid_scores_list, function(x) x["se"], numeric(1))

  if (!any(is.finite(grid_r2s))) {
    stop("All hyperparameter evaluations failed.")
  }

  ########################### GAUSSIAN PROCESS ################################

  ## GP helper for scaling a data.frame of w_h and w_v values into [0, 1]
  scale_hp <- function(z) {
    cbind(
      (z$w_h - w_h_range[1]) / diff(w_h_range),
      (z$w_v - w_v_range[1]) / diff(w_v_range)
    )
  }

  ## GP helper for converting w_h/w_v pairs into string form to find dupes
  hp_key <- function(z) {
    paste(round(z$w_h, 12), round(z$w_v, 12))
  }

  ## define dense w_h x w_v surface to use for GP estimates
  dense_grid <- expand.grid(
    w_h = seq(w_h_range[1], w_h_range[2], length.out = 100),
    w_v = seq(w_v_range[1], w_v_range[2], length.out = 100)
  )

  for (gp_i in seq_len(n_gp)) {
    tryCatch(expr = {
      ## don't include w_h/w_v pairs that
      ## a) yielded problematic r2 values as this will break GP deviance calc
      ## b) are duplicates of previous rows as this will break cholesky calc
      keep <- is.finite(grid_r2s) &
        !duplicated(hp_key(grid_cells))

      if (sum(keep) < 3) {
        stop("Too few successful hyperparameter evaluations to fit the GP.")
      }

      gp <- GP_fit(
        scale_hp(grid_cells[keep, , drop = FALSE]),
        grid_r2s[keep]
      )

      pred <- predict.GP(gp, scale_hp(dense_grid))
      mu <- as.numeric(pred$Y_hat)
      sd <- sqrt(pmax(as.numeric(pred$MSE), 0))
      imp <- mu - max(grid_r2s[keep])
      z <- imp / pmax(sd, .Machine$double.eps)

      ei <- imp * pnorm(z) + sd * dnorm(z)
      ei[sd == 0] <- pmax(imp[sd == 0], 0)
      ei[hp_key(dense_grid) %in% hp_key(grid_cells)] <- -Inf
      ei[!is.finite(ei)] <- -Inf

      if (!any(is.finite(ei))) break

      candidate <- dense_grid[which.max(ei), , drop = FALSE]
      candidate_res <- score_hp(candidate$w_h, candidate$w_v)

      r2 <- candidate_res$r2
      se <- candidate_res$se

      grid_cells <- rbind(grid_cells, candidate)
      grid_r2s <- c(grid_r2s, r2)
      grid_ses <- c(grid_ses, se)

      if (verbose) {
        cat(
          "grid search iteration", gp_i, "out of", n_gp,
          "complete; R2 is", r2, "(SE:", se, ")\n",
          "w_h is ", candidate$w_h, ", and w_v is ", candidate$w_v,
          "\n"
        )
      }
    }, error = function(e) {
      cat("Encountered an error (below) during GP iteration",
          gp_i, "of", n_gp, "\n")
      print(e)
    })
  }

  ########################### OPTIMUM DETECTION ###############################

  ## of all combinations tested, find the best pair and return
  if (large_hp_pseudoprior) {
    ## find the absolute best validation accuracy across all tested w_h/w_v
    ## pairs, as well as its standard error, if available
    max_idx <- which.max(grid_r2s)
    max_r2 <- grid_r2s[max_idx]
    max_se <- grid_ses[max_idx]

    ## define our error tolerance. preferably use standard error, but if
    ## k-folds = 1 or something went awry during the fits for the pair with the
    ## best validation accuracy, default to a tolerance of up to a 0.01 drop in
    ## validation R2
    tolerance <- if (is.na(max_se)) 0.01 else max_se

    ## identify all candidate w_h/w_v pairs that scored within tolerance of the
    ## best accuracy
    striking_distance_idx <- which(grid_r2s >= (max_r2 - tolerance))

    w_h_pool <- grid_cells$w_h[striking_distance_idx]
    w_v_pool <- grid_cells$w_v[striking_distance_idx]

    ## compute relative regularization weights (w.r.t. the range of values of
    ## w_h and w_v within the in-tolerance candidate pairs)
    range_h <- diff(range(w_h_pool))
    range_v <- diff(range(w_v_pool))
    scaled_w_h <- if (range_h > 0) (w_h_pool - min(w_h_pool)) / range_h else rep(0, length(w_h_pool)) ## nolint
    scaled_w_v <- if (range_v > 0) (w_v_pool - min(w_v_pool)) / range_v else rep(0, length(w_v_pool)) ## nolint

    ## calculate a weighted euclidean distance from the origin of the
    ## candidate pool to serve as a measure of combined hyperparameter
    ## strength. we prioritize w_h because it is more in line with the
    ## physics of catchment hydrology to expect a large degree of day-
    ## to-day smoothing than lag-to-lag smoothing
    weight_h <- 0.75
    weight_v <- 0.25
    hp_mag <- sqrt((weight_h * scaled_w_h)^2 + (weight_v * scaled_w_v)^2)

    ## select the hyperparameter pair with the largest combined hyperparameter
    ## strength that is within "striking distance" of the best performance
    best_relative_idx <- which.max(hp_mag)
    best_idx <- striking_distance_idx[best_relative_idx]

    best_w_h <- grid_cells$w_h[best_idx]
    best_w_v <- grid_cells$w_v[best_idx]

    if (verbose) {
      cat("1-SE Rule applied: Max R2 was", round(max_r2, 4),
          "| Selected R2 is", round(grid_r2s[best_idx], 4),
          "with penalties (w_h:", round(best_w_h, 2), ",
          w_v:", round(best_w_v, 2), ")\n")
    }

  } else {
    best_idx <- which.max(grid_r2s)
    best_w_h <- grid_cells$w_h[best_idx]
    best_w_v <- grid_cells$w_v[best_idx]
  }

  ################# RETURN ALL VALUES (FOR FULL VISIBILITY) ###################

  list(
    hyperparameters = list(w_h = best_w_h, w_v = best_w_v),
    search_dataframe = rbind(grid_cells, grid_r2s, grid_ses)
  )
}

### perform quantile thresholding #
##
## TODO: write a description
##
###################################
#'
#' TODO: document fully
#' TODO: finish making entirely observable post-mortem
#'
optimize_dvec_elbow <- function(beta_est,
                                m, max_lag, ars = 0,
                                x, y,
                                nq1 = 100, nq2 = 100,
                                max_q = 0.5,
                                datasplit = c(0.75, 0.25),
                                w_h, w_v,
                                grp_norm_calc_method = "new",
                                verbose = TRUE, n_cores = 1L) {

  n_cores <- max(1L, as.integer(n_cores))

  ## validate beta_est
  if (length(beta_est) != (m + 1L) * max_lag ||
        anyNA(beta_est) ||
        any(!is.finite(beta_est))) {
    stop("beta_est provided to optimize_dvec is invalid.")
  }

  ## calculate day-relative group values of the provided beta_est
  if (grp_norm_calc_method == "new") {
    stop("NEW GRP NORM CALCULATION METHOD NOT YET FINALIZED")
    ##FIXME: update this to work the way I want
    # group_values <- unlist(lapply(seq_len(m + 1), function(day_index) {
    #   beta_t <- beta_est[(((day_index - 1) * max_lag) + 1):(day_index * max_lag)]
    #   daytotal <- sum(beta_t)
    #   cumsums <- rev(cumsum(rev(beta_t)))
    #   cumsums / daytotal
    # }))
  } else {
    if (grp_norm_calc_method != "old") warning("Invalid grp_norm_calc_method. Defaulting to old.")  ## nolint
    group_values <- unlist(lapply(seq_len(m + 1), function(day_index) {
      beta_t <- beta_est[(((day_index - 1) * max_lag) + 1):(day_index * max_lag)]                   ## nolint
      rev(cumsum(rev(beta_t ** 2)))
    }))
  }
  ## helper to calculate the dvec associated with zeroing out lags where
  ## the tail makes up a fraction of less than q*100 % of total in-day
  ## coefficient value
  calculate_dvec_from_q <- function(q) {
    ## form an (m + 1) x max_lag matrix with entry=1 iff corresponding B(s, t)
    ## should be forced to zero
    mask_mtx <- matrix(group_values < q, ncol = max_lag, byrow = TRUE)
    as.numeric(apply(mask_mtx, MARGIN = 1, function(mask_row) {
      ## zero out all lags beyond the point where the tail fraction falls below
      ## q even if it later rises back above q due to non-monotonicity of group
      ## values induced by negative B(s, t) parameter estimates
      first_ineff_lag <- which(mask_row == 1)[1]
      if (is.na(first_ineff_lag)) {
        max_lag
      } else {
        first_ineff_lag - 1
      }
    }))
  }

  ## helper to find the elbow point in a (generally) increasing
  ## plot of validation performance as function of # of params
  ## if loess_span is set to NA, then no smoothing is used
  get_elbow_index <- function(x, y, min_segment_length = 3, loess_span = 0.4) {
    ## 1. Bind to original indices to prevent mapping errors after sorting
    df <- data.frame(orig_idx = seq_along(x), x = x, y = y)

    ## 2. Filter non-finite values
    df <- df[is.finite(df$x) & is.finite(df$y), ]

    ## 2.5 Remove duplicated x values to prevent LOESS fit collapse
    df <- df[!duplicated(df$x), ]

    ## Ensure enough data remains to fit LOESS model
    if (nrow(df) < (min_segment_length * 2)) {
      stop("get_elbow_index() failed: Insufficient finite pairs for piecewise regression.") ## nolint
    }

    ## 3. Sort strictly by x to correctly sequence the curve
    df <- df[order(df$x), ]

    ## 4. Extract the underlying signal from the validation noise
    ## A span of 0.4 provides strong protection against local spikes while
    ## preserving the macroscopic concavity of the curve
    if (is.na(loess_span)) {
      df$y_smooth <- y
    } else {
      df$y_smooth <- predict(loess(y ~ x, data = df, span = loess_span))
    }

    ## 5. Segmented Regression (Structural Break Detection)
    n <- nrow(df)
    total_rss <- numeric(n)
    total_rss[] <- Inf

    ## Iterate through all valid knot points, fitting a line to the "steep"
    ## parameter-gain phase and a line to the "flat" asymptotic phase.
    for (i in min_segment_length:(n - min_segment_length + 1)) {

      ## Fit the left segment (from parameter minimum to candidate elbow)
      left_fit <- lm(y_smooth ~ x, data = df[1:i, ])
      rss_left <- sum(residuals(left_fit)^2)

      ## Fit the right segment (from candidate elbow to parameter maximum)
      right_fit <- lm(y_smooth ~ x, data = df[i:n, ])
      rss_right <- sum(residuals(right_fit)^2)

      total_rss[i] <- rss_left + rss_right
    }

    ## 6. Find the knot that minimized the combined squared error of both lines
    best_sorted_idx <- which.min(total_rss)

    ## return all values needed to reconstruct this optimization post-mortem
    list(
      best_value = y[df$orig_idx[best_sorted_idx]],
      elbow_df = df,
      loess_span = loess_span,
      min_segment_length = min_segment_length
    )
  }

  ## get training and validation split (or just duplicate
  ## training data if using the old fitting routine)
  if (is.na(datasplit)) {
    x_train <- x
    y_train <- y
    x_val <- x
    y_val <- y
    first_day_index_train <- 1
    first_day_index_val <- 1
  } else {
    stopifnot(typeof(datasplit) == "numeric" && length(datasplit) == 2 && sum(datasplit) == 1)  ## nolint
    ## get training data
    n_train <- round(length(x) * datasplit[1])
    x_train <- x[1 : n_train]
    y_train <- y[1 : n_train]
    first_day_index_train <- 1

    ## get validation data
    x_val <- x[(n_train + 1):length(x)]
    y_val <- y[(n_train + 1):length(y)]
    first_day_index_val <- ((n_train) %% (m + 1)) + 1
  }

  ## test different values of neg_thresh (round 1)
  qs_1 <- seq(0.0, max_q, length.out = nq1)
  nonzero_coeffs_1 <- rep(NA_real_, nq1)
  r2s_1 <- rep(NA_real_, nq1)
  for (i in seq_len(nq1)) {
    q_i <- qs_1[i]
    dvec_i <- calculate_dvec_from_q(q_i)
    nonzero_coeffs_1[i] <- sum(dvec_i)

    r2s_1[i] <- tryCatch(expr = {
      bvec_i <- perform_fit(m = m, max_lag = max_lag, ars = ars,
                            w_h = w_h, w_v = w_v,
                            x = x_train, y = y_train,
                            first_day_index = first_day_index_train,
                            dvec = dvec_i)

      ## if ars > 0, deconvolve bvec into beta vector as we seek to
      ## validate using the underlying IRF coefficients, not the full
      ## ar-inclusive b vector as to encourage our IRF coefficients to
      ## capture as much of the long-tailed behaviour as possible
      if (ars > 0) {
        bvec_i <- get_beta_from_b_safely(bvec_i,
                                         m = m, max_lag = max_lag, ars = ars,
                                         dvec = dvec_i)
      }

      calculate_predictive_r2(bvec_i,
                              x = x_val, y = y_val,
                              first_day_index = first_day_index_val,
                              m = m, max_lag = max_lag, ars = 0)
    }, error = function(e) {
      print(e)
      NA_real_
    })
  }

  ## detect elbow point in graph of validation accuracy as a function
  ## of number of retained parameters
  elbow_check_1 <- get_elbow_index(x = nonzero_coeffs_1, y = r2s_1)
  best_q1 <- elbow_check_1$best_value
  if (verbose) {
    cat("first round of dvec optimization complete. q selected:",
        best_q1, "\n")
  }

  ## try to refine q2 estimate (but be careful of potentially
  ## overly tight search neighborhood)
  elbow_check_2 <- tryCatch(expr = {
    ## test different values of neg_thresh in neighborhood of initial estimate
    qs_2 <- seq(from = max(best_q1 * 0.75, 0),
                to = min(best_q1 * 1.25, max_q),
                length.out = nq2)
    nonzero_coeffs_2 <- rep(NA_real_, nq2)
    r2s_2 <- rep(NA_real_, nq2)
    for (i in seq_len(nq2)) {
      q_i <- qs_2[i]
      dvec_i <- calculate_dvec_from_q(q_i)
      nonzero_coeffs_2[i] <- sum(dvec_i)

      r2s_2[i] <- tryCatch(expr = {
        bvec_i <- perform_fit(m = m, max_lag = max_lag, ars = ars,
                              w_h = w_h, w_v = w_v,
                              x = x_train, y = y_train,
                              first_day_index = first_day_index_train,
                              dvec = dvec_i)

        ## if ars > 0, deconvolve bvec into beta vector as we seek to
        ## validate using the underlying IRF coefficients, not the full
        ## ar-inclusive b vector as to encourage our IRF coefficients to
        ## capture as much of the long-tailed behaviour as possible
        if (ars > 0) {
          bvec_i <- get_beta_from_b_safely(bvec_i,
                                           m = m, max_lag = max_lag, ars = ars,
                                           dvec = dvec_i)
        }

        calculate_predictive_r2(bvec_i,
                                x = x_val, y = y_val,
                                first_day_index = first_day_index_val,
                                m = m, max_lag = max_lag, ars = 0)
      }, error = function(e) {
        print(e)
        NA_real_
      })
    }

    ## detect elbow point in graph of validation accuracy as a function
    ## of number of retained parameters (again)
    elbow_check_2 <- get_elbow_index(x = nonzero_coeffs_2, y = r2s_2)
    if (verbose) {
      cat("second round of dvec optimization complete. q selected:",
          elbow_check_2$best_value, "\n")
    }
    elbow_check_2
  }, error = function(e) {
    cat("Second round of dvec optimization collapsed due to",
        "lack of sparsity variation in searched q-neighborhood...\n")
    NA
  })
  best_q_ovr <- if (!is.na(elbow_check_2)) elbow_check_2$best_value else best_q1

  ## calculate optimal delta vector using optimal q found above
  final_delta_vec <- calculate_dvec_from_q(best_q_ovr)

  list(
    delta = final_delta_vec,
    elbow_check_1 = elbow_check_1,
    elbow_check_2 = elbow_check_2
  )
}

###############################################################################
###############################################################################
###############################################################################

## outer function; fits a sequence of hflmds models with increasing degrees of
## AR error correction until residuals are sufficiently whitened, then returns
## the B(s, t) parameters estimated by the first model with sufficiently white
## residuals
#' TODO: document this fully
hflmds <- function(x, y, center = FALSE, verbose = TRUE,
                   ar_lim_stat = 0.05, ar_lim_practical = 0.15,
                   max_ar_to_check = 10,
                   mode = "inferential",
                   n_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),                                    ## nolint
                   dvec_selection_method = "elbow") {                                                               ## nolint

  ## set global variables
  m <- 364
  max_lag <- 150

  ## perform both a statistical and practical test on PACF of residuals
  ## (up to lag log_10(n)); return TRUE if either test passes
  residual_whiteness <- function(residuals) {
    n <- length(residuals)

    ## Determine max lag h using the 10 * log10(N) heuristic
    h <- as.integer(10 * log10(n))

    ## Run ACF and PACF
    ## Note: acf() includes lag 0 (which is always 1), so we drop it with [-1]
    acf_vals <- as.numeric(acf(residuals, lag.max = h, plot = FALSE, na.action = na.pass)$acf)[-1]        ## nolint
    pacf_vals <- as.numeric(pacf(residuals, lag.max = h, plot = FALSE, na.action = na.pass)$acf)          ## nolint

    ## practical test
    practical_pass_pacf <- (!any(abs(pacf_vals) > ar_lim_practical))
    practical_pass_acf <- (!any(abs(acf_vals) > ar_lim_practical))
    practical_pass <- practical_pass_acf && practical_pass_pacf

    ## Calculate the individual lag significance threshold (two tailed 95%)
    neff <- sum(is.finite(residuals))
    threshold <- 1.96 / sqrt(neff)

    ## Count the number of exceedances
    acf_exc <- sum(abs(acf_vals) > threshold)
    pacf_exc <- sum(abs(pacf_vals) > threshold)

    ## Calculate global p-values: P(X >= x)
    ## pbinom with lower.tail = FALSE computes P(X > q), so we evaluate at (x - 1)                ## nolint
    ## We return 1 if x == 0, as pbinom(-1, ...) will return 1.
    p_val_acf <- pbinom(acf_exc - 1, size = h, prob = 0.05, lower.tail = FALSE)
    p_val_pacf <- pbinom(pacf_exc - 1, size = h, prob = 0.05, lower.tail = FALSE)                 ## nolint

    ## The binomial test passes if the probability of seeing this many
    ## exceedances by chance is greater than the user-specified ARprob tolerance
    binomial_pass <- (p_val_acf >= ar_lim_stat) && (p_val_pacf >= ar_lim_stat)

    ## if verbose, inform user about this test
    if (verbose) {
      cat("residual whiteness check results:\n",
          "  > max abs(pacf) value:", max(abs(pacf_vals)), "\n",
          "  > max abs(acf) value:", max(abs(acf_vals)), "\n",
          "  > practical pass?", practical_pass, "\n",
          "  > binomial acf/pacf probs:", p_val_acf, "/", p_val_pacf, "\n",
          "  > binomial pass:", binomial_pass)
    }

    ## Return TRUE if at least one condition is satisfied
    practical_pass || binomial_pass
  }

  ## check data validity
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  stopifnot(length(x) == length(y))

  ## reserve test data if indicated
  if (mode == "predictive") {
    trainfrac <- 0.8
    n_train <- floor(trainfrac * length(x) / (m + 1)) * (m + 1)
    x_train <- x[1 : n_train]
    y_train <- y[1 : n_train]
    x_test <- x[(n_train + 1) : length(x)]
    y_test <- y[(n_train + 1) : length(y)]
    first_day_test <- (n_train %% (m + 1)) + 1
  } else if (mode == "inferential") {
    x_train <- x
    y_train <- y
    x_test <- x
    y_test <- y
    first_day_test <- 1
  } else {
    stop("invalid mode provided to hflmds()...")
  }

  ars <- -1
  residuals_whitened <- FALSE
  while ((!residuals_whitened) && (ars < max_ar_to_check)) {
    ars <- ars + 1
    if (verbose) cat("Fitting HFLMDS model with AR(", ars, ") error modeling...\n")                                 ## nolint
    mod_i <- hflmds_inner(x = x_train, y = y_train, ars = ars,
                          m = m, max_lag = max_lag,
                          center = center, verbose = verbose,
                          n_cores = n_cores,
                          dvec_selection_method = dvec_selection_method)
    residuals_whitened <- tryCatch(
      expr = {
        residual_whiteness(mod_i$b_residuals)
      }, error = function(e) {
        cat("Error (below) encountered while checking residual whiteness...\n")
        print(e)
        stop()
      }
    )
    if (length(residuals_whitened) != 1L ||
          is.na(residuals_whitened)) {
      stop("Residual whiteness could not be determined.")
    }
  }

  if (!residuals_whitened) {
    warning(sprintf(
      "Residuals were not successfully whitened within max_ar_to_check = %d.",
      max_ar_to_check
    ))
  }

  ## calculate BETA ONLY predictive r2 on RAW data
  r2_score <- calculate_predictive_r2(
    bvec = mod_i$pred_beta,
    x = x_test,
    y = y_test,
    first_day_index = first_day_test,
    m = m,
    max_lag = max_lag,
    ars = 0
  )

  ## calculate BETA ONLY predictive kge on RAW data
  kge_score <- calculate_predictive_kge(
    bvec = mod_i$pred_beta,
    x = x_test,
    y = y_test,
    first_day_index = first_day_test,
    m = m,
    max_lag = max_lag,
    ars = 0
  )

  list(
    pred_beta = mod_i$pred_beta,
    pred_ars = ars,
    scores = list(kge = kge_score, r2 = r2_score)
  )
}

## inner function; fits an hflmds model with specified degree of AR process
## used to model residuals.
#' @return list of pred_beta, the estimated $B(s, t) coefficients and b_residuals,                                  ## nolint
#' which are the residuals of X %% b - y. Note that b != B. These residuals are
#' the ones present when the convovled b coefficients are used for prediction
#' TODO: document this fully
hflmds_inner <- function(x, y, ars, center = FALSE, verbose = TRUE, m, max_lag,
                         n_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),                              ## nolint
                         dvec_selection_method) {                                                                   ## nolint

  ## set hyperparameter optimization variables
  ## #FIXME: MAKE THESE ALL PASSED PARAMS - NO MORE HARDCODING
  n_grid <- 20
  n_gp <- 30
  w_h_range <- c(-20, 30)
  w_v_range <- c(-20, 30)
  hp_datasplit <- c(0.75, 0.25)

  ## set quantile thresholding variables
  nq1 <- 200
  nq2 <- 200
  max_q <- 1
  dvec_datasplit <- c(0.75, 0.25) ## nolint

  ## perform checks on x and y
  if (length(x) != length(y)) {
    stop("x and y must have equal lengths.")
  }
  if (length(x) %% (m + 1) != 0L) {
    stop("x and y must contain complete annual cycles.")
  }

  ## if indicated, remove the seasonal signal from x and y
  if (center) {
    xmat <- matrix(x, ncol = 365, byrow = TRUE)
    x_signal <- colMeans(xmat, na.rm = TRUE)
    xmat <- sweep(xmat, 2, x_signal, "-")

    ymat <- matrix(y, ncol = 365, byrow = TRUE)
    y_signal <- colMeans(ymat, na.rm = TRUE)
    ymat <- sweep(ymat, 2, y_signal, "-")

    x <- as.vector(t(xmat))
    y <- as.vector(t(ymat))
  }

  ## optimize hyperparameters
  init_hp_fit_info <- optimize_hyperparameters(m = m, max_lag = max_lag, ars = ars,                                 ## nolint
                                               n_grid = n_grid, n_gp = n_gp,
                                               w_h_range = w_h_range, w_v_range = w_v_range,                        ## nolint
                                               x = x, y = y, dvec = NULL, datasplit = hp_datasplit,                 ## nolint
                                               verbose = verbose, n_cores = n_cores)                                ## nolint
  hyperparams_init <- init_hp_fit_info$hyperparameters

  ## fit (an) initial model
  b_init <- perform_fit(m = m, max_lag = max_lag, ars = ars,
                         w_h = hyperparams_init$w_h, w_v = hyperparams_init$w_v,                                    ## nolint
                         x = x, y = y, first_day_index = 1, dvec = NULL)                                            ## nolint
  if (ars > 0) {
    tryCatch(expr = {
      b_init <- get_beta_from_b_safely(b_init,
                                       m = m,
                                       max_lag = max_lag,
                                       ars = ars,
                                       dvec = NULL)
    }, error = function(e) {
      stop("error non-stationary estimate produced during initial fit\n")
    })
  }

  ## optimize dvec using EITHER
  ##  1. elbow method (cheaper but more aggressive)
  ##  2. "sparsest model within standard error" method (less aggressive but more computationally expensive)         ## nolint
  if (dvec_selection_method == "elbow") {
    dvec_fit_info <- optimize_dvec_elbow(beta_est = b_init, m = m, max_lag = max_lag, ars = ars,                    ## nolint
                                         x = x, y = y, w_h = hyperparams_init$w_h, w_v = hyperparams_init$w_v,      ## nolint
                                         nq1 = nq1, nq2 = nq2, max_q = max_q,
                                         datasplit = dvec_datasplit,
                                         verbose = verbose, n_cores = n_cores)
  } else {
    if (dvec_selection_method != "within_se") {
      warning(paste0(
        "Invalid dvec_selection_method (", dvec_selection_method, ")",
        ". Defaulting to \"within_se\""
      ))
    }
    dvec_fit_info <- NA
    stop("NEW OPTIMIZE_DVEC METHOD NOT IMPLEMENTED")
  }
  optimal_dvec <- dvec_fit_info$delta

  ## refit hyperparameters using optimal dvec
  final_hp_fit_info <- optimize_hyperparameters(m = m, max_lag = max_lag, ars = ars,                                ## nolint
                                                n_grid = n_grid, n_gp = n_gp,
                                                w_h_range = w_h_range, w_v_range = w_v_range,                       ## nolint
                                                x = x, y = y, dvec = optimal_dvec, datasplit = hp_datasplit,        ## nolint
                                                verbose = verbose, n_cores = n_cores)                               ## nolint
  hyperparams_final <- final_hp_fit_info$hyperparameters

  cat(
    "\n\n performing final fit with params:\n",
    "  > w_h =", hyperparams_final$w_h, "\n",
    "  > w_v =", hyperparams_final$w_v, "\n",
    "  > optimal_dvec saved as \"./optimal_dvec_check.rda\""
  )
  save(optimal_dvec, file = "./optimal_dvec_check.rda")

  ## perform final model fit
  b_final <- perform_fit(m = m, max_lag = max_lag, ars = ars,
                         w_h = hyperparams_final$w_h, w_v = hyperparams_final$w_v,                                  ## nolint
                         x = x, y = y, first_day_index = 1, dvec = optimal_dvec)                                    ## nolint

  ## calculate the residuals of the convolved b coefficients when used on
  ## full x, y series
  prediction <- calculate_yhat(bvec = b_final, x = x, y = y,
                               first_day_index = 1,
                               m = m, max_lag = max_lag, ars = ars)
  y_hat <- prediction$y_hat
  y_trunc <- prediction$y_true

  b_residuals <- y_trunc - y_hat

  if (ars > 0) {
    tryCatch(expr = {
      b_final <- get_beta_from_b_safely(b_final,
                                        m = m, max_lag = max_lag,
                                        ars = ars, dvec = optimal_dvec)
    }, error = function(e) {
      stop("error non-stationary estimate produced during final fit\n")
    })
  }

  ## FIXME: make sure to return EVERYTHING so that this code is fully observable post-mortem ## nolint
  list(pred_beta = b_final, b_residuals = b_residuals)
}