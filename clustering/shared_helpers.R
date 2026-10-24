library(lattice)
library(RColorBrewer)
library(ggplot2)
library(ggplotify)
library(grid)
library(stats)
library(Matrix)

#' @param bvec is a flattened, lag-minor indexed vector of B(s, t) coefficients
#'
#' @return a (m + 1)-length vector representing the maximum effective lag on
#' each day of the t-cycle, capped at max_lag
calculate_delta <- function(bvec, m = 364, max_lag = 150) {
  d <- as.numeric(sapply(seq_len(m + 1), function(t) {
    b_t <- bvec[(((t - 1) * (max_lag)) + 1):(t * (max_lag))]
    rev(which(b_t != 0))[1]
  }))
  d[is.na(d)] <- 0  ## any day with all 0 coeffs will get NA above
  d
}

#' TODO: document this
smooth_delta <- function(dvec, diff_pen = 1, lambda_1 = 1000) {
  if (var(dvec) == 0) {
    return(dvec)
  }

  adaptive_lambda <- lambda_1 / sqrt(var(dvec))

  ident <- Diagonal(length(dvec))

  diff_mat <- diff(ident, differences = 1)

  smoothed_delta <- solve(ident + adaptive_lambda * crossprod(diff_mat), dvec)

  as.numeric(smoothed_delta)
}

## FIXME: check that I am doing the indexing right
#' @param bvec is a flattened, lag-minor indexed vector of B(s, t) coefficients
#'
#' @return a (m + 1)-length vector representing the total transfer on day t*
#' (i.e. the sum of all B(s, t + s) coefficients with t = t* and s \in S)
calculate_total_transfer <- function(bvec, m = 364, max_lag = 150) {
  as.numeric(sapply(seq_len(m + 1), function(t) {
    sum(sapply(seq_len(max_lag), function(s) {
      day_index <- t + s - 1 ## account for the fact that s = 1 corresponds to lag 0                  ## nolint
      true_day_index <- ((day_index - 1) %% (m + 1)) + 1    ## wrap the day index around the cycle    ## nolint
      bvec[((true_day_index - 1) * max_lag) + s]
    }))
  }))
}

#' @param bvec is a flattened, lag-minor indexed vector of B(s, t) coefficients
#'
#' @return a (m + 1)-length vector representing the lag index of the maximally-
#' valued B(s, t) on day t for each t
calculate_max_eff_lag_vec <- function(bvec, m = 364, max_lag = 150) {
  ## build matrix where mat[i, j] is B(j, i)
  beta_m <- matrix(bvec, nrow = m + 1, ncol = max_lag, byrow = TRUE)

  ## find and return vector with maximal index in each row (day)
  as.numeric(apply(beta_m, 1, function(row) which.max(row)))
}

calculate_plot_info <- function(
  bvec,
  title = "",
  ncolors = 20,
  m = 364,
  max_lag = 150,
  max_delta = NULL,
  max_tt = NULL,
  max_sme = NULL) {

  d_vec <- calculate_delta(bvec, m = m, max_lag = max_lag)
  d_vec_smooth <- smooth_delta(d_vec)
  tt_vec <- calculate_total_transfer(bvec, m = m, max_lag = max_lag)
  sme_vec <- calculate_max_eff_lag_vec(bvec, m = m, max_lag = max_lag)

  list(
    delta = d_vec,
    smooth_delta = d_vec_smooth,
    total_transfer = tt_vec,
    max_effective_lag = sme_vec,

    title = title,
    ncolors = ncolors,
    m = m,
    max_lag = max_lag,

    max_delta = if (is.null(max_delta)) max(d_vec) * 1.2 else max_delta,
    max_tt = if (is.null(max_tt)) max(tt_vec) * 1.2 else max_tt,
    max_sme = if (is.null(max_sme)) max(sme_vec) * 1.2 else max_sme
  )
}
