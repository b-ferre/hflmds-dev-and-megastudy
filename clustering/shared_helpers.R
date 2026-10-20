library(lattice)
library(RColorBrewer)
library(ggplot2)
library(ggplotify)
library(grid)

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

#' @param bvec is a flattened, lag-minor indexed vector of B(s, t) coefficients
#'
#' @return a (m + 1)-length vector representing the total transfer on day t*
#' (i.e. the sum of all B(s, t) coefficients with t = t*)
calculate_total_transfer <- function(bvec, m = 364, max_lag = 150) {
  as.numeric(sapply(seq_len(m + 1), function(t) {
    sum(bvec[(((t - 1) * (max_lag)) + 1):(t * max_lag)])
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

plot_beta_and_derived_feats <- function(
  bvec,
  title = "",
  ncolors = 20,
  m = 364,
  max_lag = 150,
  max_delta = NULL,
  max_tt = NULL,
  max_sme = NULL) {

  d_vec <- calculate_delta(bvec, m = m, max_lag = max_lag)
  tt_vec <- calculate_total_transfer(bvec, m = m, max_lag = max_lag)
  sme_vec <- calculate_max_eff_lag_vec(bvec, m = m, max_lag = max_lag)

  list(
    delta = d_vec,
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

#' @param bvec is a flattened, lag-minor indexed vector of B(s, t) coefficients
#' @param title is a string to display above the levelplot as a plot title
#' @param ncolors is the number of unique non-white colors to use in the palette
#'
#' @return a levelplot() object of B(s, t) values with x = day and y = lag
plot_beta <- function(bvec, title = "", ncolors = 20, m = 364, max_lag = 150) {
  bgrid <- cbind(expand.grid(s = seq_len(max_lag), t = seq_len(m + 1)), b = bvec)           ## nolint
  dvec <- calculate_delta(bvec)

  # truncate matrix to avoid plotting tens of lags with B(s, t) = 0
  cutoff <- min(ceiling(max(dvec) / 10) * 10, max_lag)
  bgrid <- bgrid[bgrid$s <= cutoff, , drop = FALSE]

  # Generate colour key: smallest magnitudes are still distinct from white
  blue_pal <- function(n) {
    colorRampPalette(c("#d5ddfe", "#6ca4fe", "#2159e7"))(n)
  }
  red_pal <- function(n) {
    colorRampPalette(c("#f9aa98", "#d85a5e", "#67000D"))(n)
  }

  has_pos <- any(bgrid$b > 0)
  has_neg <- any(bgrid$b < 0)

  if (has_pos && !has_neg) {
    cols <- blue_pal(ncolors)
    at_vals <- seq(0, max(bgrid$b), length.out = ncolors + 1)

  } else if (!has_pos && has_neg) {
    cols <- rev(red_pal(ncolors))
    at_vals <- seq(min(bgrid$b), 0, length.out = ncolors + 1)

  } else if (has_pos && has_neg) {
    rp <- max(bgrid$b)
    rn <- abs(min(bgrid$b))

    # Ensure both signs receive at least one colour
    ncolors <- max(2L, ncolors)
    nb <- max(1L, min(ncolors - 1L, round(ncolors * rp / (rp + rn))))
    nr <- ncolors - nb

    cols <- c(rev(red_pal(nr)), blue_pal(nb))
    at_vals <- c(
      seq(-rn, 0, length.out = nr + 1),
      seq(0, rp, length.out = nb + 1)[-1]
    )

  } else {
    cols <- "#FFFFFF"
    at_vals <- c(-1, 1)
  }

  bgrid$b[bgrid$b == 0] <- NA_real_

  levelplot(b ~ t * s,
            bgrid,
            col.regions = cols,
            at = at_vals,
            xlab = list("Day of the year", cex = 1.4),
            ylab = list("Lag (Days)", cex = 1.4),
            scales = list(cex = 1.2),
            colorkey = list(axis.text = list(cex = 1.2)),
            aspect = "xy",
            main = title)
}

#' @param bvec is a flattened, lag-minor indexed vector of B(s, t) coefficients
#' @param title is a string to display above the levelplot as a plot title
#' @param ncolors is the number of unique non-white colors to use in the palette
#'
#' @return a levelplot() object of B(s, t) values with x = day and y = lag
plot_beta_and_derived_feats_OLD <- function(bvec,
                                        title = "",
                                        ncolors = 20,
                                        m = 364,
                                        max_lag = 150,
                                        max_delta = NULL,
                                        max_tt = NULL,
                                        max_sme = NULL) {
  beta_plot <- as.ggplot(plot_beta(bvec, title = "Beta Surface", ncolors = ncolors, m = m, max_lag = max_lag))      ## nolint

  day_vec <- seq_len(m + 1)

  d_vec <- calculate_delta(bvec = bvec, m = m, max_lag = max_lag)
  if (is.null(max_delta)) {
    max_delta <- max(d_vec) * 1.2
  }
  d_data <- data.frame(
    x = day_vec,
    y = d_vec
  )
  d_plot <- ggplot(d_data, aes(x = x, y = y)) +                                                                     ## nolint
    geom_line() +
    coord_cartesian(ylim = c(0, max_delta)) +
    ggtitle("Delta Vector") +
    theme(plot.title = element_text(hjust = 0.5))

  tt_vec <- calculate_total_transfer(bvec = bvec, m = m, max_lag = max_lag)
  if (is.null(max_tt)) {
    max_tt <- max(tt_vec) * 1.2
  }
  tt_data <- data.frame(
    x = day_vec,
    y = tt_vec
  )
  tt_plot <- ggplot(tt_data, aes(x = x, y = y)) +                                                                   ## nolint
    geom_line() +
    coord_cartesian(ylim = c(0, max_tt)) +
    ggtitle("Total Transfer Vector") +
    theme(plot.title = element_text(hjust = 0.5))

  sme_vec <- calculate_max_eff_lag_vec(bvec = bvec, m = m, max_lag = max_lag)
  if (is.null(max_sme)) {
    max_sme <- max(sme_vec) * 1.2
  }
  sme_data <- data.frame(
    x = day_vec,
    y = sme_vec
  )
  sme_plot <- ggplot(sme_data, aes(x = x, y = y)) +                                                                   ## nolint
    geom_line() +
    coord_cartesian(ylim = c(0, max_sme)) +
    ggtitle("Maximally Effective Lag Vector") +
    theme(plot.title = element_text(hjust = 0.5))

  list(bplot = beta_plot,
       dplot = d_plot,
       ttplot = tt_plot,
       smeplot = sme_plot,
       overall_title_string = title)
}