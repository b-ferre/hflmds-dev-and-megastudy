library(signal)
library(stats)
library(bssm)
source("./clustering/shared_helpers.R")

beta_mat <- readRDS("./clustering/precomputed_beta_mtx.rds")

example_1_str <- "4976_centered"
example_2_str <- "1126_uncentered"

beta_1 <- beta_mat[example_1_str, ]
beta_2 <- beta_mat[example_2_str, ]

delta_1 <- calculate_delta(beta_1)
delta_2 <- calculate_delta(beta_2)

smooth_delta <- function(dvec, obl_size = 7, smooth_window = 18) {
  # first, obliterate short spikes of length < 5
  dvec_no_spike <- runmed(dvec, k = obl_size)

  # smooth the resulting vector using rolling average
  window_weights <- rep(1 / smooth_window, smooth_window)
  dvec_smooth <- stats::filter(dvec_no_spike, window_weights, sides = 2)

  dvec_smooth
}

plot(delta_1, type = "l", col = "blue")
lines(smooth_delta(delta_1), col = "red")

plot(delta_2, type = "l", col = "blue")
lines(smooth_delta(delta_2), col = "red")
