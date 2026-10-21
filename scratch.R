library(signal)
library(stats)
library(bssm)
source("./clustering/shared_helpers.R")

if (!exists("beta_mat")) {
  beta_mat <- readRDS("./clustering/precomputed_beta_mtx.rds")
}

example_strings <- list("4976_centered", "1126_uncentered", "5337_uncentered")

example_betas <- lapply(example_strings, function(id) beta_mat[id, ])

example_deltas <- lapply(example_betas, function(bvec) calculate_delta(bvec))

smooth_delta_sg <- function(dvec, obl_size, smooth_window) {
  # first, obliterate short spikes of length < 5
  dvec_no_spike <- runmed(dvec, k = obl_size)

  # smooth the resulting vector using rolling average
  dvec_smooth <- sgolayfilt(dvec_no_spike, p = 3, n = smooth_window)

  dvec_smooth[is.na(dvec_smooth)] <- dvec_no_spike[is.na(dvec_smooth)]

  dvec_smooth
}

smooth_delta_sg_iter <- function(dvec, n_iter = 30, obl_size = 35, smooth_window = 71) {   ## nolint
  for (i in seq_len(n_iter)) {
    dvec <- smooth_delta_sg(dvec, obl_size = obl_size, smooth_window = smooth_window)       ## nolint
  }
  dvec
}

smooth_delta_fast_spline <- function(dvec, obl_size = 35, smooth_penalty = 0.8) { ## nolint
  # 1. Remove spikes once
  dvec_no_spike <- stats::runmed(dvec, k = obl_size)

  # 2. Fit a smoothing spline
  # spar controls smoothness (typically 0 to 1). Closer to 1 = heavier smoothing. ## nolint
  spline_fit <- smooth.spline(seq_along(dvec_no_spike), dvec_no_spike, spar = smooth_penalty)  ## nolint

  # 3. Extract the smoothed y values
  dvec_smooth <- spline_fit$y

  dvec_smooth
}

smooth_delta_fast_gaussian <- function(dvec, obl_size = 35, sigma = 20) {
  front_pad_len <- (3 * sigma) - 1
  back_pad_len <- (3 * sigma)

  dvec <- c(tail(dvec, front_pad_len), dvec, head(dvec, back_pad_len))

  # 1. Remove spikes once
  dvec_no_spike <- stats::runmed(dvec, k = obl_size)

  # 2. Create a Gaussian kernel
  # sigma controls the width of the smoothing (equivalent to more iterations)
  kernel_size <- ceiling(sigma * 6)
  x <- seq(-kernel_size / 2, kernel_size / 2, length.out = kernel_size)
  gauss_kernel <- exp(-0.5 * (x / sigma)^2)
  gauss_kernel <- gauss_kernel / sum(gauss_kernel) # Normalize to sum to 1

  # 3. Convolve the signal with the kernel
  dvec_smooth <- stats::filter(dvec_no_spike, gauss_kernel, sides = 2)

  cat("original length:", length(dvec), "\n")
  cat("sigma:", sigma, "\n")
  cat("num NAs at front:", which(!is.na(dvec_smooth))[1] - 1, "\n")
  cat("num NAs at back:", length(dvec_smooth) - rev(which(!is.na(dvec_smooth)))[1], "\n") ## nolint

  dvec_smooth <- dvec_smooth[!is.na(dvec_smooth)]
  cat("length of returning vec:", length(dvec_smooth), "\n\n")

  as.numeric(dvec_smooth)
}

for (example_index in seq_along(example_strings)) {
  beta_id <- example_strings[[example_index]]

  delta_vec <- example_deltas[[example_index]]

  pdf(paste0("smoothing_vis_id_", beta_id, ".pdf"))

  plot(delta_vec, type = "l")
  lines(smooth_delta_sg_iter(delta_vec), col = "red")
  lines(smooth_delta_fast_spline(delta_vec, smooth_penalty = 0.8), col = "blue")
  lines(smooth_delta_fast_spline(delta_vec, smooth_penalty = 0.95), col = "green") ## nolint
  lines(smooth_delta_fast_gaussian(delta_vec, sigma = 20), col = "purple") ## nolint
  lines(smooth_delta_fast_gaussian(delta_vec, sigma = 50), col = "yellow") ## nolint

  dev.off()
}