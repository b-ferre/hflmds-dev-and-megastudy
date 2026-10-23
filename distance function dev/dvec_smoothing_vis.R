tryCatch(expr = setwd("./distance function dev"), error = function(e) {})

library(signal)
library(stats)
library(bssm)
source("../clustering/shared_helpers.R")

if (!exists("beta_mat")) {
  beta_mat <- readRDS("../clustering/precomputed_beta_mtx.rds")
}

example_strings <- list(
  "4976_centered",
  "1126_uncentered",
  "5337_uncentered",
  "5016_uncentered",
  "2606_uncentered"
)
example_betas <- lapply(example_strings, function(id) beta_mat[id, ])
example_deltas <- lapply(example_betas, function(bvec) calculate_delta(bvec))


###############################################################################
iter_med_avg <- function(
  dvec,
  order = "med_avg",
  n_iter,
  med_window,
  avg_window
) {
  med_window <- round(length(dvec) * med_window)
  med_window <- med_window - (med_window %% 2) + 1

  avg_window <- round(avg_window * length(dvec))
  avg_weights <- rep(1, avg_window) / avg_window

  for (i in seq_len(n_iter)) {
    if (order == "med_avg") {
      dvec <- runmed(x = dvec, k = med_window)
      dvec <- stats::filter(
        x = dvec, filter = avg_weights, sides = 2, circular = TRUE
      )
    } else if (order == "avg_med") {
      dvec <- stats::filter(
        x = dvec, filter = avg_weights, sides = 2, circular = TRUE
      )
      dvec <- runmed(x = dvec, k = med_window)
    }
  }

  dvec
}

avg_med_avg_l <- function(
  dvec,
  avg_window_1,
  avg_window_2,
  med_window
) {
  if (avg_window_2 <= avg_window_1) {
    warning("incorrect use case for avg_med_avg_l")
  }

  med_window <- round(length(dvec) * med_window)
  med_window <- med_window - (med_window %% 2) + 1

  avg_window1 <- round(avg_window_1 * length(dvec))
  avg_weights1 <- rep(1, avg_window1) / avg_window1
  avg_window2 <- round(avg_window_2 * length(dvec))
  avg_weights2 <- rep(1, avg_window2) / avg_window2

  dvec <- stats::filter(
    x = dvec, filter = avg_weights1, sides = 2, circular = TRUE
  )
  dvec <- runmed(x = dvec, k = med_window)
  dvec <- stats::filter(
    x = dvec, filter = avg_weights2, sides = 2, circular = TRUE
  )

  dvec
}

library(Matrix)

whittaker_smooth <- function(y, lambda, d = 2) {
  n <- length(y)
  
  # Construct a sparse identity matrix
  I <- Diagonal(n)
  
  # Compute the d-th order difference matrix
  D <- diff(I, differences = d)
  
  # Solve the penalized least squares system (I + lambda * D^T D)z = y
  # crossprod(D) efficiently calculates D^T D
  z <- solve(I + lambda * crossprod(D), y)
  
  return(as.numeric(z))
}

## methods to compare:
##  - iterate (rolling median (x% window) -> rolling average (x% window))
##  - iterate (rolling median (x% window) -> rolling average (y% window)) with y > x                                            ## nolint
##  - iterate (rolling average (x% window) -> rolling median (x% window))
##  - iterate (rolling average (y% window) -> rolling median (x% window)) with y > x                                            ## nolint
##  - one round : rolling average (x% window) -> rolling median (x% window) -> rolling average (y% window) with y > x           ## nolint
##  - iterate : rolling median (x% window)

## best so far: 10 iters of 7% med -> 3% avg

method_ <- function(dvec, tuning_var) {
 whittaker_smooth(y = dvec, lambda = (tuning_var / sqrt(var(dvec))), d = 2)
}
dir_name <- "./whittaker_d2/"
if (!dir.exists(dir_name)) dir.create(dir_name)

for (example_index in seq_along(example_strings)) {

  beta_id <- example_strings[[example_index]]

  delta_vec <- example_deltas[[example_index]]
  cat("id:", beta_id, ", var:", var(delta_vec), "\n")

  pdf(paste0(dir_name, beta_id, ".pdf"))

  plot(delta_vec, type = "l", col = "blue")
  lines(method_(delta_vec, tuning_var = 10), col = " red")
  lines(method_(delta_vec, tuning_var = 100), col = "green")
  lines(method_(delta_vec, tuning_var = 1000), col = "purple")
  lines(method_(delta_vec, tuning_var = 10000), col = "orange")
  # lines(method_(delta_vec, tuning_var = 71), col = "pink")
  # lines(method_(delta_vec, tuning_var = 101), col = "black")

  legend(
    "bottomleft",
    legend = c(
      "Base Vector",
      "1",
      "5",
      "10",
      "20",
      "30",
      "50"
    ),
    col = c("blue", "red", "green", "purple", "orange", "pink", "black"),
    lty = 1
  )

  dev.off()
}