tryCatch(setwd("zzz"), error = function(e) {})

source("../model code/hflmdsINDEV.R")
source("./plotters.R")
library(stringr)
library(ggplot2)
library(lattice)

## set global params
m <- 364
max_lag <- 150

rerunfit <- FALSE
if (rerunfit) {
  candidate_gc <- "5017_centered"
  grid_key <- str_split(candidate_gc, "_")[[1]][1]
  dszn_key <- str_split(candidate_gc, "_")[[1]][2]

  ## get info about fit performed on CC
  load(paste0("../results/inf_model_", dszn_key, "_grid", grid_key, ".rda"))
  init_coeffs <- model$pred_beta
  init_scores <- model$scores
  fixed_ars <- model$pred_ars

  ## get relevant data
  x <- as.vector(t(readRDS(paste0("../data/processed/grid", grid_key, "_rain.rds"))))          ## nolint
  y <- as.vector(t(readRDS(paste0("../data/processed/grid", grid_key, "_streamflow.rds"))))    ## nolint

  ## re-perform fit using ars = fixed_ars
  new_beta_no_wh_prior <- hflmds_inner(
    x = x, y = y, ars = fixed_ars,
    center = (dszn_key == "centered"),
    verbose = TRUE, m = m, max_lag = max_lag, n_cores = 1
  )
}

viewlplots <- FALSE
if (viewlplots) {
  ## view levelplot of validation R2/KGE as a function of w_h and w_v
  metric <- "r2"
  round_ <- "2"

  hp_search_grid <- readRDS(paste0("hp_search_result_", str_to_lower(metric), "_round", round_, ".rds"))     ## nolint

  best <- hp_search_grid[which.max(hp_search_grid$score), ]

  ggplot(hp_search_grid, aes(w_h, w_v)) +
    geom_point(aes(colour = score), size = 2.5) +
    geom_point(data = best, shape = 1, size = 4.5, stroke = 1) +
    scale_colour_viridis_c() +
    coord_equal() +
    theme_minimal() +
    ggtitle(paste0(str_to_upper(metric), " ROUND ", round_))
}

optimal_dvec <- rep(150, 365)
ars <- 3    ## this was the optimally chosen ars under original model for 5017

ngrid <- 4
grid <- expand.grid(
  w_h = seq(-10, 15, length.out = ngrid),
  w_v = seq(-10, 15, length.out = ngrid)
)
n <- nrow(grid)
grid <- cbind(grid, data.frame(
  score_dszn = rep(NA, n),
  score_raw = rep(NA, n)
))

p <- txtProgressBar(min = 0, max = n, style = 3)
for (i in seq_len(n)) {
  for (dszn in c(TRUE, FALSE)) {
    hyperparams_final <- list(w_h = grid[i, "w_h"], w_v = grid[i, "w_v"])

    n_train <- 365 * (38 * 0.5)

    candidate_gc <- "5017_centered"
    grid_key <- str_split(candidate_gc, "_")[[1]][1]
    dszn_key <- str_split(candidate_gc, "_")[[1]][2]
    x <- as.vector(t(readRDS(paste0("../data/processed/grid", grid_key, "_rain.rds"))))          ## nolint
    y <- as.vector(t(readRDS(paste0("../data/processed/grid", grid_key, "_streamflow.rds"))))    ## nolint

    if (dszn) {
      x <- as.numeric(t(sweep(matrix(x, ncol = 365, byrow = TRUE), 2, colMeans(matrix(x, ncol = 365, byrow = TRUE), na.rm = TRUE), "-")))     ## nolint
      y <- as.numeric(t(sweep(matrix(y, ncol = 365, byrow = TRUE), 2, colMeans(matrix(y, ncol = 365, byrow = TRUE), na.rm = TRUE), "-")))     ## nolint
    }

    x_tr <- x[1:n_train]
    y_tr <- y[1:n_train]
    x_tst <- x[(n_train + 1):length(x)]
    y_tst <- y[(n_train + 1):length(y)]

    b_final <- perform_fit(m = m, max_lag = max_lag, ars = ars,
                            w_h = hyperparams_final$w_h, w_v = hyperparams_final$w_v,                                                         ## nolint
                            x = x_tr, y = y_tr, first_day_index = 1, dvec = optimal_dvec)                                                     ## nolint

    b_final <- get_beta_from_b_safely(b_final,
                                      m = m, max_lag = max_lag,
                                      ars = ars, dvec = optimal_dvec)

    # plot_beta(b_final)

    prediction <- calculate_yhat(bvec = b_final, x = x_tst, y = y_tst,
                                 first_day_index = 1,
                                 m = m, max_lag = max_lag, ars = 0)
    y_hat <- prediction$y_hat
    y_trunc <- prediction$y_true
    score_key <- if (dszn) "score_dszn" else "score_raw"
    grid[i, score_key] <- KGE(sim = y_hat, obs = y_trunc)
  }
  setTxtProgressBar(p, i)
}
close(p)

levelplot(
  score_dszn ~ w_h * w_v,
  grid
)
levelplot(
  score_raw ~ w_h * w_v,
  grid
)

best_dszn <- which.max(grid$score_dszn)
best_dszn_beta <- get_beta_from_b_safely(
  perform_fit(m = m, max_lag = max_lag, ars = ars,
              w_h = grid[best_dszn, "w_h"], w_v = grid[best_dszn, "w_v"],
              x = x, y = y, first_day_index = 1, dvec = optimal_dvec),
  m = m, max_lag = max_lag,
  ars = ars, dvec = optimal_dvec
)
plot_beta(best_dszn_beta)

best_raw <- which.max(grid$score_raw)
best_raw_beta <- get_beta_from_b_safely(
  perform_fit(m = m, max_lag = max_lag, ars = ars,
              w_h = grid[best_raw, "w_h"], w_v = grid[best_raw, "w_v"],
              x = x, y = y, first_day_index = 1, dvec = optimal_dvec),
  m = m, max_lag = max_lag,
  ars = ars, dvec = optimal_dvec
)
plot_beta(best_raw_beta)
