tryCatch(setwd("./zzz"), error = function(e) {})

source("./hflmdsINDEV.R")
source("./plotters.R")
library(stringr)
library(ggplot2)
library(lattice)

catchment_to_model <- "5017"
m <- 364
max_lag <- 150

## get rainfall (x) and streamflow (y) data
x <- as.vector(t(readRDS(paste0("../data/processed/grid", catchment_to_model, "_rain.rds"))))         ## nolint
y <- as.vector(t(readRDS(paste0("../data/processed/grid", catchment_to_model, "_streamflow.rds"))))   ## nolint
stopifnot((length(x) == length(y)) && ((length(x) %% (m + 1)) == 0))

## TODO: generalize
## there are 39 years for catchment 5017, use the last 6 as a test set
n_train <- ((length(x) / (m + 1)) - 6) * (m + 1)
x_trn <- x[1:(n_train)]
y_trn <- y[1:(n_train)]
x_test <- x[(n_train + 1):length(x)]
y_test <- y[(n_train + 1):length(y)]

## TODO: generalize:
## set the AR order to the one used in the CC run for catchment 5017
ars <- 3

## set global cross-run parameters
n_grid <- 20
n_gp <- 50
w_h_range <- c(-20, 20)
w_v_range <- c(-20, 20)
hp_datasplit <- c(0.75, 0.25)

nq1 <- 350
nq2 <- 100
max_q <- 1
dvec_datasplit <- c(0.75, 0.25)

sup_file_dir <- paste0("./grid", catchment_to_model, "_results_small_search/")   ## FIXME: remove "_small_search"     ## nolint

for (centering in c(TRUE, FALSE)) {
  for (opt_metric in c("r2")) {                       ## consider trying kge too
    for (gp_norm_mode in c("1", "2", "3")) {
      ## load final hyperparameters and dvec
      dvc <- readRDS(paste0(
        sup_file_dir, "ctr_", centering, "_opt_", opt_metric, "_gmode_", gp_norm_mode, "/optimal_dvec.rds"            ## nolint
      ))
      hp_grid <- readRDS(paste0(
        sup_file_dir, "ctr_", centering, "_opt_", opt_metric, "_gmode_", gp_norm_mode, "/hp_search_results.rds"       ## nolint
      ))
      best_idx <- which.max(hp_grid$score)
      best_w_h <- hp_grid$w_h[best_idx]
      best_w_v <- hp_grid$w_v[best_idx]

      ## fit final kernel using appropriate training data, hp, and dvec
      kern <- perform_fit(m = m, max_lag = max_lag, ars = ars,
                          w_h = best_w_h, w_v = best_w_v,
                          x = x_trn, y = y_trn, first_day_index = 1, dvec = dvc)                                     ## nolint
      if (ars > 0) {
        tryCatch(expr = {
          kern <- get_beta_from_b_safely(kern,
                                         m = m,
                                         max_lag = max_lag,
                                         ars = ars,
                                         dvec = dvc)
        }, error = function(e) {
          stop("error non-stationary estimate produced during initial fit\n")
        })
      }

      this_kge <- calculate_predictive_kge(
        kern, x_test, y_test, first_day_index = 1, m = m, max_lag = max_lag, ars = 0
      )

      cat(
        paste0(
          "KGE of model with settings (",
          centering, ", ",
          opt_metric, ", ",
          "w_h :", best_w_h, ", ",
          "w_v :", best_w_v, ", ",
          gp_norm_mode, ") : ", this_kge, "\n"
        )
      )
    }
  }
  load(paste0("../results/pred_model_", if (centering) "centered" else "uncentered", "_grid", catchment_to_model, ".rda"))
  cat("\n\nFor reference, model from CC run achieved test KGE of", model$scores$kge, "\n\n")
}

