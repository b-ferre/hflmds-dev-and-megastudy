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

sup_file_dir <- paste0("./grid", catchment_to_model, "_results/")

for (centering in c(TRUE, FALSE)) {
  for (opt_metric in c("r2")) {                       ## consider trying kge too
    for (gp_norm_mode in c("1", "2", "3")) {
      kern_i <- hflmds_inner(
        x = x_trn, y = y_trn, ars = ars, center = centering, verbose = TRUE,
        m = m, max_lag = max_lag, supervise = TRUE, sup_file_dir = sup_file_dir,
        n_grid = n_grid, n_gp = n_gp,
        w_h_range = w_h_range, w_v_range = w_v_range,
        hp_datasplit = hp_datasplit,
        nq1 = nq1, nq2 = nq2, max_q = 1, dvec_datasplit = dvec_datasplit,
        opt_metric = opt_metric, gp_norm_mode = gp_norm_mode
      )
    }
  }
}