source("./model code/hflmdsINDEV.R")
library(stringr)

candidate_gc <- "5017_centered"
grid_key <- str_split(candidate_gc, "_")[[1]][1]
dszn_key <- str_split(candidate_gc, "_")[[1]][2]

## get info about fit performed on CC
load(paste0("./results/inf_model_", dszn_key, "_grid", grid_key, ".rda"))
init_coeffs <- model$pred_beta
init_scores <- model$scores
fixed_ars <- model$pred_ars

## get relevant data
x <- as.vector(t(readRDS(paste0("./data/processed/grid", grid_key, "_rain.rds"))))          ## nolint
y <- as.vector(t(readRDS(paste0("./data/processed/grid", grid_key, "_streamflow.rds"))))    ## nolint

## set global params
m <- 364
max_lag <- 150

## re-perform fit using ars = fixed_ars
new_beta_no_wh_prior <- hflmds_inner(
  x = x, y = y, ars = fixed_ars,
  center = (dszn_key == "centered"),
  verbose = TRUE, m = m, max_lag = max_lag, n_cores = 1
)

## view levelplot of validation R2 as a function of w_h and w_v
hp_search_grid <- readRDS("hp_search_result.rds")

levelplot(
  score ~ w_h * w_v,
  data = hp_search_grid
)


library(ggplot2)

best <- hp_search_grid[which.max(hp_search_grid$score), ]

ggplot(hp_search_grid, aes(w_h, w_v)) +
  geom_point(aes(colour = score), size = 2.5) +
  geom_point(data = best, shape = 1, size = 4.5, stroke = 1) +
  scale_colour_viridis_c() +
  coord_equal() +
  theme_minimal()



### find and print the logs for this fit
load("./deploy_to_cc/gc_index.rda")
