tryCatch(expr = setwd("./clustering"), error = function(e) {})

library(stringr)
library(lattice)
library(RColorBrewer)                                                                                             ## nolint
library(stats)

source("./shared_helpers.R")

## define global vars
m <- 364
max_lag <- 150
reindex_to_water_year <- TRUE

## get a list of all gridcodes that we have map data for
metadata <- read.csv("../data/raw/Metadata/merged_catchments_metadata.csv")
mapped_gc <- metadata$gridcode

## get a list of all gridcodes that we have estimated B(s, t) for using both
## the centered and uncentered model and in both inference and predictive mode
beta_est_file_dir <- "../results/"
all_beta_f <- list.files(beta_est_file_dir)
get_gc <- function(fn) str_remove_all(fn, "[^\\d]")

## FIXME: uncomment below once predictive testing is done
# centered_pred_gc <- as.numeric(sapply(all_beta_f[str_detect(all_beta_f, "pred_model_centered")], get_gc))       ## nolint
# centered_inf_gc <- as.numeric(sapply(all_beta_f[str_detect(all_beta_f, "inf_model_centered")], get_gc))         ## nolint
# uncentered_pred_gc <- as.numeric(sapply(all_beta_f[str_detect(all_beta_f, "pred_model_uncentered")], get_gc))   ## nolint
# uncentered_inf_gc <- as.numeric(sapply(all_beta_f[str_detect(all_beta_f, "inf_model_uncentered")], get_gc))     ## nolint
# all_centered_gc <- intersect(centered_pred_gc, centered_inf_gc)   ## nolint
# all_uncentered_gc <- intersect(uncentered_inf_gc, uncentered_pred_gc) ## nolint
# estimated_gc <- intersect(all_centered_gc, all_uncentered_gc)   ## nolint

centered_gc <- as.numeric(sapply(all_beta_f[str_detect(all_beta_f, "_centered_")], get_gc))   ## nolint
uncentered_gc <- as.numeric(sapply(all_beta_f[str_detect(all_beta_f, "_uncentered_")], get_gc)) ## nolint
estimated_gc <- intersect(centered_gc, uncentered_gc)

## find intersection between mappable gridcodes and modeled gridcodes
grids_to_cluster <- intersect(mapped_gc, estimated_gc)

## create list of all model identifiers of the form "gridID_deseasonBool"
model_idents <- as.vector(sapply(grids_to_cluster, function(g) {
  c(paste0(g, "_centered"), paste0(g, "_uncentered"))
}))
n <- length(model_idents)
cat("Running Preprocessing for", n, "Unique Models...\n")

## get latitude and longitude df for all clusterable models
map_data <- data.frame(
  model_identifier = model_idents,
  latitude = rep(0, n),
  longitude = rep(0, n),
  hemisphere = rep(0, n)
)
cat("Building Map DataFrame...\n")
p <- txtProgressBar(min = 0, max = n, style = 3)
for (identifier in model_idents) {
  gc_only_i <- str_split(identifier, "_")[[1]][1]
  latitude_i <- metadata[metadata$gridcode == gc_only_i, "latitude"]
  longitude_i <- metadata[metadata$gridcode == gc_only_i, "longitude"]
  hemisphere_i <- if (latitude_i > 0) "N" else "S"
  map_data[map_data$model_identifier == identifier, ] <- list(
    model_identifier = identifier,
    latitude = latitude_i,
    longitude = longitude_i,
    hemisphere = hemisphere_i
  )
  setTxtProgressBar(p, which(model_idents == identifier))
}
close(p)


###############################################################################
#################### FEATURE MATRIX PRE-PRECOMPUTE ############################
###############################################################################
## pre-populate beta matrix for fast beta vector lookup
cat("Constructing Beta Matrix...\n")
p <- txtProgressBar(min = 0, max = n, style = 3)
beta_mtx <- matrix(nrow = n, ncol = max_lag * (m + 1), dimnames = list(model_idents, NULL))                         ## nolint
for (i in seq_len(nrow(beta_mtx))) {
  model_identifier_i <- rownames(beta_mtx)[i]
  gridcode_i <- str_split(model_identifier_i, "_")[[1]][1]
  dszn_index_i <- str_split(model_identifier_i, "_")[[1]][2]
  load(paste0(beta_est_file_dir, "inf_model_", dszn_index_i, "_grid", gridcode_i, ".rda"))                          ## nolint
  beta_mtx[model_identifier_i, ] <- model$pred_beta                                                                 ## nolint
  setTxtProgressBar(p, i)
}
close(p)

## re-index beta matrix to the water year (AMETSOC) if indicated
## note that this change will cascade into all other vectors/plots
## as everythinggets beta from the matrix
if (reindex_to_water_year) {
  cat("Reindexing Beta Matrix to Hydrological Year Indices...\n")
  p <- txtProgressBar(min = 0, max = n, style = 3)
  for (i in seq_len(nrow(beta_mtx))) {
    model_identifier_i <- rownames(beta_mtx)[i]
    hemisphere_i <- map_data[map_data$model_identifier == model_identifier_i, "hemisphere"]                         ## nolint
    reindex_day_mask <- if (hemisphere_i == "N") c(274:(m + 1), 1:273) else c(182:(m + 1), 1:181)                   ## nolint
    bmat <- matrix(beta_mtx[model_identifier_i, ], nrow = max_lag, ncol = m + 1)
    beta_mtx[model_identifier_i, ] <- as.numeric(bmat[, reindex_day_mask])
    setTxtProgressBar(p, i)
  }
  close(p)
}

## pre-populate delta matrix for fast distance calculations
## pre-populate beta matrix for fast beta vector lookup
cat("Constructing Delta Matrix...\n")
p <- txtProgressBar(min = 0, max = n, style = 3)
delta_mtx <- matrix(nrow = n, ncol = m + 1, dimnames = list(model_idents, NULL))                                    ## nolint
for (i in seq_len(nrow(beta_mtx))) {
  model_identifier_i <- rownames(beta_mtx)[i]
  beta_i <- beta_mtx[model_identifier_i, ]
  delta_mtx[model_identifier_i, ] <- calculate_delta(beta_i, m = m, max_lag = max_lag)                              ## nolint
  setTxtProgressBar(p, i)
}
close(p)

## pre-populate smoothed delta matrix for fast distance calculations
## pre-populate beta matrix for fast beta vector lookup
cat("Constructing Smoothed Delta Matrix...\n")
p <- txtProgressBar(min = 0, max = n, style = 3)
smooth_delta_mtx <- matrix(nrow = n, ncol = m + 1, dimnames = list(model_idents, NULL))                             ## nolint
for (i in seq_len(nrow(beta_mtx))) {
  model_identifier_i <- rownames(beta_mtx)[i]
  delta_vec_i <- delta_mtx[model_identifier_i, ]
  smooth_delta_mtx[model_identifier_i, ] <- smooth_delta(delta_vec_i)
  setTxtProgressBar(p, i)
}
close(p)

## pre-populate total transfer matrix for fast distance calculations
## pre-populate beta matrix for fast beta vector lookup
cat("Constructing Total Transfer Matrix...\n")
p <- txtProgressBar(min = 0, max = n, style = 3)
tt_mtx <- matrix(nrow = n, ncol = m + 1, dimnames = list(model_idents, NULL))                                       ## nolint
for (i in seq_len(nrow(beta_mtx))) {
  model_identifier_i <- rownames(beta_mtx)[i]
  beta_i <- beta_mtx[model_identifier_i, ]
  tt_mtx[model_identifier_i, ] <- calculate_total_transfer(beta_i, m = m, max_lag = max_lag)                        ## nolint
  setTxtProgressBar(p, i)
}
close(p)

## pre-populate maximally effective lag matrix for fast distance calculations
## pre-populate beta matrix for fast beta vector lookup
cat("Constructing SME Matrix...\n")
p <- txtProgressBar(min = 0, max = n, style = 3)
sme_mtx <- matrix(nrow = n, ncol = m + 1, dimnames = list(model_idents, NULL))                                       ## nolint
for (i in seq_len(nrow(beta_mtx))) {
  model_identifier_i <- rownames(beta_mtx)[i]
  beta_i <- beta_mtx[model_identifier_i, ]
  sme_mtx[model_identifier_i, ] <- as.numeric(sapply(seq_len(m + 1), function(t) {                                   ## nolint
    which.max(beta_i[((max_lag * (t - 1)) + 1) : (max_lag * t)])
  }))
  setTxtProgressBar(p, i)
}
close(p)


###############################################################################
###########################  DISTANCE HELPERS #################################
###############################################################################


## define function for fast calculation of the pairwise distances between
## all rows of a provided matrix by using stats::dist wherever possible
fast_distance <- function(x, method) {
  good <- rowSums(!is.finite(x)) == 0L
  if (all(good)) return(as.matrix(stats::dist(x, method = method)))
  n <- nrow(x)
  ans <- matrix(NA_real_, n, n, dimnames = list(rownames(x), rownames(x)))
  if (any(good)) {
    ans[good, good] <- as.matrix(stats::dist(x[good, , drop = FALSE], method = method)) ## nolint
  }
  for (i in which(!good)) {
    for (j in seq_len(n)) {
      d <- x[i, ] - x[j, ]
      value <- switch(method, euclidean = sqrt(sum(d * d)),
                      manhattan = sum(abs(d)), maximum = max(abs(d)))
      ans[i, j] <- ans[j, i] <- value
    }
  }
  ans
}

## even faster calculation of pairwise distance between rows of a matrix
## using the identity |x-y| = |x| + |y| - 2
fast_euclidean_distance <- function(x) {
  norm_sq <- rowSums(x * x)

  d2 <- outer(norm_sq, norm_sq, "+") - 2 * tcrossprod(x)

  # adjust for flop rounding
  d2[d2 < 0] <- 0
  diag(d2) <- 0

  as.matrix(as.dist(sqrt(d2)))
}

## helper function which takes a pairwise distance matrix of chebyshev
## distances, d, where d[i, j] is the cheby dist between rows i and j
## of the provided matrix x and returns a pairwise distance matrix of
## scaled chebyshev distances calculated as:
##                0.5 * d/mean(x_i) + 0.5 * d/mean(x_j)                                                             ## nolint
scale_chebyshev <- function(d, x) {
  inverse_mean <- 1 / rowMeans(x)
  0.5 * sweep(d, 1, inverse_mean, `*`) +
    0.5 * sweep(d, 2, inverse_mean, `*`)
}


###############################################################################
########################## DISTANCE CALCULATIONS ##############################
###############################################################################
all_feat_names <- c(
  "delta_l2",
  "delta_l1",
  "delta_cheby",
  "delta_scaled_cheby",

  "smooth_delta_l2",
  "smooth_delta_l1",
  "smooth_delta_cheby",
  "smooth_delta_scaled_cheby",

  "sme_l2",
  "sme_l1",
  "sme_cheby",
  "sme_scaled_cheby",

  "tt_l2",
  "tt_l1",
  "tt_cheby",
  "tt_scaled_cheby",

  "global_l2"
)

precomputed_feature_dist_mats <- setNames(vector("list", length(all_feat_names)),                                   ## nolint
                                          all_feat_names)


cat("Computing Pairwise Delta Distances...\n")
precomputed_feature_dist_mats$delta_l2 <- fast_distance(delta_mtx, "euclidean")                                     ## nolint
precomputed_feature_dist_mats$delta_l1 <- fast_distance(delta_mtx, "manhattan")                                     ## nolint
precomputed_feature_dist_mats$delta_cheby <- fast_distance(delta_mtx, "maximum")                                    ## nolint
precomputed_feature_dist_mats$delta_scaled_cheby <- scale_chebyshev(
  precomputed_feature_dist_mats$delta_cheby, delta_mtx
)

cat("Computing Pairwise Smooth Delta Distances...\n")
precomputed_feature_dist_mats$smooth_delta_l2 <- fast_distance(smooth_delta_mtx, "euclidean")                                     ## nolint
precomputed_feature_dist_mats$smooth_delta_l1 <- fast_distance(smooth_delta_mtx, "manhattan")                                     ## nolint
precomputed_feature_dist_mats$smooth_delta_cheby <- fast_distance(smooth_delta_mtx, "maximum")                                    ## nolint
precomputed_feature_dist_mats$smooth_delta_scaled_cheby <- scale_chebyshev(
  precomputed_feature_dist_mats$delta_cheby, smooth_delta_mtx
)

cat("Computing Pairwise SME Distances...\n")
precomputed_feature_dist_mats$sme_l2 <- fast_distance(sme_mtx, "euclidean")                                         ## nolint
precomputed_feature_dist_mats$sme_l1 <- fast_distance(sme_mtx, "manhattan")                                         ## nolint
precomputed_feature_dist_mats$sme_cheby <- fast_distance(sme_mtx, "maximum")                                        ## nolint
precomputed_feature_dist_mats$sme_scaled_cheby <- scale_chebyshev(
  precomputed_feature_dist_mats$sme_cheby, sme_mtx
)

cat("Computing Pairwise Total Transfer Distances...\n")
precomputed_feature_dist_mats$tt_l2 <- fast_distance(tt_mtx, "euclidean")                                           ## nolint
precomputed_feature_dist_mats$tt_l1 <- fast_distance(tt_mtx, "manhattan")                                           ## nolint
precomputed_feature_dist_mats$tt_cheby <- fast_distance(tt_mtx, "maximum")                                          ## nolint
precomputed_feature_dist_mats$tt_scaled_cheby <- scale_chebyshev(
  precomputed_feature_dist_mats$tt_cheby, tt_mtx
)

cat("Computing Overall Pairwise Euclidean Distances...\n")
precomputed_feature_dist_mats$global_l2 <- fast_euclidean_distance(beta_mtx)

## check that distance matrices are symmetric
stopifnot(!any(!unlist(lapply(precomputed_feature_dist_mats, isSymmetric))))                                        ## nolint

## normalize all components into [0, 1] to put tuning weights on the same scale
for (mat_index in seq_along(precomputed_feature_dist_mats)) {
  temp_mat <- precomputed_feature_dist_mats[[mat_index]]
  scaling_factor <- max(temp_mat, na.rm = TRUE)
  if (scaling_factor > 0) temp_mat <- temp_mat / scaling_factor
  precomputed_feature_dist_mats[[mat_index]] <- temp_mat
}

## check that normalization was done correctly
stopifnot(!any(unlist(lapply(precomputed_feature_dist_mats, function(mat) {
  (max(mat) > (1 + 1e-12))
}))))

## check that map_data$model_identifier is still equal to rownames/colnames,
## that rownames and colnames are still equal to each other, and that both
## are still equal to model_idents
stopifnot(!any(unlist(lapply(precomputed_feature_dist_mats, function(mat) {
  !setequal(rownames(mat), map_data$model_identifier) ||
    !setequal(colnames(mat), map_data$model_identifier) ||
    !setequal(colnames(mat), rownames(mat)) ||
    !setequal(colnames(mat), model_idents)
}))))

## check that there are no NA's, NULL's, etc. in any matrix
stopifnot(!any(unlist(lapply(precomputed_feature_dist_mats, function(m) any(!is.finite(m))))))                      ## nolint

###############################################################################
###################### PLOT & SCORE DF GENERATION #############################
###############################################################################
# pre-plot all B(s, t) plots
precomputed_beta_plots <- setNames(vector("list", length(model_idents)),
                                   as.character(model_idents))
precomputed_score_df <- data.frame(
  model_identifier = model_idents,
  # pred_kge = rep(NA, n),                                                                                          ## nolint
  # pred_r2 = rep(NA, n),                                                                                           ## nolint
  inf_kge = rep(NA, n),
  inf_r2 = rep(NA, n)
)

cat("Pre-Plotting Beta Surface Plots...\n")
p <- txtProgressBar(min = 0, max = length(model_idents), style = 3)
for (identifier in model_idents) {

  grid_i <- str_split(identifier, "_")[[1]][1]
  dszn_key_i <- str_split(identifier, "_")[[1]][2]

  ## get predictive scores FIXME: uncomment this stuff once pred results are in
  # load(paste0(beta_est_file_dir, "pred_model_", dszn_key_i, "_grid", grid_i, ".rda"))                             ## nolint
  # precomputed_score_df[precomputed_score_df$model_identifier == identifier, "pred_r2"] <- model$scores$r2         ## nolint
  # precomputed_score_df[precomputed_score_df$model_identifier == identifier, "pred_kge"] <- model$scores$kge       ## nolint
  pred_r2 <- NA_real_ #round(model$scores$r2, 3)                                                                    ## nolint
  pred_kge <- NA_real_ #round(model$scores$kge, 3)                                                                  ## nolint

  ## get inferential scores
  load(paste0(beta_est_file_dir, "inf_model_", dszn_key_i, "_grid", grid_i, ".rda"))                                ## nolint
  precomputed_score_df[precomputed_score_df$model_identifier == identifier, "inf_r2"] <- model$scores$r2            ## nolint
  precomputed_score_df[precomputed_score_df$model_identifier == identifier, "inf_kge"] <- model$scores$kge          ## nolint
  inf_r2 <- round(model$scores$r2, 3)
  inf_kge <- round(model$scores$kge, 3)

  ## build title string
  title_string <- paste0("GRID: ", grid_i, "\n",
                         "(Predictive)   KGE: ", pred_kge,
                         "   R2: ", pred_r2, "\n",
                         "(Inferential)   KGE: ", inf_kge,
                         "   R2: ", inf_r2, "\n")

  precomputed_beta_plots[[identifier]] <- calculate_plot_info(beta_mtx[identifier, ],                               ## nolint
                                                              title = title_string)                                 ## nolint
  setTxtProgressBar(p, which(model_idents == identifier))
}
close(p)

###############################################################################
############################# ASSET SAVING ####################################
###############################################################################

## map data
cat("Saving Map Data...\n")
saveRDS(map_data, file = "map_data.rds")

## precomputed plots
cat("Saving Plots...\n")
saveRDS(precomputed_beta_plots, file = "precomputed_beta_plots.rds")

## beta matrix
cat("Saving Beta Matrix...\n")
saveRDS(beta_mtx, file = "precomputed_beta_mtx.rds")

## pairwise distance matrices
cat("Saving Precomputed Feature Distance Matrices...\n")
saveRDS(precomputed_feature_dist_mats, file = "precomputed_feature_dist_mtxs.rds")                                  ## nolint

## score df
cat("Saving Score Data Frame...\n")
saveRDS(precomputed_score_df, file = "precomputed_score_df.rds")
