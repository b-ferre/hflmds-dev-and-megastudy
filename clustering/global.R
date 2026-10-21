library(shiny)
library(leaflet)
library(RColorBrewer)
library(grid)
library(lattice)
library(stringr)
#library(dbscan)
library(leafpop)
library(ggplot2)
library(later)
library(future)
library(ggplotify)
library(patchwork)
library(bslib)
library(reticulate)

source("shared_helpers.R")

###############################################################################
############################# GLOBAL RESOURCES ################################
###############################################################################
## define global variables
max_lag <- 150
m <- 364

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
  "tt_scaled_cheby"
)

all_feat_weight_names <- sapply(
  all_feat_names,
  function(s) {
    paste0(s, "_wt")
  }
)

INITCLUSTERSETTINGS <- c(                                                                                                             ## nolint
  setNames(as.list(
    c(
      0, 0, 0, 0,
      0.5, 0.2, 0.4, 1,
      0, 0, 0, 0,
      0, 0, 0, 0
    )
  ), all_feat_weight_names),                                                                                                          ## nolint
  list(
    comb_func = "test2",
    n_pts_core = 2,
    n_pts_clust = 10,
    deseasonalization_selection = "Wherever Meaningfully Better Than Uncentered",                                                     ## nolint
    dszn_thresh_selection = 0.15,
    min_kge_for_clust = 0.25
  )
)
debugmode <- TRUE

## load the map_data dataframe that contains 3 columns (lat, lon, grid)
map_data <- readRDS("map_data.rds")

grids_to_cluster <- map_data$gridcode

## load the list of precomputed pairwise feature distances
precomputed_feature_dist_mats <- readRDS("precomputed_feature_dist_mtxs.rds")

## load precomputed beta plots
precomputed_beta_plots <- readRDS("precomputed_beta_plots.rds")

## load beta matrix - each row is one estimated B(s, t) surface
## and rownames are gridcodes; rows of mtx are indexed to water year
precomputed_beta_matrix <- readRDS("precomputed_beta_mtx.rds")

## load scores dataframe
precomputed_score_df <- readRDS("precomputed_score_df.rds")

## get a list of unique GRIDCODES (not unique model_identifiers)
unique_gridcodes <- unique(
  as.vector(
    sapply(
      rownames(precomputed_beta_matrix),
      function(identifier_) {
        str_split(identifier_, "_")[[1]][1]
      }
    )
  )
)

###############################################################################
############################# HELPER FUNCTIONS ################################
###############################################################################

## convert a boolean representing whether the centered
## model should be used to a key used for identifying the
## model
to_dszn_key <- function(centered_bool) {
  dszn_data_keys <- list(
    "TRUE" = "centered",
    "FALSE" = "uncentered"
  )
  dszn_data_keys[[as.character(centered_bool)]]
}

#' @param gridcode is the gridcode of the catchment of the B(s, t) to get
#'
#' @return a flattened, lag-minor indexed vector of B(s, t) coefficients
#' with t = 1 corresponding to the first day of the hydrological year
get_beta <- function(gridcode, dszn_bool) {
  precomputed_beta_matrix[paste0(gridcode, "_", to_dszn_key(dszn_bool)), ]
}

#' @param gridcode is the gridcode of the catchment of the inf kge to get
#'
#' @return a 4 digit float representing inferential kge
get_inf_kge <- function(gridcode, dszn_bool) {
  round(precomputed_score_df$inf_kge[precomputed_score_df$model_identifier == paste0(gridcode, "_", to_dszn_key(dszn_bool))], 4)      ## nolint
}

#' @param gridcode is the gridcode of the catchment of the inf kge to get
#'
#' @return a 4 digit float representing inferential kge
get_inf_r2 <- function(gridcode, dszn_bool) {
  round(precomputed_score_df$inf_r2[precomputed_score_df$model_identifier == paste0(gridcode, "_", to_dszn_key(dszn_bool))], 4)       ## nolint
}

#' @param featwts is a vector of scalars by which each component of the overall
#' distance function should be weighted before combination
#' @param comb_func is one of "+", "*", "euc" and dictates how the weighted
#' feature-wise distance components are combined into an overall distance
#' @param dszn_bool_vec is a sequence of booleans - one for each unique
#' gridcode in the dataset - that represents whether the centered estimate
#' should be used for that gridcode
#'
#' @return a matrix with overall pairwise distances
get_pairwise_distance_mtx <- function(featwts, comb_func, dszn_bool_vec, kge_threshold = 0) {
  stopifnot(length(dszn_bool_vec) == length(unique_gridcodes))

  nfeats <- length(precomputed_feature_dist_mats)
  if (is.null(featwts)) {
    warning(paste0("NULL weight vector passed to get_pairwise_distance().",
                   "  Defaulting to a weight = 1 for each feature."))
    featwts <- rep(1, nfeats)
    names(featwts) <- names(precomputed_feature_dist_mats)                                                                            ## nolint
  }

  if (length(featwts) != nfeats) {
    warning(paste0("Length of weight vector passed to get_pairwise_distance()",
                   " not equal to the number of features to weight.",
                   "  Defaulting to a weight = 1 for each feature. (Mismatch: ",
                   length(featwts), " vs. ", nfeats, ")"))
    cat("   > weight names:", names(featwts), "\n",
        "   > mat names:", names(precomputed_feature_dist_mats), "\n")                                                                ## nolint
    featwts <- rep(1, nfeats)
    names(featwts) <- names(precomputed_feature_dist_mats)                                                                            ## nolint
  }

  ests_to_keep <- sapply(unique_gridcodes, function(gc) {
    paste0(gc, "_", to_dszn_key(dszn_bool_vec[which(unique_gridcodes == gc)]))
  })

  removed_due_to_kge <- 0
  ests_to_use <- ests_to_keep[sapply(ests_to_keep, function(ident_i) {
    inf_score_i <- precomputed_score_df$inf_kge[precomputed_score_df$model_identifier == ident_i]                              ## nolint
    retval <- inf_score_i > kge_threshold
    if (!retval) removed_due_to_kge <<- removed_due_to_kge + 1
    retval
  })]

  num_centered <- sum(as.numeric(sapply(
    ests_to_use,
    function(fn_) {
      str_split(fn_, "_")[[1]][2] == "centered"
    }
  )))
  cat("clusterable catchment list is done being trimmed. of the initial",
      length(unique_gridcodes), "available catchments,",
      removed_due_to_kge, "were removed due to low inferential KGE. of the remaining",  ## nolint
      length(ests_to_use), "models,", num_centered, "are deseasonalized\n")

  weighted_dist_components <- lapply(names(precomputed_feature_dist_mats),                                                            ## nolint
                                     function(name_j) {
                                       feat_mat_j <- precomputed_feature_dist_mats[[name_j]][ests_to_use, ests_to_use, drop = FALSE]  ## nolint
                                       featwts[[name_j]] * feat_mat_j
                                     })

  # Function to convert each entry in a matrix to its empirical CDF rank
  rank_matrix_cdf <- function(mat) {
    # Validate input
    if (!is.matrix(mat)) {
      stop("Input must be a matrix.")
    }
    # Flatten the matrix into a vector
    vals <- as.numeric(mat)
    # Handle NA values separately
    if (anyNA(vals)) {
      warning("NA values detected; they will remain NA in the output.")
    }
    # Compute empirical CDF ranks
    ranks <- rank(vals, ties.method = "average", na.last = "keep") / sum(!is.na(vals))                                                ## nolint
    # Reshape back to original matrix dimensions
    result <- matrix(ranks, nrow = nrow(mat), ncol = ncol(mat))
    result
  }

  if (comb_func == "test") {
    plus_mat <- Reduce("+", weighted_dist_components)
    new_mat <- rank_matrix_cdf(plus_mat)
    dimnames(new_mat) <- dimnames(plus_mat)
    return(new_mat) ## nolint
  } else if (comb_func == "test2") {
    weighted_dist_components <- lapply(seq_along(weighted_dist_components), function(i) {
      m_i <- weighted_dist_components[[i]]

      if (any(!is.finite(m_i))) {
        stop(paste0("ERROR: non-finite values in component", names(weighted_dist_components)[i]))
      }
      transf_m_i <- rank_matrix_cdf(m_i)
      dimnames(transf_m_i) <- dimnames(m_i)
      as.matrix(transf_m_i)
    })
    plus_mat <- Reduce("+", weighted_dist_components)
    new_mat <- rank_matrix_cdf(plus_mat)
    dimnames(new_mat) <- dimnames(plus_mat)
    return(new_mat) ## nolint"
  } else if (comb_func == "test3") {
    weighted_dist_components <- lapply(seq_along(weighted_dist_components), function(i) {
      m_i <- weighted_dist_components[[i]]

      if (any(!is.finite(m_i))) {
        stop(paste0("ERROR: non-finite values in component", names(weighted_dist_components)[i]))
      }
      transf_m_i <- rank_matrix_cdf(m_i)
      dimnames(transf_m_i) <- dimnames(m_i)
      as.matrix(transf_m_i)
    })
    new_mat <- Reduce("+", weighted_dist_components)
    return(new_mat) ## nolint"
  } else if (comb_func == "test4") {
    plus_mat <- Reduce("+", weighted_dist_components)
    new_mat <- (0.8 * rank_matrix_cdf(plus_mat)) + (0.2 * plus_mat)
    dimnames(new_mat) <- dimnames(plus_mat)
    return(new_mat) ## nolint
  } else {
    warning(paste0("invalid comb_func (\"", comb_func, "\") passed to",
                   "get_pairwise_distance(). Defaulting to \"+\"."))
    Reduce("+", weighted_dist_components)
  }
}

render_beta_and_derived_feats <- function(bvec, plot_data) {
  p <- plot_data

  beta_plot <- ggplotify::as.ggplot(
    plot_beta(    ## nolint
      bvec,
      title = "Beta Surface",
      ncolors = p$ncolors,
      m = p$m,
      max_lag = p$max_lag
    )
  )

  make_line_plot <- function(values, title, ymax) {
    # Avoid a zero-height plotting range
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1

    ggplot(
      data.frame(day = seq_along(values), value = values),
      aes(x = day, y = value) ## nolint
    ) +
      geom_line() +
      coord_cartesian(ylim = c(0, ymax)) +
      labs(title = title, x = "Day", y = NULL) +
      theme(plot.title = element_text(hjust = 0.5))
  }

  d_plot <- make_line_plot(
    p$delta, "Delta Vector", p$max_delta
  )
  smooth_d_plot <- make_line_plot(
    as.numeric(unname(p$smooth_delta)), "Smoothed Delta Vector", p$max_delta
  )
  tt_plot <- make_line_plot(
    p$total_transfer, "Total Transfer Vector", p$max_tt
  )
  sme_plot <- make_line_plot(
    p$max_effective_lag, "Maximally Effective Lag", p$max_sme
  )

  combined <- (
    beta_plot + (d_plot / smooth_d_plot / tt_plot / sme_plot)
  ) +
    patchwork::plot_layout(widths = c(2, 1))

  obj <- combined + patchwork::plot_annotation(title = p$title)
  obj
}

get_beta_plot <- function(model_identifier) {
  bvec <- precomputed_beta_matrix[as.character(model_identifier), ]

  plot_data <- NULL

  if (exists("precomputed_beta_plots")) {
    plot_data <- precomputed_beta_plots[[as.character(model_identifier)]]
  }

  if (is.null(plot_data)) {
    plot_data <- plot_beta_and_derived_feats(bvec)  ## nolint
  }

  obj <- render_beta_and_derived_feats(bvec, plot_data)
  obj
}


###############################################################################
############################## RUN CLUSTERING #################################
###############################################################################
#' TODO: document this
get_deseason_bool_vec_from_scheme <- function(scheme_key, threshold = NULL) {                                                         ## nolint
  if (debugmode) {
    threshold_str <- paste0(" and threshold: ", threshold)
    if (scheme_key %in% c("Everywhere", "Nowhere")) {
      threshold_str <- ""
    }
    cat("running get_deseason_bool_vec_from_scheme() with scheme_key:", scheme_key, threshold_str, "\n")                              ## nolint
  }

  if (scheme_key == "Everywhere") {
    ret_vec <- rep(TRUE, length(unique_gridcodes))
  } else if (scheme_key == "Nowhere") {
    ret_vec <- rep(FALSE, length(unique_gridcodes))
  } else if (scheme_key == "Wherever Meaningfully Better Than Uncentered") {
    if (is.null(threshold)) {
      threshold <- 0
      warning("no threshold provided to get_deseason_bool_vec_from_scheme()\n")
    }

    centered_indexes <- sapply(precomputed_score_df$model_identifier, function(mid) {                                                 ## nolint
      str_split(mid, "_")[[1]][2] == "centered"
    })
    uncentered_indexes <- sapply(precomputed_score_df$model_identifier, function(mid) {                                               ## nolint
      str_split(mid, "_")[[1]][2] == "uncentered"
    })

    stopifnot((sum(centered_indexes) + sum(uncentered_indexes)) == nrow(precomputed_score_df))                                        ## nolint

    centered_kges <- precomputed_score_df$inf_kge[centered_indexes]
    uncentered_kges <- precomputed_score_df$inf_kge[uncentered_indexes]
    score_diff <- centered_kges - uncentered_kges

    ret_vec <- score_diff > threshold

  } else if (scheme_key == "Wherever Not Meaningfully Worse Than Uncentered") {
    if (is.null(threshold)) {
      threshold <- 0
      warning("no threshold provided to get_deseason_bool_vec_from_scheme()\n")
    }

    centered_indexes <- sapply(precomputed_score_df$model_identifier, function(mid) {                                                 ## nolint
      str_split(mid, "_")[[1]][2] == "centered"
    })
    uncentered_indexes <- sapply(precomputed_score_df$model_identifier, function(mid) {                                               ## nolint
      str_split(mid, "_")[[1]][2] == "uncentered"
    })
    stopifnot(length(intersect(centered_indexes, uncentered_indexes)) == 0)

    centered_kges <- precomputed_score_df$inf_kge[centered_indexes]
    uncentered_kges <- precomputed_score_df$inf_kge[uncentered_indexes]
    score_diff <- uncentered_kges - centered_kges

    ret_vec <- !(score_diff > threshold)

  } else {
    stop(paste0("Invalid Deseasonalization Scheme (", scheme_key, ") Passed to get_deseason_bool_vec_from_scheme()"))                 ## nolint
  }

  if (debugmode) {
    cat("  > returning bool vec that uses", sum(ret_vec), "deseasonalized estimates out of", length(ret_vec), "\n")                   ## nolint
  }
  ret_vec
}

#' @param featwts is a vector of scalars by which each component of the overall
#' distance function should be weighted before combination
#' @param comb_func is one of "+", "*", "euc" and dictates how the weighted
#' feature-wise distance components are combined into an overall distance
#' @param n_pts_core is the value of k used for density estimation
#' @param n_pts_clust is the minimum number of points that must populate a
#' neighborhood in order for it to be labeled as a cluster
#' @param dszn_bool_vec is a sequence of booleans - one for each unique
#' gridcode in the dataset - that represents whether the centered estimate
#' should be used for that gridcode
#'
#' @return a named list of labels where the name is the gridcode and the label
#' is the cluster that the B(s, t) from the relevant gridcode belongs to, where
#' a label of 0 represents an outlier
get_hdbscan_labels <- function(featwts,
                               comb_func,
                               n_pts_core,
                               n_pts_clust,
                               dszn_bool_vec,
                               kge_thresh = 0) {
  print(head(dszn_bool_vec))

  dist_mat <- get_pairwise_distance_mtx(featwts = featwts,
                                        comb_func = comb_func,
                                        dszn_bool = dszn_bool_vec,
                                        kge_thresh = kge_thresh)

  ## check dist_mat for any issues
  stopifnot(colnames(dist_mat) == rownames(dist_mat))
  stopifnot(dim(dist_mat)[1] == dim(dist_mat)[2])

  hdbscan <- import("hdbscan")

  clusterer <- hdbscan$HDBSCAN(
    min_cluster_size = as.integer(n_pts_clust),
    min_samples = as.integer(n_pts_core),
    metric = "precomputed"
  )

  labels <- clusterer$fit_predict(dist_mat)
  labels <- labels + 1

  print(length(labels))
  print(dim(dist_mat))

  names(labels) <- colnames(dist_mat)

  if (debugmode) {
    pdf("./distmatboxplotALL.pdf")
    boxplot(dist_mat[upper.tri(dist_mat)])
    dev.off()

    nonoutlier_grids <- names(labels)[labels != 0]
    temp_mat <- dist_mat[nonoutlier_grids, nonoutlier_grids]
    pdf("./distmatboxplotNONOUTLIERS.pdf")
    boxplot(temp_mat[upper.tri(temp_mat)])
    dev.off()

    cat("clustering finished; num outliers:", sum(labels == 0), "\n")
  }

  list(labels = labels, dist_mat = dist_mat, dszn_bool_vec = dszn_bool_vec)
}


#' @param in_cluster_grids is a character vector of gridcodes that form a
#' cluster
#'
#' @return a flattened, lag-minor indexed B(s, t) surface which is the
#' point-wise median of the B(s, t) surfaces of the catchments in the cluster
calculate_median_beta <- function(in_cluster_model_ids, dszn_bool) {
  if (any(!is.character(in_cluster_model_ids))) {
    stop("Error in calculate_median_beta : in_cluster_grids must be a character vector")                                            ## nolint
  }

  apply(precomputed_beta_matrix[in_cluster_model_ids, , drop = FALSE], 2, function(column) {                  ## nolint
    median(column)
  })
}

#' @param in_cluster_grids is a character vector of gridcodes that form an
#' individual cluster
#' @param distmtx is a pairwise distance matrix of all clustered grids
#'
#' @return a flattened, lag-minor indexed B(s, t) surface which
#' is the medoid of the cluster provided
calculate_medoid_beta <- function(in_cluster_model_ids,
                                  distmtx, dszn_bool) {
  if (any(!is.character(in_cluster_model_ids))) {
    stop("Error in calculate_medoid_beta : in_cluster_grids must be a character vector")                              ## nolint
  }

  d <- distmtx[in_cluster_model_ids, in_cluster_model_ids, drop = FALSE]

  medoid_id <- rownames(d)[which.min(rowSums(d))]

  medoid_gc <- str_split(medoid_id, "_")[[1]][1]
  medoid_dszn_bool <- str_split(medoid_id, "_")[[1]][2] == "centered"

  to_ret <- get_beta(medoid_gc, medoid_dszn_bool)
  to_ret
}

#' @param n is the number of maximally distinct colors
#' to generate; must be at least 2
#'
#' @return a character vector of hex colors of length n
#' where the first color is grey (for outliers)
get_distinct_colors <- function(n) {
  floor_val <- 0.3
  top_val <- 1
  interp_pts <- seq(from = top_val, to = floor_val, length.out = ceiling(n / 6))
  color_palette_overbuilt <- as.vector(sapply(seq_len(ceiling(n / 6)), function(h) {                                  ## nolint
    c(rgb(interp_pts[h], 0, 0),
      rgb(0, interp_pts[h], 0),
      rgb(0, 0, interp_pts[h]),
      rgb(interp_pts[h], 0, interp_pts[h]),
      rgb(interp_pts[h], interp_pts[h], 0),
      rgb(0, interp_pts[h], interp_pts[h]))
  }))
  c("#6e6e6e", color_palette_overbuilt[seq_len(n - 1)])
}

#' @param scores is a set of numeric values to map to the Red -> Grey -> Green
#' color pallete based on their value
#' @param threshold_value is a single number separating "good" values (above
#' threshold) from "bad" values (below threshold)
#' @param zero_color_wb_frac is a number between 0 (black) and 1 (white) which
#' dictates how close to pure white the grey color representing the threshold
#' value should be
#'
#' @return a vector of hex color codes containing a color mapped from the
#' relevant value of
map_scores_to_rg_pal <- function(scores,
                                 threshold_value,
                                 zero_color_wb_frac = 0.92) {
  valid <- is.finite(scores)
  pos <- valid & scores >= threshold_value
  neg <- valid & scores < threshold_value

  scale_magnitude <- function(x) {
    if (!length(x)) return(x)

    smallest <- min(x)
    x <- x - smallest

    largest <- max(x)
    if (largest == 0) return(rep(0, length(x)))
    x / largest
  }

  pos_scaled <- scale_magnitude(scores[pos])
  neg_scaled <- scale_magnitude(-scores[neg])

  grey <- zero_color_wb_frac
  pos_fade <- grey * (1 - pos_scaled)
  neg_fade <- grey * (1 - neg_scaled)

  colors <- rep(NA_character_, length(scores))

  colors[pos] <- rgb(
    pos_fade,
    pos_fade + pos_scaled,
    pos_fade
  )

  colors[neg] <- rgb(
    neg_fade + neg_scaled,
    neg_fade,
    neg_fade
  )

  colors
}
















##############################################
## ARCHIVE
##############################################
#' @param gridcode is the gridcode of the catchment to load a B(s, t) plot of
#'
#' @return a levelplot() object of the B(s, t) derived from the catchment with
#' the provided gridcode
get_beta_plot_OLD <- function(gridcode, dszn_bool) {                                                                  ## nolint
  if (exists("precomputed_beta_plots")) {
    if (as.character(gridcode) %in% names(precomputed_beta_plots[[to_dszn_key(dszn_bool)]])) {                        ## nolint
      print("pulling plot from precompute...")
      all_plots <- precomputed_beta_plots[[to_dszn_key(dszn_bool)]][[as.character(gridcode)]]                         ## nolint

      ## compile these into one plot object and return it
      combined_plot <- all_plots$bplot + (all_plots$dplot / all_plots$ttplot / all_plots$smeplot) + plot_layout(widths = c(2, 1))                                   ## nolint

      final_ggplot <- combined_plot +
        plot_annotation(title = all_plots$overall_title_string)

      return(as.ggplot(final_ggplot))
    } else {
      warning("no precomputed plot was found for requested gridcode. building from scratch...")                       ## nolint
    }
  } else {
    warning("precomputed_beta_plots was not loaded. applet will experience significant latency.")                     ## nolint
  }

  plot_beta_and_derived_feats(get_beta(gridcode, dszn_bool))                                                                            ## nolint
}


#' @param featwts is a vector of scalars by which each component of the overall
#' distance function should be weighted before combination
#' @param comb_func is one of "+", "*", "euc" and dictates how the weighted
#' feature-wise distance components are combined into an overall distance
#' @param minpts is a scalar which dictates how many members must exist within
#' a given neighborhood in order to be considered a cluster
#'
#' @return a named list of labels where the name is the gridcode and the label
#' is the cluster that the B(s, t) from the relevant gridcode belongs to, where
#' a label of 0 represents an outlier
get_hdbscan_labels_OLD <- function(featwts, comb_func, minpts, dszn_bool) {   ## nolint
  dist_mat <- get_pairwise_distance_mtx(featwts = featwts,
                                        comb_func = comb_func,
                                        dszn_bool = dszn_bool)

  ## check dist_mat for any issues
  stopifnot(colnames(dist_mat) == rownames(dist_mat))
  stopifnot(dim(dist_mat)[1] == dim(dist_mat)[2])

  ## run hdbscan
  fit <- dbscan::hdbscan(as.dist(dist_mat), minPts = minpts)

  ## extract and name labels
  labels <- fit$cluster
  names(labels) <- colnames(dist_mat)

  if (debugmode) {
    pdf("./distmatboxplotALL.pdf")
    boxplot(dist_mat[upper.tri(dist_mat)])
    dev.off()

    nonoutlier_grids <- names(labels)[labels != 0]
    temp_mat <- dist_mat[nonoutlier_grids, nonoutlier_grids]
    pdf("./distmatboxplotNONOUTLIERS.pdf")
    boxplot(temp_mat[upper.tri(temp_mat)])
    dev.off()
  }

  list(labels = labels, dist_mat = dist_mat, dszn_bool = dszn_bool)
}