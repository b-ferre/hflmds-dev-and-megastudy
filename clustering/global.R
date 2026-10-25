library(shiny)
library(leaflet)
library(RColorBrewer)
library(grid)
library(lattice)
library(stringr)
library(leafpop)
library(ggplot2)
library(later)
library(future)
library(ggplotify)
library(patchwork)
library(bslib)
library(reticulate)
library(sf)
library(cluster)
library(dplyr)

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
  "tt_scaled_cheby",
  "global_l2"
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
      0.4, 0, 0.1, 0.2,
      0.4, 0.25, 0.2, 0.3,
      0, 0, 0, 0,
      0.1, 0.15, 0.05, 0,
      0
    )
  ), all_feat_weight_names),                                                                                                          ## nolint
  list(
    comb_func = "eCDF(+)",
    clust_algo = "K-Medoids",
    n_pts_core = 2,
    n_pts_clust = 10,
    linkage_method = "ward.D2",
    n_clusts_hclust = 8,
    n_clusts_kmedoids = 8,
    deseasonalization_selection = "Wherever Meaningfully Better Than Uncentered",                                                     ## nolint
    dszn_thresh_selection = 0,
    min_kge_for_clust = 0.4,
    remove_negative_beta_ests = TRUE
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
######################## RANDOM HELPER FUNCTIONS ##############################
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

###############################################################################
######################### PLOTTING HELPERS ####################################
###############################################################################
#' @param bvec is a flattened, lag-minor indexed vector of B(s, t) coefficients
#' @param title is a string to display above the levelplot as a plot title
#' @param ncolors is the number of unique non-white colors to use in the palette
#'
#' @return a levelplot() object of B(s, t) values with x = day and y = lag and
#' optionally, the vectors for smooth_delta and delta overlayed
plot_beta <- function(
  bvec, title = "", ncolors = 20, m = 364, max_lag = 150,
  delta = NULL, smooth_delta = NULL, simplified = FALSE
) {
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

  line_cols <- c("Raw delta" = "#E69F00", "Smoothed delta" = "#222222")
  line_types <- c(2, 1)

  # Strip ts/matrix attributes before plotting.
  if (!is.null(delta)) delta <- as.numeric(delta)
  if (!is.null(smooth_delta)) smooth_delta <- as.numeric(smooth_delta)

  xlab_ <- if (simplified) NULL else list("Day of the year", cex = 1.4)
  ylab_ <- if (simplified) NULL else list("Lag (Days)", cex = 1.4)
  scale_cex <- if (simplified) 0.6 else 1.2
  colorkey_width <- if (simplified) 0.6 else 1

  lattice::levelplot(
    b ~ t * s,
    data = bgrid,
    col.regions = cols,
    at = at_vals,
    xlab = xlab_,
    ylab = ylab_,
    scales = list(cex = scale_cex),
    colorkey = list(width = colorkey_width, axis.text = list(cex = scale_cex)),
    aspect = "xy",
    main = list(title, cex = scale_cex),

    panel = function(...) {
      # Draw the surface first, then overlay the lines.
      lattice::panel.levelplot(...)

      if (!is.null(delta)) {
        delta <- delta + 0.5
        lattice::panel.lines(
          x = seq_along(delta),
          y = as.numeric(delta),
          col = line_cols[1],
          lty = line_types[1],
          lwd = 2
        )
      }

      if (!is.null(smooth_delta)) {
        smooth_delta <- smooth_delta + 0.5
        lattice::panel.lines(
          x = seq_along(smooth_delta),
          y = as.numeric(smooth_delta),
          col = line_cols[2],
          lty = line_types[2],
          lwd = 2
        )
      }
    },

    key = if (!is.null(delta) || !is.null(smooth_delta)) {
      present <- c(!is.null(delta), !is.null(smooth_delta))

      list(
        space = "bottom",
        columns = sum(present),
        text = list(names(line_cols)[present]),
        lines = list(
          col = unname(line_cols[present]),
          lty = line_types[present],
          lwd = 2
        )
      )
    } else {
      NULL
    }
  )
}

render_beta_and_derived_feats <- function(bvec, plot_data) {
  p <- plot_data

  beta_plot <- plot_beta(
    bvec,
    title = "Beta Surface and Delta Vectors",
    ncolors = p$ncolors,
    m = p$m,
    max_lag = p$max_lag,
    delta = p$delta,
    smooth_delta = p$smooth_delta
  )

  make_line_plot <- function(values, title, ymax) {
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1

    ggplot(
      data.frame(
        day = seq_along(values),
        value = as.numeric(values)
      ),
      aes(x = day, y = value)
    ) +
      geom_line() +
      coord_cartesian(ylim = c(0, ymax)) +
      labs(title = title, x = "Day", y = NULL) +
      theme(plot.title = element_text(hjust = 0.5))
  }

  tt_plot <- make_line_plot(
    p$total_transfer, "Total Transfer Vector", p$max_tt
  )

  sme_plot <- make_line_plot(
    p$max_effective_lag, "Maximally Effective Lag", p$max_sme
  )

  combined <- gridExtra::arrangeGrob(
    beta_plot,
    tt_plot,
    sme_plot,
    layout_matrix = rbind(
      c(1, 2),
      c(1, 3)
    ),
    widths = c(2, 1),
    top = p$title
  )

  combined
}

get_rendered_beta_plot <- function(model_identifier) {
  bvec <- precomputed_beta_matrix[as.character(model_identifier), ]

  plot_data <- NULL

  if (exists("precomputed_beta_plots")) {
    plot_data <- precomputed_beta_plots[[as.character(model_identifier)]]
  }

  if (is.null(plot_data)) {
    plot_data <- calculate_plot_info(bvec)  ## nolint
  }

  obj <- render_beta_and_derived_feats(bvec, plot_data)
  obj
}

get_simple_rendered_beta_plot <- function(model_identifier) {
  bvec <- precomputed_beta_matrix[as.character(model_identifier), ]

  obj <- plot_beta(bvec, title = model_identifier, simplified = TRUE)
  obj
}

###############################################################################
###################### MAIN CLUSTERING ENTRY-POINT ############################
###############################################################################
get_clustering_result <- function(
  clustering_algorithm,
  distance_feature_weights,
  distance_feature_comb_func,
  deseasonalization_scheme,
  deseasonalization_cutoff,
  min_kge_for_inclusion,
  ignore_negative_beta_ests,
  clustering_algorithm_params,
  outlier_thresh,
  outlier_neighbor_index
) {
  ## first, calculate deseasonalization boolean for each unique gridcode
  if (debugmode) {
    cat("  > calculating dszn_bool_vec using scheme:",
        deseasonalization_scheme,
        "and cutoff:", deseasonalization_cutoff, "\n")
  }
  dszn_bool_vec <- get_deseason_bool_vec_from_scheme(
    scheme_key = deseasonalization_scheme,
    threshold = deseasonalization_cutoff
  )
  if (debugmode) {
    cat("  > done.", sum(dszn_bool_vec), "of",
        length(dszn_bool_vec), "models are deseasonalized.\n")
  }

  ## create a list of model identifiers using bool vec above such that there
  ## is one entry for each unique gridcode we have which is either of the form
  ## gc_centered or gc_uncentered depending on whether or not deseasonalized
  ## estimate should be used for that gridcode
  models_to_consider <- as.vector(sapply(seq_along(dszn_bool_vec), function(idx_) { ## nolint
    gc_i <- names(dszn_bool_vec)[idx_]
    dszn_tag_i <- if (dszn_bool_vec[idx_]) "_centered" else "_uncentered"
    paste0(gc_i, dszn_tag_i)
  }))

  ## trim our list of model identifiers using min_kge_for_inclusion such that
  ## only models with an inferential kge above this cutoff are included in
  ## clustering
  to_keep <- as.vector(sapply(models_to_consider, function(model_id_) {
    precomputed_score_df$inf_kge[
      precomputed_score_df$model_identifier == model_id_
    ] >= min_kge_for_inclusion
  }))
  models_to_use <- models_to_consider[to_keep]
  if (debugmode) {
    cat("  > trimming low KGE models. removed",
        length(to_keep) - sum(to_keep), "of", length(to_keep),
        "models, leaving", sum(to_keep), "\n")
  }

  ## further trim our list of models to use such that we ignore any model that
  ## has more than 5% of parameters negative or any large (>10% the magnitude
  ## of the largest positive parameter) negative parameter
  if (ignore_negative_beta_ests) {
    to_keep <- sapply(models_to_use, function(model_id_) {
      beta_i_ <- precomputed_beta_matrix[model_id_, ]
      neg_frac <- sum(beta_i_ < 0) / sum(beta_i_ != 0)
      largest_neg_rel <- max(c(
        max(c(abs(beta_i_[beta_i_ < 0]), 0)) / max(c(abs(beta_i_[beta_i_ > 0]), 0)),  ## nolint
        0
      ), na.rm = TRUE)

      (neg_frac < 0.05) && (largest_neg_rel < 0.1)
    })
    models_to_use <- models_to_use[to_keep]
    if (debugmode) {
      cat("  > trimming \'mis-estimated\' models. removed",
          length(to_keep) - sum(to_keep), "of", length(to_keep),
          "models, leaving", sum(to_keep), "\n")
    }
  }

  ## feed our list of model identifiers along with our distance function params
  ## (i.e. feature weights and combination function) to get_pairwise_dist_mat
  if (debugmode) cat("  > calculating distance matrix...\n")
  dist_mat <- get_pairwise_distance_mtx(
    model_ids = models_to_use,
    featwts = distance_feature_weights,
    comb_func = distance_feature_comb_func
  )
  if (debugmode) cat("  > done.\n")

  ## if indicated, pre-trim any outliers with k-neighbor distances
  ## in the upper outlier_thresh quantile before feeding to clusterer
  if (outlier_thresh > 0) {
    if (debugmode) {
      cat(
        "  > trimming catchments with",
        paste0(outlier_neighbor_index, "-neighbor"), "distance",
        "in the top", paste0((outlier_thresh) * 100, "% of all"),
        paste0(outlier_neighbor_index, "-neighbor"), "distances...\n"
      )
    }

    neighbor_dists <- apply(X = dist_mat, MARGIN = 1, FUN = function(d_i) {
      sort(d_i, partial = outlier_neighbor_index)[outlier_neighbor_index]
    })

    cutoff <- quantile(neighbor_dists, 1 - outlier_thresh)
    keep <- neighbor_dists < cutoff

    n_keep <- sum(keep)
    if (debugmode) {
      cat("  > done. removed", length(keep) - n_keep, "catchments\n")
    }

    dist_mat <- dist_mat[keep, keep, drop = FALSE]
  }

  ## feed our pairwise distance matrix into the clustering algorithm indicated
  ## by clustering_algorithm parameter and then get a list of labels named by
  ## the corresponding model identifier
  labels <- switch(
    clustering_algorithm,
    "HDBSCAN" = {
      if (debugmode) cat("  > running hdbscan clustering...\n")
      get_hdbscan_labels(
        dist_mat = dist_mat,
        n_pts_core = clustering_algorithm_params$n_pts_core,
        n_pts_clust = clustering_algorithm_params$n_pts_clust
      )
    },
    "Hierarchical" = {
      if (debugmode) cat("  > running hierarchical clustering...\n")
      get_hclust_labels(
        dist_mat = dist_mat,
        linkage_method = clustering_algorithm_params$linkage_method,
        n_clusts = clustering_algorithm_params$n_clusts_hclust
      )
    },
    "K-Medoids" = {
      if (debugmode) cat("  > running k-medoids clustering...\n")
      get_kmedoids_labels(
        dist_mat = dist_mat,
        n_clusts = clustering_algorithm_params$n_clusts_kmedoids
      )
    }
  )
  if (debugmode) {
    cluster_sizes <- as.vector(
      sapply(sort(unique(labels[labels != 0])), function(lab_) {
        paste0(lab_, ":", sum(labels == lab_))
      })
    )
    cat("  > clustering finished. \n     num clusters:",
        length(unique(labels[labels != 0])),
        "\n     num outliers:", sum(labels == 0),
        "\n     cluster sizes: (", cluster_sizes, ")\n")
  }

  ## generate a color key based on whether or not our selected clustering
  ## algorithm returns labels of 0 to represent outliers or is incapable of
  ## automatic outlier detection
  num_non_outlier_clusts <- sum(unique(labels) != 0)
  cluster_colors <- get_distinct_colors(num_non_outlier_clusts)
  if (sum(labels == 0) > 0) {
    cluster_colors <- c("#777777", cluster_colors)
  }
  names(cluster_colors) <- as.character(sort(unique(labels)))

  ## return everything
  list(
    labels = labels,
    dist_mat = dist_mat,
    clust_colors = cluster_colors,
    all_catchments_considered = models_to_consider
  )
}

###############################################################################
####################### CLUSTERING HELPERS ####################################
###############################################################################\

#' @param scheme_key is a string explaining when/where to use deseasonalized
#' model results as opposed to "uncentered" model results. should be one of:
#' "Everywhere", "Nowhere", "Wherever Meaningfully Better than Uncentered",
#' or "Wherever not Meaningfully Worse than Uncentered"
#' @param threshold is a value used for determining what a "meaningful" gap in
#' performance (inferential KGE) between a centered and uncentered model is
get_deseason_bool_vec_from_scheme <- function(scheme_key, threshold) {                                                         ## nolint
  if (scheme_key == "Everywhere") {
    return(unlist(setNames(rep(TRUE, length(unique_gridcodes)),
                           as.character(unique_gridcodes))))
  }
  if (scheme_key == "Nowhere") {
    return(unlist(setNames(rep(FALSE, length(unique_gridcodes)),
                           as.character(unique_gridcodes))))
  }
  if (scheme_key == "Wherever Meaningfully Better Than Uncentered") {
    return(unlist(setNames(lapply(unique_gridcodes, function(gc_i) {
      uncent_kge_i <- precomputed_score_df$inf_kge[
        precomputed_score_df$model_identifier == paste0(gc_i, "_uncentered")
      ]
      cent_kge_i <- precomputed_score_df$inf_kge[
        precomputed_score_df$model_identifier == paste0(gc_i, "_centered")
      ]
      cent_kge_i - uncent_kge_i > threshold
    }), as.character(unique_gridcodes))))
  }
  if (scheme_key == "Wherever Not Meaningfully Worse Than Uncentered") {
    return(unlist(setNames(lapply(unique_gridcodes, function(gc_i) {
      uncent_kge_i <- precomputed_score_df$inf_kge[
        precomputed_score_df$model_identifier == paste0(gc_i, "_uncentered")
      ]
      cent_kge_i <- precomputed_score_df$inf_kge[
        precomputed_score_df$model_identifier == paste0(gc_i, "_centered")
      ]
      uncent_kge_i - cent_kge_i < threshold
    }), as.character(unique_gridcodes))))
  }
  stop(paste0("Invalid Deseasonalization Scheme (", scheme_key, ") Passed to get_deseason_bool_vec_from_scheme()"))                 ## nolint
}

#' @param featwts is a vector of scalars by which each component of the overall
#' distance function should be weighted before combination
#' @param comb_func is one of "+", "*", "euc" and dictates how the weighted
#' feature-wise distance components are combined into an overall distance
#' @param model_ids is a list of model identifiers (names of the cols/rows of
#' the precomputed distance matrices) to use
#'
#' @return a matrix with overall pairwise distances
get_pairwise_distance_mtx <- function(model_ids, featwts, comb_func) {
  nfeats <- length(precomputed_feature_dist_mats)
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

  weighted_dist_components <- lapply(names(precomputed_feature_dist_mats),                                                            ## nolint
                                     function(name_j) {
                                       feat_mat_j <- precomputed_feature_dist_mats[[name_j]][model_ids, model_ids, drop = FALSE]  ## nolint
                                       featwts[[name_j]] * feat_mat_j
                                     })

  # helper function which converts every entry in a matrix to its
  # empirical CDF value across all matrix values
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

  if (comb_func == "+") {
    plus_mat <- Reduce("+", weighted_dist_components)
    return(plus_mat) ## nolint
  } else if (comb_func == "eCDF(+)") {
    plus_mat <- Reduce("+", weighted_dist_components)
    new_mat <- rank_matrix_cdf(plus_mat)
    dimnames(new_mat) <- dimnames(plus_mat)
    return(new_mat) ## nolint
  } else if (comb_func == "eCDF(+(eCDF))") {
    weighted_dist_components <- lapply(seq_along(weighted_dist_components), function(i) {   ## nolint
      m_i <- weighted_dist_components[[i]]

      if (any(!is.finite(m_i))) {
        stop(paste0("ERROR: non-finite values in component", names(weighted_dist_components)[i]))   ## nolint
      }
      transf_m_i <- rank_matrix_cdf(m_i)
      dimnames(transf_m_i) <- dimnames(m_i)
      as.matrix(transf_m_i)
    })
    plus_mat <- Reduce("+", weighted_dist_components)
    new_mat <- rank_matrix_cdf(plus_mat)
    dimnames(new_mat) <- dimnames(plus_mat)
    return(new_mat) ## nolint"
  } else if (comb_func == "+(eCDF)") {
    weighted_dist_components <- lapply(seq_along(weighted_dist_components), function(i) {   ## nolint
      m_i <- weighted_dist_components[[i]]

      if (any(!is.finite(m_i))) {
        stop(paste0("ERROR: non-finite values in component", names(weighted_dist_components)[i]))   ## nolint
      }
      transf_m_i <- rank_matrix_cdf(m_i)
      dimnames(transf_m_i) <- dimnames(m_i)
      as.matrix(transf_m_i)
    })
    new_mat <- Reduce("+", weighted_dist_components)
    return(new_mat) ## nolint"
  } else if (comb_func == "(0.8)eCDF(+) + (0.2)(+)") {
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

#' @param dist_mat is a pairwise distance matrix with colum/row names based
#' on model identifiers
#' @param n_pts_core is the value of k used for density estimation
#' @param n_pts_clust is the minimum number of points that must populate a
#' neighborhood in order for it to be labeled as a cluster
#'
#' @return a named list of labels where the name is the gridcode and the label
#' is the cluster that the B(s, t) from the relevant gridcode belongs to, where
#' a label of 0 represents an outlier
get_hdbscan_labels <- function(dist_mat,
                               n_pts_core,
                               n_pts_clust) {
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
  labels <- labels + 1 ## adjust for the fact that python uses -1 for outliers
  names(labels) <- colnames(dist_mat)


  labels
}

#' @param dist_mat is a pairwise distance matrix with colum/row names based
#' on model identifiers
#' @param linkage_method is a linkage method to be passed to stats::hclust
#' @param n_clusts is the number of clusters to find
#'
#' @return a named list of labels where the name is the gridcode and the label
#' is the cluster that the B(s, t) from the relevant gridcode belongs to
get_hclust_labels <- function(dist_mat,
                              linkage_method,
                              n_clusts) {
  ## check dist_mat for any issues
  stopifnot(colnames(dist_mat) == rownames(dist_mat))
  stopifnot(dim(dist_mat)[1] == dim(dist_mat)[2])

  ## check linkage_method for any issues
  stopifnot(linkage_method %in% c(
    "complete",
    "ward.D",
    "ward.D2",
    "average",
    "single",
    "mcquitty",
    "median",
    "centroid"
  ))

  clust_tree <- hclust(d = as.dist(dist_mat), method = linkage_method)
  labels <- cutree(tree = clust_tree, k = n_clusts)
  names(labels) <- colnames(dist_mat)

  labels
}


#' @param dist_mat is a pairwise distance matrix with colum/row names based
#' on model identifiers
#' @param n_clusts is the number of clusters to find
#'
#' @return a named list of labels where the name is the gridcode and the label
#' is the cluster that the B(s, t) from the relevant gridcode belongs to
get_kmedoids_labels <- function(dist_mat,
                                n_clusts) {
  ## check dist_mat for any issues
  stopifnot(colnames(dist_mat) == rownames(dist_mat))
  stopifnot(dim(dist_mat)[1] == dim(dist_mat)[2])

  clust_res <- cluster::pam(x = dist_mat, k = n_clusts, diss = TRUE)
  labels <- clust_res$clustering

  labels
}
###############################################################################
####################### REPRESENTATIVE CALCULATIONS ###########################
###############################################################################

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

###############################################################################
######################### COLOR-CODING HELPERS ################################
###############################################################################

#' @param n is the number of maximally distinct colors
#' to generate; must be at least 2
#'
#' @return a character vector of hex colors of length n
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
  color_palette_overbuilt[seq_len(n)]
}

#' @param scores is a set of numeric values to map to the Red -> Grey -> Green
#' color pallete based on their value (increasing)
#' @param threshold_value is a single number separating "good" values (above
#' threshold) from "bad" values (below threshold)
#' @param zero_color_wb_frac is a number between 0 (black) and 1 (white) which
#' dictates how close to pure white the grey color representing the threshold
#' value should be
#'
#' @return a vector of hex color codes in the same order as scores where each
#' entry is the color mapping of the corresponding value in scores
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

#' TODO: document this
#' heatmap grid and preset radius in km
generate_heatmap_layer <- function(catchments_sf,
                                   heatmap_grid_size_km,
                                   preset_radius_km = 500) {

  if (debugmode) cat("HEATMAP: heatmap_grid_size_km:", heatmap_grid_size_km, "\n")                  ## nolint
  if (debugmode) cat("HEATMAP: heatmap_trim_radius_km:", preset_radius_km, "\n")                    ## nolint

  heatmap_grid_size <- heatmap_grid_size_km * 1000
  preset_radius <- preset_radius_km * 1000

  MAX_OPACITY <- 0.7  ## corresponds to areas that are 100% dominated by one cluster                ## nolint
  MIN_OPACITY <- 0.3  ## corresponds to areas that are equally populated by two or more clusters    ## nolint

  ## project coords to mercator for calculations
  catchments_sf <- st_transform(catchments_sf, 3857)

  ## mask out areas where there are insufficient number of catchments
  boundary_mask <- st_union(st_buffer(catchments_sf, dist = preset_radius))

  if (heatmap_grid_size == 0) {
    # when heatmap_grid_size is set to 0, it should just be a clipped
    # vornoi diagram centered around our catchments
    points_merged <- do.call(c, st_geometry(catchments_sf))
    voronoi_polys <- st_collection_extract(st_voronoi(points_merged))

    poly_data <- st_sf(geometry = voronoi_polys) |>
      st_join(catchments_sf) |>
      mutate(
        majority_cluster = cluster,
        opacity = MAX_OPACITY
      )
  } else {
    # generate the grid of cells to look at using hex grid (looks better)
    grid_polys <- st_make_grid(
      catchments_sf,
      cellsize = heatmap_grid_size,
      square = FALSE
    )
    grid_sf <- st_sf(poly_id = seq_along(grid_polys), geometry = grid_polys)

    # attach spatial data to grid polygons
    grid_pts <- st_join(catchments_sf, grid_sf)

    # for each grid point, calculate which cluster dominates and by how much
    poly_stats <- grid_pts |>
      st_drop_geometry() |>
      group_by(poly_id, cluster) |>                                                 ## nolint
      summarise(n = n(), .groups = "drop") |>
      group_by(poly_id) |>
      arrange(desc(n)) |>
      summarise(
        total_catchments = sum(n),
        majority_cluster = first(cluster),
        count_1st = first(n),
        count_2nd = nth(n, 2, default = 0)
      ) |>
      mutate(
        margin_ratio = (count_1st - count_2nd) / total_catchments,                  ## nolint
        opacity = MIN_OPACITY + (margin_ratio * (MAX_OPACITY - MIN_OPACITY))        ## nolint
      )

    poly_data <- grid_sf |>
      left_join(poly_stats, by = "poly_id") |>
      filter(!is.na(majority_cluster))                                              ## nolint
  }

  # clip polygons to boundaries
  final_heatmap <- st_intersection(poly_data, boundary_mask)

  ## convert final object back to WSG84 for leaflet rendering
  final_heatmap <- st_transform(final_heatmap, 4326)

  final_heatmap
}