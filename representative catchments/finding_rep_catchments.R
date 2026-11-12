# Select a nested pilot panel without fitting any hydrological models.
# Usage:
# Rscript select_catchments.R /path/to/hflm-ds_model /path/to/output
#
# Only base R and the recommended package 'cluster' are required. The
# applet's large RDS inputs require several GB of available memory.
# Inputs are read only. Outputs and an input-validated cache go into output.
# All ranking ties are resolved by ascending numeric gridcode.
root <- getwd()
out <- "./selected_catchments"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(requireNamespace("cluster", quietly = TRUE))
set.seed(20261111)

# Hash the data and source defining this selection. A changed input forces
# reconstruction of the feature cache, rather than silently reusing it.
files <- c(
  "deploy_to_cc/gc_index.rda",
  "data/raw/Metadata/merged_catchments_metadata.csv",
  "clustering/global.R", "clustering/precomputed_beta_mtx.rds",
  "clustering/precomputed_beta_plots.rds",
  "clustering/precomputed_feature_dist_mtxs.rds",
  "clustering/precomputed_score_df.rds",
  "clustering/precomputed_kernel_goodness_df.rds"
)
hashes <- tools::md5sum(file.path(root, files))
names(hashes) <- files
stopifnot(!anyNA(hashes))
result_files <- sort(list.files(file.path(root, "results")))
cache_path <- file.path(out, "selection_features_cache.rds")
cache_version <- 1L
cached <- if (file.exists(cache_path)) readRDS(cache_path) else NULL
reuse <- !is.null(cached) && identical(cached$hashes, hashes) &&
  identical(cached$version, cache_version) &&
  identical(cached$result_files, result_files)

if (!reuse) {
  message("Reading metadata, diagnostics, and deployment population...")
  index_env <- new.env()
  load(file.path(root, files[1]), envir = index_env)
  ids <- sort(unique(as.integer(index_env$gc_index)))
  meta <- read.csv(file.path(root, files[2]))
  stopifnot(!anyDuplicated(meta$gridcode))
  d <- meta[match(ids, meta$gridcode), , drop = FALSE]
  stopifnot(identical(as.integer(d$gridcode), ids))
  d$deployment_index <- match(ids, index_env$gc_index)
  read_app <- function(name) {
    readRDS(file.path(root, "clustering", paste0(name, ".rds")))
  }
  scores <- read_app("precomputed_score_df")
  quality <- read_app("precomputed_kernel_goodness_df")
  plots <- read_app("precomputed_beta_plots")
  beta <- read_app("precomputed_beta_mtx")
  stopifnot(!anyDuplicated(rownames(beta)))
  stopifnot(!anyDuplicated(scores$model_identifier))
  stopifnot(!anyDuplicated(quality$model_identifier))
  stopifnot(ncol(beta) == 150L * 365L)
  app_ids <- unique(as.integer(sub("_.*", "", rownames(beta))))
  stopifnot(all(app_ids %in% ids))
  d$in_applet <- ids %in% app_ids

  # Keep centered and uncentered diagnostics as separate dimensions. A
  # strong result in one regime must not hide a failure in the other.
  for (regime in c("centered", "uncentered")) {
    mids <- paste0(ids, "_", regime)
    for (nm in c("inf_kge", "pred_kge", "benchmark_kge")) {
      d[[paste0(nm, "_", regime)]] <-
        scores[[nm]][match(mids, scores$model_identifier)]
    }
    for (nm in setdiff(names(quality), "model_identifier")) {
      d[[paste0(nm, "_", regime)]] <-
        quality[[nm]][match(mids, quality$model_identifier)]
    }
    lag_features <- t(vapply(mids, function(mid) {
      p <- plots[[mid]]
      if (is.null(p)) return(rep(NA_real_, 6L))
      # The applet's delta is a support LENGTH: index 1 means lag 0.
      # Subtract one for lag in days; all-zero days remain zero.
      lag <- pmax(0, p$delta - 1)
      c(median(lag), max(lag), max(lag) - median(lag),
        mean(lag >= 100), mean(p$total_transfer),
        sd(p$total_transfer))
    }, numeric(6)))
    colnames(lag_features) <- paste0(c(
      "median_lag", "max_lag", "lag_excursion", "long_lag_fraction",
      "transfer_mean", "transfer_sd"
    ), "_", regime)
    d <- cbind(d, lag_features)
  }

  # Record available result combinations, including incomplete runs that
  # the applet preprocessing omitted. Absence does not establish a cause.
  for (mode in c("inf", "pred")) {
    for (regime in c("centered", "uncentered")) {
      nm <- paste("has", mode, regime, sep = "_")
      fn <- paste0(mode, "_model_", regime, "_grid", ids, ".rda")
      d[[nm]] <- file.exists(file.path(root, "results", fn))
    }
  }

  # Reuse the applet's exact functions and settings without sourcing its
  # Shiny UI, Python initialization, or unrelated large resources.
  # Evaluate only function definitions and the three settings assignments.
  app <- new.env(parent = globalenv())
  expressions <- parse(file.path(root, "clustering/global.R"))
  for (expr in expressions) {
    if (!is.call(expr) || !identical(expr[[1]], as.name("<-"))) next
    lhs <- as.character(expr[[2]])
    rhs <- expr[[3]]
    is_function <- is.call(rhs) &&
      identical(rhs[[1]], as.name("function"))
    if (is_function || lhs %in% c(
      "all_feat_names", "all_feat_weight_names", "INITCLUSTERSETTINGS"
    )) eval(expr, envir = app)
  }
  app$debugmode <- FALSE
  app$unique_gridcodes <- as.character(app_ids)
  app$precomputed_score_df <- scores
  app$precomputed_kernel_goodness_df <- quality
  app$precomputed_beta_matrix <- beta
  message("Reading feature distances and reproducing default clusters...")
  app$precomputed_feature_dist_mats <-
    read_app("precomputed_feature_dist_mtxs")
  settings <- app$INITCLUSTERSETTINGS
  weights <- unlist(settings[app$all_feat_weight_names])
  names(weights) <- app$all_feat_names
  reference <- app$get_clustering_result(
    clustering_algorithm = settings$clust_algo,
    distance_feature_weights = weights,
    distance_feature_comb_func = settings$comb_func,
    deseasonalization_scheme = settings$deseasonalization_selection,
    deseasonalization_cutoff = settings$dszn_thresh_selection,
    min_kge_for_inclusion = settings$min_kge_for_clust,
    min_goodness_for_inclusion = settings$min_goodness_for_clust,
    ignore_negative_beta_ests = settings$remove_negative_beta_ests,
    clustering_algorithm_params = settings,
    outlier_thresh = settings$pre_clust_outlier_trim,
    outlier_neighbor_index = settings$pre_clust_outlier_neighbor_index
  )
  reference_ids <- as.integer(sub("_.*", "", names(reference$labels)))
  d$default_cluster <- as.character(
    reference$labels[match(ids, reference_ids)]
  )
  d$default_cluster[is.na(d$default_cluster)] <- "excluded"
  default_mid <- reference$all_catchments_considered
  d$default_model <- default_mid[match(
    ids, as.integer(sub("_.*", "", default_mid))
  )]

  # Clusters are only one lens. Also build 12-cluster sensitivity views
  # separately in each regime, with NO score, sign, or outlier exclusion.
  # This uses the same distance definition, but changes the partition and
  # population. These labels are not claimed to be applet defaults.
  for (regime in c("centered", "uncentered")) {
    message("Clustering the unfiltered ", regime, " kernels...")
    mids <- paste0(sort(app_ids), "_", regime)
    dm <- app$get_pairwise_distance_mtx(
      mids, weights, settings$comb_func
    )
    stopifnot(all(is.finite(dm)))
    labels <- cluster::pam(dm, k = 12, diss = TRUE)$clustering
    d[[paste0("cluster12_", regime)]] <- as.character(
      labels[match(ids, sort(app_ids))]
    )
    rm(dm)
    gc(verbose = FALSE)
  }
  cached <- list(
    version = cache_version, hashes = hashes, data = d,
    settings = settings, weights = weights, result_files = result_files
  )
  saveRDS(cached, cache_path)
  rm(app, beta, plots, reference)
  gc(verbose = FALSE)
} else {
  message("Using feature cache validated against current input hashes.")
}
d <- cached$data

# Define transparent geographic proxies, not claims about watershed
# drainage direction or coastal boundaries. The northeast-European proxy
# is limited to the countries actually present in this deployment.
d$west_us <- with(d,
  country == "US" & longitude <= -120 & latitude >= 32 & latitude <= 50
)
d$west_gb <- with(d, country == "GB" & longitude <= -3)
d$ne_europe <- with(d,
  country %in% c("DE", "DK", "AT") & latitude >= 52 & longitude >= 8
)
d$west_us_spike <- with(d,
  west_us & max_lag_uncentered >= 100 & lag_excursion_uncentered >= 60
)
d$west_us_spike[is.na(d$west_us_spike)] <- FALSE
d$kge_difference <- d$inf_kge_centered - d$inf_kge_uncentered

# Bands are population terciles, separately for each regime. They mean
# relatively low, medium, and high, not universal hydrological cutoffs.
# Missing diagnostics are retained as a distinct band.
band <- function(x) {
  ans <- rep("missing", length(x))
  ok <- is.finite(x)
  breaks <- quantile(x[ok], c(1 / 3, 2 / 3), names = FALSE)
  ans[ok] <- c("low", "medium", "high")[
    1L + (x[ok] > breaks[1]) + (x[ok] > breaks[2])
  ]
  ans
}
for (nm in c("inf_kge", "goodness", "median_lag")) {
  for (regime in c("centered", "uncentered")) {
    key <- paste(nm, regime, sep = "_")
    d[[paste0(key, "_band")]] <- band(d[[key]])
  }
}

# Make a balanced distance across six conceptual blocks. Each numeric
# variable is converted to a population midrank in [0, 1]. This prevents
# units and extreme KGE magnitudes from dominating the geometry. Missing
# values use the middle rank, with separate missingness indicators.
# Within a block, variables have equal weight; blocks also have equal
# weight, regardless of how many correlated diagnostics they contain.
rank01 <- function(x) {
  ok <- is.finite(x)
  ans <- rep(0.5, length(x))
  if (sum(ok) > 1L) {
    ans[ok] <- (rank(x[ok], ties.method = "average") - 0.5) / sum(ok)
  }
  ans
}
numeric_distance <- function(columns) {
  x <- vapply(d[columns], rank01, numeric(nrow(d)))
  missing <- vapply(d[columns], function(x) !is.finite(x),
                    logical(nrow(d)))
  base <- as.matrix(dist(x, method = "manhattan")) / ncol(x)
  miss <- as.matrix(dist(missing, method = "manhattan")) / ncol(x)
  0.8 * base + 0.2 * miss
}
blocks <- list(
  attributes = c(
    "area_km2", "elevation_mean_m", "mean_slope_degree",
    "Median_DepthToBedrock_cm", "Aridity", "SF", "Porosity",
    "Seasonality_of_Moisture_Index", "low_high_ratio",
    "wet_days_ratio_1mm", "high_prec_freq", "low_prec_dur"
  ),
  performance = c("inf_kge_centered", "inf_kge_uncentered",
                  "kge_difference"),
  quality = grep(
    "^(negative|lag_spike|temporal|tail|patches)_comp_",
    names(d), value = TRUE
  ),
  lag_transfer = grep(
    "^(median_lag|max_lag|lag_excursion|transfer_mean|transfer_sd)_",
    names(d), value = TRUE
  )
)
blocks$lag_transfer <- blocks$lag_transfer[
  !grepl("_band$", blocks$lag_transfer)
]
message("Constructing balanced catchment distances...")
distance <- matrix(0, nrow(d), nrow(d))
for (columns in blocks) distance <- distance + numeric_distance(columns)

# Spherical coordinates avoid treating degrees of longitude as equal
# distances at every latitude. Chord distance is scaled into [0, 1].
lat <- d$latitude * pi / 180
lon <- d$longitude * pi / 180
xyz <- cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
stopifnot(all(is.finite(xyz)))
distance <- distance + as.matrix(dist(xyz)) / 2
cluster_keys <- c("default_cluster", "cluster12_centered",
                  "cluster12_uncentered")
for (key in cluster_keys) {
  x <- d[[key]]
  x[is.na(x)] <- "missing"
  distance <- distance + outer(x, x, "!=") / length(cluster_keys)
}
distance <- distance / 6
stopifnot(all(is.finite(distance)), max(abs(diag(distance))) == 0)

# Minimum coverage requirements deliberately overlap: a single catchment
# can cover a country, a cluster, a lag band, and a performance band.
# They are pilot-design safeguards, not proportional sampling weights.
membership <- list()
targets <- numeric()
add_requirement <- function(name, eligible, target) {
  eligible[is.na(eligible)] <- FALSE
  if (!any(eligible)) return(invisible(NULL))
  membership[[name]] <<- eligible
  targets[name] <<- min(target, sum(eligible))
}
for (country in sort(unique(d$country))) {
  add_requirement(paste0("country_", country), d$country == country, 1)
}
for (key in cluster_keys) {
  for (level in sort(unique(na.omit(d[[key]])))) {
    if (level == "excluded") next
    minimum <- if (key == "default_cluster") 3 else 2
    add_requirement(paste(key, level, sep = "_"),
                    d[[key]] == level, minimum)
  }
}
for (key in grep("_band$", names(d), value = TRUE)) {
  for (level in c("low", "medium", "high")) {
    add_requirement(paste(key, level, sep = "_"),
                    d[[key]] == level, 12)
  }
}
add_requirement("west_US", d$west_us, 10)
add_requirement("west_GB", d$west_gb, 10)
add_requirement("northeast_Europe", d$ne_europe, 8)
add_requirement("west_US_lag_spikes", d$west_us_spike, 5)
add_requirement("outside_applet", !d$in_applet, 10)

# Protect the contrasting lag behavior within the specifically requested
# European regions as well as across the whole population.
for (region in c("west_gb", "ne_europe")) {
  for (regime in c("centered", "uncentered")) {
    local_band <- band(ifelse(
      d[[region]], d[[paste0("median_lag_", regime)]], NA_real_
    ))
    for (level in c("low", "high")) {
      add_requirement(paste(region, regime, "lag", level, sep = "_"),
                      d[[region]] & local_band == level, 2)
    }
  }
}

# Also cover all three performance bands within each west-coast proxy.
# These are within-region terciles, so the UK and US comparisons are not
# forced to have identical absolute KGE distributions.
for (region in c("west_us", "west_gb")) {
  for (regime in c("centered", "uncentered")) {
    local_band <- band(ifelse(
      d[[region]], d[[paste0("inf_kge_", regime)]], NA_real_
    ))
    for (level in c("low", "medium", "high")) {
      add_requirement(paste(region, regime, level, sep = "_"),
                      d[[region]] & local_band == level, 2)
    }
  }
}
membership <- do.call(cbind, membership)
storage.mode(membership) <- "double"

# Anchor the user's known case plus a small set of diagnostic extremes.
# These cases test failure modes; they are not a random evaluation sample.
selected <- integer()
reasons <- character()
append_case <- function(i, reason) {
  if (i %in% selected) return(invisible(NULL))
  selected <<- c(selected, i)
  reasons <<- c(reasons, reason)
}
append_case(match(5017L, d$gridcode), "Known roughness case: 5017")
add_extreme <- function(key, eligible, number, label, high = TRUE) {
  eligible[is.na(eligible)] <- FALSE
  candidates <- which(eligible & is.finite(d[[key]]))
  direction <- if (high) -1 else 1
  ordered <- candidates[order(direction * d[[key]][candidates],
                              d$gridcode[candidates])]
  for (i in head(ordered, number)) append_case(i, label)
}
add_extreme("lag_excursion_uncentered", d$west_us_spike, 3,
            "West US uncentered lag-spike stress case")
for (regime in c("centered", "uncentered")) {
  for (component in c("temporal_comp", "tail_comp", "negative_comp")) {
    add_extreme(paste0(component, "_", regime), rep(TRUE, nrow(d)),
                1, paste("Extreme", regime, component))
  }
  add_extreme(paste0("goodness_", regime), rep(TRUE, nrow(d)), 1,
              paste("High-goodness", regime, "control"))
  add_extreme(paste0("median_lag_", regime), rep(TRUE, nrow(d)), 1,
              paste("Longest typical lag", regime, "stress case"))
}
stopifnot(!anyNA(selected))
anchor_count <- length(selected)

# Fill requirements first. Among candidates with equal requirement gain,
# prefer the largest reduction in population distance to the panel.
# Once requirements are met, combine mean-coverage improvement (70%) and
# distance from the existing panel (30%), each scaled by its current max.
# The latter prevents dense populations from consuming the entire panel.
nearest <- apply(distance[, selected, drop = FALSE], 1, min)
for (step in seq.int(length(selected) + 1L, 200L)) {
  counts <- colSums(membership[selected, , drop = FALSE])
  deficits <- pmax(0, targets - counts)
  quota_gain <- as.vector(membership %*% (deficits > 0))
  gain <- vapply(seq_len(nrow(d)), function(j) {
    mean(pmax(0, nearest - distance[, j]))
  }, numeric(1))
  scale_max <- function(x) if (max(x) > 0) x / max(x) else x
  coverage_score <- 0.7 * scale_max(gain) + 0.3 * scale_max(nearest)
  available <- setdiff(seq_len(nrow(d)), selected)
  order_i <- order(-quota_gain[available], -coverage_score[available],
                    d$gridcode[available])
  i <- available[order_i[1]]
  covered <- names(targets)[deficits > 0 & membership[i, ] > 0]
  reason <- if (length(covered)) {
    paste("Coverage:", paste(covered, collapse = "; "))
  } else {
    "Population coverage and diversity extension"
  }
  append_case(i, reason)
  nearest <- pmin(nearest, distance[, i])
}

# Export a nested sequence: the first 100 are contained in 125, 150, 200.
# Export both gridcodes and the ORIGINAL deployment indices. The existing
# cluster launcher expects an index into gc_index, not a gridcode.
sizes <- c(100L, 125L, 150L, 200L)
audit <- list()
coverage <- list()
for (size in sizes) {
  ix <- selected[seq_len(size)]
  panel <- d[ix, , drop = FALSE]
  panel$selection_order <- seq_len(size)
  panel$selection_reason <- reasons[seq_len(size)]
  panel$is_anchor <- seq_len(size) <= anchor_count
  write.csv(panel, file.path(out, paste0("catchments_", size, ".csv")),
            row.names = FALSE, na = "")
  writeLines(as.character(panel$gridcode),
             file.path(out, paste0("gridcodes_", size, ".txt")))
  writeLines(as.character(panel$deployment_index),
             file.path(out, paste0("deployment_indices_", size, ".txt")))
  counts <- colSums(membership[ix, , drop = FALSE])
  audit[[as.character(size)]] <- data.frame(
    size = size, requirement = names(targets), minimum = targets,
    population = colSums(membership), selected = counts,
    satisfied = counts >= targets
  )
  near <- apply(distance[, ix, drop = FALSE], 1, min)
  # Exclude selected cases from the coverage metric so that increasing N
  # cannot improve it just by adding more trivial zero self-distances.
  unsampled <- setdiff(seq_len(nrow(d)), ix)
  coverage[[as.character(size)]] <- data.frame(
    size = size, mean_distance = mean(near[unsampled]),
    p95_distance = unname(quantile(near[unsampled], 0.95)),
    max_distance = max(near[unsampled]),
    requirements_met = sum(counts >= targets),
    requirements_total = length(targets),
    outside_applet = sum(!d$in_applet[ix])
  )
}
write.csv(do.call(rbind, audit), file.path(out, "coverage_requirements.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, coverage), file.path(out, "size_comparison.csv"),
          row.names = FALSE)
write.csv(d, file.path(out, "population_features.csv"), row.names = FALSE,
          na = "")
write.csv(data.frame(file = names(hashes), md5 = unname(hashes)),
          file.path(out, "input_manifest.csv"), row.names = FALSE)
writeLines(result_files, file.path(out, "result_file_inventory.txt"))
saveRDS(list(settings = cached$settings, weights = cached$weights,
             feature_blocks = blocks, anchors = d$gridcode[
               selected[seq_len(anchor_count)]
             ], seed = 20261111), file.path(out, "selection_settings.rds"))
capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
stopifnot(length(unique(selected)) == 200L)
stopifnot(all(d$gridcode[selected] %in% cached$data$gridcode))

# Verify that both forcing and response input files exist for every panel
# member. File existence is not a claim about record quality or fit success.
for (kind in c("rain", "streamflow")) {
  paths <- file.path(root, "deploy_to_cc/data", paste0(
    "grid", d$gridcode[selected], "_", kind, ".rds"
  ))
  stopifnot(all(file.exists(paths)))
}

# Report distributions as well as coverage. The design intentionally
# enriches rare countries and failure cases, so an unweighted panel mean
# must not be presented as an unbiased population estimate.
summaries <- list()
metrics <- c("inf_kge_centered", "inf_kge_uncentered",
             "goodness_centered", "goodness_uncentered",
             "median_lag_centered", "median_lag_uncentered")
for (size in c(0L, sizes)) {
  ix <- if (size == 0L) seq_len(nrow(d)) else selected[seq_len(size)]
  for (key in metrics) {
    x <- d[[key]][ix]
    q <- quantile(x, c(0, 0.1, 0.5, 0.9, 1), na.rm = TRUE)
    summaries[[length(summaries) + 1L]] <- data.frame(
      size = size, metric = key, missing = sum(!is.finite(x)),
      minimum = q[1], p10 = q[2], median = q[3], p90 = q[4],
      maximum = q[5]
    )
  }
}
write.csv(do.call(rbind, summaries),
          file.path(out, "distribution_comparison.csv"), row.names = FALSE)

# A random-panel benchmark contextualizes the distance metric. It is not
# a significance test and is not used to choose a favorable random seed.
# This diagnostic is separate from the deterministic selection itself.
set.seed(20261111)
benchmark <- lapply(seq_len(50L), function(replicate) {
  ix <- sample.int(nrow(d), 100L)
  unsampled <- setdiff(seq_len(nrow(d)), ix)
  near <- apply(distance[unsampled, ix, drop = FALSE], 1, min)
  counts <- colSums(membership[ix, , drop = FALSE])
  data.frame(
    replicate = replicate, mean_distance = mean(near),
    p95_distance = unname(quantile(near, 0.95)),
    requirements_met = sum(counts >= targets)
  )
})
write.csv(do.call(rbind, benchmark),
          file.path(out, "random_100_benchmark.csv"), row.names = FALSE)

# Export a compact scientific overview using base graphics. Coordinates
# are plotted directly; the map does not imply detailed coast boundaries.
png(file.path(out, "coverage_overview.png"), width = 1800, height = 1200,
    res = 160)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
ix <- selected[seq_len(100L)]
plot(d$longitude, d$latitude, pch = 16, col = "#D5D9DD", cex = 0.6,
     xlab = "Longitude", ylab = "Latitude", main = "100-catchment panel")
points(d$longitude[ix], d$latitude[ix], pch = 21, bg = "#007C91",
       col = "white", cex = 0.9)
legend("bottomleft", c("Population", "Selected"), pch = 16,
       col = c("#D5D9DD", "#007C91"), bty = "n", cex = 0.8)
for (key in c("inf_kge_centered", "inf_kge_uncentered",
              "goodness_centered", "median_lag_uncentered")) {
  plot(ecdf(d[[key]]), main = gsub("_", " ", key),
       xlab = "Value", ylab = "Cumulative fraction", col = "#8B939B",
       verticals = TRUE, do.points = FALSE, lwd = 2)
  lines(ecdf(d[[key]][ix]), col = "#007C91", lwd = 2,
        verticals = TRUE, do.points = FALSE)
}
curve <- do.call(rbind, coverage)
plot(curve$size, curve$mean_distance, type = "b", pch = 16,
     col = "#007C91", ylim = range(c(curve$mean_distance,
                                      curve$p95_distance)),
     xlab = "Number of catchments", ylab = "Nearest-panel distance",
     main = "Coverage improves gradually")
lines(curve$size, curve$p95_distance, type = "b", pch = 17,
      col = "#BD5B39")
legend("topright", c("Mean", "95th percentile"), pch = c(16, 17),
       col = c("#007C91", "#BD5B39"), bty = "n", cex = 0.8)
invisible(dev.off())
print(do.call(rbind, coverage), row.names = FALSE)
message("Finished. Review coverage requirements before choosing panel size.")