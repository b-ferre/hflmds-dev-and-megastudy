library(leaflet)
library(stringr)
library(leafpop)

## helper for color-coding
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

precomputed_beta_matrix <- readRDS("./clustering/precomputed_beta_mtx.rds")
map_data <- readRDS("./clustering/map_data.rds")
unique_gridcodes <- unique(unname(sapply(rownames(precomputed_beta_matrix), function(mn) {   ## nolint
  str_split(mn, "_")[[1]][1]
})))

source("./zzz/plotters.R")

## Kernel quality diagnostics. Set TRUE to penalize abrupt terminal drops.
source("./clustering/kernel_quality.R")
include_tail_decay <- FALSE
i <- 0
p <- txtProgressBar(min = 0, max = nrow(precomputed_beta_matrix), style = 3)
qualities <- lapply(map_data$model_identifier, function(id_) {
  i <<- i + 1
  setTxtProgressBar(p, i)
  b <- precomputed_beta_matrix[id_, ]
  kernel_quality(b, include_tail = include_tail_decay)
})
close(p)

map_data$k_goodness <- 100 * exp(-map_data$k_ruff)

## Draw larger markers first so both models at a site can be clicked.
map_data$color <- map_scores_to_rg_pal(-1 * map_data$k_ruff, -1 * median(map_data$k_ruff), zero_color_wb_frac = 0.3)
map_data$radius <- ifelse(endsWith(map_data$model_identifier, "_centered"), 5, 11)
map_data <- map_data[order(map_data$radius, decreasing = TRUE), ]

## Render each trellis immediately so only one full plot is held at a time.
ids <- map_data$model_identifier
stopifnot(all(ids %in% rownames(precomputed_beta_matrix)))

recalculate_plots <- FALSE
if (recalculate_plots) {
  popups <- setNames(character(length(ids)), ids)
  p <- txtProgressBar(min = 0, max = length(ids), style = 3)
  for (i in seq_along(ids)) {
    id <- ids[i]
    plot <- plot_beta(precomputed_beta_matrix[id, ], title = paste0(id, "; roughness =", map_data$k_ruff[map_data$model_identifier == as.character(id)]))
    popups[[as.character(id)]] <- leafpop::popupGraph(list(plot), type = "png",
                                     width = 480, height = 320)[[1]]
    setTxtProgressBar(p, i)
  }
  close(p)
} else {
  popups <- readRDS("./popups.rds")
}

## Recover named entries from older caches that contain empty preallocated slots.
if (is.null(names(popups)) || !all(ids %in% names(popups))) {
  stop("Popup cache lacks model identifiers. Set recalculate_plots <- TRUE.")
}
popups <- popups[ids]
if (anyNA(popups) || any(!nzchar(popups))) {
  stop("Popup cache contains missing or empty plots. Set recalculate_plots <- TRUE.")
}

map <- leaflet(map_data) |>
  addTiles() |>
  addCircleMarkers(
    lng = ~longitude,
    lat = ~latitude,
    layerId = ~model_identifier,
    radius = ~radius,
    color = ~color,
    popup = unname(popups[map_data$model_identifier]),
    popupOptions = popupOptions(maxWidth = 500)
  )

saveRDS(popups, "./popups.rds")
map

cutoff <- 0.02
sum(map_data$k_ruff < cutoff) / nrow(map_data)
eligible <- which(map_data$k_ruff > cutoff)
closest <- eligible[which.min(map_data$k_ruff[eligible])]
to_view <- map_data$model_identifier[closest]
plot_beta(precomputed_beta_matrix[to_view, ])

map_data$color <- map_scores_to_rg_pal(-1 * map_data$k_ruff, -1 * cutoff, zero_color_wb_frac = 0.3)
map_data2 <- map_data[which(map_data$k_ruff < cutoff), ]
map2 <- leaflet(map_data2) |>
  addTiles() |>
  addCircleMarkers(
    lng = ~longitude,
    lat = ~latitude,
    layerId = ~model_identifier,
    radius = ~radius,
    color = ~color,
    popup = unname(popups[map_data2$model_identifier]),
    popupOptions = popupOptions(maxWidth = 500)
  )
map2

quantile(map_data$k_ruff, 0.67)

kernel_roughness_df <- data.frame(
  model_identifier = map_data$model_identifier,
  roughness = map_data$k_ruff
)

saveRDS(kernel_roughness_df, file = "kernel_roughness_df.rds")
