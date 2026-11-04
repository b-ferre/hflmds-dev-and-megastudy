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

### calculate the roughness or "squiggly-ness" of the kernels in a way
## that matches my qualitative sense
## TODO: implement
calculate_kernel_roughness <- function(kernel) {
  0
}


i <- 0
p <- txtProgressBar(min = 0, max = nrow(precomputed_beta_matrix), style = 3)
map_data$k_ruff <- sapply(map_data$model_identifier, function(id_) {
  i <<- i + 1
  setTxtProgressBar(p, i)
  b <- precomputed_beta_matrix[id_, ]
  calculate_kernel_roughness(b)
})
close(p)

## Draw larger markers first so both models at a site can be clicked.
map_data$color <- map_scores_to_rg_pal(rank(map_data$k_ruff), length(map_data$k_ruff) / 2)
map_data$radius <- ifelse(endsWith(map_data$model_identifier, "_centered"), 4, 7)
map_data <- map_data[order(map_data$radius, decreasing = TRUE), ]

## Render each trellis immediately so only one full plot is held at a time.
ids <- map_data$model_identifier
stopifnot(all(ids %in% rownames(precomputed_beta_matrix)))

recalculate_plots <- FALSE
if (recalculate_plots) {
  popups <- character(length(ids))
  p <- txtProgressBar(min = 0, max = length(ids), style = 3)
  for (i in seq_along(ids)) {
    id <- ids[i]
    plot <- plot_beta(precomputed_beta_matrix[id, ], title = id)
    popups[i] <- leafpop::popupGraph(list(plot), type = "png",
                                    width = 480, height = 320)[[1]]
    setTxtProgressBar(p, i)
  }
  close(p)
} else {
  popups <- readRDS("./popups.rds")
}


map <- leaflet(map_data) |>
  addTiles() |>
  addCircleMarkers(
    lng = ~longitude,
    lat = ~latitude,
    layerId = ~model_identifier,
    radius = ~radius,
    color = ~color,
    popup = popups,
    popupOptions = popupOptions(maxWidth = 500)
  )

saveRDS(popups, "./popups.rds")
map
