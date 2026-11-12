library(leaflet)
library(stringr)
library(ggplot2)

map_data <- readRDS("./clustering/map_data.rds")
scores <- readRDS("./clustering/precomputed_score_df.rds")

map <- leaflet(map_data) |>
  addTiles() |>
  addCircleMarkers(
    lng = ~longitude,
    lat = ~latitude,
    popup = ~model_identifier
  )



scores$group <- sapply(scores$model_identifier, function(id_) {
  if (str_detect(id_, "_centered")) "centered" else "uncentered"
})
ggplot(data  = scores, aes(x = group, y = inf_kge, fill = group)) +
  geom_violin() +
  geom_boxplot(width = 0.2)


centered_ids <- map_data$model_identifier[
  sapply(map_data$model_identifier, function(id_) {
    str_detect(id_, "_centered")
  })
]
centered_scores <- scores[scores$model_identifier %in% centered_ids, ]


uncentered_ids <- map_data$model_identifier[
  sapply(map_data$model_identifier, function(id_) {
    str_detect(id_, "_uncentered")
  })
]
uncentered_scores <- scores[scores$model_identifier %in% uncentered_ids, ]