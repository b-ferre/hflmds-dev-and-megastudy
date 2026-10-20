setwd("./to deploy")
library(stringr)

extract_gc <- function(fn) {
  str_remove_all(fn, "[^\\d]")
}

data_avail_rain <- as.numeric(sapply(list.files("./data", pattern = "rain"), extract_gc))         ## nolint
data_avail_sflow <- as.numeric(sapply(list.files("./data", pattern = "streamflow"), extract_gc))  ## nolint

gc_index <- as.numeric(intersect(data_avail_rain, data_avail_sflow))
