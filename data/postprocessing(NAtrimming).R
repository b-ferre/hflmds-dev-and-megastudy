library(stringr)

all_data_files <- list.files("./data/processed")
rain_files <- all_data_files[str_detect(all_data_files, "rain")]
sflow_files <- all_data_files[str_detect(all_data_files, "streamflow")]

rain_grids <- as.numeric(sapply(rain_files, function(fn) str_remove_all(fn, "[^\\d]")))     ## nolint
sflow_grids <- as.numeric(sapply(sflow_files, function(fn) str_remove_all(fn, "[^\\d]")))   ## nolint

stopifnot((sum(rain_grids == sflow_grids) == length(rain_grids)) &&
            (length(rain_grids) == length(sflow_grids)))

grids <- rain_grids


## things I already know
#   - all rainfall data is finite
#   - every rain file is paired with sflow file

num_missing_data <- 0
to_delete <- c()
p <- txtProgressBar(min = 0, max = length(grids), initial = 0, style = 3)
for (grid in grids) {
  sflow <- readRDS(paste0("./data/processed/grid", grid, "_streamflow.rds"))

  if ((any(!is.finite(as.vector(sflow))))) {
    corrupted_day_idxs <- which(!is.finite(sflow), arr.ind = TRUE)
    corrupted_year_idxs <- unique(corrupted_day_idxs[, 1])
    checked <- c()
    corruption_df <- data.frame(
      yearid = corrupted_year_idxs,
      num_corr_days = unlist(sapply(corrupted_year_idxs, function(year_idx) {
        sum(corrupted_day_idxs[, 1] == year_idx)
      }))
    )
    num_unrecoverable_years <- sum(corruption_df$num_corr_days > 0)
    frac_unrecoverable_years <- num_unrecoverable_years / nrow(sflow)
    if (frac_unrecoverable_years > 0.30) {
      to_delete <- c(to_delete, grid)
    }
  }

  setTxtProgressBar(p, which(grids == grid)[1])
}

# for (grid in to_delete) {
#   file.remove(paste0("./data/processed/grid", grid, "_streamflow.rds"))
#   file.remove(paste0("./data/processed/grid", grid, "_rain.rds"))
# }