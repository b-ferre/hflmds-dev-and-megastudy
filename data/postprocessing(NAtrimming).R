library(stringr)
library(lattice)
library(ggplot2)
library(tidyr)

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
corr_df_all <- data.frame(
  grid = grids,
  tot_days = rep(NA, length(grids)),
  tot_na = rep(NA, length(grids)),
  frac_na_days = rep(NA, length(grids)),
  num_years_w_na = rep(NA, length(grids)),
  num_years_total = rep(NA, length(grids)),
  frac_years_w_na = rep(NA, length(grids))
)
p <- txtProgressBar(min = 0, max = length(grids), initial = 0, style = 3)
for (grid in grids) {

  sflow <- readRDS(paste0("./data/processed/grid", grid, "_streamflow.rds"))

  na_indexes <- which(!is.finite(sflow), arr.ind = TRUE)
  na_df <- data.frame(
    year = unique(na_indexes[, 1]),
    num_corr_days = sapply(unique(na_indexes[, 1]), function(yr_) {
      sum(na_indexes[, 1] == yr_)
    })
  )

  row_idx <- which(corr_df_all$grid == grid)

  corr_df_all[row_idx, "tot_days"] <- length(as.vector(sflow))
  corr_df_all[row_idx, "tot_na"] <- sum(!is.finite(as.vector(sflow)))
  corr_df_all[row_idx, "frac_na_days"] <- corr_df_all[row_idx, "tot_na"] / corr_df_all[row_idx, "tot_days"]                     ## nolint

  corr_df_all[row_idx, "num_years_total"] <- corr_df_all[row_idx, "tot_days"] / 365                                             ## nolint
  corr_df_all[row_idx, "num_years_w_na"] <- length(na_df$year)
  corr_df_all[row_idx, "frac_years_w_na"] <- corr_df_all[row_idx, "num_years_w_na"] / corr_df_all[row_idx, "num_years_total"]   ## nolint

  setTxtProgressBar(p, which(grids == grid)[1])
}
close(p)

## violin plots for all corr_df_all characteristics
if (FALSE) {
  corr_df_all |>
    pivot_longer(
      cols = -grid,
      names_to = "metric",
      values_to = "value"
    ) |>
    ggplot(aes(x = "", y = value)) +
    geom_violin(fill = "steelblue", alpha = 0.6) +
    facet_wrap(~metric, nrow = 1, scales = "free_y") +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

if (FALSE) {
  for (day_cutoff in c(0.01, 0.05, 0.1, 0.15, 0.2)) {
    num_disqualified_catchments <- sum(corr_df_all$frac_na_days > day_cutoff)

    year_thresholds <- c(0.05, 0.1, 0.25, 0.5)
    p <- txtProgressBar()
    num_bad_years_avg <- sapply(year_thresholds, function(yr_thresh) {
      mean(sapply(corr_df_all$grid[corr_df_all$frac_na_days > day_cutoff], function(grd_) {                      ## nolint

        sflow <- readRDS(paste0("./data/processed/grid", grd_, "_streamflow.rds"))

        na_indexes <- which(!is.finite(sflow), arr.ind = TRUE)
        na_df <- data.frame(
          year = unique(na_indexes[, 1]),
          num_corr_days = sapply(unique(na_indexes[, 1]), function(yr_) {
            sum(na_indexes[, 1] == yr_)
          })
        )
        na_df$frac_corr_days <- na_df$num_corr_days / 365

        sum(na_df$frac_corr_days > yr_thresh)
      }))
    })

    yr_cutoff_info_strs <- sapply(year_thresholds, function(thrsh_) {
      paste0(thrsh_, "              ",
            round(num_bad_years_avg[which(year_thresholds == thrsh_)], 3),
            "\n")
    })

    cat(paste0("With overall cutoff :", day_cutoff,
              " we disqualify ", num_disqualified_catchments, "\n"),
        paste0("YEAR CUTOFF     NUM BAD YEARS\n"),
        paste0(yr_cutoff_info_strs, collapse = ""))
  }
}

## we want catchments that have at least
check_ovr_cut <- 1

grids_to_consider <- corr_df_all$grid[corr_df_all$frac_na_days < check_ovr_cut]
cat("Of", length(corr_df_all$grid), "catchments,", length(grids_to_consider), "have less than", paste0(check_ovr_cut * 100, "%"), "of days missing.\n")                                             ## nolint

p <- txtProgressBar(min = 0, max = length(grids_to_consider), initial = 0, style = 3)         ## nolint
info_df <- data.frame(
  grid = grids_to_consider,
  num_salvageable_years_1p = rep(0, length(grids_to_consider)),
  num_salvageable_years_5p = rep(0, length(grids_to_consider)),
  num_salvageable_years_10p = rep(0, length(grids_to_consider)),
  num_salvageable_years_25p = rep(0, length(grids_to_consider)),
  num_salvageable_years_50p = rep(0, length(grids_to_consider))
)
for (grid in grids_to_consider) {
  sflow <- readRDS(paste0("./data/processed/grid", grid, "_streamflow.rds"))

  na_indexes <- which(!is.finite(sflow), arr.ind = TRUE)
  na_df <- data.frame(
    year = seq_len(nrow(sflow)),
    num_corr_days = sapply(seq_len(nrow(sflow)), function(yr_) {
      if (yr_ %in% unique(na_indexes[, 1])) {
        sum(na_indexes[, 1] == yr_)
      } else {
        0
      }
    })
  )
  na_df$frac_corr_days <- na_df$num_corr_days / 365

  row_ <- which(info_df$grid == grid)
  info_df$num_salvageable_years_1p[row_] <- sum(na_df$frac_corr_days < 0.01)
  info_df$num_salvageable_years_5p[row_] <- sum(na_df$frac_corr_days < 0.05)
  info_df$num_salvageable_years_10p[row_] <- sum(na_df$frac_corr_days < 0.1)
  info_df$num_salvageable_years_25p[row_] <- sum(na_df$frac_corr_days < 0.25)
  info_df$num_salvageable_years_50p[row_] <- sum(na_df$frac_corr_days < 0.50)

  setTxtProgressBar(p, which(grids == grid)[1])
}
close(p)

cat(
  "Of those remaining", length(grids_to_consider), "catchments there are on average...\n",     ## nolint
  "  >", mean(info_df$num_salvageable_years_1p), "years with less than 1% corruption\n",       ## nolint
  "  >", mean(info_df$num_salvageable_years_5p), "years with less than 5% corruption\n",       ## nolint
  "  >", mean(info_df$num_salvageable_years_10p), "years with less than 10% corruption\n",     ## nolint
  "  >", mean(info_df$num_salvageable_years_25p), "years with less than 25% corruption\n",     ## nolint
  "  >", mean(info_df$num_salvageable_years_50p), "years with less than 50% corruption\n"     ## nolint
)

num_years_minimum <- 32
year_salvageability_thresh <- 1
num_useable_catchments <- sum(info_df[, paste0("num_salvageable_years_", year_salvageability_thresh, "p")] >= num_years_minimum)    ## nolint

cat(
  "If we only accept catchments with...\n",
  "  - less than", paste0(check_ovr_cut * 100, "%"), "total corruption\n",
  "  - at least", num_years_minimum, "years each with less than", paste0(year_salvageability_thresh, "%"), "corruption\n",          ## nolint
  "Then we can use", num_useable_catchments, "catchments total.\n"
)
