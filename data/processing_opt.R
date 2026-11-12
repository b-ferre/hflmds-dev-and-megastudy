library(stringr)
library(data.table)
read_func <- fread

# get file names for all available catchments
all_files <- list.files("./data/raw/Daily_data", pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)                                      ## nolint
m <- length(all_files)
cat("total # of catchments with data:", m, "\n")

# set up tracking df
goodness_df <- data.frame(
  fname_ = all_files,
  data_len_flg = rep(NA, m),
  closure_flg_1 = rep(NA, m),
  closure_flg_2 = rep(NA, m),
  snow_flag = rep(NA, m)
)

p <- txtProgressBar(min = 0, max = length(all_files), style = 3)
# in one pass, assess satisfaction of all criteria at once
for (i in seq_along(all_files)) {
  filepath <- all_files[i]
  goodness_df_index <- which(goodness_df$fname_ == filepath)

  # read file
  data <- read_func(filepath)

  # discard any catchment with < 30 sufficiently uncorrupted years of data
  # meaning years that are less than 1% corrupted, each
  stop("IMPLEMENT LINE 30")
  goodness_df$data_len_flg[goodness_df_index] <- TRUE

  # discard any catchment that is not approximately closed
  avg_p <- mean(data$precipitation_mmd, na.rm = TRUE)
  avg_pet <- mean(data$pet_mmd, na.rm = TRUE)
  avg_aet <- mean(data$aet_mm, na.rm = TRUE)
  avg_q <- mean(data$streamflow_mmd, na.rm = TRUE)

  # flag catchments that have net water loss long term
  # which suggests sginificant IGF out of catchment
  ## mass balance : P ~ Q + AET <= Q + PET
  ## technically in violation : P > Q + PET
  ## violation with "faff" pad :  (P)(0.95) > Q + PET
  ## plain english : (overest. of) volume out is less than 95% of volume in
  closure_flag_1 <- (avg_p * 0.95) > (avg_pet + avg_q)
  goodness_df$closure_flg_1[goodness_df_index] <- closure_flag_1

  # remove catchments that have net water gain long term
  # which suggests significant IGF into catchment
  ## mass balance : P ~ Q + AET
  ## technically in violation : P < Q + AET
  ## violation with "faff" pad :  P < (0.8)(Q + AET)
  ## plain english : volume in is less than 80% of volume out
  ## (larger faff pad this direction because AET estimates
  ## are highly uncertain)
  closure_flag_2 <- avg_p < (0.8 * (avg_q + avg_aet))
  goodness_df$closure_flg_2[goodness_df_index] <- closure_flag_2


  # discard any catchment that is < 95% rain dominated
  total_precip <- sum(data$precipitation_mmd, na.rm = TRUE)
  total_snow <- sum(data$precipitation_mmd[data$temperature_C < 0], na.rm = TRUE)                                                               ## nolint
  snow_flag <- (total_snow / total_precip) >= 0.05
  goodness_df$snow_flag[goodness_df_index] <- snow_flag

  setTxtProgressBar(p, i)
}
close(p)

# trim list of all files to include only those that passed checks
bad_data <- sum(goodness_df$data_len_flg)
bad_cls_1 <- sum(goodness_df$closure_flg_1)
bad_cls_2 <- sum(goodness_df$closure_flg_2)
bad_snow <- sum(goodness_df$snow_flag)
discard <- goodness_df$data_len_flg | goodness_df$closure_flg_1 | goodness_df$closure_flg_2 | goodness_df$snow_flag                                ## nolint
keep <- !discard
cat(
  "the following number of catchments did not pass the following checks:\n",
  "   > data length :", bad_data, "\n",
  "   > closure 1 (net loss) :", bad_cls_1, "\n",
  "   > closure 2 (net gain) :", bad_cls_2, "\n",
  "   > rain domination :", bad_snow, "\n",
  "in total,", sum(discard), "catchments were disqualified, leaving",
  sum(keep), "catchments to use. \n"
)

passed_files <- goodness_df$fname_[keep]

# extract grid names from valid files for saving purposes
gs <- str_remove(basename(passed_files), "\\.csv$")

# reformat, clean, and save data for all usable catchments
p <- txtProgressBar(min = 0, max = length(passed_files), style = 3)
for (j in seq_along(passed_files)) {
  filepath <- passed_files[j]
  grid_id <- gs[j]

  data <- read_func(filepath)

  # fast leap day removal
  leapdayidxs <- grepl("-02-29$", data$date)
  if (any(leapdayidxs)) {
    data <- data[!leapdayidxs, ]
  }

  # fast start date check
  if (!grepl("-01-01$", data$date[1])) {
    stop(paste("ERROR: First date is not Jan 1st for grid", grid_id))
  }

  # Vectorized matrix creation (replaces the nested for-loop)
  nyears <- as.integer(floor(nrow(data) / 365))
  valid_length <- nyears * 365

  # reshape x and y into a matrix with dimension years x days
  x <- matrix(data$precipitation_mmd[seq_len(valid_length)], nrow = nyears, ncol = 365, byrow = TRUE)                                           ## nolint
  y <- matrix(data$streamflow_mmd[seq_len(valid_length)], nrow = nyears, ncol = 365, byrow = TRUE)                                              ## nolint

  # delete any year that is >10% corrupted
  stop("IMPLEMENT LINE 120")

  # save files
  saveRDS(x, file = paste0("./data/processed/grid", grid_id, "_rain.rds"))
  saveRDS(y, file = paste0("./data/processed/grid", grid_id, "_streamflow.rds"))

  setTxtProgressBar(p, j)
}
close(p)
