library(stringr)

gs <- c()
for (i in seq_len(9)) {
  gs <- c(gs, sapply(list.files(paste0("./data/raw/Daily_data/cluster_", i, "/")), function(f) {str_remove(f, ".csv")}))      ## nolint
}
cat("total # of catchments with data : ", length(gs), "\n")

###############################################################################
### remove catchments that (likely) violate mass balance ######################
if (TRUE) {
  p <- txtProgressBar(min = 0, max = length(gs), style = 3)
  bad_c <- gs[sapply(gs, function(grid) {
    ## find the cluster that contains this file
    cluster_id <- which(sapply(seq_len(9), function(f) {
      paste0(grid, ".csv") %in% list.files(paste0("./data/raw/Daily_data/cluster_", f, "/"))          ## nolint
    }))[1]
    data <- read.csv(paste0("./data/raw/Daily_data/cluster_", cluster_id, "/", grid, ".csv"))         ## nolint
    average_precip_mmd <- mean(data$precipitation_mmd, na.rm = TRUE)
    average_pet_mmd <- mean(data$pet_mmd, na.rm = TRUE)
    average_aet_mmd <- mean(data$aet_mm, na.rm = TRUE)
    average_sflow_mmd <- mean(data$streamflow_mmd, na.rm = TRUE)

    ### flag catchments that have net water loss long term
    ### which suggests sginificant IGF out of catchment
    ## mass balance : P ~ Q + AET <= Q + PET
    ## technically in violation : P > Q + PET
    ## violation with "faff" pad :  (P)(0.95) > Q + PET
    ## plain english : (overest. of) volume out is less than 95% of volume in
    flag_1 <- (average_precip_mmd) * (0.95) > (average_pet_mmd + average_sflow_mmd)                   ## nolint

    ### remove catchments that have net water gain long term
    ### which suggests significant IGF into catchment
    ## mass balance : P ~ Q + AET
    ## technically in violation : P < Q + AET
    ## violation with "faff" pad :  P < (0.8)(Q + AET)
    ## plain english : volume in is less than 80% of volume out
    ## (larger faff pad this direction because AET estimates
    ## are highly uncertain)
    flag_2 <- average_precip_mmd < (0.8 * (average_sflow_mmd + average_aet_mmd))

    setTxtProgressBar(p, which(gs == grid)[1])

    flag_1 || flag_2
  })]
  gs <- gs[!(gs %in% bad_c)]
  cat("# of catchments remaining after mass balance filtering : ", length(gs), "\n")             ## nolint
}
###############################################################################

###############################################################################
## remove catchments where <95% of precipitation volume is rainfall ###########
if (TRUE) {
  #' helper function to filter out snow-dominated catchments
  #'
  #' @param gridcodes are the gridcodes of catchments to filter in numeric form
  #' @param tol is the fraction of total precipitation volume that
  #' must consist of snowfall for a catchment to be identified
  #' as snow dominated
  #'
  #' @return a numeric vector containing the gridcodes of only the rain
  #' dominated catchments
  remove_snow_dominated <- function(gridcodes, tol) {
    p <- txtProgressBar(min = 1, max = length(gridcodes), initial = 0, style = 3)                       ## nolint
    gridcodes[sapply(gridcodes, function(grid) {
      ## find the cluster that contains this file
      cluster_id <- which(sapply(seq_len(9), function(f) {
        paste0(grid, ".csv") %in% list.files(paste0("./data/raw/Daily_data/cluster_", f, "/"))          ## nolint
      }))[1]
      data <- read.csv(paste0("./data/raw/Daily_data/cluster_", cluster_id, "/", grid, ".csv"))         ## nolint
      total_precip <- sum(data$precipitation_mmd, na.rm = TRUE)
      total_snow <- sum(data$precipitation_mmd[data$temperature_C < 0], na.rm = TRUE)                   ## nolint
      ## update progress bar
      setTxtProgressBar(p, which(gridcodes == grid)[1])
      ## keep if total snow fraction is less than tol
      (total_snow / total_precip) < tol
    })]
  }

  gs <- remove_snow_dominated(gs, tol = 0.05)
  cat("# of catchments remaining after snow domination filtering : ", length(gs), "\n")             ## nolint
}
###############################################################################

###############################################################################
## remove catchments that have less than 20 years of useable data #############
if (TRUE) {
  p <- txtProgressBar(min = 0, max = length(gs), initial = 0, style = 3)
  gs <- gs[sapply(gs, function(g) {
    ## find the cluster that contains this file
    setTxtProgressBar(p, which(gs == g)[1])
    cluster_id <- which(sapply(seq_len(9), function(f) {
      paste0(g, ".csv") %in% list.files(paste0("./data/raw/Daily_data/cluster_", f, "/"))          ## nolint
    }))[1]
    data <- read.csv(paste0("./data/raw/Daily_data/cluster_", cluster_id, "/", g, ".csv"))         ## nolint
    if (length(data$precipitation_mmd) < (365 * 20)) {
      FALSE
    } else {
      TRUE
    }
  })]
  close(p)
}
###############################################################################

###############################################################################
## clean, format, and save data for all usable catchments######################
if (TRUE) {
  p <- txtProgressBar(min = 0, max = length(gs), style = 3)
  for (g in gs) {
    cluster_id <- which(sapply(seq_len(9), function(f) {
      paste0(g, ".csv") %in% list.files(paste0("./data/raw/Daily_data/cluster_", f, "/"))            ## nolint
    }))[1]
    data <- read.csv(paste0("./data/raw/Daily_data/cluster_", cluster_id, "/", g, ".csv"))           ## nolint
    leapdayidxs <- sapply(data$date, function(date) {
      splitdta <- tail(str_split(date, "-")[[1]], 2)
      paste0(splitdta[1], splitdta[2]) == "0229"
    })
    data <- data[!leapdayidxs, ]
    first_date_split_data <- tail(str_split(data$date[1], "-")[[1]], 2)
    if (paste0(first_date_split_data[1], first_date_split_data[2]) != "0101") {
      print()
      stop("ERROR")
    }
    nyears <- as.integer(floor(length(data$precipitation_mmd) / 365))   ## ignore partial years at the end of data  ## nolint
    x <- matrix(0, nrow = nyears, ncol = 365)
    y <- matrix(0, nrow = nyears, ncol = 365)
    for (year_idx in seq_len(nyears)) {
      start <- 1 + ((year_idx - 1) * 365)
      end <- year_idx * 365
      x[year_idx, ] <- data$precipitation_mmd[start : end]
      y[year_idx, ] <- data$streamflow_mmd[start : end]
    }
    saveRDS(x, file = paste0("./data/processed/grid", g, "_rain.rds"))  ## nolint
    saveRDS(y, file = paste0("./data/processed/grid", g, "_streamflow.rds"))  ## nolint
    setTxtProgressBar(p, which(gs == g)[1])
  }
}
###############################################################################