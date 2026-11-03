library(stringr)

log_dir <- "./hflmds_inf_logs/"

indices_to_check <- list.files(log_dir)    ## nolint
n <- length(indices_to_check)

check <- data.frame(
  fname = indices_to_check,
  gc_index = sapply(indices_to_check, function(i_) {
    str_split(str_split(i_, ".out")[[1]][1], "_")[[1]][3]
  }),
  max_unconstr = rep(NA, n),
  max_constr = rep(NA, n),
  abs_diff = rep(NA, n),
  rel_diff = rep(NA, n),
  q_selection = rep(NA, n)
)

p <- txtProgressBar(min = 0, max = n, style = 3)
for (i_ in seq_len(n)) {
  fname <- check$fname[i_]

  ## 51 lines from end of file to beginning of final fit
  lines <- tail(readLines(paste0("./hflmds_inf_logs/", fname)), 51)

  ## lines 12-19 are r2s from unconstrained fit GP
  max_unconstrained_r2 <- max(sapply(seq(from = 12, to = 19, by = 1), function(j) {     ## nolint
    as.numeric(str_split(str_split(lines[j], "R2 is ")[[1]][2], " , w_h is")[[1]][1])   ## nolint
  }))

  ## line 21 is the chosen q
  q_sel <- as.numeric(str_split(lines[21], "q selected: ")[[1]][2])

  ## lines 33-40 are r2s from thresholded fit GP
  max_constrained_r2 <- max(sapply(seq(from = 33, to = 40, by = 1), function(j) {       ## nolint
    as.numeric(str_split(str_split(lines[j], "R2 is ")[[1]][2], " , w_h is")[[1]][1])   ## nolint
  }))

  ## calculate the absolute and fractional difference and save it
  abs_diff <- max_unconstrained_r2 - max_constrained_r2
  rel_diff <- abs_diff / max_unconstrained_r2

  ## save all data
  check$max_unconstr[i_] <- max_unconstrained_r2
  check$max_constr[i_] <- max_constrained_r2
  check$abs_diff[i_] <- abs_diff
  check$rel_diff[i_] <- rel_diff
  check$q_selection[i_] <- q_sel
  check$gc_index[i_] <- str_split(str_split(fname, ".out")[[1]][1], "_")[[1]][3]

  ## update progress bar
  setTxtProgressBar(p, i_)
}
close(p)