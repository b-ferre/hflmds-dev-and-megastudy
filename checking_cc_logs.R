library(stringr)

log_dir <- "./hflmds_inf_logs/"

indices_to_check <- sapply(list.files(log_dir), function(fn) str_remove_all(fn, "[^\\d]"))    ## nolint