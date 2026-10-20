source("./hflmds.R")
load("./gc_index.rda")

args <- commandArgs(trailingOnly = TRUE)
index_ <- as.integer(args[1])

if (length(args) != 1L ||
      is.na(index_) ||
      index_ < 1L ||
      index_ > length(gc_index)) {
  stop("Invalid catchment index.")
}

grid <- gc_index[index_]

## load and flatten data vectors
x <- readRDS(paste0("./data/grid", grid, "_rain.rds"))
x <- as.vector(t(x))
y <- readRDS(paste0("./data/grid", grid, "_streamflow.rds"))
y <- as.vector(t(y))

t0 <- Sys.time()

## run and save centered model
model <- hflmds(x, y, center = TRUE, verbose = TRUE, max_ar_to_check = 5, mode = "predictive")       ## nolint
save(model, file = paste0("./results/pred_model_centered_grid", grid, ".rda"))

tc <- Sys.time()
cat("\n\n ##################################################\n",
    "Finished fitting centered model. Total time elapsed:",
    difftime(tc, t0, units = "hours"),
    "\n##################################################\n\n")

## run and save uncentered model
model <- hflmds(x, y, center = FALSE, verbose = TRUE, max_ar_to_check = 5, mode = "predictive")      ## nolint
save(model, file = paste0("./results/pred_model_uncentered_grid", grid, ".rda"))

tc <- Sys.time()
cat("\n\n ##################################################\n",
    "Finished fitting uncentered model. Total time elapsed:",
    difftime(tc, t0, units = "hours"),
    "\n##################################################\n\n")
