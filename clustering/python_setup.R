# find python on whatever machine is actually running this app
## and require reticulate >= 1.41.
library(reticulate)
if (packageVersion("reticulate") < "1.41") {
  stop("This app requires reticulate >= 1.41 for managed Python environments.")
}
Sys.setenv(RETICULATE_PYTHON = "managed")
reticulate::py_require(
  packages = readLines("requirements.txt", warn = FALSE),
  python_version = "==3.12.10"
)
hdbscan <- reticulate::import("hdbscan")
message("HDBSCAN initialized with Python ", reticulate::py_config()$version_string)   ## nolint
