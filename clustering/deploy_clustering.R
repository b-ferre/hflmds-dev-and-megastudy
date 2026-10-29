# script to deploy the app to rsconnect
deploy_clustering <- function(
  app_dir = "C:/Users/bfeen/OneDrive/Documents/hydro_ML/hflm-ds_model/clustering",    ## nolint
  account = "bfeenrre", app_name = "clustering"
) {
  stopifnot(requireNamespace("rsconnect", quietly = TRUE))
  old <- options(rsconnect.python.enabled = FALSE,
                 rsconnect.max.bundle.size = 5 * 1024^3)
  on.exit(options(old), add = TRUE)
  files <- c("global.R", "server.R", "ui.R", "shared_helpers.R",
             "python_setup.R", "requirements.txt", "map_data.rds",
             "precomputed_beta_mtx.rds", "precomputed_beta_plots.rds",
             "precomputed_feature_dist_mtxs.rds", "precomputed_score_df.rds")
  stopifnot(all(file.exists(file.path(app_dir, files))))
  # set app settings ahead of time to avoid issues on startup
  info <- rsconnect::accountInfo(account, "shinyapps.io")
  client <- rsconnect:::clientForAccount(info)
  apps <- rsconnect::applications(account = account, server = "shinyapps.io")
  hit <- apps[apps$name == app_name, , drop = FALSE]
  if (nrow(hit) == 0L) {
    app <- client$createApplication(app_name, "HFLMDS Clustering Applet",
                                    "shiny", info$accountId, "shiny")
    id <- app$id
  } else {
    stopifnot(nrow(hit) == 1L)
    id <- hit$id[[1]]
  }
  client$setApplicationProperty(id, "application.instances.template", "xxxlarge")       ## nolint
  client$setApplicationProperty(id, "application.shiny.scheduler.max.processes", 1)     ## nolint
  client$setApplicationProperty(id, "application.shiny.timeout.init", 300)
  rsconnect::deployApp(appDir = app_dir, appFiles = files, appName = app_name,
                       appId = id, account = account, server = "shinyapps.io",
                       launch.browser = FALSE)
}
