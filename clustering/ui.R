shinyUI(fluidPage(
  titlePanel("HFLMDS Clustering Results"),
  theme = bslib::bs_theme(version = 5),

  ## distance function tuning + map display
  sidebarLayout(
    ## sidebar panel containing clustering tuning params
    sidebarPanel = sidebarPanel(
      accordion(
        accordion_panel(
          title = "Delta Comparisons",
          sliderInput(
            inputId = "delta_l2_wt",
            label = "L2 Distance",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$delta_l2_wt, step = 0.05
          ),
          sliderInput(
            inputId = "delta_l1_wt",
            label = "L1 Total Transfer Similarity",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$delta_l1_wt, step = 0.05
          ),
          sliderInput(
            inputId = "delta_cheby_wt",
            label = "Chebyshev Distance",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$delta_cheby_wt, step = 0.05
          ),
          sliderInput(
            inputId = "delta_scaled_cheby_wt",
            label = "Scaled Chebyshev Distance",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$delta_scaled_cheby_wt, step = 0.05
          )
        ),
        accordion_panel(
          title = "Delta Comparisons",
          sliderInput(
            inputId = "smooth_delta_l2_wt",
            label = "L2 Distance",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$delta_l2_wt, step = 0.05
          ),
          sliderInput(
            inputId = "smooth_delta_l1_wt",
            label = "L1 Total Transfer Similarity",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$delta_l1_wt, step = 0.05
          ),
          sliderInput(
            inputId = "smooth_delta_cheby_wt",
            label = "Chebyshev Distance",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$delta_cheby_wt, step = 0.05
          ),
          sliderInput(
            inputId = "smooth_delta_scaled_cheby_wt",
            label = "Scaled Chebyshev Distance",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$delta_scaled_cheby_wt, step = 0.05
          )
        ),
        accordion_panel(
          title = "Maximally Effective Lag Vector Comparisons",
          sliderInput(
            inputId = "sme_l2_wt",
            label = "L2 Distance",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$sme_l2_wt, step = 0.05
          ),
          sliderInput(
            inputId = "sme_l1_wt",
            label = "L1 Total Transfer Similarity",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$sme_l1_wt, step = 0.05
          ),
          sliderInput(
            inputId = "sme_cheby_wt",
            label = "Chebyshev Distance",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$sme_cheby_wt, step = 0.05
          ),
          sliderInput(
            inputId = "sme_scaled_cheby_wt",
            label = "Scaled Chebyshev Distance",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$sme_scaled_cheby_wt, step = 0.05
          )
        ),
        accordion_panel(
          title = "Total Transfer Comparisons",
          sliderInput(
            inputId = "tt_l2_wt",
            label = "L2 Total Transfer Similarity",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$tt_l2_wt, step = 0.05
          ),
          sliderInput(
            inputId = "tt_l1_wt",
            label = "L1 Total Transfer Similarity",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$tt_l1_wt, step = 0.05
          ),
          sliderInput(
            inputId = "tt_cheby_wt",
            label = "Chebyshev Distance",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$tt_cheby_wt, step = 0.05
          ),
          sliderInput(
            inputId = "tt_scaled_cheby_wt",
            label = "Scaled Chebyshev Distance",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$tt_scaled_cheby_wt, step = 0.05
          )
        ),
        accordion_panel(
          title = "Other Parameters",
          sliderInput(
            inputId = "n_pts_core",
            label = "Core Neighborhood Size",
            min = 2, max = 20,
            value = INITCLUSTERSETTINGS$n_pts_core, step = 1
          ),
          sliderInput(
            inputId = "n_pts_clust",
            label = "Minimum Cluster Size",
            min = 2, max = 20,
            value = INITCLUSTERSETTINGS$n_pts_clust, step = 1
          ),
          selectInput(
            inputId = "comb_func",
            label = "Combination Function",
            choices = c("test", "test2", "test3", "test4"),
            selected = INITCLUSTERSETTINGS$comb_func
          ),
          selectInput(
            inputId = "deseasonalization_selection",
            label = "Select Degree of Use of Centered Estimates",
            choices = c("Everywhere",
                        "Nowhere",
                        "Wherever Meaningfully Better Than Uncentered",
                        "Wherever Not Meaningfully Worse Than Uncentered"),
            selected = INITCLUSTERSETTINGS$deseasonalization_selection
          ),
          sliderInput(
            inputId = "dszn_thresh_selection",
            label = "Definition of \"Meaningfully Better\" (Above)",
            min = 0, max = 0.5, value = INITCLUSTERSETTINGS$dszn_thresh_selection, step = 0.01
          ),
          sliderInput(
            inputId = "min_kge_for_clust",
            label = "Min Inferential KGE to Use in Clustering",
            min = -1, max = 1, value = INITCLUSTERSETTINGS$min_kge_for_clust, step = 0.01
          )
        ),
        accordion_panel(
          "Display Options",
          fluidRow(
            checkboxInput(
              inputId = "hide_outliers_toggle",
              label = "Hide Outliers",
              value = TRUE
            )
          ),
          fluidRow(
            selectInput(
              inputId = "coloring_regime",
              label = "Coloring Scheme",
              choices = c(
                "Cluster Membership",
                "Inferential KGE",
                "Uncentered KGE - Centered KGE"
                #"Inferential KGE - Seasonal Mean KGE"
              ),
              selected = "Cluster Membership"
            )
          ),
          fluidRow(
            sliderInput(
              inputId = "kge_goodness_thresh",
              label = "\"Good\" KGE Cutoff",
              min = 0,
              max = 1,
              value = 0.2,
              step = 0.05
            )
          )
        ),
        multiple = TRUE,
        open = FALSE
      ),
      fluidRow(
        column(actionButton(
          inputId = "recluster_btn",
          label = "Rerun Clustering",
          icon = icon("sync")
        ), width = 12)
      )
    ),
    ## main panel containing map of catchments and representatives
    mainPanel = mainPanel(
      fluidRow(leafletOutput("clustering_map", height = 500)),
      fluidRow(
        column(uiOutput("cluster_rep_area"), width = 12)
      )
    )
  )
))
