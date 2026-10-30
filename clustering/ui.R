shinyUI(fluidPage(
  titlePanel("HFLMDS Clustering Applet"),
  theme = bslib::bs_theme(version = 5),

  tags$head(tags$style(".no-pointer-events { pointer-events: none; }")),

  ## distance function tuning + map display
  page_navbar(
    ## panel containing clustering tuning params
    nav_panel(
      title = "Clustering and Mapping Parameters",
      accordion(
        accordion_panel(
          title = "Dissimilarity Component Weights",
          ## delta distances
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
              label = "L1 Distance",
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
          ## smoothed delta distances
          accordion_panel(
            title = "Smoothed Delta Comparisons",
            sliderInput(
              inputId = "smooth_delta_l2_wt",
              label = "L2 Distance",
              min = 0, max = 1,
              value = INITCLUSTERSETTINGS$smooth_delta_l2_wt, step = 0.05
            ),
            sliderInput(
              inputId = "smooth_delta_l1_wt",
              label = "L1 Distance",
              min = 0, max = 1,
              value = INITCLUSTERSETTINGS$smooth_delta_l1_wt, step = 0.05
            ),
            sliderInput(
              inputId = "smooth_delta_cheby_wt",
              label = "Chebyshev Distance",
              min = 0, max = 1,
              value = INITCLUSTERSETTINGS$smooth_delta_cheby_wt, step = 0.05
            ),
            sliderInput(
              inputId = "smooth_delta_scaled_cheby_wt",
              label = "Scaled Chebyshev Distance",
              min = 0, max = 1,
              value = INITCLUSTERSETTINGS$smooth_delta_scaled_cheby_wt,
              step = 0.05
            )
          ),
          ## sme vector distances
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
              label = "L1 Distance",
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
          ## tt vector distances
          accordion_panel(
            title = "Total Transfer Comparisons",
            sliderInput(
              inputId = "tt_l2_wt",
              label = "L2 Distance",
              min = 0, max = 1,
              value = INITCLUSTERSETTINGS$tt_l2_wt, step = 0.05
            ),
            sliderInput(
              inputId = "tt_l1_wt",
              label = "L1 Distance",
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
          ## global l2 distance
          accordion_panel(
            title = "Overall L2 Distance",
            sliderInput(
              inputId = "global_l2_wt",
              label = "L2 Distance",
              min = 0, max = 1,
              value = INITCLUSTERSETTINGS$global_l2_wt, step = 0.05
            )
          ),
        ),
        ## clustering algorithm selection & hyperparams
        accordion_panel(
          title = "Clustering Algorithm and Hyperparameters",
          accordion_panel(
            title = "HDBSCAN Parameters",
            sliderInput(
              inputId = "n_pts_core",
              label = "Core Neighborhood Size (HDBSCAN)",
              min = 2, max = 20,
              value = INITCLUSTERSETTINGS$n_pts_core, step = 1
            ),
            sliderInput(
              inputId = "n_pts_clust",
              label = "Minimum Cluster Size (HDBSCAN)",
              min = 2, max = 20,
              value = INITCLUSTERSETTINGS$n_pts_clust, step = 1
            )
          ),
          accordion_panel(
            title = "Hierarchical Clustering Parameters",
            selectInput(
              inputId = "linkage_method",
              label = "Linkage Method (Hierarchical)",
              choices = c("complete",
                          "ward.D",
                          "ward.D2",
                          "single",
                          "average",
                          "mcquitty",
                          "median",
                          "centroid"),
              selected = INITCLUSTERSETTINGS$linkage_method
            ),
            sliderInput(
              inputId = "n_clusts_hclust",
              label = "Number of Clusters (Hierarchical)",
              min = 1, max = 50, step = 1,
              value = INITCLUSTERSETTINGS$n_clusts_hclust
            )
          ),
          accordion_panel(
            title = "K-Medoids Clustering Parameters",
            sliderInput(
              inputId = "n_clusts_kmedoids",
              label = "Number of Clusters (K-Medoids)",
              min = 1, max = 50, step = 1,
              value = INITCLUSTERSETTINGS$n_clusts_kmedoids
            )
          ),
          selectInput(
            inputId = "clust_algo",
            label = "Clustering Algorithm",
            choices = c("HDBSCAN", "Hierarchical", "K-Medoids"),
            selected = INITCLUSTERSETTINGS$clust_algo
          ),
          sliderInput(
            inputId = "pre_clust_outlier_trim",
            label = "Pre-Clustering Outlier Filtering Quantile",
            min = 0, max = 1,
            value = INITCLUSTERSETTINGS$pre_clust_outlier_trim,
            step = 0.01
          ),
          sliderInput(
            inputId = "pre_clust_outlier_neighbor_index",
            label = "Pre-Clustering Outlier Filtering Neighbor Index",
            min = 2, max = 10,
            value = INITCLUSTERSETTINGS$pre_clust_outlier_neighbor_index,
            step = 1
          )
        ),
        ## other params (dist comb, kge threshing, and dszn regime)
        accordion_panel(
          title = "Other Parameters",
          selectInput(
            inputId = "comb_func",
            label = "Combination Function for Distance Components",
            choices = c("+",
                        "eCDF(+)",
                        "+(eCDF)"),
            selected = INITCLUSTERSETTINGS$comb_func
          ),
          sliderInput(
            inputId = "min_kge_for_clust",
            label = "Minimum Inferential KGE to Use in Clustering",
            min = -1, max = 1,
            value = INITCLUSTERSETTINGS$min_kge_for_clust, step = 0.01
          ),
          checkboxInput(
            inputId = "remove_negative_beta_ests",
            label = "Remove Catchments With Many Negative B(s, t) Estimates",
            value = INITCLUSTERSETTINGS$remove_negative_beta_ests
          ),
          ## deseasonalization regime
          accordion_panel(
            title = "Deseasonalization Regime",
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
              min = 0, max = 0.5,
              value = INITCLUSTERSETTINGS$dszn_thresh_selection, step = 0.01
            )
          )
        ),
        ## display options
        accordion_panel(
          "Display Options",
          fluidRow(
            selectInput(
              inputId = "coloring_regime",
              label = "Coloring Scheme",
              choices = c(
                "Cluster Membership",
                "Inferential KGE",
                "Predictive KGE",
                "Predictive KGE - Inferential KGE",
                "Uncentered KGE - Centered KGE",
                "Inferential KGE - Benchmark KGE"
              ),
              selected = INITCLUSTERSETTINGS$coloring_regime
            )
          ),
          fluidRow(
            uiOutput("coloring_regime_tuner_UI")
          )
        ),
        ## presets
        accordion_panel(
          title = "Preset Selection",
          fluidRow(
            actionButton(
              inputId = "default_settings",
              label = "Default Clustering Settings"
            ),
            actionButton(
              inputId = "preset_1",
              label = "Preset 1 (Equal Weights, Hierarchical)"
            ),
            actionButton(
              inputId = "preset_2",
              label = "Preset 2 (Global L2 Only, Hierarchical)"
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
      ),
    ),
    ## panel containing map of catchments and representatives
    nav_panel(
      title = "Results and Visualizations",
      tabsetPanel(
        tabPanel(
          title = "Catchment-Wise Cluster Visualizer",
          div(
            style = "position: relative;",

            leafletOutput("clustering_map", height = 500),

            div(
              style = paste(
                "position: absolute;",
                "top: 12px;",
                "right: 12px;",
                "width: 22%;",
                "height: 170px;",
                "z-index: 1000;",
                "background: rgba(255, 255, 255, 0.95);",
                "border-radius: 8px;",
                "box-shadow: 0 2px 10px rgba(0, 0, 0, 0.2);",
                "overflow: hidden;",
                "pointer-events: none;"
              ),
              plotOutput("plot_vis_area", height = "170px")
            )
          )
        ),
        tabPanel(
          title = "Cluster Heatmap Visualizer",
          fluidRow(
            column(sliderInput(
              inputId = "heatmap_grid_size",
              label = "Heatmap Grid Size (km, Log Base 10 Scale)",
              min = 1, max = 3,
              value = INITCLUSTERSETTINGS$heatmap_grid_size,
              step = 0.01,
              width = "100%"
            ), width = 12)
          ),
          fluidRow(
            leafletOutput("clustering_heatmap", height = 500)
          )
        )
      ),
      fluidRow(
        column(uiOutput("cluster_rep_area"), width = 12)
      )
    )
  )
))
