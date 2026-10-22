shinyServer(function(input, output, session) {
  print("hello world")

  ## define and initialize reactive values used throughout the app
  clustering_result <- reactiveVal()
  cluster_rep_cache <- reactiveVal()
  selected_representative <- reactiveVal(list(cluster = "", type = ""))

  ## render a map based on the most recent clustering result
  observeEvent(list(
    clustering_result(),
    input$hide_outliers_toggle,
    input$coloring_regime,
    input$kge_goodness_thresh
  ), {
    cat("SERVER: rebuilding map...\n")

    clust_res <- clustering_result()
    req(clust_res)

    coloring_regime <- input$coloring_regime
    req(coloring_regime)

    labels <- clust_res$labels
    clust_colors <- clust_res$clust_colors
    considered_catchments <- clust_res$all_catchments_considered
    clustered_catchments <- names(labels)

    ## first, get a df of all considered catchments (even those with
    ## insufficient performance to be used for clustering)
    map_df <- map_data
    map_df <- map_data[map_df$model_identifier %in% considered_catchments, ]

    ## for all catchments that were actually clustered, set their cluster
    ## value in the df to their label in {0, 1, ..., n clusters} where a
    ## label of 0 indicates an outlier. if a catchment was not clustered,
    ## it will have an -1 in this column
    map_df$cluster <- rep("-1", nrow(map_df))
    map_df$cluster[
      map_df$model_identifier %in% clustered_catchments
    ] <- as.numeric(unname(
      labels[as.character(map_df$model_identifier[
        map_df$model_identifier %in% clustered_catchments
      ])]
    ))

    stopifnot(!anyNA(map_df$cluster))
    stopifnot(!any(duplicated(as.vector(sapply(map_df$model_identifier, function(id_) { # nolint
      str_split(id_, "_")[[1]][1]
    })))))

    if (input$coloring_regime == "Cluster Membership") {
      ## if viewing clustering results, disregard all unclustered catchments as
      ## to avoid obscuring map with a bunch of random stuff
      map_df <- map_df[map_df$model_identifier %in% clustered_catchments, ]

      map_df$color <- clust_colors[as.character(map_df$cluster)]
      map_df$opacity <- ifelse(map_df$cluster == 0, 0.2, 0.85)
      map_df$radius <- ifelse(map_df$cluster == 0, 2, 3)

      hide_outliers_toggle <- input$hide_outliers_toggle
      if (isTRUE(hide_outliers_toggle)) {
        map_df <- map_df[map_df$cluster != 0L, , drop = FALSE]
      }

    } else if (input$coloring_regime == "Inferential KGE") {
      kge_goodness_thresh <- input$kge_goodness_thresh
      req(kge_goodness_thresh)

      scores <- precomputed_score_df$inf_kge[
        precomputed_score_df$model_identifier %in% considered_catchments
      ]

      map_df$color <- map_scores_to_rg_pal(scores, threshold = kge_goodness_thresh)                                             ## nolint
      map_df$opacity <- rep(1, nrow(map_df))
      map_df$radius <- rep(3, nrow(map_df))

    } else if (input$coloring_regime == "Uncentered KGE - Centered KGE") {

      cent_scores <- precomputed_score_df$inf_kge[
        as.vector(unlist(lapply(str_split(map_df$model_identifier, "_"), function(spl_i) {    ## nolint
          spl_i[2] == "centered"
        })))
      ]
      uncent_scores <- precomputed_score_df$inf_kge[
        as.vector(unlist(lapply(str_split(map_df$model_identifier, "_"), function(spl_i) {    ## nolint
          spl_i[2] == "uncentered"
        })))
      ]
      scores <- uncent_scores - cent_scores

      map_df$color <- map_scores_to_rg_pal(scores, threshold = 0)
      map_df$opacity <- rep(1, nrow(map_df))
      map_df$radius <- rep(3, nrow(map_df))
    }

    output$clustering_map <- renderLeaflet({
      leaflet() |>
        addTiles() |>
        addCircleMarkers(
          data = map_df,
          lng = ~longitude,
          lat = ~latitude,
          radius = ~radius,
          opacity = ~opacity,
          color = ~color,
          layerId = ~model_identifier
        )
    })

    cat("SERVER: map construction finished...\n")
  }, ignoreInit = TRUE, priority = 0)

  ## handle map marker clicks and popup graph displays
  observeEvent(input$clustering_map_marker_click, {

    req(input$clustering_map_marker_click$id)
    clicked_gridcode <- as.character(input$clustering_map_marker_click$id)

    beta_plot <- get_beta_plot(clicked_gridcode)

    output$beta_popup_plot <- renderPlot({
      beta_plot
    })

    showModal(modalDialog(
      title = paste("Catchment", clicked_gridcode),
      plotOutput("beta_popup_plot", height = "500px"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  }, ignoreInit = TRUE)

  ## render (& reactively re-render) representative explorer UI
  ## based on clustering results
  observeEvent(clustering_result(), {
    output$cluster_rep_area <- renderUI({
      if (debugmode) cat("SERVER: running cluster rep UI handling...\n")

      res_i <- req(clustering_result())
      cluster_ids <- sort(unique(res_i$labels[res_i$labels > 0]))
      cluster_colors <- res_i$clust_colors
      nclusters <- length(cluster_ids)

      if (nclusters > 0) {
        named_clust_id_list <- as.character(cluster_ids)
        fluidRow(
          ## build HTML tags that set the background color of any inputs with
          ## the name "cluster_i_rep" for some integer i to the corresponding
          ## cluster color. has the purpose of color-coding the "Select Cluster
          ## to Inspect" dropdown that is created below
          tags$head(
            tags$style(HTML(
              str_flatten(as.vector(sapply(seq_len(nclusters), function(index) {                                  ## nolint
                colorval <- cluster_colors[index]
                selected_id <- named_clust_id_list[index]
                c(
                  ".option[data-value='",
                  selected_id, "']",
                  " { background-color: ", colorval, " !important; }"
                )
              })))
            ))
          ),
          ## actually build the dropdowns with the correct values
          div(style = "padding:30px"),
          fluidRow(
            column(selectInput(
              inputId = "selected_rep_type",
              label = "Select Representative Type",
              choices = c("", "Medoid", "Median"),
              selected = ""
            ), width = 6),
            column(selectInput(
              inputId = "selected_rep_cluster",
              label = "Select Cluster to Inspect",
              choices = c("", named_clust_id_list),
              selected = ""
            ), width = 6)
          ),
          fluidRow(
            plotOutput("representative_plot_display")
          )
        ) |> tagList()
      } else {
        div(paste0("There are no cluster representatives to display as no",
                  " clusters were detected. Consider lowering minpts."),
            style = "color:black; padding-top:30px; padding-left:55px")
      }
    })
  }, ignoreInit = TRUE)

  ## render (& reactively re-render) color scheme selection UI
  ## high priority so that this runs before any map recoloring that
  ## relies on kge_goodness_thresh
  observeEvent(input$coloring_regime, {
    if (debugmode) cat("SERVER: running color scheme tuning UI update with coloring_regime:", input$coloring_regime, "\n")  ## nolint

    if (input$coloring_regime != "Cluster Membership") {
      output$coloring_regime_tuner_UI <- renderUI({
        fluidRow(
          sliderInput(
            inputId = "kge_goodness_thresh",
            label = "Break Point Between Bad (Red) and Good (Green) Values of Selected Coloring Feature", ## nolint
            min = 0,
            max = 1,
            value = 0.2,
            step = 0.05
          )
        )
      })
    } else {
      output$coloring_regime_tuner_UI <- renderUI({
        fluidRow(
          checkboxInput(
            inputId = "hide_outliers_toggle",
            label = "Hide Outliers",
            value = TRUE
          )
        )
      })
    }
  }, ignoreInit = FALSE, priority = 100)

  ## observe selected_rep_type (built inside of UI above)
  ## and keep selected_representative up to date
  observeEvent(input$selected_rep_type, {
    if (debugmode) cat("SERVER: updating selected_representative to have type:",
                       input$selected_rep_type, "...\n")
    isolate({
      selected_representative(list(
        cluster = selected_representative()$cluster,
        type = input$selected_rep_type
      ))
    })
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  ## observe selected_rep_cluster (built inside of UI above)
  ## and keep selected_representative up to date
  observeEvent(input$selected_rep_cluster, {
    if (debugmode) cat("SERVER: updating selected_representative to view cluster:", ## nolint
                       input$selected_rep_cluster, "...\n")
    isolate({
      selected_representative(list(
        cluster = as.character(input$selected_rep_cluster),
        type = selected_representative()$type
      ))
    })
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  ## whenever a new clustering result returns, run this first (high priority)
  ## and clear previous representative cache + reset rep type and rep cluster
  ## values
  observeEvent(clustering_result(), {
    freezeReactiveValue(input, "selected_rep_type")
    freezeReactiveValue(input, "selected_rep_cluster")

    cluster_rep_cache(NULL)
    selected_representative(list(cluster = "", type = ""))

    if (debugmode) {
      cat("SERVER: Cleared representative info. Selected rep:",
          "(cl:", selected_representative()$cluster,
          ", type:", selected_representative()$type, ")\n")
    }
  }, ignoreInit = TRUE, priority = 10)

  ## helper function to render "placeholder" representative plot
  clear_rep_plot <- function(cluster_missing, type_missing) {
    output$representative_plot_display <- renderPlot(
      ggplot(data = data.frame(), aes(x = 1, y = 1)) +
        xlim(0, 2) +
        ylim(0, 2) +
        annotate("text",
                 x = 1, y = 1,
                 label = paste0(
                   rep(paste0(
                     "Select a ",
                     rep("cluster", cluster_missing),
                     rep(" and a ", cluster_missing && type_missing),
                     rep("representative type", type_missing),
                     " to view using the dropdowns above."
                   ), cluster_missing || type_missing),
                   rep(paste0(
                     "loading requested representative plot..."
                   ), !(cluster_missing || type_missing))
                 ),
                 size = 7,
                 fontface = "italic",
                 color = "black") +
        theme_void()
    )
  }

  ## continually display the appropriate cluster representative
  ## given the value of selected_representative and clustering result
  observeEvent(selected_representative(), {
    rep_info <- selected_representative()
    res_i <- clustering_result()
    req(rep_info)
    req(res_i)

    clear_rep_plot(
      cluster_missing = (rep_info$cluster == ""),
      type_missing = (rep_info$type == "")
    )

    if (rep_info$cluster != "" && rep_info$type != "") {
      ## run this after a small delay so that clear_rep_plot runs first
      later::later(function() {
        rep_info <- isolate(selected_representative())
        res_i <- isolate(clustering_result())
        current_rep_cache <- isolate(cluster_rep_cache())
        rep_cache_key <- paste0(rep_info$type, rep_info$cluster)
        if (rep_cache_key %in% names(current_rep_cache)) {                                                                           ## nolint
          if (debugmode) cat("SERVER: loading rep with key \"", rep_cache_key, "\" from cache...\n")                            ## nolint
          output$representative_plot_display <- renderPlot(
            render_beta_and_derived_feats(
              current_rep_cache[[paste0(rep_info$type, rep_info$cluster)]],
              plot_beta_and_derived_feats(
                current_rep_cache[[paste0(rep_info$type, rep_info$cluster)]]
              )
            )
          )
        } else {
          in_cluster_grids_i <- as.vector(names(res_i$labels)[res_i$labels == rep_info$cluster])                        ## nolint
          if (rep_info$type == "Medoid") {
            new_rep <- calculate_medoid_beta(
              in_cluster_model_ids = in_cluster_grids_i,
              distmtx = res_i$dist_mat,
              dszn_bool = res_i$dszn_bool_vec
            )
            current_rep_cache[[paste0(rep_info$type, rep_info$cluster)]] <- new_rep                                     ## nolint
            cluster_rep_cache(current_rep_cache)
            output$representative_plot_display <- renderPlot(
              render_beta_and_derived_feats(
                new_rep,
                plot_beta_and_derived_feats(
                  new_rep
                )
              )
            )
          } else if (rep_info$type == "Median") {
            new_rep <- calculate_median_beta(
              in_cluster_model_ids = in_cluster_grids_i,
              dszn_bool = res_i$dszn_bool_vec
            )
            current_rep_cache[[paste0(rep_info$type, rep_info$cluster)]] <- new_rep                                     ## nolint
            cluster_rep_cache(current_rep_cache)
            output$representative_plot_display <- renderPlot(
              render_beta_and_derived_feats(
                new_rep,
                plot_beta_and_derived_feats(
                  new_rep
                )
              )
            )
          } else {
            stop("invalid representative type (",
                 rep_info$type,
                 ") discovered while attempting to plot")
          }
        }
      }, delay = 0.005)
    }

  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  ## perform re-clustering and update clustering_result
  observeEvent(list(input$recluster_btn), {
    if (debugmode) cat("SERVER: running clustering...\n")

    featwts <- setNames(
      sapply(all_feat_weight_names, function(nm) {
        # cat("about to req", nm, "\n") ## nolint
        # print(input[[nm]]) ## nolint
        req(input[[nm]])
        # cat("done req-ing", nm, "\n") ## nolint
        input[[nm]]
      }),
      all_feat_names
    )
    req(featwts)

    comb_func <- input$comb_func
    req(comb_func)

    clust_algo <- input$clust_algo
    req(clust_algo)

    n_pts_core <- input$n_pts_core
    req(n_pts_core)

    n_pts_clust <- input$n_pts_clust
    req(n_pts_clust)

    linkage_method <- input$linkage_method
    req(linkage_method)

    n_clusts_hclust <- input$n_clusts_hclust
    req(n_clusts_hclust)

    algo_hyperparams <- list(
      n_pts_clust = n_pts_clust,
      n_pts_core = n_pts_core,
      linkage_method = linkage_method,
      n_clusts_hclust = n_clusts_hclust
    )

    dszn_scheme <- input$deseasonalization_selection
    req(dszn_scheme)

    dszn_cutoff <- input$dszn_thresh_selection
    req(dszn_cutoff)

    min_kge_for_clust <- input$min_kge_for_clust
    req(min_kge_for_clust)

    ignore_negative_betas <- input$remove_negative_beta_ests

    new_clustering_res <- get_clustering_result(
      clustering_algorithm = clust_algo,
      distance_feature_weights = featwts,
      distance_feature_comb_func = comb_func,
      deseasonalization_scheme = dszn_scheme,
      deseasonalization_cutoff = dszn_cutoff,
      min_kge_for_inclusion = min_kge_for_clust,
      ignore_negative_beta_ests = ignore_negative_betas,
      clustering_algorithm_params = algo_hyperparams
    )

    if (debugmode) cat("SERVER: clustering finished\n")
    clustering_result(new_clustering_res)
  }, ignoreInit = FALSE)

  ## preset button one (all equal feature weights + hclust)
  observeEvent(input$preset_1, {
    updateSliderInput(session, "delta_l2_wt", value = 1)
    updateSliderInput(session, "delta_l1_wt", value = 1)
    updateSliderInput(session, "delta_cheby_wt", value = 1)
    updateSliderInput(session, "delta_scaled_cheby_wt", value = 1)

    updateSliderInput(session, "smooth_delta_l2_wt", value = 1)
    updateSliderInput(session, "smooth_delta_l1_wt", value = 1)
    updateSliderInput(session, "smooth_delta_cheby_wt", value = 1)
    updateSliderInput(session, "smooth_delta_scaled_cheby_wt", value = 1)

    updateSliderInput(session, "sme_l2_wt", value = 1)
    updateSliderInput(session, "sme_l1_wt", value = 1)
    updateSliderInput(session, "sme_cheby_wt", value = 1)
    updateSliderInput(session, "sme_cheby_wt", value = 1)

    updateSliderInput(session, "tt_l2_wt", value = 1)
    updateSliderInput(session, "tt_l1_wt", value = 1)
    updateSliderInput(session, "tt_cheby_wt", value = 1)
    updateSliderInput(session, "tt_cheby_wt", value = 1)

    updateSliderInput(session, "global_l2_wt", value = 1)

    updateSliderInput(session, "n_pts_core", value = 3)
    updateSliderInput(session, "n_pts_clust", value = 9)

    updateSelectInput(session, "linkage_method", selected = "ward.D2")
    updateSliderInput(session, "n_clusts_hclust", value = 8)

    updateSelectInput(session, "clust_algo", selected = "Hierarchical")

    updateSelectInput(session, "comb_func", selected = "eCDF(+)")
    updateSliderInput(session, "min_kge_for_clust", value = 0.4)
    updateCheckboxInput(session, "remove_negative_beta_ests", value = TRUE)

    updateSelectInput(session, "deseasonalization_selection",
                      selected = "Wherever Meaningfully Better Than Uncentered")
    updateSliderInput(session, "dszn_thresh_selection", value = 0.15)

    updateSelectInput(session, "coloring_regime", selected = "Cluster Membership")
  }, ignoreInit = TRUE)

  ## preset button two (just use global l2 + hclust)
  observeEvent(input$preset_2, {
    updateSliderInput(session, "delta_l2_wt", value = 0)
    updateSliderInput(session, "delta_l1_wt", value = 0)
    updateSliderInput(session, "delta_cheby_wt", value = 0)
    updateSliderInput(session, "delta_scaled_cheby_wt", value = 0)

    updateSliderInput(session, "smooth_delta_l2_wt", value = 0)
    updateSliderInput(session, "smooth_delta_l1_wt", value = 0)
    updateSliderInput(session, "smooth_delta_cheby_wt", value = 0)
    updateSliderInput(session, "smooth_delta_scaled_cheby_wt", value = 0)

    updateSliderInput(session, "sme_l2_wt", value = 0)
    updateSliderInput(session, "sme_l1_wt", value = 0)
    updateSliderInput(session, "sme_cheby_wt", value = 0)
    updateSliderInput(session, "sme_cheby_wt", value = 0)

    updateSliderInput(session, "tt_l2_wt", value = 0)
    updateSliderInput(session, "tt_l1_wt", value = 0)
    updateSliderInput(session, "tt_cheby_wt", value = 0)
    updateSliderInput(session, "tt_cheby_wt", value = 0)

    updateSliderInput(session, "global_l2_wt", value = 1)

    updateSliderInput(session, "n_pts_core", value = 3)
    updateSliderInput(session, "n_pts_clust", value = 9)

    updateSelectInput(session, "linkage_method", selected = "ward.D2")
    updateSliderInput(session, "n_clusts_hclust", value = 8)

    updateSelectInput(session, "clust_algo", selected = "Hierarchical")

    updateSelectInput(session, "comb_func", selected = "+")
    updateSliderInput(session, "min_kge_for_clust", value = -1)
    updateCheckboxInput(session, "remove_negative_beta_ests", value = FALSE)

    updateSelectInput(session, "deseasonalization_selection",
                      selected = "Wherever Not Meaningfully Worse Than Uncentered")   ## nolint
    updateSliderInput(session, "dszn_thresh_selection", value = 0.05)

    updateSelectInput(session, "coloring_regime", selected = "Cluster Membership")    ## nolint
  }, ignoreInit = TRUE)
})