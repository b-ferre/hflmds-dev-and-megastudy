shinyServer(function(input, output, session) {
  print("hello world")

  ## define and initialize reactive values used throughout the app
  clustering_result <- reactiveVal(
    get_hdbscan_labels(
      featwts = setNames(unlist(INITCLUSTERSETTINGS[all_feat_weight_names]), all_feat_names),   ## nolint
      comb_func = INITCLUSTERSETTINGS$comb_func,
      n_pts_core = INITCLUSTERSETTINGS$n_pts_core,
      n_pts_clust = INITCLUSTERSETTINGS$n_pts_clust,
      dszn_bool = get_deseason_bool_vec_from_scheme(INITCLUSTERSETTINGS$deseasonalization_selection, INITCLUSTERSETTINGS$dszn_thresh_selection)     ## nolint
    )
  )
  cluster_rep_cache <- reactiveVal()
  selected_representative <- reactiveVal(list(cluster = "", type = ""))

  ## perform re-clustering and update clustering_result
  observeEvent(list(input$recluster_btn, input$deseasonalization_selection), {
    cat("rerunning clustering...\n")

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

    n_pts_core <- input$n_pts_core
    req(n_pts_core)

    n_pts_clust <- input$n_pts_clust
    req(n_pts_clust)

    comb_func <- input$comb_func
    req(comb_func)

    new_clustering_res <- get_hdbscan_labels(
      featwts = featwts,
      comb_func = comb_func,
      n_pts_core = n_pts_core,
      n_pts_clust = n_pts_clust,
      dszn_bool = get_deseason_bool_vec_from_scheme(input$deseasonalization_selection, input$dszn_thresh_selection)     ## nolint
    )
    cat("clustering finished\n")
    clustering_result(new_clustering_res)
  })

  ## render a map based on the most recent clustering result
  observeEvent(list(
    clustering_result(),
    input$hide_outliers_toggle,
    input$coloring_regime,
    input$deseasonalization_selection,
    input$kge_goodness_thresh
  ), {
    cat("rebuilding map...\n")

    clust_res <- clustering_result()
    req(clust_res)

    clust_res <- clustering_result()$labels

    colorkey <- get_distinct_colors(max(c(0L, clust_res)) + 1L)

    map_df <- map_data
    map_df <- map_data[map_df$model_identifier %in% names(clust_res), ]
    map_df$cluster <- unname(
      clust_res[as.character(map_df$model_identifier)]
    )

    stopifnot(!anyNA(map_df$cluster))

    if (input$coloring_regime == "Cluster Membership") {

      map_df$color <- colorkey[map_df$cluster + 1L]
      map_df$opacity <- ifelse(map_df$cluster == 0L, 0.5, 0.9)

    } else if (input$coloring_regime == "Inferential KGE") {
      scores <- precomputed_score_df$inf_kge[precomputed_score_df$model_identifier %in% names(clust_res)]                             ## nolint

      map_df$color <- map_scores_to_rg_pal(scores, threshold = input$kge_goodness_thresh)                                             ## nolint
      map_df$opacity <- rep(1, nrow(map_df))

    } else if (input$coloring_regime == "Uncentered KGE - Centered KGE") {
      cent_scores <- as.numeric(sapply(map_df$model_identifier, function(q) get_inf_kge(q, TRUE)))                                            ## nolint
      uncent_scores <- as.numeric(sapply(map_df$model_identifier, function(q) get_inf_kge(q, FALSE)))                                         ## nolint
      scores <- uncent_scores - cent_scores

      map_df$color <- map_scores_to_rg_pal(scores, threshold = 0)
      map_df$opacity <- rep(1, nrow(map_df))

    }

    ## FIXME: fix precompute to also compute seasonal signal predictive KGEs and then reactivate below
    # else if (input$coloring_regime == "Inferential KGE - Seasonal Mean KGE") {
    #   mod_scores <- precomputed_score_df$inf_kge[precomputed_score_df$model_identifier %in% names(clust_res)]                         ## nolint                  ## nolint
    #   seasonal_score_df <- precomputed_score_dfs[["seasonal"]]
    #   mean_scores <- as.numeric(sapply(map_df$model_identifier, function(q) {
    #     seasonal_score_df$inf_kge[seasonal_score_df$gridcode == q]
    #   }))
    #   scores <- mod_scores - mean_scores

    #   map_df$color <- map_scores_to_rg_pal(scores, threshold = 0)
    #   map_df$opacity <- rep(1, nrow(map_df))
    # }

    if (isTRUE(input$hide_outliers_toggle)) {
      map_df <- map_df[map_df$cluster != 0L, , drop = FALSE]
    }

    output$clustering_map <- renderLeaflet({
      leaflet() |>
        addTiles() |>
        addCircleMarkers(
          data = map_df,
          lng = ~longitude,
          lat = ~latitude,
          radius = 3,
          opacity = ~opacity,
          color = ~color,
          layerId = ~model_identifier
        )
    })

    cat("map construction finished...\n")
  }, priority = 0)

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
  })

  ## render (& reactively re-render) representative explorer UI
  ## based on clustering results
  output$cluster_rep_area <- renderUI({
    if (debugmode) cat("running cluster rep UI handling...\n")

    res_i <- req(clustering_result())
    cluster_ids <- sort(unique(res_i$labels[res_i$labels > 0]))
    nclusters <- length(cluster_ids)

    if (nclusters > 0) {
      named_clust_id_list <- as.character(cluster_ids)
      colors <- get_distinct_colors(nclusters + 1)[2:(nclusters + 1)]
      fluidRow(
        ## build HTML tags that set the background color of any inputs with the
        ## name "cluster_i_rep" for some integer i to the corresponding cluster
        ## color. has the purpose of color-coding the "Select Cluster to
        ## Inspect" dropdown that is created below
        tags$head(
          tags$style(HTML(
            str_flatten(as.vector(sapply(seq_len(nclusters), function(index) {                                  ## nolint
              colorval <- colors[index]
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

  ## observe selected_rep_type (built inside of UI above)
  ## and keep selected_representative up to date
  observeEvent(input$selected_rep_type, {
    if (debugmode) cat("updating selected_representative to have type:",
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
    if (debugmode) cat("updating selected_representative to view cluster:",
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
      cat("Just cleared representative info. Selected rep:\n")
      print(selected_representative())
    }
  }, priority = 10)

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
          if (debugmode) cat("loading rep with key \"", rep_cache_key, "\" from cache...\n")                            ## nolint
          output$representative_plot_display <- renderPlot(
            plot_beta(current_rep_cache[[paste0(rep_info$type, rep_info$cluster)]])                                     ## nolint
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
              plot_beta(new_rep)
            )
          } else if (rep_info$type == "Median") {
            new_rep <- calculate_median_beta(
              in_cluster_model_ids = in_cluster_grids_i,
              dszn_bool = res_i$dszn_bool_vec
            )
            current_rep_cache[[paste0(rep_info$type, rep_info$cluster)]] <- new_rep                                     ## nolint
            cluster_rep_cache(current_rep_cache)
            output$representative_plot_display <- renderPlot(
              plot_beta(new_rep)
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
})