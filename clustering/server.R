shinyServer(function(input, output, session) {
  print("hello world")

  ## define and initialize reactive values used throughout the app
  clustering_result <- reactiveVal()
  cluster_rep_cache <- reactiveVal()
  selected_representative <- reactiveVal(list(cluster = "", type = ""))

  ## render a catchment-scale map based on the most recent clustering result
  ## this should be the second thing that runs after a new clust result returns
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

    ## first, get a df of all considered catchments regardless of thresholding
    ## (all available gridcodes are included once, with either centerered or
    ## uncentered tag)
    map_df <- map_data
    map_df <- map_data[map_df$model_identifier %in% considered_catchments, ]
    stopifnot(nrow(map_df) == length(unique_gridcodes))

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

    } else if (input$coloring_regime == "Predictive KGE") {
      kge_goodness_thresh <- input$kge_goodness_thresh
      req(kge_goodness_thresh)

      scores <- precomputed_score_df$pred_kge[
        precomputed_score_df$model_identifier %in% considered_catchments
      ]

      map_df$color <- map_scores_to_rg_pal(scores, threshold = kge_goodness_thresh)                                             ## nolint
      map_df$opacity <- rep(1, nrow(map_df))
      map_df$radius <- rep(3, nrow(map_df))

    } else if (input$coloring_regime == "Uncentered KGE - Centered KGE") {

      cent_scores <- as.numeric(sapply(unique_gridcodes, function(gc_i_) {
        precomputed_score_df$inf_kge[
          precomputed_score_df$model_identifier == paste0(gc_i_, "_centered")
        ]
      }))
      uncent_scores <- as.numeric(sapply(unique_gridcodes, function(gc_i_) {
        precomputed_score_df$inf_kge[
          precomputed_score_df$model_identifier == paste0(gc_i_, "_uncentered")
        ]
      }))
      scores <- uncent_scores - cent_scores

      map_df$color <- map_scores_to_rg_pal(scores, threshold = 0)
      map_df$opacity <- rep(1, nrow(map_df))
      map_df$radius <- rep(3, nrow(map_df))
    } else if (input$coloring_regime == "Predictive KGE - Inferential KGE") {
      pred_scores <- precomputed_score_df$pred_kge[
        precomputed_score_df$model_identifier %in% considered_catchments
      ]
      inf_scores <- precomputed_score_df$inf_kge[
        precomputed_score_df$model_identifier %in% considered_catchments
      ]
      scores <- pred_scores - inf_scores

      map_df$color <- map_scores_to_rg_pal(scores, threshold = 0)                                             ## nolint
      map_df$opacity <- rep(1, nrow(map_df))
      map_df$radius <- rep(3, nrow(map_df))
    } else {
      stop(paste0("invalid color_regime (",
                  input$coloring_regime,
                  ") passed to map generation in server.\n"))
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
  }, ignoreInit = TRUE, priority = 9)

  ## render a heatmap based on the most recent clustering result
  ## this should be the third thing that runs after a new clust result returns
  heatmap_click_data <- reactiveVal(NULL)
  heatmap_revision <- 0
  observeEvent(list(
    clustering_result(),
    input$hide_outliers_toggle,
    input$heatmap_grid_size
  ), {
    cat("SERVER: rebuilding heatmap...\n")

    clust_res <- clustering_result()
    req(clust_res)

    labels <- clust_res$labels
    clust_colors <- clust_res$clust_colors
    clustered_catchments <- names(labels)

    ## first, get a df of all considered catchments (even those with
    ## insufficient performance to be used for clustering)
    map_df <- map_data
    map_df <- map_data[map_df$model_identifier %in% clustered_catchments, ]

    map_df$cluster <- as.numeric(unname(
      labels[as.character(map_df$model_identifier)]
    ))

    ## trim outliers if indicated
    map_df <- map_df[map_df$cluster != "0", ]

    stopifnot(!anyNA(map_df$cluster))
    stopifnot(!any(duplicated(as.vector(sapply(map_df$model_identifier, function(id_) { # nolint
      str_split(id_, "_")[[1]][1]
    })))))

    if (debugmode) cat("SERVER: about to call generate_heatmap_layer()...\n")
    heatmap_sf_object <- generate_heatmap_layer(
      ## the CRS I use here corresponds to WGS84 coords. ## TODO: make sure this aligns with metadata ## nolint
      catchments_sf = st_as_sf(map_df, coords = c("longitude", "latitude"), crs = 4326),              ## nolint
      heatmap_grid_size_km = (10 ** input$heatmap_grid_size)
    )
    if (debugmode) cat("SERVER: done.\n")

    heatmap_revision <- heatmap_revision + 1L
    heatmap_sf_object$click_id <- paste0(
      heatmap_revision, "_", seq_len(nrow(heatmap_sf_object))
    )
    # Use the same projected CRS used to construct the heatmap.
    polygon_members <- sf::st_intersects(
      sf::st_transform(heatmap_sf_object, 3857),
      sf::st_transform(
        sf::st_as_sf(
          map_df,
          coords = c("longitude", "latitude"),
          crs = 4326
        ),
        3857
      )
    )

    model_ids_by_polygon <- lapply(polygon_members, function(rows) {
      unique(as.character(map_df$model_identifier[rows]))
    })
    cluster_ids_by_polygon <- lapply(polygon_members, function(rows) {
      unique(as.character(map_df$cluster[rows]))
    })

    names(model_ids_by_polygon) <- heatmap_sf_object$click_id
    names(cluster_ids_by_polygon) <- heatmap_sf_object$click_id
    heatmap_click_data(list(
      model_ids = model_ids_by_polygon,
      cluster_ids = cluster_ids_by_polygon,
      dszn_bool = clust_res$dszn_bool_vec
    ))

    # # Recalculated when clustering or heatmap geometry changes.
    # cat("about to calculate polygon representatives...\n")
    # representatives_by_polygon <- lapply(
    #   model_ids_by_polygon,
    #   function(model_ids) {
    #     if (length(model_ids) == 0L) return(NULL)

    #     betas <- precomputed_beta_matrix[
    #       model_ids, , drop = FALSE
    #     ]

    #     distances <- clust_res$dist_mat[
    #       model_ids, model_ids, drop = FALSE
    #     ]

    #     medoid_id <- model_ids[which.min(rowSums(distances))]

    #     list(
    #       n = length(model_ids),
    #       medoid_id = medoid_id,
    #       betas = list(
    #         Mean = colMeans(betas),
    #         Median = apply(betas, 2, stats::median),
    #         Medoid = precomputed_beta_matrix[medoid_id, ]
    #       )
    #     )
    #   }
    # )
    # cat("done.\n")

    # heatmap_click_data(representatives_by_polygon)

    heatmap_pal <- colorFactor(
      palette = unname(clust_colors),
      domain = names(clust_colors)
    )

    output$clustering_heatmap <- renderLeaflet({
      cat("HEATMAP: renderer started\n")

      hm <- heatmap_sf_object
      hm <- hm[!sf::st_is_empty(hm), ]

      validate(need(nrow(hm) > 0, "Heatmap contains no nonempty geometries."))

      # Use the existing named cluster colours directly.
      hm$fill_colour <- unname(
        clust_colors[as.character(hm$majority_cluster)]
      )

      validate(need(
        !anyNA(hm$fill_colour),
        "Some heatmap clusters have no matching colour."
      ))

      bb <- sf::st_bbox(hm)

      result <- leaflet(hm) |>
        addTiles() |>
        # Explicit panes keep catchment dots beneath the grid polygons.
        addMapPane("catchmentDots", zIndex = 410) |>
        addMapPane("heatmapPolygons", zIndex = 420) |>
        addCircleMarkers(
          data = map_df,
          lng = ~longitude,
          lat = ~latitude,
          radius = 1,
          stroke = FALSE,
          fillColor = "#000000",
          fillOpacity = 1,
          group = "CatchmentDots",
          options = pathOptions(
            pane = "catchmentDots",
            interactive = FALSE
          )
        ) |>
        addPolygons(
          layerId = ~click_id,
          fillColor = ~fill_colour,
          fillOpacity = ~opacity,
          stroke = TRUE,
          color = "#333333",
          weight = 1,
          group = "HeatmapLayer",
          options = pathOptions(pane = "heatmapPolygons")
        ) |>
        fitBounds(
          lng1 = unname(bb["xmin"]),
          lat1 = unname(bb["ymin"]),
          lng2 = unname(bb["xmax"]),
          lat2 = unname(bb["ymax"])
        )

      cat("HEATMAP: widget constructed\n")
      result
    })

    cat("SERVER: heatmap construction finished...\n")
  }, ignoreInit = TRUE, priority = 8)

  ## handle clicks on polygons in heatmap
  # observeEvent(input$clustering_heatmap_shape_click, {
  #   click <- input$clustering_heatmap_shape_click
  #   state <- heatmap_click_data()
  #   req(click$id, state)

  #   model_ids <- state$model_ids[[as.character(click$id)]]
  #   cluster_ids <- state$cluster_ids[[as.character(click$id)]]

  #   # ignore clicks from an earlier rendering of the heatmap.
  #   req(!is.null(model_ids))
  #   req(!is.null(cluster_ids))

  #   if (length(model_ids) == 0L) {
  #     showNotification(
  #       "No catchments are located in this polygon.",
  #       type = "message"
  #     )
  #     return()
  #   }

  #   mean_beta <- calculate_mean_beta(
  #     in_cluster_model_ids = model_ids,
  #     dszn_bool = state$dszn_bool
  #   )

  #   mean_plot <- render_beta_and_derived_feats(
  #     mean_beta,
  #     calculate_plot_info(mean_beta)
  #   )

  #   output$heatmap_median_plot <- renderPlot({
  #     grid::grid.newpage()
  #     grid::grid.draw(mean_plot)
  #   })

  #   showModal(modalDialog(
  #     title = paste0(
  #       "Pointwise Mean — ",
  #       length(model_ids),
  #       " Catchment",
  #       if (length(model_ids) == 1) "" else "s",
  #       ", ", length(unique(cluster_ids)),
  #       " Cluster",
  #       if (length(unique(cluster_ids)) == 1) "" else "s"
  #     ),
  #     plotOutput("heatmap_median_plot", height = "500px"),
  #     size = "l",
  #     easyClose = TRUE,
  #     footer = modalButton("Close")
  #   ))
  # }, ignoreInit = TRUE)

  summary_types <- c("Mean", "Median", "Medoid")

  summary_loading_ui <- function() {
    div(
      style = paste(
        "height: 500px;",
        "display: flex;",
        "align-items: center;",
        "justify-content: center;",
        "gap: 12px;",
        "font-size: 20px;",
        "color: #666;"
      ),
      icon("spinner", class = "fa-spin"),
      "Loading..."
    )
  }


  summary_popup <- reactiveValues(
    model_ids = NULL,
    clustering = NULL,
    index = 1L,
    status = "loading",
    error = NULL,
    plot = NULL,
    medoid_id = NULL
  )

  summary_plot_cache <- list()
  summary_request <- 0L

  output$heatmap_summary_label <- renderText({
    ids <- req(summary_popup$model_ids)
    type <- summary_types[summary_popup$index]

    label <- paste0(
      type, " — ", length(ids), " catchment",
      if (length(ids) == 1L) "" else "s"
    )

    if (type == "Medoid" && !is.null(summary_popup$medoid_id)) {
      label <- paste0(label, " — ", summary_popup$medoid_id)
    }

    label
  })

  output$heatmap_summary_plot <- renderPlot({
    req(summary_popup$status == "ready")
    req(summary_popup$plot)

    grid::grid.newpage()
    grid::grid.draw(summary_popup$plot)
  })

  queue_summary_plot <- function() {
    ids <- isolate(summary_popup$model_ids)
    clust <- isolate(summary_popup$clustering)
    type <- summary_types[isolate(summary_popup$index)]
    req(ids, clust)

    summary_request <<- summary_request + 1   ## nolint
    request_id <- summary_request

    summary_popup$status <- "loading"
    summary_popup$error <- NULL
    summary_popup$plot <- NULL
    summary_popup$medoid_id <- NULL
    set_summary_body(summary_loading_ui())

    # Send the loading UI before starting the calculation.
    session$onFlushed(function() {
      later::later(function() {
        if (session$isClosed()) return()
        if (request_id != summary_request) return()

        tryCatch({
          cached <- summary_plot_cache[[type]]

          if (is.null(cached)) {
            medoid_id <- NULL

            if (type == "Medoid") {
              distances <- clust$dist_mat[ids, ids, drop = FALSE]

              medoid_id <- ids[which.min(rowSums(distances))]
              beta <- precomputed_beta_matrix[medoid_id, ]
            } else {
              betas <- precomputed_beta_matrix[
                ids, , drop = FALSE
              ]

              beta <- switch(
                type,
                Mean = colMeans(betas),
                Median = apply(betas, 2, stats::median)
              )
            }

            cached <- list(
              plot = render_beta_and_derived_feats(
                beta,
                calculate_plot_info(beta)
              ),
              medoid_id = medoid_id
            )

            summary_plot_cache[[type]] <<- cached
          }

          if (request_id != summary_request) return()

          summary_popup$plot <- cached$plot
          summary_popup$medoid_id <- cached$medoid_id
          summary_popup$status <- "ready"
          set_summary_body(
            plotOutput("heatmap_summary_plot", height = "500px")
          )
        }, error = function(e) {
          if (request_id != summary_request) return()

          summary_popup$error <- paste(
            "Could not calculate this summary:",
            conditionMessage(e)
          )
          summary_popup$status <- "error"
          set_summary_body(
            div(
              style = "padding: 30px; color: #b00020;",
              summary_popup$error
            )
          )
        })
      }, delay = 0.05)
    }, once = TRUE)
  }

  set_summary_body <- function(ui) {
    removeUI(
      selector = "#heatmap_summary_body > *",
      multiple = TRUE,
      immediate = TRUE,
      session = session
    )

    insertUI(
      selector = "#heatmap_summary_body",
      where = "beforeEnd",
      ui = ui,
      immediate = TRUE,
      session = session
    )
  }

  observeEvent(input$clustering_heatmap_shape_click, {
    click <- input$clustering_heatmap_shape_click
    state <- heatmap_click_data()
    req(click$id, state)

    ids <- state$model_ids[[as.character(click$id)]]

    # Ignore clicks from an older heatmap rendering.
    req(!is.null(ids))

    if (length(ids) == 0L) {
      showNotification(
        "No catchments are located in this polygon.",
        type = "message"
      )
      return()
    }

    clust <- req(clustering_result())

    # Start a fresh cache for this polygon.
    summary_plot_cache <<- list()

    summary_popup$model_ids <- ids
    summary_popup$clustering <- clust
    summary_popup$index <- 2L  # Open on Median.
    summary_popup$status <- "loading"
    summary_popup$error <- NULL
    summary_popup$plot <- NULL
    summary_popup$medoid_id <- NULL

    showModal(modalDialog(
      title = "Catchment summaries",

      div(
        style = paste(
          "display: flex;",
          "align-items: center;",
          "justify-content: space-between;",
          "gap: 16px;",
          "margin-bottom: 12px;"
        ),

        actionButton(
          "heatmap_summary_previous",
          label = NULL,
          icon = icon("chevron-left"),
          title = "Previous summary"
        ),

        div(
          style = "flex: 1; text-align: center; font-weight: 600;",
          textOutput("heatmap_summary_label")
        ),

        actionButton(
          "heatmap_summary_next",
          label = NULL,
          icon = icon("chevron-right"),
          title = "Next summary"
        )
      ),

      div(
        id = "heatmap_summary_body",
        summary_loading_ui()
      ),

      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))

    queue_summary_plot()
  }, ignoreInit = TRUE)


  observeEvent(input$heatmap_summary_previous, {
    req(summary_popup$model_ids)

    summary_popup$index <- (
      (summary_popup$index - 2L) %% length(summary_types)
    ) + 1L

    queue_summary_plot()
  }, ignoreInit = TRUE)


  observeEvent(input$heatmap_summary_next, {
    req(summary_popup$model_ids)

    summary_popup$index <- (
      summary_popup$index %% length(summary_types)
    ) + 1L

    queue_summary_plot()
  }, ignoreInit = TRUE)


  ## handle map marker clicks and popup graph displays for beta surface
  ## and derived features
  observeEvent(input$clustering_map_marker_click, {

    req(input$clustering_map_marker_click$id)
    clicked_gridcode <- as.character(input$clustering_map_marker_click$id)

    rendered_beta_plot <- get_rendered_beta_plot(clicked_gridcode)

    output$beta_popup_plot <- renderPlot({
      grid::grid.newpage()
      grid::grid.draw(rendered_beta_plot)
    })

    showModal(modalDialog(
      title = paste("Catchment", clicked_gridcode),
      plotOutput("beta_popup_plot", height = "500px"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  }, ignoreInit = TRUE)

  ## handle map marker mouseovers and plot_vis_area display
  ## for JUST B(s, t) surface
  # define a placeholder plot so I don't have to regenerate it many times
  placeholder_plot <- grid::textGrob(
    paste(
      "Mouse over a catchment",
      "to view its corresponding",
      "\u03b2(s, t) plot.",
      sep = "\n"
    ),
    gp = grid::gpar(
      fontsize = 14,
      col = "#555555",
      lineheight = 1.4
    )
  )

  current_mouseover_map_id <- reactiveVal(NULL)

  # whenever a new marker is moused over, update our reactive
  # value that tracks this value
  observeEvent(input$clustering_map_marker_mouseover, {
    id <- input$clustering_map_marker_mouseover$id
    req(id)

    current_mouseover_map_id(as.character(id))
  }, ignoreInit = TRUE)

  # whenever we mouse out from a marker, update our reactive
  # value to NULL if the moused out marker matches the currently
  # mouse over marker
  observeEvent(input$clustering_map_marker_mouseout, {
    id <- input$clustering_map_marker_mouseout$id

    if (!is.null(id) &&
          identical(as.character(id), current_mouseover_map_id())) {
      current_mouseover_map_id(NULL)
    }
  }, ignoreInit = TRUE)

  # to avoid issues with "queueing" tons of updates to the beta
  # visualizer, create a debounced version of the current_mouseover_id
  # for the plotter to track
  settled_mouseover_map_id <- shiny::debounce(
    reactive(current_mouseover_map_id()),
    millis = function() {
      # wait for a while before returning to the placeholder in order to permit
      # quick "flips" between catchments in the same cluster for user ease of
      # comparison, and add a much shorter debounce to transfers to new
      # model ids so that accidentally mousing over a bunch of them while
      # moving cursor across the screen doesn't break anything
      if (is.null(current_mouseover_map_id())) {
        750
      } else {
        30  # Brief pause to filter accidental hovers
      }
    }
  )

  # plot a beta plot for whatever the current value of
  # settled_mouseover_map_id iss
  output$plot_vis_area <- renderPlot({
    id <- settled_mouseover_map_id()

    if (is.null(id)) {
      grid::grid.newpage()
      grid::grid.draw(placeholder_plot)
    } else {
      print(get_simple_rendered_beta_plot(id))
    }
  })

  ## render (& reactively re-render) color scheme selection UI
  output$coloring_regime_tuner_UI <- renderUI({
    if (debugmode) cat("SERVER: running color scheme tuning UI update with coloring_regime:", input$coloring_regime, "\n")  ## nolint

    if (input$coloring_regime != "Cluster Membership") {
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
    } else {
      fluidRow(
        checkboxInput(
          inputId = "hide_outliers_toggle",
          label = "Hide Outliers",
          value = TRUE
        )
      )
    }
  })
  ## make sure that opening the display option accordion panel doesn't
  ## trigger an update for the dynamic UI and cause random map rerendering
  outputOptions(output, "coloring_regime_tuner_UI", suspendWhenHidden = FALSE)

  ## render (& reactively re-render) representative explorer UI
  ## based on clustering results
  observeEvent(clustering_result(), {
    output$cluster_rep_area <- renderUI({
      if (debugmode) cat("SERVER: running cluster rep UI handling...\n")

      res_i <- req(clustering_result())
      cluster_ids <- sort(unique(res_i$labels[res_i$labels > 0]))
      cluster_colors <- res_i$clust_colors[as.character(cluster_ids)]
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
              choices = c("", "Medoid", "Median", "Mean"),
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
          plt <- render_beta_and_derived_feats(
            current_rep_cache[[paste0(rep_info$type, rep_info$cluster)]],
            calculate_plot_info(
              current_rep_cache[[paste0(rep_info$type, rep_info$cluster)]]
            )
          )
          output$representative_plot_display <- renderPlot({
            grid::grid.newpage()
            grid::grid.draw(plt)
          })
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
            plt <- render_beta_and_derived_feats(
              new_rep,
              calculate_plot_info(
                new_rep
              )
            )
            output$representative_plot_display <- renderPlot({
              grid::grid.newpage()
              grid::grid.draw(plt)
            })
          } else if (rep_info$type == "Median") {
            new_rep <- calculate_median_beta(
              in_cluster_model_ids = in_cluster_grids_i,
              dszn_bool = res_i$dszn_bool_vec
            )
            current_rep_cache[[paste0(rep_info$type, rep_info$cluster)]] <- new_rep                                     ## nolint
            cluster_rep_cache(current_rep_cache)
            plt <- render_beta_and_derived_feats(
              new_rep,
              calculate_plot_info(
                new_rep
              )
            )
            output$representative_plot_display <- renderPlot({
              grid::grid.newpage()
              grid::grid.draw(plt)
            })
          } else if (rep_info$type == "Mean") {
            new_rep <- calculate_mean_beta(
              in_cluster_model_ids = in_cluster_grids_i,
              dszn_bool = res_i$dszn_bool_vec
            )
            current_rep_cache[[paste0(rep_info$type, rep_info$cluster)]] <- new_rep                                     ## nolint
            cluster_rep_cache(current_rep_cache)
            plt <- render_beta_and_derived_feats(
              new_rep,
              calculate_plot_info(
                new_rep
              )
            )
            output$representative_plot_display <- renderPlot({
              grid::grid.newpage()
              grid::grid.draw(plt)
            })
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

    n_clusts_kmedoids <- input$n_clusts_kmedoids
    req(n_clusts_kmedoids)

    algo_hyperparams <- list(
      n_pts_clust = n_pts_clust,
      n_pts_core = n_pts_core,
      linkage_method = linkage_method,
      n_clusts_hclust = n_clusts_hclust,
      n_clusts_kmedoids = n_clusts_kmedoids
    )

    dszn_scheme <- input$deseasonalization_selection
    req(dszn_scheme)

    dszn_cutoff <- input$dszn_thresh_selection
    req(dszn_cutoff)

    min_kge_for_clust <- input$min_kge_for_clust
    req(min_kge_for_clust)

    ignore_negative_betas <- input$remove_negative_beta_ests

    outlier_thresh <- input$pre_clust_outlier_trim
    req(outlier_thresh)

    outlier_neighbor_index <- input$pre_clust_outlier_neighbor_index
    req(outlier_neighbor_index)

    new_clustering_res <- get_clustering_result(
      clustering_algorithm = clust_algo,
      distance_feature_weights = featwts,
      distance_feature_comb_func = comb_func,
      deseasonalization_scheme = dszn_scheme,
      deseasonalization_cutoff = dszn_cutoff,
      min_kge_for_inclusion = min_kge_for_clust,
      ignore_negative_beta_ests = ignore_negative_betas,
      clustering_algorithm_params = algo_hyperparams,
      outlier_thresh = outlier_thresh,
      outlier_neighbor_index = outlier_neighbor_index
    )

    if (debugmode) cat("SERVER: clustering finished\n")
    clustering_result(new_clustering_res)
  }, ignoreInit = FALSE)

  ## helper function for presets that takes a list of named
  ## values for the applet and sets all relevant inputs to them
  make_settings <- function(settings_list) {
    ## distance function component weights
    updateSliderInput(session, "delta_l2_wt",
                      value = settings_list$delta_l2_wt)
    updateSliderInput(session, "delta_l1_wt",
                      value = settings_list$delta_l1_wt)
    updateSliderInput(session,
                      "delta_cheby_wt",
                      value = settings_list$delta_cheby_wt)
    updateSliderInput(session,
                      "delta_scaled_cheby_wt",
                      value = settings_list$delta_scaled_cheby_wt)
    updateSliderInput(session,
                      "smooth_delta_l2_wt",
                      value = settings_list$smooth_delta_l2_wt)
    updateSliderInput(session,
                      "smooth_delta_l1_wt",
                      value = settings_list$smooth_delta_l1_wt)
    updateSliderInput(session,
                      "smooth_delta_cheby_wt",
                      value = settings_list$smooth_delta_cheby_wt)
    updateSliderInput(session,
                      "smooth_delta_scaled_cheby_wt",
                      value = settings_list$smooth_delta_scaled_cheby_wt)
    updateSliderInput(session,
                      "sme_l2_wt",
                      value = settings_list$sme_l2_wt)
    updateSliderInput(session,
                      "sme_l1_wt",
                      value = settings_list$sme_l1_wt)
    updateSliderInput(session,
                      "sme_cheby_wt",
                      value = settings_list$sme_cheby_wt)
    updateSliderInput(session,
                      "sme_scaled_cheby_wt",
                      value = settings_list$sme_scaled_cheby_wt)
    updateSliderInput(session,
                      "tt_l2_wt",
                      value = settings_list$tt_l2_wt)
    updateSliderInput(session,
                      "tt_l1_wt",
                      value = settings_list$tt_l1_wt)
    updateSliderInput(session,
                      "tt_cheby_wt",
                      value = settings_list$tt_cheby_wt)
    updateSliderInput(session,
                      "tt_scaled_cheby_wt",
                      value = settings_list$tt_scaled_cheby_wt)
    updateSliderInput(session,
                      "global_l2_wt",
                      value = settings_list$global_l2_wt)

    ## clustering algorithm hyperparams
    updateSliderInput(session,
                      "n_pts_core",
                      value = settings_list$n_pts_core)
    updateSliderInput(session,
                      "n_pts_clust",
                      value = settings_list$n_pts_clust)
    updateSelectInput(session,
                      "linkage_method",
                      selected = settings_list$linkage_method)
    updateSliderInput(session,
                      "n_clusts_hclust",
                      value = settings_list$n_clusts_hclust)
    updateSliderInput(session,
                      "n_clusts_kmedoids",
                     value = settings_list$n_clusts_kmedoids)
    updateSelectInput(session,
                      "clust_algo",
                      selected = settings_list$clust_algo)

    ## other clustering-relevant inputs
    updateSelectInput(session,
                      "comb_func",
                      selected = settings_list$comb_func)
    updateSliderInput(session,
                      "min_kge_for_clust",
                      value = settings_list$min_kge_for_clust)
    updateCheckboxInput(session,
                        "remove_negative_beta_ests",
                        value = settings_list$remove_negative_beta_ests)
    updateSliderInput(session,
                      "pre_clust_outlier_trim",
                      value = settings_list$pre_clust_outlier_trim)
    updateSliderInput(session,
                      "pre_clust_outlier_neighbor_index",
                      value = settings_list$pre_clust_outlier_neighbor_index)   ## nolint

    ## deseasonalization regime
    updateSelectInput(session, "deseasonalization_selection",
                      selected = settings_list$deseasonalization_selection)     ## nolint
    updateSliderInput(session, "dszn_thresh_selection",
                      value = settings_list$dszn_thresh_selection)

    ## display parameters
    updateSelectInput(session, "coloring_regime",
                      selected = settings_list$coloring_regime)
    updateSliderInput(session, "heatmap_grid_size",
                      value = settings_list$heatmap_grid_size)
  }

  ## default button (resets all values to those in INITCLUSTERSETTINGS)
  observeEvent(input$default_settings, {
    make_settings(INITCLUSTERSETTINGS)
  }, ignoreInit = TRUE)

  ## preset button one:
  ## all features weighted equally and combined with "+"
  ## all catchments > 0.3 KGE considered
  ## remove negative B(s, t) estimates
  ## 5% outlier trimming with k = 2
  ## use whichever (centered, uncentered) has higher inf kge
  ## default heatmap grid size
  ## color according to cluster membership
  ## cluster using hclust with 12 clusters
  observeEvent(input$preset_1, {
    settings_to_use <- c(
      setNames(as.list(
        rep(1, 17)
      ), all_feat_weight_names),
      list(
        comb_func = "+",
        clust_algo = "Hierarchical",
        n_pts_core = 2,
        n_pts_clust = 10,
        linkage_method = "ward.D2",
        n_clusts_hclust = 12,
        n_clusts_kmedoids = 12,
        deseasonalization_selection = "Wherever Meaningfully Better Than Uncentered",                                                     ## nolint
        dszn_thresh_selection = 0,
        min_kge_for_clust = 0.3,
        pre_clust_outlier_trim = 0.05,
        pre_clust_outlier_neighbor_index = 2,
        remove_negative_beta_ests = TRUE,
        heatmap_grid_size = INITCLUSTERSETTINGS$heatmap_grid_size,
        coloring_regime = "Cluster Membership"
      )
    )
    make_settings(settings_to_use)
  }, ignoreInit = TRUE)

  ## preset button two:
  ## only use global L2
  ## all catchments > 0 KGE considered
  ## remove negative B(s, t) estimates
  ## 5% outlier trimming with k = 2
  ## use whichever (centered, uncentered) has higher inf kge
  ## default heatmap grid size
  ## color according to cluster membership
  ## cluster using hclust with 12 clusters
  observeEvent(input$preset_2, {
    settings_to_use <- c(
      setNames(as.list(
        c(rep(0, 16), 1)
      ), all_feat_weight_names),
      list(
        comb_func = "+",
        clust_algo = "Hierarchical",
        n_pts_core = 2,
        n_pts_clust = 10,
        linkage_method = "ward.D2",
        n_clusts_hclust = 12,
        n_clusts_kmedoids = 12,
        deseasonalization_selection = "Wherever Meaningfully Better Than Uncentered",                                                     ## nolint
        dszn_thresh_selection = 0,
        min_kge_for_clust = 0.3,
        pre_clust_outlier_trim = 0.05,
        pre_clust_outlier_neighbor_index = 2,
        remove_negative_beta_ests = TRUE,
        heatmap_grid_size = INITCLUSTERSETTINGS$heatmap_grid_size,
        coloring_regime = "Cluster Membership"
      )
    )
    make_settings(settings_to_use)
  }, ignoreInit = TRUE)
})