# library(stringr)

# scores <- readRDS("./clustering/precomputed_score_df.rds")

# unique_gcs <- unique(sapply(scores$model_identifier, function(fn) str_split(fn, "_")[[1]][1]))  ## nolint

# centered_scores <- data.frame(
#   grid = unique_gcs,
#   inf_kge = sapply(unique_gcs, function(g) scores$inf_kge[scores$model_identifier == paste0(g, "_centered")]),    ## nolint
#   pred_kge = sapply(unique_gcs, function(g) scores$pred_kge[scores$model_identifier == paste0(g, "_centered")]),   ## nolint
#   inf_r2 = sapply(unique_gcs, function(g) scores$inf_r2[scores$model_identifier == paste0(g, "_centered")]),    ## nolint
#   pred_r2 = sapply(unique_gcs, function(g) scores$pred_r2[scores$model_identifier == paste0(g, "_centered")])    ## nolint
# )

# uncentered_scores <- data.frame(
#   grid = unique_gcs,
#   inf_kge = sapply(unique_gcs, function(g) scores$inf_kge[scores$model_identifier == paste0(g, "_uncentered")]),    ## nolint
#   pred_kge = sapply(unique_gcs, function(g) scores$pred_kge[scores$model_identifier == paste0(g, "_uncentered")]),   ## nolint
#   inf_r2 = sapply(unique_gcs, function(g) scores$inf_r2[scores$model_identifier == paste0(g, "_uncentered")]),    ## nolint
#   pred_r2 = sapply(unique_gcs, function(g) scores$pred_r2[scores$model_identifier == paste0(g, "_uncentered")])    ## nolint
# )

# cat("num grids with uncentered better according to metric (inf kge, pred kge, inf r2, pred r2) : \n",   ## nolint
#     "(", sum(uncentered_scores$inf_kge > centered_scores$inf_kge, na.rm = TRUE),
#     ",", sum(uncentered_scores$pred_kge > centered_scores$pred_kge, na.rm = TRUE),  ## nolint
#     ",", sum(uncentered_scores$inf_r2 > centered_scores$inf_r2, na.rm = TRUE),
#     ",", sum(uncentered_scores$pred_r2 > centered_scores$pred_r2, na.rm = TRUE), ") \n")  ## nolint


# ## average (KGE, R2) win margin for centered scores in inferential models
# diff_inf_kge <- centered_scores$inf_kge - uncentered_scores$inf_kge
# median(diff_inf_kge[diff_inf_kge > 0], na.rm = TRUE)
# median(diff_inf_kge, na.rm = TRUE)

# diff_inf_r2 <- centered_scores$inf_r2 - uncentered_scores$inf_r2
# median(diff_inf_r2[diff_inf_r2 > 0], na.rm = TRUE)
# median(diff_inf_r2, na.rm = TRUE)


# ## average (KGE, R2) win margin for uncentered scores in predictive models
# diff_pred_kge <- uncentered_scores$pred_kge - centered_scores$pred_kge
# median(diff_pred_kge[diff_pred_kge > 0], na.rm = TRUE)
# median(diff_pred_kge, na.rm = TRUE)

# diff_pred_r2 <- uncentered_scores$pred_r2 - centered_scores$pred_r2
# median(diff_pred_r2[diff_pred_r2 > 0], na.rm = TRUE)
# median(diff_pred_r2, na.rm = TRUE)


# pred_kge_thresh_test <- 0.4
# num_left <- sum(scores$pred_kge >= pred_kge_thresh_test, na.rm = TRUE)
# unique_gcs_left <- unique(sapply(
#   scores[scores$pred_kge > pred_kge_thresh_test, "model_identifier"],
#   function(f) {
#     str_split(f, "_")[[1]][1]
#   }
# ))
# length(unique_gcs_left)
# print(length(unique_gcs_left) / length(unique_gcs))
# frac_left <- num_left / nrow(scores)
# frac_left


#### checking to see if using pred kge for dszn/seasonal selection changes 
#### how many negative models we have to disregard down the line
library(stringr)
beta_mtx <- readRDS("./clustering/precomputed_beta_mtx.rds")
precomputed_score_df <- readRDS("./clustering/precomputed_score_df.rds")

cent_models <- precomputed_score_df$model_identifier[sapply(
  precomputed_score_df$model_identifier,
  function(s) str_detect(s, "_centered")
)]
uncent_models <- precomputed_score_df$model_identifier[sapply(
  precomputed_score_df$model_identifier,
  function(s) str_detect(s, "_uncentered")
)]

n_neg_cent_models <- sum(apply(beta_mtx[cent_models, ], 1, function(row) {
  any(row < 0)
}))
n_neg_uncent_models <- sum(apply(beta_mtx[uncent_models, ], 1, function(row) {
  any(row < 0)
}))

get_gc <- function(fn) str_remove_all(fn, "[^\\d]")

all_grids <- unname(sapply(cent_models, get_gc))
stopifnot(setequal(unname(sapply(uncent_models, get_gc)), all_grids))

models_we_use1 <- unname(sapply(all_grids, function(gc) {
  if (
    precomputed_score_df$inf_kge[precomputed_score_df$model_identifier == paste0(gc, "_centered")] > precomputed_score_df$inf_kge[precomputed_score_df$model_identifier == paste0(gc, "_uncentered")]
  ) {
    paste0(gc, "_centered")
  } else {
    paste0(gc, "_uncentered")
  }
}))
models_we_use2 <- unname(sapply(all_grids, function(gc) {
    cent_score <- precomputed_score_df$pred_kge[precomputed_score_df$model_identifier == paste0(gc, "_centered")]
    uncent_score <- precomputed_score_df$pred_kge[precomputed_score_df$model_identifier == paste0(gc, "_uncentered")]
    if (is.na(cent_score)) cent_score <- -Inf
    if (is.na(uncent_score)) uncent_score <- -Inf
  if (
    cent_score > uncent_score
  ) {
    paste0(gc, "_centered")
  } else {
    paste0(gc, "_uncentered")
  }
}))

num_centered_in_scheme_1 <- sum(as.numeric(sapply(models_we_use1, function(q) str_split(q, "_")[[1]][2] == "centered")))    ## nolint
num_centered_in_scheme_2 <- sum(as.numeric(sapply(models_we_use2, function(q) str_split(q, "_")[[1]][2] == "centered")))    ## nolint

dq_on_score1 <- sapply(models_we_use1, function(id_i) {
  pred_score <- precomputed_score_df$pred_kge[precomputed_score_df$model_identifier == id_i]
  bench_score <- precomputed_score_df$benchmark_kge[precomputed_score_df$model_identifier == id_i]
  if (!is.finite(bench_score) && is.finite(pred_score)) {
    FALSE
  } else if (!is.finite(pred_score)) {
    TRUE
  } else {
    stopifnot(is.finite(pred_score) && is.finite(bench_score))
    bench_score >= pred_score
  }
})
dq_on_score2 <- sapply(models_we_use2, function(id_i) {
  pred_score <- precomputed_score_df$pred_kge[precomputed_score_df$model_identifier == id_i]
  bench_score <- precomputed_score_df$benchmark_kge[precomputed_score_df$model_identifier == id_i]
  if (!is.finite(bench_score) && is.finite(pred_score)) {
    FALSE
  } else if (!is.finite(pred_score)) {
    TRUE
  } else {
    stopifnot(is.finite(pred_score) && is.finite(bench_score))
    bench_score >= pred_score
  }
})

num_models_dq_on_score1 <- sum(dq_on_score1)
num_models_dq_on_score2 <- sum(dq_on_score2)

remaining_models1 <- models_we_use1[!dq_on_score1]
remaining_models2 <- models_we_use2[!dq_on_score2]

## of the remaining models, how many have at least 1 negative parameter
check1 <- as.numeric(apply(beta_mtx[remaining_models1, ], 1, function(row) any(row < 0)))
sum(check1)
check2 <- as.numeric(apply(beta_mtx[remaining_models2, ], 1, function(row) any(row < 0)))
sum(check2)







###################
do_pre_cluster_filtering <- function(
  gridcodes_, selection_mode = "predictive",
  neg_trim_val = 0
) {
  ## first, select which model (centered/uncentered) to use for each gridcode
  ## based on either predictive or inferential KGE
  if (selection_mode == "predictive") {
    score_key <- "pred_kge"
  } else {
    score_key <- "inf_kge"
  }
  model_ids_to_filter <- sapply(gridcodes_, function(gc_i) {
    cent_score <- precomputed_score_df[precomputed_score_df$model_identifier == paste0(gc_i, "_centered"), score_key]       ## nolint
    uncent_score <- precomputed_score_df[precomputed_score_df$model_identifier == paste0(gc_i, "_uncentered"), score_key]   ## nolint
    if (is.na(cent_score)) {
      cent_score <- -Inf
    }
    if (is.na(uncent_score)) {
      uncent_score <- -Inf
    }
    if (cent_score > uncent_score) paste0(gc_i, "_centered") else paste0(gc_i, "_uncentered")                               ## nolint
  })
  n_cent_models_selected <- sum(sapply(model_ids_to_filter, function(ii) {
    str_detect(ii, "_centered")
  }))

  cat("done choosing uncentered/centered model for", length(gridcodes_),
      "gridcodes based on", selection_mode, "KGE. selected",
      n_cent_models_selected, "centered models and",
      length(gridcodes_) - n_cent_models_selected, "uncentered models.\n")

  ## then, filter model_ids based on whether their pred KGE is above benchmark
  model_ids_good_scores <- model_ids_to_filter[sapply(model_ids_to_filter, function(id_i) {                                 ## nolint
    pred_score <- precomputed_score_df$pred_kge[precomputed_score_df$model_identifier == id_i]                              ## nolint
    bench_score <- precomputed_score_df$benchmark_kge[precomputed_score_df$model_identifier == id_i]                        ## nolint

    if (is.na(pred_score)) return(FALSE)
    if (is.na(bench_score)) return(TRUE)
    bench_score < pred_score
  })]
  cat("finished filtering all candidate models based on their performance",
      "compared to benchmark KGE.", length(model_ids_good_scores), "models",
      "remain\n")

  ## finally, trim out any models that have negative beta(s, t) ordinates
  model_ids_to_use <- model_ids_good_scores[!sapply(model_ids_good_scores, function(id_i_) {                                 ## nolint
    b <- beta_mtx[id_i_, ]

    num_neg_params <- sum(b < 0)
    num_nonzero_params <- sum(b != 0)
    frac_neg <- num_neg_params / num_nonzero_params
    # cat("num_neg_params:", num_neg_params,
    #     "num_nonzero_params:", num_nonzero_params,
    #     "frac_neg:", frac_neg, "\n")
    many_neg <- frac_neg > neg_trim_val

    largest_neg_val <- max(c(abs(b[b < 0]), 0))
    largest_pos_val <- max(b[b > 0])
    largest_neg_prop <- largest_neg_val / largest_pos_val

    # cat("largest_neg_val:", largest_neg_val,
    #     "largest_pos_val:", largest_pos_val,
    #     "largest_neg_prop:", largest_neg_prop, "\n")

    big_neg <- largest_neg_prop > neg_trim_val

    # cat("returning: many_neg (", many_neg,
    #     ") || big_neg (", big_neg, ") ==",
    #     many_neg || big_neg, "\n")
    many_neg || big_neg
  })]
  cat("finished filtering any models with negative b(s, t) ordinates",
      "with a threshold of", neg_trim_val, "manner.",
      length(model_ids_to_use), "models remain\n")

  model_ids_to_use
}

## predictive model selection
check <- do_pre_cluster_filtering(all_grids, neg_trim_val = 0)
check <- do_pre_cluster_filtering(all_grids, neg_trim_val = 0.05)
check <- do_pre_cluster_filtering(all_grids, neg_trim_val = 0.1)

## inferential model selection
check <- do_pre_cluster_filtering(all_grids, selection_mode = "inferential", neg_trim_val = 0)
check <- do_pre_cluster_filtering(all_grids, selection_mode = "inferential", neg_trim_val = 0.05)
check <- do_pre_cluster_filtering(all_grids, selection_mode = "inferential", neg_trim_val = 0.1)
