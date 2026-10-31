library(stringr)

scores <- readRDS("./clustering/precomputed_score_df.rds")

unique_gcs <- unique(sapply(scores$model_identifier, function(fn) str_split(fn, "_")[[1]][1]))  ## nolint

centered_scores <- data.frame(
  grid = unique_gcs,
  inf_kge = sapply(unique_gcs, function(g) scores$inf_kge[scores$model_identifier == paste0(g, "_centered")]),    ## nolint
  pred_kge = sapply(unique_gcs, function(g) scores$pred_kge[scores$model_identifier == paste0(g, "_centered")]),   ## nolint
  inf_r2 = sapply(unique_gcs, function(g) scores$inf_r2[scores$model_identifier == paste0(g, "_centered")]),    ## nolint
  pred_r2 = sapply(unique_gcs, function(g) scores$pred_r2[scores$model_identifier == paste0(g, "_centered")])    ## nolint
)

uncentered_scores <- data.frame(
  grid = unique_gcs,
  inf_kge = sapply(unique_gcs, function(g) scores$inf_kge[scores$model_identifier == paste0(g, "_uncentered")]),    ## nolint
  pred_kge = sapply(unique_gcs, function(g) scores$pred_kge[scores$model_identifier == paste0(g, "_uncentered")]),   ## nolint
  inf_r2 = sapply(unique_gcs, function(g) scores$inf_r2[scores$model_identifier == paste0(g, "_uncentered")]),    ## nolint
  pred_r2 = sapply(unique_gcs, function(g) scores$pred_r2[scores$model_identifier == paste0(g, "_uncentered")])    ## nolint
)

cat("num grids with uncentered better according to metric (inf kge, pred kge, inf r2, pred r2) : \n",   ## nolint
    "(", sum(uncentered_scores$inf_kge > centered_scores$inf_kge, na.rm = TRUE),
    ",", sum(uncentered_scores$pred_kge > centered_scores$pred_kge, na.rm = TRUE),  ## nolint
    ",", sum(uncentered_scores$inf_r2 > centered_scores$inf_r2, na.rm = TRUE),
    ",", sum(uncentered_scores$pred_r2 > centered_scores$pred_r2, na.rm = TRUE), ") \n")  ## nolint


## average (KGE, R2) win margin for centered scores in inferential models
diff_inf_kge <- centered_scores$inf_kge - uncentered_scores$inf_kge
median(diff_inf_kge[diff_inf_kge > 0], na.rm = TRUE)
median(diff_inf_kge, na.rm = TRUE)

diff_inf_r2 <- centered_scores$inf_r2 - uncentered_scores$inf_r2
median(diff_inf_r2[diff_inf_r2 > 0], na.rm = TRUE)
median(diff_inf_r2, na.rm = TRUE)


## average (KGE, R2) win margin for uncentered scores in predictive models
diff_pred_kge <- uncentered_scores$pred_kge - centered_scores$pred_kge
median(diff_pred_kge[diff_pred_kge > 0], na.rm = TRUE)
median(diff_pred_kge, na.rm = TRUE)

diff_pred_r2 <- uncentered_scores$pred_r2 - centered_scores$pred_r2
median(diff_pred_r2[diff_pred_r2 > 0], na.rm = TRUE)
median(diff_pred_r2, na.rm = TRUE)


pred_kge_thresh_test <- 0.4
num_left <- sum(scores$pred_kge >= pred_kge_thresh_test, na.rm = TRUE)
unique_gcs_left <- unique(sapply(
  scores[scores$pred_kge > pred_kge_thresh_test, "model_identifier"],
  function(f) {
    str_split(f, "_")[[1]][1]
  }
))
length(unique_gcs_left)
print(length(unique_gcs_left) / length(unique_gcs))
frac_left <- num_left / nrow(scores)
frac_left

cent_models <- precomputed_score_df$model_identifier[sapply(
  precomputed_score_df$model_identifier,
  function(s) str_detect(s, "_centered")
)]
uncent_models <- precomputed_score_df$model_identifier[sapply(
  precomputed_score_df$model_identifier,
  function(s) str_detect(s, "_uncentered")
)]

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

remaining_models <- models_we_use1[!dq_on_score1]








## of the remaining models, how many have at least 1 negative parameter
check1 <- as.numeric(apply(beta_mtx[remaining_models, ], 1, function(row) any(row < 0)))
sum(check1)
check2 <- as.numeric(apply(beta_mtx[remaining_models, ], 1, function(row) {
  many_neg <- (sum(row < 0) / sum(row != 0)) > 0.05
  big_neg <- ((abs(min(row)) / abs(max(row))) > 0.1) && (min(row) < 0)
  many_neg || big_neg
}))
