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
