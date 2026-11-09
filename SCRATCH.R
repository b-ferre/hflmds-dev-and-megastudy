library(ggplot2)

hp_grid <- readRDS("zzz/grid5017_results/ctr_TRUE_opt_r2_gmode_1/hp_search_results_round2.rds")   ## nolint

ggplot(hp_grid, aes(w_h, w_v, colour = score)) +
  geom_point(size = 3) +
  geom_point(
    data = dplyr::slice_max(hp_grid, score, n = 1, with_ties = FALSE),
    shape = 1,
    size = 6,
    stroke = 1.2,
    colour = "black"
  )
boxplot(hp_grid$score)
quantile(hp_grid$score, 0.75)

hp_grid$score2 <- pmax(hp_grid$score, rep(quantile(hp_grid$score, 0.75), nrow(hp_grid)))          ## nolint

ggplot(hp_grid, aes(w_h, w_v, colour = score2)) +
  geom_point(size = 3) +
  geom_point(
    data = dplyr::slice_max(hp_grid, score2, n = 1, with_ties = FALSE),
    shape = 1,
    size = 6,
    stroke = 1.2,
    colour = "black"
  )
boxplot(hp_grid$score2)
quantile(hp_grid$score2, 0.75)


nbuckets <- 5
bottom_bucket <- 0.75
buckets <- seq(from = bottom_bucket, to = 1, length.out = nbuckets)
bucket_lims <- sapply(buckets, function(z) quantile(hp_grid$score, z))
hp_grid$bucket_scores <- sapply(hp_grid$score, function(s_i) {
  bucket_lims[which(bucket_lims > s_i)[1]]
})

ggplot(hp_grid, aes(w_h, w_v, colour = bucket_scores)) +
  geom_point(size = 3) +
  geom_point(
    data = dplyr::slice_max(hp_grid, score, n = 1, with_ties = FALSE),
    shape = 1,
    size = 6,
    stroke = 1.2,
    colour = "black"
  )
