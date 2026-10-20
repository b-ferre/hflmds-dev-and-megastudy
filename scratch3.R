set.seed(1)

get_delta_1 <- function() {
  x <- seq_len(365)
  y <- rep(NA_real_, 365)

  i <- x <= 50
  y[i] <- x[i]^2 / 75 + 5

  i <- x > 50 & x <= 100
  y[i] <- (-(x[i] - 90)^2 + 50^2) / 8^2 + 25

  i <- x > 100 & x <= 110
  y[i] <- 62.5 + 0.25 * (x[i] - 100)^3 / 1000 +
    (x[i] - 100)^2 / 20

  i <- x > 110 & x <= 130
  y[i] <- 70 - (x[i] - 125)^2 / 100

  i <- x > 130 & x <= 170
  y[i] <- 69.75 - (x[i] - 130)^2 / 8^2

  i <- x > 170 & x <= 180
  y[i] <- 1.903 + 44.75 - (x[i] - 165)^2.4 / 5^2

  i <- x > 180 & x <= 240
  y[i] <- ((x[i] - 195)^2 + 200) / 20

  i <- x > 240 & x <= 275
  y[i] <- -(x[i] - 250)^2 / 10 + 121.25

  i <- x > 275 & x <= 290
  y[i] <- (x[i] - 290)^2 / 10 - 22.5 + 58.75

  i <- x > 290
  y[i] <- (x[i] - 365)^2 / 180 + 5

  y
}
get_delta_2 <- function() {
  x <- seq_len(365)
  y_bell <- seq(-182, 182)

  peak_day <- 0
  decay_width <- 80  # Larger = slower decay; smaller = faster decay

  i <- x >= peak_day
  y_bell[i] <- exp(-0.5 * ((x[i] - peak_day) / decay_width)^2)

  y_bell
}

base_delta <- c(seq(from = 5, to = 7, length.out = 170),
                seq(from = 7, to = 6.2, length.out = 60),
                seq(from = 6.2, to = 8.3, length.out = 50),
                seq(from = 8.3, to = 5, length.out = 85))
base_delta <- base_delta / max(base_delta)
diff_delta <- get_delta_1()
diff_delta <- diff_delta / max(diff_delta)

cheby <- function(vec1, vec2) {
  max(abs(vec1 - vec2))
}
cheby_scaled_a <- function(vec1, vec2) {
  cheby_standard <- cheby(vec1, vec2)
  (0.5 * (cheby_standard / mean(vec1))) + (0.5 * (cheby_standard / mean(vec2)))
}
cheby_scaled_b <- function(vec1, vec2) {
  cheby_standard <- cheby(vec1, vec2)
  cheby_standard / mean((vec1 + vec2) / 2)
}
cheby_scaled_c <- function(vec1, vec2) {
  cheby(
    vec1 / sqrt(sum(vec1 ** 2)),
    vec2 / sqrt(sum(vec2 ** 2))
  )
}

scale_1 <- 10
scale_2 <- 50
base_scale_1 <- 10
base_scale_2 <- 50
diff_scale <- 5

d1a <- scale_1 + (base_scale_1 * base_delta)
d1b <- d1a + (diff_scale * diff_delta)

d2a <- scale_2 + (base_scale_2 * base_delta)
d2b <- d2a + (diff_scale * diff_delta)

scores_1 <- data.frame(
  cheby = cheby(d1a, d1b),
  cheby_s_a = cheby_scaled_a(d1a, d1b),
  cheby_s_b = cheby_scaled_b(d1a, d1b),
  cheby_s_c = cheby_scaled_c(d1a, d1b)
)

scores_2 <- data.frame(
  cheby = cheby(d2a, d2b),
  cheby_s_a = cheby_scaled_a(d2a, d2b),
  cheby_s_b = cheby_scaled_b(d2a, d2b),
  cheby_s_c = cheby_scaled_c(d2a, d2b)
)

plauto <- TRUE
if (plauto) {
  plot(d1b, ylim = c(0, 150))
  points(d1a, col = "blue")
  points(d2b, col = "green")
  points(d2a, col = "red")

  print(scores_1)
  print(scores_2)
}