n_days <- 365L
max_lag <- 149L

## Return a, periodic, zero-mean pattern with range one which is L/R shifted randomly                 ## nolint
shape_pattern <- function(shape, n = n_days) {
  day <- 0:(n - 1L)

  dayshift <- runif(1, min = 1, max = 364)

  raw <- switch(
    shape,
    constant = rep(0, n),
    semi_flat = sin(2 * pi * (day + dayshift) / n),
    jagged = 8 * sin(pi * (2.5 * (day + dayshift) / n)) + 5 * sin(1.3 * pi * (day + dayshift) / n),   ## nolint
    stop("Unknown shape: ", shape)
  )

  if (shape == "constant") return(raw)
  raw <- raw - mean(raw)
  raw / diff(range(raw))
}

## relative_range is (max(x) - min(x)) / mean(x).
make_delta <- function(shape, magnitude) {
  relative_range <- switch(shape,
                           constant = 0,
                           semi_flat = 0.30,
                           jagged = 0.90)

  delta <- round(magnitude *
                   (1 + relative_range * shape_pattern(shape)))
  pmax(1L, pmin(max_lag, as.integer(delta)))
}

eg_delta_0 <- make_delta(shape = "semi_flat", magnitude = 30)
eg_delta_1 <- eg_delta_0 + make_delta(shape = "jagged", magnitude = 10)
plot(eg_delta_1, ylim = c(0, 150), col = "red", main = "mag30 raw")
points(eg_delta_0, col = "blue")

eg_delta_2 <- make_delta(shape = "semi_flat", magnitude = 90)
eg_delta_3 <- eg_delta_2 + make_delta(shape = "jagged", magnitude = 10)
plot(eg_delta_3, ylim = c(0, 150), col = "red", main = "mag90 raw")
points(eg_delta_2, col = "blue")

plot(eg_delta_0 / eg_delta_1, main = "mag30 ratio")
plot(eg_delta_3 / eg_delta_2, main = "mag90 ratio")

cor(eg_delta_0, eg_delta_1)
