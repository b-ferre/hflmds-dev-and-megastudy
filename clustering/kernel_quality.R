# heuristic measure of the "realism" of a given B(s, t) surface
kernel_qualities <- function(kernel,
                             max_lag = 150L,
                             support_tolerance = 0) {

  ## validate support tolerance and weight values
  if (support_tolerance < 0 || support_tolerance >= 1) {
    stop("support_tolerance must be between zero and one.")
  }

  ## transform kernel into a matrix for easier computation of roughnesses
  if (is.matrix(kernel)) {
    z <- kernel
  } else {
    if (max_lag < 2 || length(kernel) %% max_lag != 0) {
      stop("max_lag must be an integer >= 2 dividing the vector length.")
    }
    z <- matrix(kernel, ncol = max_lag, byrow = TRUE)
  }

  ## normalize kernel matrix so that our roughness measure is scale invariant
  z <- z / max(abs(z))

  ## pre-calculate a bunch of useful features
  ##  > magnitude is just |B(s, t)|
  ##  > mass is the sum of all magnitudes
  ##  > daily peak is a (m+1)-length vector of the maximum in-day magnitude
  ##  > duration is just delta_vector when support_tolerance == 0, otherwise
  ## it is the index of the last significant in-day coefficient
  magnitude <- abs(z)
  mass <- sum(magnitude)
  n <- nrow(z)
  previous <- c(n, seq_len(n - 1L))
  following <- c(seq.int(2L, n), 1L)
  daily_peak <- apply(magnitude, 1, max)
  active <- magnitude > support_tolerance * daily_peak
  duration <- max.col(active, ties.method = "last")
  duration[daily_peak == 0] <- 0

  # calculate temporal roughness by comparing each dvec(t) to a (circular)
  # median across two window sizes (~25% and ~50% of total cycle length)
  # note: we make use of a factor of 4 tolerance to reduce penalization of
  # normal/physically realistic increases in hydrological memory
  excursions <- list(wdw_91 = rep(0, n), wdw_183 = rep(0, n))
  for (requested_window in c(91, 183)) {
    window <- min(requested_window, n - as.integer(n %% 2 == 0L))
    half <- (window - 1) / 2
    extended <- c(tail(duration, half), duration, head(duration, half))
    local_duration <- runmed(extended, k = window, endrule = "constant")[half + seq_len(n)] ## nolint

    # given local medians, calculate each dvec(t)'s "excursion" from the
    # median as under log transform
    raw_exc <- log(pmax(duration, 1) / local_duration) - log(4)
    ## anywhere where dvec(t) is less than 4 times median dvec is ignored
    bounded_exc <- pmax(0, raw_exc)
    excursions[[paste0("wdw_", requested_window)]] <- bounded_exc
  }
  ## finally, report the larger excursion value across the two windows
  excursion <- pmax(excursions[[1]], excursions[[2]])


  ## to capture secondary peaks in lag, we compare to a unimodal variation
  ## calculate lag variation as B(s + 1, t) - B(s, t)
  lag_variation <- rowSums(abs(z[, -1, drop = FALSE] -
                                 z[, -ncol(z), drop = FALSE]))
  ## unimodal var is (peak - B(1, t)) + (peak - B(150, t))
  unimodal_variation <- 2 * daily_peak - magnitude[, 1L] - magnitude[, ncol(z)]
  excess_variation <- pmax(0, lag_variation - unimodal_variation) /
    pmax(2 * daily_peak, .Machine$double.eps)
  # permit modest secondary shoulders and only charge excess variation beyond
  # one quarter of the main peak
  lag_spikes <- sqrt(mean(excursion^2)) + mean(pmax(0, excess_variation - 0.25))


  # calculate (normalized) first and second differences along t axis in order
  # to detect "squiggles" a la Ali or "striping" a la me
  first <- sum(abs(z[following, , drop = FALSE] - z)) / mass                                        ## nolint
  second <- sum(abs(z[following, , drop = FALSE] - (2 * z) + z[previous, , drop = FALSE])) / mass   ## nolint
  # Include deviations from weekly and monthly local averages, so smooth-looking
  # but rapidly repeated seasonal oscillations are also visible to the score.
  residuals <- vapply(c(7, 31), function(width) {
    radius <- (width - 1) / 2
    smoothed <- z * 0
    for (offset in seq.int(-radius, radius)) {
      smoothed <- smoothed + z[((seq_len(n) - 1L + offset) %% n) + 1L, , drop = FALSE] / width      ## nolint
    }
    sum(abs(z - smoothed)) / mass
  }, numeric(1))
  temporal <- mean(c(first, second, residuals))


  ## count and penalize the total negative parameter weight
  negative <- sum(pmax(-z, 0)) / mass
  negative_peak <- max(pmax(-z, 0))
  negative <- sqrt(negative) + negative_peak


  # penalize "patchy" estimates where coefficients oscillate between positive
  # and negative values across either the lag or time axis
  first_active <- max.col(active, ties.method = "first")
  span <- pmax(duration - first_active + 1L, 1L)
  gap_fraction <- pmax(0, span - rowSums(active)) / span
  # for the time axis, we allow >3 day patches but penalize shorter ones
  short_patch_fraction <- vapply(which(colSums(active) > 0), function(j) {
    occupied <- active[, j]
    if (all(occupied)) return(0)
    start <- which(!occupied)[1L]
    circular <- occupied[((seq_len(n) + start - 1L) %% n) + 1L]
    runs <- rle(circular)
    lengths <- runs$lengths[runs$values]
    if (length(lengths) <= 1L) return(0)
    sum(lengths[lengths <= 3L]) / sum(lengths)
  }, numeric(1))
  patches <- mean(gap_fraction[daily_peak > 0]) + mean(short_patch_fraction)


  # use the final two non-zero lags on each day to estimate B((d(t) + 1), t)
  # and penalize the magnitude of this implied coefficient in order to catch
  # kernels that are not tapering to 0
  days <- which(daily_peak > 0 & duration >= 2L)
  endpoint <- magnitude[cbind(days, duration[days])]
  penultimate <- magnitude[cbind(days, duration[days] - 1L)]
  projected <- pmax(0, endpoint - pmax(0, penultimate - endpoint))
  tail_penalty <- if (length(days)) mean(projected / daily_peak[days]) else 0


  ## combine all component weights into final kernel goodness measure
  list(negative = negative, lag_spikes = lag_spikes,
       temporal = temporal, tail = tail_penalty, patches = patches)
}


kernel_goodness <- function(
  kernel_qualities,
  weights = c(negative = 12, lag_spikes = 6,
              temporal = 10, tail = 1,
              patches = 2),
  include_tail = FALSE
) {

  required <- c("negative", "lag_spikes", "patches", "temporal", "tail")
  if (!setequal(names(weights), required) || any(weights < 0) || !setequal(names(kernel_qualities), required)) {   ## nolint
    stop("weights must be nonnegative and named for the five components.")
  }

  contributions <- as.numeric(kernel_qualities) * weights[names(kernel_qualities)]  ## nolint
  if (!include_tail) contributions["tail"] <- 0
  roughness <- sum(contributions)
  goodness <- 100 * exp(-roughness)

  list(
    roughness = roughness,
    goodness = goodness
  )
}