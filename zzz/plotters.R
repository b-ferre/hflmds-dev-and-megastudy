library(RColorBrewer)
library(lattice)

calculate_delta <- function(bvec, m = 364, max_lag = 150) {
  d <- as.numeric(sapply(seq_len(m + 1), function(t) {
    b_t <- bvec[(((t - 1) * (max_lag)) + 1):(t * (max_lag))]
    rev(which(b_t != 0))[1]
  }))
  d[is.na(d)] <- 0  ## any day with all 0 coeffs will get NA above
  d
}

plot_beta <- function(
  bvec, title = "", ncolors = 20, m = 364, max_lag = 150,
  delta = NULL, smooth_delta = NULL, simplified = FALSE
) {
  bgrid <- cbind(expand.grid(s = seq_len(max_lag), t = seq_len(m + 1)), b = bvec)           ## nolint
  dvec <- calculate_delta(bvec)

  # truncate matrix to avoid plotting tens of lags with B(s, t) = 0
  cutoff <- min(ceiling(max(dvec) / 10) * 10, max_lag)
  bgrid <- bgrid[bgrid$s <= cutoff, , drop = FALSE]

  # Generate colour key: smallest magnitudes are still distinct from white
  blue_pal <- function(n) {
    colorRampPalette(c("#d5ddfe", "#6ca4fe", "#2159e7"))(n)
  }
  red_pal <- function(n) {
    colorRampPalette(c("#f9aa98", "#d85a5e", "#67000D"))(n)
  }

  has_pos <- any(bgrid$b > 0)
  has_neg <- any(bgrid$b < 0)

  if (has_pos && !has_neg) {
    cols <- blue_pal(ncolors)
    at_vals <- seq(0, max(bgrid$b), length.out = ncolors + 1)

  } else if (!has_pos && has_neg) {
    cols <- rev(red_pal(ncolors))
    at_vals <- seq(min(bgrid$b), 0, length.out = ncolors + 1)

  } else if (has_pos && has_neg) {
    rp <- max(bgrid$b)
    rn <- abs(min(bgrid$b))

    # Ensure both signs receive at least one colour
    ncolors <- max(2L, ncolors)
    nb <- max(1L, min(ncolors - 1L, round(ncolors * rp / (rp + rn))))
    nr <- ncolors - nb

    cols <- c(rev(red_pal(nr)), blue_pal(nb))
    at_vals <- c(
      seq(-rn, 0, length.out = nr + 1),
      seq(0, rp, length.out = nb + 1)[-1]
    )

  } else {
    cols <- "#FFFFFF"
    at_vals <- c(-1, 1)
  }

  bgrid$b[bgrid$b == 0] <- NA_real_

  line_cols <- c("Raw delta" = "#E69F00", "Smoothed delta" = "#222222")
  line_types <- c(2, 1)

  # Strip ts/matrix attributes before plotting.
  if (!is.null(delta)) delta <- as.numeric(delta)
  if (!is.null(smooth_delta)) smooth_delta <- as.numeric(smooth_delta)

  xlab_ <- if (simplified) NULL else list("Day of the year", cex = 1.4)
  ylab_ <- if (simplified) NULL else list("Lag (Days)", cex = 1.4)
  scale_cex <- if (simplified) 0.6 else 1.2
  colorkey_width <- if (simplified) 0.6 else 1

  lattice::levelplot(
    b ~ t * s,
    data = bgrid,
    col.regions = cols,
    at = at_vals,
    xlab = xlab_,
    ylab = ylab_,
    scales = list(cex = scale_cex),
    colorkey = list(width = colorkey_width, axis.text = list(cex = scale_cex)),
    aspect = "xy",
    main = list(title, cex = scale_cex),

    panel = function(...) {
      # Draw the surface first, then overlay the lines.
      lattice::panel.levelplot(...)

      if (!is.null(delta)) {
        delta <- delta + 0.5
        lattice::panel.lines(
          x = seq_along(delta),
          y = as.numeric(delta),
          col = line_cols[1],
          lty = line_types[1],
          lwd = 2
        )
      }

      if (!is.null(smooth_delta)) {
        smooth_delta <- smooth_delta + 0.5
        lattice::panel.lines(
          x = seq_along(smooth_delta),
          y = as.numeric(smooth_delta),
          col = line_cols[2],
          lty = line_types[2],
          lwd = 2
        )
      }
    },

    key = if (!is.null(delta) || !is.null(smooth_delta)) {
      present <- c(!is.null(delta), !is.null(smooth_delta))

      list(
        space = "bottom",
        columns = sum(present),
        text = list(names(line_cols)[present]),
        lines = list(
          col = unname(line_cols[present]),
          lty = line_types[present],
          lwd = 2
        )
      )
    } else {
      NULL
    }
  )
}