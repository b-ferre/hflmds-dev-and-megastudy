m <- 2
max_lag <- 150

x <- seq(0, 1, length.out = max_lag)

alpha <- 1.2
beta <- 10
basekern <- dbeta(x, shape1 = alpha, shape2 = beta)
basekern <- basekern / max(basekern)

day1kern <- basekern * 0.95
day2kern <- basekern * 0.75
day3kern <- basekern * 0.3
ovr_kern <- c(day1kern, day2kern, day3kern)

plot(ovr_kern)

## Joe's group norm measure
group_norms_1 <- function(kernel_vals) {
  group_norms <- c()
  for (dayindex in seq_len(m + 1)) {
    day_slice <- kernel_vals[(((dayindex - 1) * max_lag) + 1):(dayindex * max_lag)]   ## nolint
    day_norms <- rev(cumsum(rev(day_slice ** 2)))
    group_norms <- c(group_norms, day_norms)
  }
  group_norms
}

## My group norm measure
group_norms_2 <- function(beta_est) {
  unlist(lapply(seq_len(m + 1), function(day_index) {
      beta_t <- beta_est[(((day_index - 1) * max_lag) + 1):(day_index * max_lag)]     ## nolint
      daytotal <- sum(beta_t)
      cumsums <- rev(cumsum(rev(beta_t)))
      cumsums / daytotal
    }))
}

q <- 0.7

### JOE'S METHOD ####################################################
group_norms_1 <- group_norms_1(ovr_kern)
thresh_val1 <- quantile(group_norms_1, q)
to_zero <- group_norms_1 < thresh_val1

nonzero_kern <- ovr_kern
nonzero_kern[to_zero] <- NA

zeroed_kern <- ovr_kern
zeroed_kern[!to_zero] <- NA

day_i_kern_nz <- nonzero_kern[1:max_lag]
day_f_kern_nz <- nonzero_kern[(2 * max_lag + 1) : (3 * max_lag)]
day_i_kern_z <- zeroed_kern[1:max_lag]
day_f_kern_z <- zeroed_kern[(2 * max_lag + 1):(3 * max_lag)]

## plot day 1
plot(day_i_kern_nz)
points(day_i_kern_z, col = "red")

## plot day 2
plot(day_f_kern_nz)
points(day_f_kern_z, col = "red")
#####################################################################


### MY METHOD ####################################################
group_norms_2 <- group_norms_2(ovr_kern)
thresh_val1 <- 0.05
to_zero <- group_norms_2 < thresh_val1

nonzero_kern <- ovr_kern
nonzero_kern[to_zero] <- NA

zeroed_kern <- ovr_kern
zeroed_kern[!to_zero] <- NA

day_i_kern_nz <- nonzero_kern[1:max_lag]
day_f_kern_nz <- nonzero_kern[(2 * max_lag + 1) : (3 * max_lag)]
day_i_kern_z <- zeroed_kern[1:max_lag]
day_f_kern_z <- zeroed_kern[(2 * max_lag + 1):(3 * max_lag)]

## plot day 1
plot(day_i_kern_nz)
points(day_i_kern_z, col = "red")

## plot day 2
plot(day_f_kern_nz)
points(day_f_kern_z, col = "red")
#####################################################################