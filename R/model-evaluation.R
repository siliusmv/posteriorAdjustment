
# We can probably remove a lot of this...

#' @export
stwcrps_norm = function(y, μ, σ, y0 = -Inf) {
  S = mapply(expected_twcrps_norm, μ_truth = μ, σ_truth = σ, y0 = y0)
  twcrps_norm(y, μ, σ, y0) / abs(S) + log(abs(S))
}

#' @export
twcrps_norm = function(y, μ, σ, y0 = -Inf) {
  I = function(a) ifelse(a, 1, 0)
  F = function(x, μ, σ) pnorm((x - μ) / σ)
  ϕ = function(x, μ, σ) dnorm((x - μ) / σ)
  ifelse(is.infinite(y0), 0, - y0 * (F(y0, μ, σ) - I(y <= y0)) ^ 2) +
    y * I(y >= y0) * (F(y, μ, σ) * (1 + I(y > y0)) - 1) -
    μ * (1 - F(y0, μ, σ) ^ 2) -
    2 * σ * ( (1 - F(y0, μ, σ / sqrt(2))) / (2 * sqrt(pi)) + ϕ(y0, μ, σ) * F(y0, μ, σ)) +
    2 * μ * (1 - F(pmax(y0, y), μ, σ)) +
    2 * σ * ϕ(pmax(y0, y), μ, σ)
}

#' @export
permutation_test = function(x, y, fun = mean, n = 1000) {
  stopifnot(length(x) == length(y))
  true_val = fun(x) - fun(y)
  perm_vals = rep(0, n)
  for (i in seq_len(n)) {
    indices = sample.int(2, length(x), replace = TRUE)
    perm_x = c(x[indices == 1], y[indices == 2])
    perm_y = c(y[indices == 1], x[indices == 2])
    perm_vals[i] = fun(perm_x) - fun(perm_y)
  }
  list(prob_perms_larger_than_truth = mean(perm_vals >= true_val),
       prob_abs_perms_larger_than_abs_truth = mean(abs(perm_vals) >= abs(true_val)))
}

#' @export
diebold_mariano = function(x, y) {
  stopifnot(length(x) == length(y))
  if (all(x - y == 0)) return(.5)
  t = mean(x - y) / sd(x - y) * sqrt(length(x))
  1 - pt(t, length(x) - 1)
}

#' @export
rank_of_y_in_x = function(y, x) {
  # This function is based on the paper by Thorarinsdottir et al. (2016) Doi: 10.1080/10618600.2014.977447
  stopifnot(length(y) == nrow(x))
  n = ncol(x) # number of ensemble members
  x = cbind(y, x)

  # for each dimension in x, compute the rank of all n ensemble members and y in that dimension
  ranks = matrixStats::rowRanks(x)

  # Compute the preranks for all n ensemble members and y
  preranks = (n + 1 - ranks) * (ranks - 1) + n
  preranks = matrixStats::colMeans2(preranks)

  average_ranks = matrixStats::colMeans2(ranks)

  # Standardise the preranks to take values 1, 2, ..., n + 1
  preranks = rank(preranks, ties.method = "random")
  average_ranks = rank(average_ranks, ties.method = "random")

  c(preranks[1], average_ranks[1])
}

#' @export
variogram_score = function(y, x, weights = 1, p = 0.5) {
  stopifnot(length(y) == nrow(x))
  y_diff = dist(y, method = "manhattan") ^ p
  x_diff = 0
  for (j in seq_len(ncol(x))) {
    x_diff = x_diff + dist(x[, j], method = "manhattan") ^ p
  }
  x_diff = x_diff / ncol(x)
  mean(weights * as.matrix((y_diff - x_diff) ^ 2), na.rm = TRUE)
}

#' Estimate the extremal correlation where we group together all observations
#' that have similar distances from s0
#' @export
extremal_correlation_in_dist_groups = function(y, y0, dist_to_s0, dist_groups, p) {

  # Place all the distances into their respective dist_groups
  dist_groups = sort(dist_groups)
  if (is.finite(tail(dist_groups, 1))) dist_groups = c(dist_groups, Inf)
  dist_group_index = vector("list", length(dist_groups) - 1)
  for (j in seq_along(dist_groups[-1])) {
    dist_group_index[[j]] = which(dist_to_s0 >= dist_groups[j] & dist_to_s0 < dist_groups[j + 1])
  }

  # Compute the total number of large observations and the total number of observations
  # in each of the dist_groups, for each of the probabilities
  n_big_in_group = matrix(0, nrow = length(dist_groups) - 1, ncol = length(p))
  n_obs_in_group = matrix(0, nrow = length(dist_groups) - 1, ncol = length(p))
  for (i in seq_along(y0)) {
    for (j in seq_along(p)) {
      if (y0[i] > qlaplace(p[j])) {
        for (k in seq_along(dist_groups[-1])) {
          n_big_in_group[k, j] = n_big_in_group[k, j] + sum(y[i, dist_group_index[[k]]] > qlaplace(p[j]))
          n_obs_in_group[k, j] = n_obs_in_group[k, j] + length(dist_group_index[[k]])
        }
      }
    }
  }

  # Compute the percentages of large observations in each dist_group and each probability
  res = n_big_in_group / n_obs_in_group
  colnames(res) = p

  # Remove NaNs
  if (any(n_obs_in_group == 0)) {
    res[which(n_obs_in_group == 0)] = 0
  }

  res
}


#' Calculate chi and chibar at the quantiles u for two random variables
#' This one is heavily based on the implementation in
#' https://github.com/harrysouthworth/texmex/blob/master/R/chi.R
#' @export
extremal_correlation = function(p, x, y, truncate = TRUE) {
  n = length(x)
  if (n != length(y)) stop("The data have different length")
  x_quantiles = rank(x, ties.method = "first") / (n + 1)
  y_quantiles = rank(y, ties.method = "first") / (n + 1)
  quantiles = cbind(x_quantiles, y_quantiles)

  rowmax = apply(quantiles, 1, max)
  rowmin = apply(quantiles, 1, min)
  if (min(p) < min(rowmax)) stop("lower quantile limit is too low")
  if (max(p) > max(rowmax)) stop("upper quantile limit is too high")

  # Find probability of both components being larger than p (C) and of
  # both components being less than p (C_bar)
  C = sapply(p, function(p) mean(rowmax < p))
  C_bar = sapply(p, function(p) mean(rowmin > p))

  # Get chi and chibar
  chi = 2 - log(C) / log(p) # 3.2 of Coles, Heffernan, Tawn
  chi_bar = 2 * log(1 - p) / log(C_bar) - 1 # Page 348 of Coles, Heffernan, Tawn

  # Estimate standard deviation using delta method
  σ_chi = sqrt((1 / log(p)^2) / C * (1 - C) / n)
  σ_chi_bar = sqrt((((4 * log(1 - p)^2) / (log(C_bar)^4 * C_bar^2)) *
                      C_bar * (1 - C_bar)) / n)

  # Truncate values of chi or chibar that are outside theoretical bounds
  if (truncate) {
    chivals = truncate_chi(chi, chi_bar, p)
    chi = chivals$chi
    chi_bar = chivals$chi_bar
  }

  data.frame(
    val = c(chi, chi_bar),
    sd = c(σ_chi, σ_chi_bar),
    par = rep(c("chi", "chibar"), each = length(p)),
    p = rep(p, 2))
}

truncate_chi = function(chi, chi_bar, p) {
  chi_lb = 2 - log(pmax(2 * p - 1, 0)) / log(p)
  chi_bar_lb = 2 * log(1 - p) / log(1 - 2 * p + pmax(2 * p - 1, 0)) - 1
  chi[chi > 1] = 1
  chi_bar[chi_bar > 1] = 1
  chi = sapply(seq_along(chi), function(i) pmax(chi[i], chi_lb[i]))
  chi_bar = sapply(seq_along(chi_bar), function(i) pmax(chi_bar[i], chi_bar_lb[i]))
  list(chi = chi, chi_bar = chi_bar)
}

expected_stwcrps_norm = function(μ_truth, σ_truth, μ = μ_truth, σ = σ_truth, y0 = -Inf) {
  S = expected_twcrps_norm(μ, σ, y0 = y0)
  f = function(x) dnorm(x, μ_truth, σ_truth) * twcrps_norm(x, μ, σ, y0)
  lower = max(y0, qnorm(1e-6, μ, σ))
  upper = qnorm(1 - 1e-6, μ, σ)
  integral = integrate(f, lower, upper, subdivisions = 500L)$value
  integral / abs(S) + log(abs(S))
}

expected_twcrps_norm = function(μ_truth, σ_truth, μ = μ_truth, σ = σ_truth, y0 = -Inf) {
  f = function(x) dnorm(x, μ_truth, σ_truth) * twcrps_norm(x, μ, σ, y0)
  lower = max(y0, qnorm(1e-6, μ, σ))
  upper = qnorm(1 - 1e-6, μ, σ)
  integrate(f, lower, upper, subdivisions = 500L)$value
}
