
#' @export
marginal_distribution = function(x, q) {
  zero_prob = mean(x == 0, na.rm = TRUE)
  positive_index = which(x > 0)
  F_empirical = ecdf(x[positive_index])
  max_possible_prob = 1 - 1 / length(positive_index)
  if (q > max_possible_prob) {
    warning("Too large threshold! No observations are considered extreme")
    rm(x, positive_index)
    F_tail = F_empirical
    u = Inf
  } else {
    u = unname(quantile(x[positive_index], q, na.rm = TRUE))
    extreme_index = which(x > u)
    pars = ismev::gpd.fit(x[extreme_index], threshold = u, show = FALSE)$mle
    if (pars[2] < 0) warning("Fit gave negative tail parameter")
    rm(x, positive_index, extreme_index)
    F_tail = function(x) q + (1 - q) * pgpd(x, u = u, σ = pars[1], ξ = pars[2])
  }

  distribution = function(x, jitter = FALSE, remove_zeros = FALSE) {
    res = rep(NA_real_, length(x))
    x_below_u = which(x <= u)
    x_above_u = which(x > u)
    x_is_zero = which(x == 0)
    if (any(x_below_u)) res[x_below_u] = zero_prob + (1 - zero_prob) * F_empirical(x[x_below_u])
    if (any(x_above_u)) res[x_above_u] = zero_prob + (1 - zero_prob) * F_tail(x[x_above_u])
    if (any(x_is_zero)) {
      if (remove_zeros) {
        res[x_is_zero] = NA_real_
      } else if (jitter) {
        res[x_is_zero] = runif(length(x_is_zero), 0, zero_prob)
      } else {
        res[x_is_zero] = zero_prob
      }
    }
    res
  }

  distribution
}

#' @export
inverse_distribution = function(F) {
  u = environment(F)$u
  q = environment(F)$q
  zero_prob = environment(F)$zero_prob
  F_empirical = environment(F)$F_empirical
  gpd_pars = environment(F)$pars
  rm(F)
  Q = function(p) {
    res = rep(NA_real_, length(p))
    # Standardise p from F(x) to F(x | x > 0)
    p = (p - zero_prob) / (1 - zero_prob)
    p_below_zero_prob = which(p <= 0)
    p_below_q = which(0 < p & p <= q)
    p_above_q = which(q < p & p <= 1)
    if (any(p_below_zero_prob)) res[p_below_zero_prob] = 0
    if (any(p_below_q)) {
      res[p_below_q] = unname(quantile(F_empirical, p[p_below_q], type = 1))
    }
    if (any(p_above_q)) {
      res[p_above_q] = qgpd(
        p = (p[p_above_q] - q) / (1 - q),
        u = u,
        σ = gpd_pars[1],
        ξ = gpd_pars[2])
    }
    res
  }
  Q
}

#' @export
pgpd = function(x, u, σ, ξ) {
  if (length(σ) > 1 || length(ξ) > 1) stop("pgpd() only accepts parameters of length 1")
  z = pmax((x - u) / σ, 0)
  if (ξ == 0) {
    res = 1 - exp(-z)
  } else {
    res = 1 - pmax(1 + ξ * z, 0) ^ (-1 / ξ)
  }
  res
}

#' @export
qgpd = function(p, u, σ, ξ) {
  if (length(σ) > 1 || length(ξ) > 1) stop("qgpd() only accepts parameters of length 1")
  if (ξ == 0) {
    res = u - σ * log(1 - p)
  } else {
    res = u + σ * ((1 - p)^-ξ - 1) / ξ
  }
  res
}

#' @export
dgpd = function(x, u, σ, ξ) {
  if (length(σ) > 1 || length(ξ) > 1) stop("pgpd() only accepts parameters of length 1")
  res = rep(0, length(x))
  z = (x - u) / σ
  if (ξ == 0) {
    res[z >= 0] = exp(-z[z >= 0]) / σ
  } else if (ξ > 0){
    res[z >= 0] = (1 / σ) * (1 + ξ * z[z >= 0]) ^ (-(1 / ξ + 1))
  } else if (ξ < 0){
    res[z >= 0 & z < -1 / ξ] = (1 / σ) * (1 + ξ * z[z >= 0 & z < -1 / ξ]) ^ (-(1 / ξ + 1))
  }
  res
}


#' @export
plaplace = function(x) {
  res = x
  below_zero = which(x <= 0)
  above_zero = which(x > 0)
  if (any(below_zero)) res[below_zero] = exp(x[below_zero]) / 2
  if (any(above_zero)) res[above_zero] = 1 - exp(-x[above_zero]) / 2
  res
}

#' @export
qlaplace = function(p) {
  res = p
  below_median = which(p <= .5 & p >= 0)
  above_median = which(p > .5 & p <= 1)
  if (any(below_median)) res[below_median] = log(2 * p[below_median])
  if (any(above_median)) res[above_median] = -log(2 * (1 - p[above_median]))
  res
}

#' @export
dlaplace = function(x) {
  res = x
  below_zero = which(x <= 0)
  above_zero = which(x > 0)
  if (any(below_zero)) res[below_zero] = exp(x[below_zero]) / 2
  if (any(above_zero)) res[above_zero] = exp(-x[above_zero]) / 2
  res
}
