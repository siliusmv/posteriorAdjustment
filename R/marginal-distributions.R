
#' Compute the empirical marginal distribution of the vector x and return
#' a function F(x) that is equal to the marginal distribution when x > 0
#' and returns NA when x <= 0.
#' @export
marginal_distribution = function(x) {
  positive_index = which(x > 0)
  stopifnot(any(positive_index))
  F = ecdf(x[positive_index])
  rm(x, positive_index)

  distribution = function(x) {
    res = F(x)
    x_not_positive = which(x <= 0)
    if (any(x_not_positive)) res[x_not_positive] = NA_real_
    res
  }

  distribution
}

#' The quantile function of the Laplace distribution
#' @export
qlaplace = function(p) {
  res = p
  below_median = which(p <= .5 & p >= 0)
  above_median = which(p > .5 & p <= 1)
  if (any(below_median)) res[below_median] = log(2 * p[below_median])
  if (any(above_median)) res[above_median] = -log(2 * (1 - p[above_median]))
  res
}
