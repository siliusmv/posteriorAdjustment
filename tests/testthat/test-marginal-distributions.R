test_that("plaplace() is bounded between 0 and 1", {
  F_vals = plaplace(seq(-2, 2, length.out = 100))
  expect_true(all(F_vals >= 0 & F_vals <= 1))
})

test_that("plaplace() is monotonic increasing", {
  F_vals = plaplace(seq(-2, 2, length.out = 100))
  expect_true(all(F_vals[-1] - F_vals[-100] >= 0))
})

test_that("plaplace() is the integral of dlaplace", {
  p_low = runif(10, 0, 1)
  x_low = qlaplace(p_low)
  p_high = runif(10, p_low, 1)
  x_high = qlaplace(p_high)
  integrals = sapply(seq_along(p_low), function(i) integrate(dlaplace, x_low[i], x_high[i])$value)
  expect_equal(plaplace(x_high) - plaplace(x_low), integrals, tolerance = 1e-4)
})

test_that("qlaplace() gives the inverse of plaplace()", {
  p = runif(1000)
  expect_equal(p, plaplace(qlaplace(p)))
})

test_that("dpgd() integrates to pgpd()", {
  σ = runif(20)
  ξ = c(rep(0, 10), rnorm(10))
  u = rnorm(20)
  p1 = runif(20)
  p2 = runif(20, p1, 1)
  p = p2 - p1
  lower = mapply(qgpd, p = p1, u = u, σ = σ, ξ = ξ)
  upper = mapply(qgpd, p = p2, u = u, σ = σ, ξ = ξ)
  int = NULL
  for (i in 1:20) {
    int[i] = integrate(dgpd, lower = lower[i], upper = upper[i], σ = σ[i], ξ = ξ[i], u = u[i])$value
  }
  expect_equal(p, int, tolerance = 1e-3)
})

test_that("pgpd() is the inverse of qgpd()", {
  σ = runif(20)
  ξ = c(rep(0, 10), rnorm(10))
  p = runif(100)
  u = rnorm(100)
  p1 = t(matrix(rep(p, each = 20), ncol = length(p)))
  p2 = mapply(function(p, u, σ, ξ) pgpd(qgpd(p, u, σ, ξ), u, σ, ξ), σ = σ, ξ = ξ, MoreArgs = list(p = p, u = u))
  expect_equal(p1, p2, tolerance = 1e-5)
})

test_that("marginal_distribution() works with too large q", {
  x = runif(100)
  q = .9999
  expect_warning(marginal_distribution(x, q))
  F = suppressWarnings(marginal_distribution(x, q))
  F2 = ecdf(x)
  y = runif(1000)
  expect_equal(F(y), F2(y))
})

test_that("marginal_distribution() handles NA by ignoring it", {
  x = runif(1000)
  q = .8
  F1 = suppressWarnings(marginal_distribution(x, q))
  F2 = suppressWarnings(marginal_distribution(c(x, rep(NA, 300)), q))
  expect_equal(F1(x), F2(x))
})

test_that("marginal_distribution() is bounded between 0 and 1", {
  q = .8
  F = suppressWarnings(marginal_distribution(runif(1000), q))
  F_vals = F(seq(-2, 2, length.out = 100))
  expect_true(all(F_vals >= 0 & F_vals <= 1))
})

test_that("marginal_distribution() is monotonic increasing", {
  q = .8
  F = suppressWarnings(marginal_distribution(runif(1000), q))
  F_vals = F(seq(0, 1, length.out = 5000))
  expect_true(all(F_vals[-1] - F_vals[-100] >= 0))
})

test_that("inverse_distribution() gives the inverse of F when q is too large", {
  x = runif(1000)
  q = 1
  F = suppressWarnings(marginal_distribution(x, q))
  Q = inverse_distribution(F)
  expect_equal(x, Q(F(x)))
  p = runif(100)
  expect_equal(p, F(Q(p)), tolerance = 1e-2)
})

test_that("inverse_distribution() gives the inverse of F", {
  q = .8
  x = runif(1000)
  F = suppressWarnings(marginal_distribution(x, q))
  Q = inverse_distribution(F)
  expect_equal(x, Q(F(x)))
  p = runif(100)
  expect_equal(p, F(Q(p)), tolerance = 1e-2)
})

test_that("marginal_distribution() handles NA by returning NA", {
  q = .8
  x = runif(1000)
  F = suppressWarnings(marginal_distribution(x, q))
  F1 = F(x)
  F2 = F(c(x, rep(NA, 300)))
  expect_equal(F2, c(F1, rep(NA, 300)))
})

test_that("inverse_distribution() handles NA by returning NA", {
  q = .8
  x = runif(1000)
  F = suppressWarnings(marginal_distribution(x, q))
  Q = inverse_distribution(F)
  p = runif(1000)
  Q1 = Q(p)
  Q2 = Q(c(p, rep(NA, 300)))
  expect_equal(Q2, c(Q1, rep(NA, 300)))
})

test_that("plaplace() handles NA by returning NA", {
  x = rnorm(1000)
  F1 = plaplace(x)
  F2 = plaplace(c(x, rep(NA, 300)))
  expect_equal(F2, c(F1, rep(NA, 300)))
})

test_that("qlaplace() handles NA by returning NA", {
  p = runif(1000)
  Q1 = qlaplace(p)
  Q2 = qlaplace(c(p, rep(NA, 300)))
  expect_equal(Q2, c(Q1, rep(NA, 300)))
})
