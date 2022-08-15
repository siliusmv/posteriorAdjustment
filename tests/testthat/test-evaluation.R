
# # This is probably useless...
# 
# test_that("twcrps_norm is equal to scoringRuls::crps_norm with y0 = -Inf", {
#   n = 100
#   μ = rnorm(n)
#   σ = exp(rnorm(n))
#   y = rnorm(n, μ, σ)
#   S1 = scoringRules::crps_norm(y, μ, σ)
#   S2 = twcrps_norm(y, μ, σ, y0 = -Inf)
#   expect_equal(S1, S2)
# })
# 
# test_that("twcrps_norm is equal to what we get by numerical integration", {
#   n = 100
#   μ = rnorm(n)
#   σ = exp(rnorm(n))
#   y = rnorm(n, μ, σ)
#   y0 = y + c(runif(n / 2), -runif(n / 2))
#   S = twcrps_norm(y, μ, σ, y0)
#   f = function(z, y, μ, σ) (pnorm(z, μ, σ) - ifelse(y <= z, 1, 0)) ^ 2
#   int = function(y, μ, σ, y0) integrate(f, lower = y0, upper = qnorm(.99999, μ, σ), μ = μ, σ = σ, y = y)$value
#   integrals = rep(0, n)
#   for (i in seq_len(n)) integrals[i] = int(y[i], μ[i], σ[i], y0[i])
#   expect_equal(S, integrals, tolerance = 1e-4)
# })
