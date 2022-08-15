library(pracma)
library(ggplot2)

# ===============================================================
# Define the two bivariate sampling algorithms for the conditional
# extremes model where
# x and y have exponential marginals,
# a(x) = μ + α * x,
# b(x) = σ * x ^ β,
# and the residual distribution is standardised Gaussian
# ===============================================================

keef_sampling_bivariate = function(n, α, β, μ, σ, t) {
  s0_index = sample.int(2, n, replace = TRUE)
  res = matrix(NA_real_, nrow = n, ncol = 2)
  not_finished = rep(TRUE, n)
  while (any(not_finished)) {
    first = (s0_index == 1) & not_finished
    second = (s0_index == 2) & not_finished
    if (any(first)) {
      res[first, 1] = rexp(sum(first)) + t
      res[first, 2] = rnorm(sum(first), μ + α * res[first, 1], σ * res[first, 1] ^ β)
    }
    if (any(second)) {
      res[second, 2] = rexp(sum(second)) + t
      res[second, 1] = rnorm(sum(second), μ + α * res[second, 2], σ * res[second, 2] ^ β)
    }
    not_finished = (first & res[, 2] > res[, 1]) | (second & res[, 1] > res[, 2])
  }
  res
}

wadsworth_sampling_bivariate = function(n, α, β, μ, σ, t, n_thin = 20) {
  s0_index = sample(c(1, 2), n_thin * n, replace = TRUE)
  res = matrix(NA_real_, nrow = n_thin * n, ncol = 2)
  first = (s0_index == 1)
  second = (s0_index == 2)
  res[first, 1] = rexp(sum(first)) + t
  res[first, 2] = rnorm(sum(first), μ + α * res[first, 1], σ * res[first, 1] ^ β)
  res[second, 2] = rexp(sum(second)) + t
  res[second, 1] = rnorm(sum(second), μ + α * res[second, 2], σ * res[second, 2] ^ β)
  importance_weights = 1 / ((res[, 1] > t) + (res[, 2] > t))
  importance_sampling_index = sample.int(n_thin * n, n, prob = importance_weights)
  res[importance_sampling_index, ]
}

# ===============================================================
# show that the two integrals in the appendix differ,
# but that they agree with the two different simulation algorithms
# ===============================================================

α = .9
β = .8
σ = 1
μ = 0
t = 4
upper = 50
n_samples = 1e4
n_thin = 20

f = function(y0, y) dnorm(y, μ + α * y0, σ * y0^β) * dexp(y0 - t)
keef_prob = local({
  a = pracma::integral2(
  fun = f,
  xmin = t,
  xmax = upper,
  ymin = t,
  ymax = function(x) pmin(x, upper))$Q
  b = pracma::integral2(
  fun = f,
  xmin = t,
  xmax = upper,
  ymin = -20,
  ymax = function(x) pmin(x, upper))$Q
  a / b
})
wadsworth_prob = local({
  a = pracma::integral2(
    fun = f,
    xmin = t,
    xmax = upper,
    ymin = t,
    ymax = upper)$Q
  b = pracma::integral2(
    fun = f,
    xmin = t,
    xmax = 100,
    ymin = -20,
    ymax = t)$Q
  c = pracma::integral2(
    fun = f,
    xmin = t,
    xmax = 100,
    ymin = t,
    ymax = 100)$Q
  a / (2 * b + c)
})
c(keef_prob, wadsworth_prob)

y_keef = keef_sampling_bivariate(n_samples, α, β, μ, σ, t)
y_wadsworth = wadsworth_sampling_bivariate(n_samples, α, β, μ, σ, t, n_thin)

keef_prob_mc = mean(t < y_keef[, 1] & y_keef[, 1] < upper & t < y_keef[, 2] & y_keef[, 2] < upper)
wadsworth_prob_mc = mean(t < y_wadsworth[, 1] & y_wadsworth[, 1] < upper & t < y_wadsworth[, 2] & y_wadsworth[, 2] < upper)

c(keef_prob_mc, wadsworth_prob_mc)
c(keef_prob, wadsworth_prob)

# =====================================================
# Show how this affects theta^*, which becomes different from theta
# =====================================================

ll = function(θ, y, t) {
  α = θ[1]
  β = θ[2]
  if (α < 0 || α > 1 || β < 0 || β > 1) return(-1e199)
  index = which(y[, 1] > t)
  res = sum(dnorm(y[index, 2], α * y[index, 1], y[index, 1] ^ β, log = TRUE))
  res
}

n_repl = 100

set.seed(1)
thetahat = list()
for (i in 1:n_repl) {
  y = keef_sampling_bivariate(n_samples, α, β, μ, σ, t)
  keef_pars = optim(
    par = c(α, β),
    fn = ll,
    y = y,
    t = t,
    control = list(fnscale = -1))$par
  y = wadsworth_sampling_bivariate(n_samples, α, β, μ, σ, t)
  wadsworth_pars = optim(
    par = c(α, β),
    fn = ll,
    y = y,
    t = t,
    control = list(fnscale = -1))$par
  thetahat[[i]] = data.frame(
    value = c(keef_pars, wadsworth_pars),
    par = c("α", "β", "α", "β"),
    name = c("Keef", "Keef", "Wadsworth", "Wadsworth"))
  if (i %% 10 == 0) message("Progress: ", i, " / ", n_repl)
}
thetahat = do.call(rbind, thetahat)

ggplot(thetahat) +
  geom_density(aes(x = value)) +
  geom_vline(
    data = data.frame(par = c("α", "β"), value = c(α, β)),
    mapping = aes(xintercept = value)) +
  facet_grid(name ~ par, scales = "free")
