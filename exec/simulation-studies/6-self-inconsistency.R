library(pracma)
library(ggplot2)
library(tidyr)
library(tidyselect)

devtools::load_all()

# ===============================================================
# Define the two bivariate sampling algorithms for the conditional
# extremes model where
# x = (x_1, x_2), and both x_1 and x_2 have exponential marginals,
# a(x) = mu + alpha * x,
# b(x) = sigma * x ^ beta,
# and the residual distribution is standardised Gaussian
# ===============================================================

# Sample from the global conditional extremes distribution of [(x, y) | max(x, t) > t],
# using the simulation algorithm of Keef et al. (2013).
#
# The input variables are:
# n: How many realisations of [x in R^2 | max(x) > t] should we sample?
# alpha: A parameter in the function a(x)
# beta: A parameter in the function b(x)
# mu: A parameter in the function a(x)
# sigma: A parameter in the function b(x)
# t: The threshold used for defining the conditional extremes model
keef_sampling_bivariate = function(n, alpha, beta, mu, sigma, t) {

  # Preallocate the output
  res = matrix(NA_real_, nrow = n, ncol = 2)

  # Draw which location should be used as conditioning site for each of the
  # n samples
  s0_index = sample.int(2, n, replace = TRUE)

  # A boolean vector. Element nr. i describes if we are finished yet at tuple nr. i,
  # i.e. if the conditioning site is larger than the other random variable
  not_finished = rep(TRUE, n)

  # Keep simulating until all conditioning sites have the max
  while (any(not_finished)) {

    # first and second are indices for where we need to sample more using x_1 or x_2
    # as conditioning sites, respectively. With x = (x_1, x_2)
    first = (s0_index == 1) & not_finished
    second = (s0_index == 2) & not_finished

    # Simulate tuples with x_1 as the conditioning variable
    if (any(first)) {
      res[first, 1] = rexp(sum(first)) + t
      res[first, 2] = rnorm(
        n = sum(first),
        mean = mu + alpha * res[first, 1],
        sd = sigma * res[first, 1] ^ beta)
    }

    # Simulate tuples with x_2 as the conditioning variable
    if (any(second)) {
      res[second, 2] = rexp(sum(second)) + t
      res[second, 1] = rnorm(
        n = sum(second),
        mean = mu + alpha * res[second, 2],
        sd = sigma * res[second, 2] ^ beta)
    }

    # Keep track of which iterations the conditioning site is still not the maximum
    not_finished = (first & res[, 2] > res[, 1]) | (second & res[, 1] > res[, 2])
  }

  # Return the result
  res
}

# Sample from the global conditional extremes distribution of [(x, y) | max(x, t) > t],
# using the simulation algorithm of Wadsworth et al. (2019).
#
# The input variables are:
# n: How many realisations of [x in R^2 | max(x) > t] should we sample?
# alpha: A parameter in the function a(x)
# beta: A parameter in the function b(x)
# mu: A parameter in the function a(x)
# sigma: A parameter in the function b(x)
# t: The threshold used for defining the conditional extremes model
# n_thin: The algorithm is based on an importance sampling algorithm. First, we sample
#   n_1 random variables, then all n_1 random variables get different weights, and
#   we select n < n_1 of the original samples using weighted sampling without
#   replacement. We set n_1 = n_thin * n
wadsworth_sampling_bivariate = function(n, alpha, beta, mu, sigma, t, n_thin = 20) {

  # Preallocate a matrix for all n_1 = n * n_thin random tuples
  res = matrix(NA_real_, nrow = n_thin * n, ncol = 2)

  # Draw which location should be used as conditioning site for each of the
  # n samples
  s0_index = sample(c(1, 2), n_thin * n, replace = TRUE)

  # first and second are indices for where we use x_1 or x_2 as
  # conditioning sites, respectively
  first = (s0_index == 1)
  second = (s0_index == 2)

  # Simulate tuples with x_1 as the conditioning variable
  res[first, 1] = rexp(sum(first)) + t
  res[first, 2] = rnorm(
    n = sum(first),
    mean = mu + alpha * res[first, 1],
    sd = sigma * res[first, 1] ^ beta)

  # Simulate tuples with x_2 as the conditioning variable
  res[second, 2] = rexp(sum(second)) + t
  res[second, 1] = rnorm(
    n = sum(second),
    mean = mu + alpha * res[second, 2],
    sd = sigma * res[second, 2] ^ beta)

  # Compute the importance weights
  importance_weights = 1 / ((res[, 1] > t) + (res[, 2] > t))

  # Draw the n indices we are going to return
  importance_sampling_index = sample.int(n_thin * n, n, prob = importance_weights)

  # Return the n selected tuples
  res[importance_sampling_index, ]
}

# ===============================================================
# show that the two integrals in the appendix differ from each other,
# but that they agree with the two different simulation algorithms.
# Both integrals compute P(x_1 > t, x_2 > t | max(x) > t)
# ===============================================================

# Parameters of the conditional extremes model
alpha = .9
beta = .8
sigma = 1
mu = 0
t = 3 # Threshold for the conditional extremes model

n_samples = 1e4 # How many tuples should we sample?
n_thin = 20 # Thinning parameter used for the Wadsworth sampling
# Tecnicality: We cannot integrate from t to infinity, so we instead integrate from
# t to u
u = 50

# Function used for integrating over the distribution of [x | max(x) > t]
f = function(y0, y) dnorm(y, mu + alpha * y0, sigma * y0^beta) * dexp(y0 - t)

# Compute P(t < x_1 < u, t < x_2 < u | max(x) > t) using both integrals in the appendix
keef_prob = local({
  a = pracma::integral2(
  fun = f,
  xmin = t,
  xmax = u,
  ymin = t,
  ymax = function(x) pmin(x, u))$Q
  b = pracma::integral2(
  fun = f,
  xmin = t,
  xmax = u,
  ymin = -20,
  ymax = function(x) pmin(x, u))$Q
  a / b
})
wadsworth_prob = local({
  a = pracma::integral2(
    fun = f,
    xmin = t,
    xmax = u,
    ymin = t,
    ymax = u)$Q
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


# Estimate P(t < x_1 < u, t < x_2 < u | max(x) > t) by sampling from the two
# simulation algorithms

# First, simulate some data
x_keef = keef_sampling_bivariate(n_samples, alpha, beta, mu, sigma, t)
x_wadsworth = wadsworth_sampling_bivariate(n_samples, alpha, beta, mu, sigma, t, n_thin)

# Then, compute the Monte Carlo estimates
keef_prob_mc = mean(
  t < x_keef[, 1] & x_keef[, 1] < u & t < x_keef[, 2] & x_keef[, 2] < u)
wadsworth_prob_mc = mean(
  t < x_wadsworth[, 1] & x_wadsworth[, 1] < u & t < x_wadsworth[, 2] & x_wadsworth[, 2] < u)

# Compare the integrals with the Monte Carlo estimates
c(keef_prob_mc, wadsworth_prob_mc)
c(keef_prob, wadsworth_prob)

# To ensure that we are correct, repeat the Monte Carlo estimation several times
# to evaluate uncertainties
n_repl = 100
pb = progress_bar(n_repl)
probs = list()
set.seed(1)
for (i in 1:n_repl) {
  x_keef = keef_sampling_bivariate(n_samples, alpha, beta, mu, sigma, t)
  x_wadsworth = wadsworth_sampling_bivariate(n_samples, alpha, beta, mu, sigma, t, n_thin)
  keef_prob_mc = mean(
    t < x_keef[, 1] & x_keef[, 1] < u & t < x_keef[, 2] & x_keef[, 2] < u)
  wadsworth_prob_mc = mean(
    t < x_wadsworth[, 1] & x_wadsworth[, 1] < u & t < x_wadsworth[, 2] & x_wadsworth[, 2] < u)
  probs[[i]] = c(keef_prob_mc, wadsworth_prob_mc)
  pb$tick()
}
probs = as.data.frame(do.call(rbind, probs))
names(probs) = c("Keef", "Wadsworth")
pb$terminate()

# Create histograms of the Monte Carlo estimates together with the numerical
# solutions of the integrals
probs |>
  tidyr::pivot_longer(tidyselect::everything()) |>
  ggplot() +
  geom_histogram(aes(x = value, y = ..density..)) +
  geom_vline(
    data = data.frame(name = c("Keef", "Wadsworth"), x = c(keef_prob, wadsworth_prob)),
    aes(xintercept = x)) +
  facet_wrap(~name, scales = "free", nrow = 1)

# =================================================================
# Show how this affects θ*, which becomes different from θ
# =================================================================

# Compute the log-likelihood for the global model of [(x_1, x_2) | max(x_1, x_2) > t]
loglik = function(theta, x, t) {
  alpha = theta[1]
  beta = theta[2]
  if (alpha < 0 || alpha > 1 || beta < 0 || beta > 1) return(-1e199)
  index = which(x[, 1] > t)
  res = sum(dnorm(x[index, 2], alpha * x[index, 1], x[index, 1] ^ beta, log = TRUE))
  res
}


# Simulate large amounts of data using the two simulation algorithms, then
# estimate the KLD minimiser θ* and show that it differs considerably
# between the two simulation algorithms, and also that it differs considerably
# from the parameters θ = (α, β), that we use for simulating the data
n_repl = 100
thetahat = list()
set.seed(1)
pb = progress_bar(n_repl)
for (i in 1:n_repl) {

  # Simulate from the Keef algorithm, then estimate θ*
  x = keef_sampling_bivariate(n_samples, alpha, beta, mu, sigma, t)
  keef_pars = optim(
    par = c(alpha, beta),
    fn = loglik,
    x = x,
    t = t,
    control = list(fnscale = -1))$par

  # Simulate from the Wadsworth algorithm, then estimate θ*
  x = wadsworth_sampling_bivariate(n_samples, alpha, beta, mu, sigma, t)
  wadsworth_pars = optim(
    par = c(alpha, beta),
    fn = loglik,
    x = x,
    t = t,
    control = list(fnscale = -1))$par

  # Save the results in a data.frame
  thetahat[[i]] = data.frame(
    value = c(keef_pars, wadsworth_pars),
    par = c("alpha", "beta", "alpha", "beta"),
    name = c("Keef", "Keef", "Wadsworth", "Wadsworth"))

  pb$tick()
}
thetahat = do.call(rbind, thetahat)
pb$terminate()

# Plot the estimated densities of the MLEs using the two different
# simulation algorithms, and compare them with the original θ values
ggplot(thetahat) +
  geom_density(aes(x = value)) +
  geom_vline(
    data = data.frame(par = c("alpha", "beta"), value = c(alpha, beta)),
    mapping = aes(xintercept = value)) +
  facet_grid(name ~ par, scales = "free")
