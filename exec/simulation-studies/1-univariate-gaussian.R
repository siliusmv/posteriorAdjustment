devtools::load_all()
library(INLA)
library(parallel)
library(numDeriv)
library(purrr)

sigma = 1 # Standard deviance of the misspecified Gaussian
n_train = 100 # Number of training data
n_repl = 1000 # Number of replications of the experiment
n_cores = 15 # Run code in parallel
n_posterior_samples = 1e5 # Number of samples from posterior for estimating credible interval
probs = c(.9, .95, .99) # Credible interval probabilities
filename = file.path(results_dir(), "univariate-gaussian.rds") # Filename for results

# Log-likelihood for the misspecified Gaussian
loglik = function(mu, y, sigma, return_vec = FALSE) {
  res = dnorm(y, mean = mu, sd = sigma, log = TRUE)
  if (!return_vec) res = sum(res)
  res
}

# Run all the n_repl experiments
set.seed(1, kind = "L'Ecuyer-CMRG")
quantiles = parallel::mclapply(
  X = seq_len(n_repl),
  mc.cores = n_cores,
  FUN = function(i) {

    # Simulate data from the true distribution
    y = rt(n_train, 1)

    # Fit the misspecified distribution
    res = inla(
      formula = y ~ 1,
      data = data.frame(y = y),
      control.compute = list(config = TRUE),
      control.fixed = list(prec.intercept = 0.001),
      control.family = list(hyper = list(prec = list(
        initial = log(sigma^-2),
        fixed = TRUE))),
      control.inla = list(int.strategy = "grid", dz = .05, diff.logdens = 10),
      num.threads = 1)

    # Estimate θ*
    mode = inla.mmarginal(res$marginals.fixed[[1]])

    # Estimate J, H and C
    H = - numDeriv::hessian(loglik, x = mode, y = y, sigma = sigma)
    grads = numDeriv::jacobian(loglik, x = mode, y = y, sigma = sigma, return_vec = TRUE)
    J = sum(grads^2)
    C = get_C(H, J)

    # Sample from the posterior
    samples = inla.posterior.sample(
      n = n_posterior_samples,
      result = res,
      selection = list("(Intercept)" = 1),
      seed = 1L)
    samples = sapply(samples, `[[`, "latent")

    # Adjust the samples to improve frequency properties
    samples = list(
      unadjusted = samples,
      adjusted = mode + (samples - mode) * as.numeric(C))

    # Estimate credible intervals
    quantiles = list()
    for (name in names(samples)) {
      quantiles[[name]] = quantile(
        samples[[name]],
        probs = c((1 - probs) / 2, 1 - (1 - probs) / 2))
    }

    if (i %% 10 == 0) message("Progress: ", i, " / ", n_repl)
    quantiles
  })

# Restructure and save the results
quantiles = purrr::transpose(quantiles)
for (name in names(quantiles)) {
  quantiles[[name]] = do.call(rbind, quantiles[[name]])
}
saveRDS(quantiles, filename)

# Load the results
quantiles = readRDS(filename)

# Compute coverage percentages
coverage = list()
for (name in names(quantiles)) {
  zero_above_lower = quantiles[[name]][, seq_along(probs)] < 0
  zero_below_upper = quantiles[[name]][, length(probs) + seq_along(probs)] > 0
  stopifnot(all(dim(zero_above_lower) == dim(zero_below_upper)))
  inside_interval = zero_below_upper & zero_above_lower
  coverage[[name]] = apply(inside_interval, 2, mean)
  names(coverage[[name]]) = paste(100 * probs, "%")
}
coverage

# ==============================================
# Reformat the results into latex tabular format
# ==============================================
tmp = coverage |>
  do.call(what = cbind) |>
  (\(x) cbind(probs, x))() |>
  (\(x) round(x * 100, digits = 0))()
colnames(tmp) = c("Aim", "Unadjusted", "Adjusted")
rownames(tmp) = NULL

table = paste(paste(colnames(tmp), collapse = " & "), "\\\\")
table[2] = "\\hline"
for (i in 1:nrow(tmp)) {
  table[[length(table) + 1]] = paste(paste0("$", tmp[i, ], "\\%$", collapse = " & "), "\\\\")
}
table = paste(paste(table, collapse = "\n"), "\n")

cat(table)
