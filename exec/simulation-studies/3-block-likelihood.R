devtools::load_all()
library(INLA)
library(mvtnorm)
library(ggplot2)

domain_size = 25 # Spatial domain is quadratic with lengths equal to domain_size
block_size = 5 # Size of the quadratic blocks used in the block likelihood
n_blocks = (domain_size / block_size)^2 # Number of blocks
n_loc = 400 # Number of locations containing data
n_train = 100 # Number of replications of the Gaussian field in the training data
n_theta_star = 1e5 # Number of replications of the Gaussian field for estimating theta*
n_cores = 7 # Run code in parallel
n_repl = 300 # Number of replications of the experiment
n_posterior_samples = 1e5 # Number of samples from posterior for estimating credible interval
probs = c(.9, .95, .99) # Credible interval probabilities
filename = file.path(results_dir(), "block-likelihood.rds") # Filename for results
sigma = 1 # Standard deviance of the Gaussian random field
rho = 12 # Range of the Gaussian random field
tau = 100 # Precision of the nugget effect

# Draw all the locations randomly
set.seed(1)
loc = matrix(runif(n_loc * 2) * domain_size, ncol = 2)

# Create the mesh for both the true and the composite distributions.
# This is a really coarse mesh, chosen in order to speed up inference as much as possible
mesh = inla.mesh.2d(
  loc.domain = loc,
  boundary = list(inla.nonconvex.hull(loc, convex = -.3),
                  inla.nonconvex.hull(loc, convex = -1)),
  cutoff = 6,
  max.edge = c(9, 50))
#plot(mesh); points(loc)

# Define the SPDE and the projection matrix for the mesh
spde = inla.spde2.pcmatern(
  mesh,
  prior.range = c(rho, 0.5),
  prior.sigma = c(sigma, 0.5))
A = inla.spde.make.A(mesh, loc)

# Compute the precision matrix of the true distribution
Q = inla.spde2.precision(spde, theta = c(log(rho), log(sigma)))

# Log-likelihood for the block composite likelihood.
# block-index is a list of vectors of indices, representing different blocks,
# where each block is considered independent of each other
block_loglik = function(theta,
                        y,
                        block_index,
                        spde,
                        A,
                        return_mat = FALSE) {
  tau = exp(theta[1])
  rho = exp(theta[2])
  sigma = exp(theta[3])
  Q = INLA::inla.spde2.precision(spde, theta = c(log(rho), log(sigma)))
  cov_mat = A %*% Matrix::solve(Q, Matrix::t(A))
  cov_mat = as.matrix(cov_mat) + diag(1 / tau, nrow(cov_mat))
  res = lapply(
    X = block_index,
    FUN = function(x) matrix(NA_real_, nrow = length(x), ncol = ncol(y)))
  for (i in seq_along(block_index)) {
    if (any(block_index[[i]])) {
      res[[i]] = mvtnorm::dmvnorm(
        t(y[block_index[[i]], ]),
        sigma = cov_mat[block_index[[i]], block_index[[i]]],
        log = TRUE)
    }
  }
  res = do.call(rbind, res)
  if (!return_mat) res = sum(res)
  res
}

# Create the block indices for the block composite likelihood
block_index = list()
for (j in 1:(domain_size / block_size)) {
  for (i in 1:(domain_size / block_size)) {
    y_lims = c(i - 1, i) * block_size
    x_lims = c(j - 1, j) * block_size
    block_index[[length(block_index) + 1]] = which(
      x_lims[1] < loc[, 1] & loc[, 1] <= x_lims[2] &
      y_lims[1] < loc[, 2] & loc[, 2] <= y_lims[2])
  }
}

# Run all the n_repl experiments
set.seed(1)
model_fits = list()
attributes(model_fits)$block_index = block_index
attributes(model_fits)$theta_hat = log(c(tau, rho, sigma))
pb = progress_bar(n_repl)
for (i in 1:n_repl) {

  # Simulate data from the true distribution
  z = rnorm_spde(n_train, Q)
  y = as.matrix(A %*% z)
  y = y + rnorm(length(y), sd = tau^-.5)

  # Create the necessary objects for running R-INLA
  effects = list(
    spatial = inla.spde.make.index("spatial", n.spde = mesh$n, n.repl = n_train * n_blocks))
  y_inla = as.numeric(y[unlist(block_index), ])
  A_inla = list(spatial = inla.spde.make.A(
    mesh = mesh,
    loc = loc,
    index = rep(unlist(block_index), n_train),
    repl = rep(1:(n_blocks * n_train), rep(sapply(block_index, length), n_train))))
  stack = inla.stack(data = list(y = y_inla), A = A_inla, effects = effects)
  formula = y ~ -1 + f(spatial, model = spde, replicate = spatial.repl, nrep = n_train * n_blocks)

  # Fit the model
  fit = inla(
    formula = formula,
    data = inla.stack.data(stack),
    family = "gaussian",
    num.threads = n_cores,
    only.hyperparam = TRUE,
    control.inla = list(control.vb = list(enable = FALSE)),
    inla.mode = "experimental",
    control.predictor = list(A = inla.stack.A(stack)))

  # Estimate θ*
  theta = fit$mode$theta

  # Estimate H, J and C
  H = solve(fit$misc$cov.intern)
  grads = numDeriv::jacobian(
    function(x) block_loglik(x, y, block_index, spde, A, return_mat = TRUE),
    x = theta)
  J = 0
  for (j in 1:n_train) {
    for (k in 1:n_blocks) {
      for (l in 1:n_blocks) {
        J = J + grads[(j - 1) * n_blocks + k, ] %*%
          grads[(j - 1) * n_blocks + l, , drop = FALSE]
      }
    }
  }
  C = get_C(H, J)

  # Keep only the relevant parts of the model fit, to reduce requirements on file storage
  fit = list(
    misc = fit$misc,
    internal.marginals.hyperpar = fit$internal.marginals.hyperpar,
    model = list(theta = fit$mode$theta))
  fit$misc$reordering = NULL
  fit$misc$family = NULL
  fit$misc$linkfunctions = NULL
  class(fit) = "inla"

  # Save the temporary results
  model_fits[[i]] = list(fit = fit, C = C)
  saveRDS(model_fits, filename)
  pb$tick()
}
pb$terminate()

# Load the results
res = readRDS(filename)
block_index = attributes(res)$block_index
theta_hat = attributes(res)$theta_hat
n_theta = length(theta_hat)
theta_names = c("log_precision", "log_rho", "log_sigma")

# Estimate credible intervals
pb = progress_bar(length(res))
intervals = list()
for (i in seq_along(res)) {
  intervals[[i]] = list()
  theta_uncorrected = inla.hyperpar.sample(n_posterior_samples, res[[i]]$fit, intern = TRUE)
  for (k in 1:2) {
    if (k == 1) {
      theta_corrected = theta_uncorrected
      label = "Unadjusted"
    } else {
      theta_corrected = matrix(rep(res[[i]]$fit$mode$theta, each = n_posterior_samples), ncol = n_theta)
      theta_corrected = theta_corrected + (theta_uncorrected - theta_corrected) %*% t(res[[i]]$C)
      label = "Adjusted"
    }
    intervals[[i]][[length(intervals[[i]]) + 1]] = rbind(
      as.data.frame(apply(theta_corrected, 2, quantile, probs = (1 - probs) / 2)) |>
        dplyr::mutate(label = label, prob = probs, lim = "lower"),
      as.data.frame(apply(theta_corrected, 2, quantile, probs = 1 - (1 - probs) / 2)) |>
        dplyr::mutate(label = label, prob = probs, lim = "upper"))
    row.names(intervals[[i]][[length(intervals[[i]])]]) = NULL
    names(intervals[[i]][[length(intervals[[i]])]])[1:n_theta] = theta_names
  }
  intervals[[i]] = do.call(rbind, intervals[[i]]) |>
    dplyr::mutate(iter = i)
  pb$tick()
}
intervals = do.call(rbind, intervals)
pb$terminate()

# Find if θ* is inside the credible intervals or not
is_inside_interval = local({
  tmp1 = dplyr::filter(intervals, lim == "lower")
  tmp2 = dplyr::filter(intervals, lim == "upper")
  for (i in seq_along(theta_names)) {
    tmp1[[i]] = ifelse(tmp1[[i]] < theta_hat[i], TRUE, FALSE)
    tmp2[[i]] = ifelse(tmp2[[i]] > theta_hat[i], TRUE, FALSE)
  }
  tmp = (tmp1[, 1:n_theta] & tmp2[, 1:n_theta]) |>
    as.data.frame() |>
    dplyr::mutate(label = tmp1$label, prob = tmp1$prob, iter = tmp1$iter)
  tmp
})

# Compute coverage percentages
tmp = is_inside_interval |>
  tidyr::pivot_longer(all_of(theta_names)) |>
  dplyr::group_by(prob, label, name) |>
  dplyr::summarise(coverage = mean(value)) |>
  dplyr::mutate(
    name = factor(
      x = name,
      levels = theta_names,
      labels = paste0("$\\", c("tau", "rho", "sigma"))),
    label = factor(
      x = label,
      levels = c("Unadjusted", "Adjusted"),
      labels = c("$", "_\\text{adj}$"))) |>
  tidyr::pivot_wider(names_from = c(name, label), values_from = coverage)
names(tmp) = sub("_", "", names(tmp))
tmp = tmp[, c(1, 5, 2, 6, 3, 7, 4)]
print(tmp)

# ==============================================
# Reformat the results into latex tabular format
# ==============================================
table = paste(paste(c("Aim", names(tmp)[-1]), collapse = " & "), "\\\\")
table[2] = "\\hline"
j = 3
for (i in 1:nrow(tmp)) {
  table[j] = paste(
    c(paste0("$", round(100 * tmp$prob[i], digits = 0), "\\%$"),
      paste0("$", round(100 * tmp[i, -1], digits = 0), "\\%$")),
    collapse = " & ")
  table[j] = paste(table[j], "\\\\")
  j = j + 1
}
table = paste(paste(table, collapse = "\n"), "\n")

cat(table)
