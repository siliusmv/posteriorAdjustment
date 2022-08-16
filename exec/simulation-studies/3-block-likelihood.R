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
# The block index list for the full likelihood is the list with only one index vector,
# which contains all n_loc indices
full_index = list(1:n_loc)

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

  # Create the necessary objects for running R-INLA with both the
  # full and the composite likelihoods
  effects_blocks = list(
    spatial = inla.spde.make.index("spatial", n.spde = mesh$n, n.repl = n_train * n_blocks))
  effects_full = list(
    spatial = inla.spde.make.index("spatial", n.spde = mesh$n, n.repl = n_train))
  y_inla_blocks = as.numeric(y[unlist(block_index), ])
  y_inla_full = as.numeric(y[unlist(full_index), ])
  A_blocks = list(spatial = inla.spde.make.A(
    mesh = mesh,
    loc = loc,
    index = rep(unlist(block_index), n_train),
    repl = rep(1:(n_blocks * n_train), rep(sapply(block_index, length), n_train))))
  A_full = list(spatial = inla.spde.make.A(
    mesh = mesh,
    loc = loc,
    index = rep(unlist(full_index), n_train),
    repl = rep(1:n_train, rep(sapply(full_index, length), n_train))))
  stack_blocks = inla.stack(data = list(y = y_inla_blocks), A = A_blocks, effects = effects_blocks)
  stack_full = inla.stack(data = list(y = y_inla_full), A = A_full, effects = effects_full)
  formula_blocks = y ~ -1 + f(spatial, model = spde, replicate = spatial.repl, nrep = n_train * n_blocks)
  formula_full = y ~ -1 + f(spatial, model = spde, replicate = spatial.repl, nrep = n_train)

  # Fit the two models
  block_fit = inla(
    formula = formula_blocks,
    data = inla.stack.data(stack_blocks),
    family = "gaussian",
    num.threads = n_cores,
    only.hyperparam = TRUE,
    control.inla = list(control.vb = list(enable = FALSE)),
    inla.mode = "experimental",
    control.predictor = list(A = inla.stack.A(stack_blocks)))
  full_fit = inla(
    formula = formula_full,
    data = inla.stack.data(stack_full),
    family = "gaussian",
    num.threads = n_cores,
    only.hyperparam = TRUE,
    control.inla = list(control.vb = list(enable = FALSE)),
    inla.mode = "experimental",
    control.predictor = list(A = inla.stack.A(stack_full)))

  # Estimate θ* for both models, using only the training data
  theta_block = block_fit$mode$theta
  theta_full = full_fit$mode$theta

  # Estimate H, J and C
  H_block = solve(block_fit$misc$cov.intern)
  H_full = solve(full_fit$misc$cov.intern)
  grads_block = numDeriv::jacobian(
    function(x) block_loglik(x, y, block_index, spde, A, return_mat = TRUE),
    x = theta_block)
  grads_full = numDeriv::jacobian(
    function(x) block_loglik(x, y, full_index, spde, A, return_mat = TRUE),
    x = theta_full)
  J_block = 0
  J_full = 0
  for (j in 1:n_train) {
    for (k in 1:n_blocks) {
      for (l in 1:n_blocks) {
        J_block = J_block + grads_block[(j - 1) * n_blocks + k, ] %*%
          grads_block[(j - 1) * n_blocks + l, , drop = FALSE]
      }
    }
    J_full = J_full + grads_full[j, ] %*% grads_full[j, , drop = FALSE]
  }
  C_block = get_C(H_block, J_block)
  C_full = get_C(H_full, J_full)

  # Keep only the relevant parts of the model fits, to reduce requirements on file storage
  block_fit = list(
    misc = block_fit$misc,
    internal.marginals.hyperpar = block_fit$internal.marginals.hyperpar,
    model = list(theta = block_fit$mode$theta))
  block_fit$misc$reordering = NULL
  block_fit$misc$family = NULL
  block_fit$misc$linkfunctions = NULL
  class(block_fit) = "inla"

  full_fit = list(
    misc = full_fit$misc,
    internal.marginals.hyperpar = full_fit$internal.marginals.hyperpar,
    model = list(theta = full_fit$mode$theta))
  full_fit$misc$reordering = NULL
  full_fit$misc$family = NULL
  full_fit$misc$linkfunctions = NULL
  class(full_fit) = "inla"

  # Save the temporary results
  model_fits[[i]] = list(
    block_fit = block_fit,
    full_fit = full_fit,
    C_block = C_block,
    C_full = C_full)
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
  for (j in 1:2) {
    if (j == 1) {
      fit = res[[i]]$block_fit
      C = res[[i]]$C_block
      model = "Block"
    } else {
      fit = res[[i]]$full_fit
      C = res[[i]]$C_full
      model = "Full"
    }
    theta_uncorrected = inla.hyperpar.sample(n_posterior_samples, fit, intern = TRUE)
    for (k in 1:2) {
      if (k == 1) {
        theta_corrected = theta_uncorrected
        label = "Unadjusted"
      } else {
        theta_corrected = matrix(rep(fit$mode$theta, each = n_posterior_samples), ncol = n_theta)
        theta_corrected = theta_corrected + (theta_uncorrected - theta_corrected) %*% t(C)
        label = "Adjusted"
      }
      intervals[[i]][[length(intervals[[i]]) + 1]] = rbind(
        as.data.frame(apply(theta_corrected, 2, quantile, probs = (1 - probs) / 2)) |>
          dplyr::mutate(model = model, label = label, prob = probs, lim = "lower"),
        as.data.frame(apply(theta_corrected, 2, quantile, probs = 1 - (1 - probs) / 2)) |>
          dplyr::mutate(model = model, label = label, prob = probs, lim = "upper"))
      row.names(intervals[[i]][[length(intervals[[i]])]]) = NULL
      names(intervals[[i]][[length(intervals[[i]])]])[1:n_theta] = theta_names
    }
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
    dplyr::mutate(model = tmp1$model, label = tmp1$label,
                  prob = tmp1$prob, iter = tmp1$iter)
  tmp
})

# Compute coverage percentages
tmp = is_inside_interval |>
  tidyr::pivot_longer(all_of(theta_names)) |>
  dplyr::group_by(prob, label, model, name) |>
  dplyr::summarise(coverage = mean(value)) |>
  dplyr::mutate(name = factor(
    x = name,
    levels = theta_names,
    labels = paste0("$\\", c("tau", "rho", "sigma"), "$"))) |>
  tidyr::pivot_wider(names_from = name, values_from = coverage)
tmp = tmp[1:nrow(tmp) + rep(c(2, 2, -2, -2), 3), ]
print(tmp)

# ==============================================
# Reformat the results into latex tabular format
# ==============================================
table = paste(paste(c("\\(1 - \\alpha\\)", "Method", names(tmp)[-(1:2)]), collapse = " & "), "\\\\")
table = paste(paste(
  c("", "", "\\multicolumn{3}{c}{Block likelihood}", "\\multicolumn{3}{c}{Full likelihood}"),
  collapse = " & "), "\\\\")
table[2] = "\\cmidrule(r){3-5} \\cmidrule(l){6-8}"
table[3] = paste(paste(c("\\(1 - \\alpha\\)", "Method", rep(names(tmp)[-(1:3)], 2)), collapse = " & "), "\\\\")
table[4] = "\\midrule"
j = 5
for (i in seq(1, nrow(tmp) - 1, by = 2)) {
  table[j] = paste(
    c(tmp$label[i],
      paste0("$", round(100 * tmp[i, -(1:3)], digits = 0), "\\%$"),
      paste0("$", round(100 * tmp[i + 1, -(1:3)], digits = 0), "\\%$")),
    collapse = " & ")
  if (i %% 4 == 1) {
    table[j] = paste(
      c(paste0("$", round(100 * tmp$prob[i], digits = 0), "\\%$"), table[j]), collapse = " & ")
  } else {
    table[j] = paste(c("", table[j]), collapse = " & ")
  }
  table[j] = paste(table[j], "\\\\")
  j = j + 1
  if (i %% 4 == 3 && i != (nrow(tmp) - 1)) {
    table[j] = "\\midrule"
    j = j + 1
  }
}
table = paste(paste(table, collapse = "\n"), "\n")

cat(table)
