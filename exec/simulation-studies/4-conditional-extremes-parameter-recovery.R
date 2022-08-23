devtools::load_all()
library(INLA)
library(ggplot2)
library(sf)
library(dplyr)
library(tidyr)

# Compile and link all the cgeneric scripts, if this has not already been done
make_cgeneric("all")

n_repl = 100 # Number of replications of the experiments
n_train = 100 # Number of replications of the conditional extremes model in the training data
n_cores = 15 # Run code in parallel (this might be too many cores for you)
threshold = qlaplace(.95) # The threshold t for defining the conditional extremes model
filename = file.path(results_dir(), "parameter-recovery.rds") # Where to store results
case_study_filename = file.path(results_dir(), "final-modelling.rds") # Results from the case study

# ==============================================================================
# Load the case study data to create a similar setting in the simulation
# as for the case study
# ==============================================================================
radar = readRDS(file.path(downloads_dir(), "radar.rds"))
coords = st_coordinates(radar$coords)
n_loc = nrow(coords)
rm(radar)

# Create the mesh and the SPDE
mesh = inla.mesh.2d(
  loc = coords,
  boundary = list(inla.nonconvex.hull(coords, convex = -.1),
                  inla.nonconvex.hull(coords, convex = -1)),
  max.edge = c(5, 30))
spde = inla.spde2.pcmatern(mesh, prior.range = c(60, .95), prior.sigma = c(4, .05))

# Create the projection matrix for the mesh and SPDE
A = inla.spde.make.A(mesh, coords)

# Choose one conditioning site to use for the experiments
s0 = matrix(c(281, 7092), nrow = 1)
s0_index = which(coords[, 1] == s0[1] & coords[, 2] == s0[2])

# Compute distances to s0 from all the locations and from the mesh
dist_to_s0 = dist_euclid(coords, s0)
dist_to_s0_from_mesh = dist_euclid(mesh$loc[, 1:2], s0)

# Get all the parameter values from the global model fits in the case study
tmp = readRDS(case_study_filename)
ss = tmp$fit[[1]]$summary.hyperpar
tau = ss[1, 1]
rho = exp(ss[2, 1])
sigma = exp(ss[3, 1])
lambda = exp(ss[5, 1])
kappa = exp(ss[6, 1])
rho_b = exp(tmp$log_rho_b)
rm(tmp, ss)

# Create functions for a() and b()
a_func = function(y0, dist_to_s0, lambda, kappa) {
  y0 * exp(- (dist_to_s0 / lambda) ^ kappa)
}
b_func = function(dist_to_s0, rho_b) {
  tmp = dist_to_s0 / rho_b
  tmp[tmp < 1e-9] = 1e-9
  sqrt(1 - exp(-2 * tmp))
}

# Create the precision matrix of the SPDE approximation
Q = INLA::inla.spde2.precision(spde, log(c(rho, sigma)))

# Define the composite log-likelihood for the global conditional
# extremes model. This is necessary for computing the MLE.
# The function below is a wrapper for the function loglik_conditional
#
# Input variables:
# theta: A vector of parameters of length 5 or 6, depending on the value of rho_b.
# y, y0, dist_to_s0, dist_to_s0_from_mesh, A, n_cores: See ?loglik_conditional for more info
# rho_b: Either NULL or a double. If is.null(rho_b), then theta has length 6, and
#   we try to estimate rho_b. Else, theta has length 5, and rho_b is fixed and equal
#   to the given value.
# sum_terms: A boolean describing wheter we should compute the sum of the log-likelihood,
#   or if we should remove one value for each threshold exceedance.
loglik = function(theta,
                  y,
                  y0,
                  dist_to_s0,
                  dist_to_s0_from_mesh,
                  A,
                  rho_b = NULL,
                  sum_terms = TRUE,
                  n_cores = 1) {
  if (is.null(rho_b)) {
    stopifnot(length(theta) == 6)
    rho_b = exp(theta[3])
    log_rho = theta[4]
    log_sigma = theta[5]
    tau = exp(theta[6])
  } else {
    stopifnot(length(theta) == 5)
    log_rho = theta[3]
    log_sigma = theta[4]
    tau = exp(theta[5])
  }
  lambda = exp(theta[1])
  kappa = exp(theta[2])
  Q = INLA::inla.spde2.precision(spde, c(log_rho, log_sigma))
  cov_mat = as.matrix(Matrix::solve(Q))
  a_func = function(y, dist) {
    alpha = exp(- (dist / lambda)^kappa)
    matrix(rep(y, each = length(alpha)) * rep(alpha, length(y)),
           nrow = length(dist), ncol = length(y))
  }
  b_func = function(y, dist) {
    tmp = dist / rho_b
    tmp[tmp < 1e-9] = 1e-9
    b = sqrt(1 - exp(-2 * tmp))
    matrix(rep(b, length(y)), nrow = length(dist), ncol = length(y))
  }
  res = loglik_conditional(
    y = y,
    y0 = y0,
    a_func = a_func,
    b_func = b_func,
    sigma = cov_mat,
    tau = tau,
    dist_to_s0 = dist_to_s0,
    dist_to_s0_from_mesh = dist_to_s0_from_mesh,
    A = A,
    n_cores = n_cores)
  if (sum_terms) res = sum(res)
  res
}

# =======================================================================
# Simulate from the model y = a(y0, |s - s0|) + b(y0, |s - s0|) * Z + ε,
# and recover parameters using R-INLA with the cgeneric models
# =======================================================================

cl = parallel::makeForkCluster(n_cores, mc.set.seed = TRUE)
parallel::clusterSetRNGStream(cl, 1)
params = pbapply::pblapply(
  X = 1:n_repl,
  cl = cl,
  FUN = function(i) {

    # Sample y0
    y0 = threshold + rexp(n_train)

    # Create y
    z = rnorm_spde(n = n_train, Q = Q)
    y = as.matrix(A %*% z)
    for (j in seq_along(y0)) {
      y[, j] = y[, j] * b_func(dist_to_s0, rho_b) +
        a_func(y0[j], dist_to_s0, lambda, kappa)
    }
    y = y + rnorm(length(y), sd = tau^-.5)

    # Compute the maximum likelihood estimator. We use good initial values to
    # speed things up a bit
    est = optim(
      par = log(c(lambda, kappa, rho_b, rho, sigma, tau)),
      fn = loglik,
      y = list(y),
      y0 = list(y0),
      dist_to_s0 = list(dist_to_s0),
      dist_to_s0_from_mesh = list(dist_to_s0_from_mesh),
      A = list(A),
      control = list(fnscale = -1, maxit = 300))

    # R-INLA requires the observations y to be on a vector format
    y_inla = as.numeric(y)

    # Our implementation of the cgeneric model for a requires y0 and dist_to_s0 as input.
    # However, it requires one value of y0 and dist_to_s0 for each of the observations y.
    # We therefore replicate y0 and dist_to_s0 in the correct way so we get to vectors with
    # equal length to y_inla, where y_inla[i] was observed when y(s0) was equal to
    # y0_inla[i], and the distance between them was dist_to_s0_inla[i].
    dist_to_s0_inla = rep(dist_to_s0, n_train)
    y0_inla = rep(y0, each = n_loc)

    # Define the cgeneric model for a
    a_priors = list(lambda = c(1, 3, 3), kappa = c(1, -.4, 3))
    a_model = a_generic_model(
      y0 = y0_inla,
      dist_to_s0 = dist_to_s0_inla,
      init = est$par[1:2],
      priors = a_priors)

    # Define the cgeneric model for Z_b
    spde_priors = list(
      rho = c(1, 60, .95),
      sigma = c(1, 5, .05),
      rho_b = c(1, log(6), 2))
    spde_model = spde_generic_model_with_b_func(
      spde = spde,
      n = n_train,
      init = est$par[c(4, 5, 3)],
      priors = spde_priors,
      dist_to_s0 = dist_to_s0_from_mesh)

    # Create the remaining necessary objects for running R-INLA with the conditional extremes model
    formula = y ~ -1 +
      f(idx, model = a_model) +
      f(spatial, model = spde_model)

    effects = list(
      spatial = seq_len(n_train * mesh$n),
      idx = seq_along(y_inla))

    A_inla = inla.spde.make.A(
      mesh = mesh,
      loc = coords,
      index = rep(1:n_loc, n_train),
      repl = rep(1:n_train, each = n_loc))

    stack = inla.stack(
      data = list(y = y_inla),
      A = list(spatial = A_inla, 1),
      effects = effects)

    # Estimate parameters
    fit = tryCatch({
      inla(
        formula,
        data = inla.stack.data(stack),
        control.predictor = list(A = inla.stack.A(stack)),
        num.threads = 1,
        only.hyperparam = TRUE,
        control.mode = list(theta = est$par[c(6, 1:2, 4:5, 3)], restart = TRUE),
        #verbose = TRUE,
        silent = 2,
        control.inla = list(control.vb = list(enable = FALSE)),
        inla.mode = "experimental")},
      error = function(e) NULL)
    if (is.null(fit)) {
      message("Error occured for iter ", i)
      return(list(y = y_inla, y0 = y0_inla))
    }

    # Sample parameters from the posterior and compute summary statistics
    hyper = inla.hyperpar.sample(1e5, fit) |>
      (\(x) {x[is.infinite(x)] = NA; x})() |>
      apply(2, \(x) c(mean(x, na.rm = TRUE), quantile(x, c(.025, .975), na.rm = TRUE)))

    # Create a data.frame with relevant information
    res = data.frame(
      value = hyper[1, ] |> (\(x) c(1 / x[1], exp(x[2:3]), exp(x[4:6])))(),
      lower = hyper[2, ] |> (\(x) c(1 / x[1], exp(x[2:3]), exp(x[4:6])))(),
      upper = hyper[3, ] |> (\(x) c(1 / x[1], exp(x[2:3]), exp(x[4:6])))(),
      name = c("nugget", "lambda", "kappa", "rho", "sigma", "rho_b"),
      truth = c(1 / tau, lambda, kappa, rho, sigma, rho_b),
      cpu = fit$cpu[4],
      mlik = fit$mlik[1],
      i = i)
    res[1, 2:3] = res[1, 3:2]
    row.names(res) = NULL

    res
  })
parallel::stopCluster(cl)

# Save the results
saveRDS(params, filename)

# =============================================================
# Read the data and evaluate the results
# =============================================================
params = readRDS(filename)

# If we had any bad runs, remove them
bad_runs = which(sapply(params, length) == 2)
if (any(bad_runs)) params = params[-bad_runs]

# Examine MSE, MAE, coverage percentage for the 95% credible interval,
# and other summary statistics
params |>
  do.call(what = rbind) |>
  dplyr::mutate(err = truth - value, rel_err = err / truth) |>
  dplyr::mutate(is_included = truth > lower & truth < upper) |>
  dplyr::group_by(name) |>
  dplyr::summarise(mean = mean(err), mean_rel = mean(rel_err),
                   sd = sd(err), max = max(abs(err)),
                   mse = mean(err^2), mae = mean(abs(err)),
                   coverage = mean(is_included), mlik = mean(mlik))

# Plot the densities of the posterior means of θ, together with the true values
params |>
  do.call(what = rbind) |>
  dplyr::mutate(
    name = factor(
      name,
      levels = c("nugget", "rho", "sigma", "rho_b", "lambda", "kappa"),
      labels = c("nugget", "rho", "sigma", "rho_b", "lambda", "kappa"))) |>
  ggplot() +
  geom_density(aes(x = value)) +
  facet_wrap(~name, scales = "free") +
  geom_vline(aes(xintercept = truth))
