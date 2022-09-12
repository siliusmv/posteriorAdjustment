devtools::load_all()
library(INLA)
library(ggplot2)
library(dplyr)
library(Matrix)
library(sf)
library(pbapply)
library(parallel)

# Compile and link all the cgeneric scripts, if this has not already been done
make_cgeneric("all")

n_repl = 300 # Number of replications of the experiments
n_train = 10 # Number of replications of the conditional extremes model in the training data
n_thin = 50 # Thinning variable used as input for the Wadsworth sampling algorithm
n_cores = 15 # Run code in parallel (this might be too many cores for you)
n_posterior_samples = 1e5 # Number of samples from posterior for estimating credible interval
threshold = qlaplace(.95) # The threshold t for defining the conditional extremes model
probs = c(.9, .95, .99) # Credible interval probabilities
filename = file.path(results_dir(), "conditional-adjustment.rds") # Where to store results
theta_star_filename = file.path(results_dir(), "conditional-theta-star.rds") # Where to store θ*
case_study_filename = file.path(results_dir(), "final-modelling.rds") # Results from the case study

if (!file.exists(filename)) saveRDS(list(), filename)

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

# Create functions for a(y0, d; θ) and b(y0, d; θ)
get_a_func = function(lambda, kappa) {
  function(y0, dist_to_s0) {
    alpha = exp(- (dist_to_s0 / lambda)^kappa)
    matrix(rep(y0, each = length(alpha)) * rep(alpha, length(y0)),
           nrow = length(dist_to_s0),
           ncol = length(y0))
  }
}
get_b_func = function(rho_b) {
  function(y0, dist_to_s0) {
    tmp = dist_to_s0 / rho_b
    tmp[tmp < 1e-9] = 1e-9
    b = sqrt(1 - exp(-2 * tmp))
    matrix(rep(b, length(y0)), nrow = length(dist_to_s0), ncol = length(y0))
  }
}
a_func = get_a_func(lambda, kappa)
b_func = get_b_func(rho_b)

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
  res = loglik_conditional(
    y = y,
    y0 = y0,
    a_func = get_a_func(lambda, kappa),
    b_func = get_b_func(rho_b),
    sigma = cov_mat,
    tau = tau,
    dist_to_s0 = dist_to_s0,
    dist_to_s0_from_mesh = dist_to_s0_from_mesh,
    A = A,
    n_cores = n_cores)
  if (sum_terms) res = sum(res)
  res
}

# Compute the gradient of the log-likelihood using numerical derivation
loglik_grad = function(theta, sum_terms = TRUE, ...) {
  res = numDeriv::jacobian(loglik, theta, ..., sum_terms = FALSE)
  if (sum_terms) res = apply(res, 2, sum)
  res
}

# ==============================================================================
# Choose all conditioning sites for the global model
# ==============================================================================

# Create a grid of conditioning sites with resolution (delta_s0 x delta_s0)
delta_s0 = 4
x_coords = coords[, 1] |> unique() |> sort()
y_coords = coords[, 2] |> unique() |> sort()
s0_locs = expand.grid(
  x = x_coords[seq(delta_s0, length(x_coords), by = delta_s0)],
  y = y_coords[seq(delta_s0, length(y_coords), by = delta_s0)]) |>
  as.matrix()

# Find out which indices the chosen conditioning sites correspond to
s0_index = lapply(
  X = 1:nrow(s0_locs),
  FUN = function(i) which(coords[, 1] == s0_locs[i, 1] & coords[, 2] == s0_locs[i, 2]))
s0_index = s0_index[sapply(s0_index, length) > 0] |>
  unlist() |>
  unname()

# ==============================================================================
# Simulate a large amount of samples from the distribution of [y | max(y) > t],
# using the algorithm of Keef et al. (2013).
# This is done in order to estimate the KLD minimiser θ*
# ==============================================================================

A_mats = lapply(1:n_loc, \(i) A[-i, ])
dist_to_s0 = lapply(1:n_loc, \(i) dist_euclid(coords[-i, ], coords[i, , drop = FALSE]))
dist_to_s0_from_mesh = lapply(1:n_loc, \(i) dist_euclid(mesh$loc[, 1:2], coords[i, , drop = FALSE]))

# Simulate all the observations.
# We split the simulations into multiple steps, to avoid memory problems
obs = list()
pb = progress_bar(50)
set.seed(1)
for (i in seq_len(50)) {
  obs[[i]] = wadsworth_sampling(
    n = 1e3,
    n_thin = n_thin,
    a_func = a_func,
    b_func = b_func,
    Q = Q,
    tau = tau,
    threshold = threshold,
    dist_to_s0 = dist_to_s0,
    dist_to_s0_from_mesh = dist_to_s0_from_mesh,
    A = A_mats,
    replace = FALSE)
  pb$tick()
}
obs = do.call(cbind, obs)
pb$terminate()

# Extract observations in a nonregular grid around conditioning sites that
# exceed the threshold
thinning = c(1, 2, 4, 6, 8, 16, 32)
data = extract_extreme_fields(
  data = t(obs),
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n = thinning,
  r = cumsum(4 * thinning))

# Free up some memory
rm(obs)

# Compute the distances to the conditioning sites from the mesh nodes,
# and add the information to `data`
data$dist_to_s0_from_mesh = list()
for (i in seq_along(data$s0)) {
  data$dist_to_s0_from_mesh[[i]] = dist_euclid(mesh$loc[, 1:2], data$s0[[i]])
}

# Compute the maximum likelihood estimator, which is approximately equal to θ*.
# Use θ as initial values
est = list(par = log(c(lambda, kappa, rho_b, rho, sigma, tau)), convergence = 1)
while (est$convergence != 0) {
  est = optim(
    par = est$par,
    fn = loglik,
    y = data$y,
    y0 = data$y0,
    dist_to_s0 = data$dist_to_s0,
    dist_to_s0_from_mesh = data$dist_to_s0_from_mesh,
    A = lapply(data$obs_index, function(x) A[x, ]),
    n_cores = n_cores,
    control = list(fnscale = -1, maxit = 300, trace = 6, REPORT = 1))
}
theta_star = est$par

# Save the estimated value of θ*
saveRDS(theta_star, theta_star_filename)

# ===============================================================================
# Sample smaller amounts of data and fit the conditional extremes model
# many times, in order to estimate coverage percentages
# ===============================================================================

# Save the estimated value of θ*
theta_star = readRDS(theta_star_filename)

# Run all the n_repl experiments
cl = parallel::makeForkCluster(n_cores, mc.set.seed = TRUE)
parallel::clusterSetRNGStream(cl, 1)
pbapply::pblapply(
  X = 1:n_repl,
  cl = cl,
  FUN = function(i) {

    # Store the state of the random seed before we start simulating random data
    seed = .Random.seed

    # Simulate data using the Wadsworth et al. (2019) algorithm
    obs = wadsworth_sampling(
      n = n_train,
      n_thin = n_thin,
      a_func = a_func,
      b_func = b_func,
      Q = Q,
      tau = tau,
      threshold = threshold,
      dist_to_s0 = dist_to_s0,
      dist_to_s0_from_mesh = dist_to_s0_from_mesh,
      A = A_mats)

    # Extract observations in a nonregular grid around conditioning sites that
    # exceed the threshold
    thinning = c(1, 2, 4, 6, 8, 16, 32)
    data = extract_extreme_fields(
      data = t(obs),
      coords = coords,
      s0_index = s0_index,
      threshold = threshold,
      n = thinning,
      r = cumsum(4 * thinning))

    # Compute distances from the mesh to the conditioning sites
    data$dist_to_s0_from_mesh = list()
    for (j in seq_along(data$s0)) {
      data$dist_to_s0_from_mesh[[j]] = dist_euclid(mesh$loc[, 1:2], data$s0[[j]])
    }

    # Estimate the maximum likelihood estimator to get good initial values for the
    # model fitting. We only allow maxit = 120 because this gives initial values
    # that are good enough, whitout spending too much time locating the actual maximum
    est = optim(
      par = theta_star,
      fn = loglik,
      y = data$y,
      y0 = data$y0,
      dist_to_s0 = data$dist_to_s0,
      dist_to_s0_from_mesh = data$dist_to_s0_from_mesh,
      A = lapply(data$obs_index, function(x) A[x, ]),
      n_cores = 1,
      control = list(fnscale = -1, maxit = 120))

    # R-INLA requires the observations y to be on a vector format
    y_inla = unlist(data$y)

    # Our implementation of the cgeneric model for a requires y0 and dist_to_s0 as input.
    # However, it requires one value of y0 and dist_to_s0 for each of the observations y.
    # We therefore replicate y0 and dist_to_s0 in the correct way so we get to vectors with
    # equal length to y_inla, where y_inla[i] was observed when y(s0) was equal to
    # y0_inla[i], and the distance between them was dist_to_s0_inla[i].
    dist_to_s0_inla = data$dist_to_s0[rep(seq_along(data$n), data$n)]
    y0_inla = rep(unlist(data$y0), sapply(dist_to_s0_inla, length))
    dist_to_s0_inla = unlist(dist_to_s0_inla)

    # Define the cgeneric model for a
    a_priors = list(lambda = c(1, 3, 3), kappa = c(1, -.4, 3))
    a_model = a_generic_model(
      y0 = y0_inla,
      dist_to_s0 = dist_to_s0_inla,
      init = theta_star[1:2],
      priors = a_priors)

    # Define the cgeneric model for Z_b
    spde_priors = list(
      rho = c(1, 60, .95),
      sigma = c(1, 5, .05),
      rho_b = c(1, log(6), 2))
    spde_model = spde_generic_model_with_b_func(
      spde = spde,
      n = data$n,
      init = theta_star[c(4, 5, 3)],
      priors = spde_priors,
      dist_to_s0 = do.call(rbind, data$dist_to_s0_from_mesh))

    # Create the necessary objects for running R-INLA with the conditional extremes model
    formula = y ~ -1 +
      f(spatial, model = spde_model) +
      f(idx, model = a_model)

    effects = list(
      spatial = seq_len(sum(data$n) * mesh$n),
      idx = seq_along(y_inla))

    A_inla = local({
      all_location_indices = rep(data$obs_index, data$n)
      all_repl_indices = rep(seq_along(all_location_indices), sapply(all_location_indices, length))
      all_location_indices = unlist(all_location_indices)
      A = inla.spde.make.A(
        mesh,
        loc = coords,
        index = all_location_indices,
        repl = all_repl_indices)
      A
    })

    # Build the stack
    stack = inla.stack(
      data = list(y = y_inla),
      A = list(spatial = A_inla, 1),
      effects = effects)

    # Fit the model
    fit = tryCatch({
      inla(
        formula = formula,
        data = inla.stack.data(stack),
        control.predictor = list(A = inla.stack.A(stack)),
        only.hyperparam = TRUE,
        control.mode = list(theta = est$par[c(6, 4:5, 3, 1:2)], restart = TRUE),
        control.inla = list(control.vb = list(enable = FALSE)),
        #verbose = TRUE,
        num.threads = 1,
        inla.mode = "experimental")
    }, error = function(e) NULL)
    if (is.null(fit)) next

    # Which way should the parameters of the model fits be reordered to
    # give the correct input to the log-likelihood functions?
    theta_reordering = c(5, 6, 4, 2, 3, 1)

    # Estimate H
    H = tryCatch(solve(fit$misc$cov.intern), error = \(e) NULL)
    if (is.null(H)) {
      warning("Inverting cov.intern failed for i = ", i)
      next
    }
    H = H[theta_reordering, theta_reordering]

    # Compute all terms of the gradient of the log-likelihood, so we can estimate J
    grads = loglik_grad(
      theta = fit$mode$theta[theta_reordering],
      y = data$y,
      y0 = data$y0,
      dist_to_s0 = data$dist_to_s0,
      dist_to_s0_from_mesh = data$dist_to_s0_from_mesh,
      A = lapply(data$obs_index, function(x) A[x, ]),
      n_cores = 1,
      sum_terms = FALSE)

    # Estimate J. We know that observations from two different times are completely independent
    times = unlist(data$time_index)
    J = 0
    for (k in seq_along(times)) {
      index = which(times == times[k])
      for (j in index) {
        J = J + grads[k, ] %*% grads[j, , drop = FALSE]
      }
    }

    # Compute the estimate for C
    C = get_C(H, J)

    # Keep only the relevant parts of the model fit, to reduce requirements on file storage
    fit = list(misc = fit$misc, internal.marginals.hyperpar = fit$internal.marginals.hyperpar,
               mode = list(theta = fit$mode$theta))
    fit$misc$reordering = NULL
    fit$misc$family = NULL
    fit$misc$linkfunctions = NULL
    class(fit) = "inla"

    # Create a list containing all the data of interest
    model_fits = list(
      fit = fit,
      C = C,
      n_terms = nrow(grads),
      seed = seed,
      theta_reordering = theta_reordering,
      time_index = data$time_index,
      n = sum(data$n))

    # Try to save the data. This may fail, since we are running several simulations in parallel.
    # If it fails, then wait a little while and try again
    while (TRUE) {
      success = tryCatch({
        tmp = readRDS(filename)
        tmp[[length(tmp) + 1]] = model_fits
        saveRDS(tmp, filename)
        rm(tmp)
        TRUE
      }, error = function(e) FALSE)
      if (success) break

      # Wait a random time between 1 and 10 seconds if we failed to save the data
      Sys.sleep(runif(1, 1, 10))
    }
  })
parallel::stopCluster(cl)

# ====================================================================
# Load the data and evaluate coverage properties
# ====================================================================
res = readRDS(filename)
theta_star = readRDS(theta_star_filename)
theta_reordering = res[[1]]$theta_reordering
n_theta = length(theta_reordering)

# Parameter names used for creating a latex-friendly table of the results
theta_names = c("log_lambda", "log_kappa", "log_rho_b", "log_rho", "log_sigma", "log_precision")
theta_tex_names = c("lambda", "kappa", "rho_b", "rho", "sigma", "tau")

# If we had any bad runs, remove them
bad_index = which(sapply(res, length) == 0)
if (any(bad_index)) res = res[-bad_index]

# Estimate credible intervals
pb = progress_bar(length(res))
intervals = list()
for (i in seq_along(res)) {
  intervals[[i]] = list()
  theta_uncorrected = inla.hyperpar.sample(n_posterior_samples, res[[i]]$fit, intern = TRUE)
  theta_uncorrected = theta_uncorrected[, theta_reordering]
  for (k in 1:2) {
    if (k == 1) {
      label = "Unadjusted"
      theta_corrected = theta_uncorrected
    } else {
      label = "Adjusted"
      theta_corrected = matrix(
        data = rep(res[[i]]$fit$mode$theta[theta_reordering], each = n_posterior_samples),
        ncol = n_theta)
      theta_corrected = theta_corrected + (theta_uncorrected - theta_corrected) %*% t(res[[i]]$C)
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
    tmp1[[i]] = ifelse(tmp1[[i]] < theta_star[i], TRUE, FALSE)
    tmp2[[i]] = ifelse(tmp2[[i]] > theta_star[i], TRUE, FALSE)
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
      labels = paste0("$\\", theta_tex_names)),
    label = factor(
      x = label,
      levels = c("Unadjusted", "Adjusted"),
      labels = c("$", "_\\text{adj}$"))) |>
  tidyr::pivot_wider(names_from = c(name, label), values_from = coverage)
names(tmp) = sub("_\\$", "$", names(tmp))
names(tmp) = sub("_\\\\", "\\\\", names(tmp))
names(tmp) = sub("rho_b_", "{rho_b}_", names(tmp))
tmp = tmp[, c(1, 9, 3, 8, 2, 12, 6, 11, 5, 13, 7, 10, 4)]
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
