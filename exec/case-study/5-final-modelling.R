devtools::load_all()
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(INLA)
library(inlabru)

# Compile and link all the cgeneric scripts, if this has not already been done
make_cgeneric("all")

# ==============================================================================
# Decide necessary model parameters:
# ==============================================================================
threshold = qlaplace(.9975) # The threshold t for defining the conditional extremes model
n_cores = 15 # Run code in parallel. Check yourself that this is an ok number
r = 5 # Radius used for computing aggregated empirical distribution functions

filename = file.path(results_dir(), "final-modelling.rds")

# ==============================================================================
# Load the radar data
# ==============================================================================
radar = readRDS(file.path(downloads_dir(), "radar.rds"))
coords = st_coordinates(radar$coords)
n_loc = nrow(coords)

# ==============================================================================
# Transform the data to Laplace margins
# ==============================================================================
cl = parallel::makeForkCluster(n_cores)
transformed_data = pbapply::pblapply(
  X = 1:n_loc,
  cl = cl,
  FUN = function(i) {
    F = aggregated_ecdf(radar$data, coords, coords[i, ], radius = r)
    u = F(radar$data[, i])
    # Ensure that we don't get infinities
    if (any(u == 1, na.rm = TRUE)) {
      u[which(u == 1)] = (1 + max(u[which(u != 1)])) / 2
    }
    qlaplace(u)
  })
parallel::stopCluster(cl)
transformed_data = do.call(cbind, transformed_data)

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
# Extract observations in a nonregular grid around conditioning sites that
# exceed the threshold
# ==============================================================================
thinning = c(1, 2, 4, 6, 8, 16, 32)
data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n = thinning,
  r = cumsum(4 * thinning))

# Number of threshold exceedances
sum(data$n)
# Summary of the number of observations used for inference per conditioning site
data$obs_index |> sapply(length) |> summary()
# Summary of the number of conditioning site that simultaneously exceed the threshold
# for all the unique times a threshold exceedance occurs.
data$time_index |> unlist() |> table() |> as.numeric() |> summary()

# ==============================================================================
# Create the mesh and spde
# ==============================================================================
mesh = inla.mesh.2d(
  loc = coords,
  boundary = list(inla.nonconvex.hull(coords, convex = -.1),
                  inla.nonconvex.hull(coords, convex = -1)),
  max.edge = c(5, 30))
mesh$n
# plot(mesh)
# points(coords)

spde = inla.spde2.pcmatern(mesh, prior.range = c(60, .95), prior.sigma = c(5, .05))

# Compute the distances to the conditioning sites from the mesh nodes,
# and add the information to `data`
data$dist_to_s0_from_mesh = list()
for (i in seq_along(data$s0)) {
  data$dist_to_s0_from_mesh[[i]] = dist_euclid(mesh$loc[, 1:2], data$s0[[i]])
}

# Plot an example of the locations used for performing inference for a given
# conditioning site, together with a display of the triangular mesh
plot = local({
  i = 12
  df = as.data.frame(coords[data$obs_index[[i]], ]) |>
    dplyr::mutate(tag = "Used")
  df2 = as.data.frame(coords[-c(data$obs_index[[i]], data$s0_index[[i]]), ]) |>
    dplyr::mutate(tag = "Not used")
  s0_df = as.data.frame(coords[data$s0_index[[i]], , drop = FALSE]) |>
    dplyr::mutate(tag = "$\\bm s_0$")
  rbind(df, df2, s0_df) |>
    dplyr::mutate(tag = factor(tag, levels = c("$\\bm s_0$", "Used", "Not used"))) |>
    ggplot() +
    inlabru::gg(mesh, int.color = "grey", ext.color = "grey") +
    geom_point(aes(x = X, y = Y, shape = tag, size = tag, color = tag)) +
    scale_color_manual(values = c("red", "black", "black")) +
    scale_shape_manual(values = c(17, 19, 19)) +
    scale_size_manual(values = c(3, 2.4, .5)) +
    theme_light() +
    guides(shape = "none", size = "none", color = "none") +
    labs(x = "Easting", y = "Northing") +
    theme(axis.ticks = element_blank(), axis.text = element_blank(),
          text = element_text(size = 23),
          panel.grid = element_blank())
})
tikz_plot(file.path(image_dir(), "design-of-experiment.pdf"),
          plot, width = 12, height = 12)

# ==============================================================================
# Compute the MLE and use it as initial values for R-INLA
# ==============================================================================

# First, define the composite log-likelihood for the global conditional
# extremes model. This is necessary for computing the MLE, but also for
# estimating J(θ*) and computing log-scores later in the code.
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

# Compute the maximum likelihood estimator. We repeat the optimisation process
# multiple times in case it does not manage to converge withing the first
# maxit = 500 number of iterations
est = list(par = c(3.0, -.4, -2, 2.6, .6, 3.1), convergence = 1)
while (est$convergence != 0) {
  est = optim(
    par = est$par,
    fn = loglik,
    y = data$y,
    y0 = data$y0,
    dist_to_s0 = data$dist_to_s0,
    dist_to_s0_from_mesh = data$dist_to_s0_from_mesh,
    A = lapply(data$obs_index, function(x) inla.spde.make.A(mesh, coords[x, ])),
    n_cores = n_cores,
    control = list(fnscale = -1, maxit = 300, trace = 6, REPORT = 1))
}

# ==============================================================================
# Create all the model components for fitting the global conditional extremes
# model to data with R-INLA
# ==============================================================================

# R-INLA requires the observations y to be on a vector format
y_inla = unlist(data$y)
# Keep track of all indices where y_inla is NA, so we can remove these
# before feeding the data into R-INLA, for reduced requirements on CPU and RAM.
na_index = which(is.na(y_inla))

# Our implementation of the cgeneric model for a requires y0 and dist_to_s0 as input.
# However, it requires one value of y0 and dist_to_s0 for each of the observations y.
# We therefore replicate y0 and dist_to_s0 in the correct way so we get to vectors with
# equal length to y_inla, where y_inla[i] was observed when y(s0) was equal to
# y0_inla[i], and the distance between them was dist_to_s0_inla[i].
dist_to_s0_inla = data$dist_to_s0[rep(seq_along(data$n), data$n)]
y0_inla = rep(unlist(data$y0), sapply(dist_to_s0_inla, length))
dist_to_s0_inla = unlist(dist_to_s0_inla)

# Remove data at all fo the NA indices
y_inla = y_inla[-na_index]
y0_inla = y0_inla[-na_index]
dist_to_s0_inla = dist_to_s0_inla[-na_index]

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
  n = data$n,
  init = est$par[c(4, 5, 3)],
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
  # Creation of the projection matrix A is a bit more complex.
  # First we must find the location index for each observation in y_inla.
  # Then we must find which time index (all_repl_indices) the observations come from.
  # Then we can finally build the matrix using inla.spde.make.A()
  all_location_indices = rep(data$obs_index, data$n)
  all_repl_indices = rep(seq_along(all_location_indices), sapply(all_location_indices, length))
  all_location_indices = unlist(all_location_indices)
  A = inla.spde.make.A(
    mesh,
    loc = coords,
    index = all_location_indices,
    repl = all_repl_indices)
  A = A[-na_index, ] # Remove the rows of A that correspond to NA observations
  A
})

# Build the stack
stack = inla.stack(
  data = list(y = y_inla),
  A = list(spatial = A_inla, 1),
  effects = effects)

# ==============================================================================
# Fit the model
# ==============================================================================
fit = inla(
  formula = formula,
  data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack)),
  control.mode = list(theta = est$par[c(6, 4:5, 3, 1:2)], restart = TRUE),
  only.hyperparam = TRUE,
  control.inla = list(control.vb = list(enable = FALSE)),
  verbose = TRUE,
  num.threads = 4,
  inla.mode = "experimental")

# Remove unnecessary objects from the environment and model fit that use a lot of RAM
fit$model.matrix = NULL
fit$offset.linear.predictor = NULL
rm(effects, A_inla)

# ==============================================================================
# Fit models using only one conditioning site
# ==============================================================================

# Select ten conditioning sites. The first site is equal to the conditioning
# site from the training data that has the most threshold exceedances. The other
# 9 conditioning sites are drawn randomly from both the training data and the
# test data
set.seed(123)
single_site_s0s = c(data$s0_index[which.max(data$n)], sample.int(nrow(coords), 9))

# Plot the locations of the randomly selected conditioning sites
local({
  coords = as.data.frame(coords)
  coords$tag = FALSE
  coords$tag[single_site_s0s] = TRUE
  ggplot(coords) +
    geom_point(aes(x = X, y = Y, col = tag))
})

# Extract observations in nonregular grids around the ten conditioning sites
single_site_data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = single_site_s0s,
  threshold = threshold,
  n = thinning,
  r = cumsum(4 * thinning))

# Compute dist_to_s0_from_mesh for all single sites
single_site_data$dist_to_s0_from_mesh = list()
for (i in seq_along(single_site_data$s0)) {
  single_site_data$dist_to_s0_from_mesh[[i]] = dist_euclid(
    mesh$loc[, 1:2],
    single_site_data$s0[[i]])
}

# Perform inference for all of the single site fits
single_site_fits = vector("list", length(single_site_s0s))
pb = progress_bar(length(single_site_s0s))
for (i in seq_along(single_site_s0s)) {

  # Repeat what we did for the global model fit, for preparing data on a format
  # that R-INLA likes
  y_mini = as.numeric(single_site_data$y[[i]])
  na_index = which(is.na(y_mini))
  mini_n = single_site_data$n[i]
  dist_to_s0_mini = rep(single_site_data$dist_to_s0[[i]], mini_n)
  y0_mini = rep(
    single_site_data$y0[[i]],
    each = length(single_site_data$dist_to_s0[[i]]))
  y_mini = y_mini[-na_index]
  y0_mini = y0_mini[-na_index]
  dist_to_s0_mini = dist_to_s0_mini[-na_index]

  # Create the model for a and Z_b
  a_model_mini = a_generic_model(
    y0 = y0_mini,
    dist_to_s0 = dist_to_s0_mini,
    init = est$par[1:2],
    priors = a_priors)
  spde_model_mini = spde_generic_model_with_b_func(
    spde = spde,
    n = mini_n,
    init = est$par[c(4, 5, 3)],
    priors = spde_priors,
    dist_to_s0 = matrix(single_site_data$dist_to_s0_from_mesh[[i]], nrow = 1))

  # Create the stack
  effects_mini = list(
    spatial = seq_len(sum(mini_n) * mesh$n),
    idx = seq_along(y_mini))
  A_inla_mini = local({
    all_location_indices = rep(single_site_data$obs_index[i], mini_n)
    all_repl_indices = rep(
      seq_along(all_location_indices),
      sapply(all_location_indices, length))
    all_location_indices = unlist(all_location_indices)
    A = inla.spde.make.A(
      mesh,
      loc = coords,
      index = all_location_indices,
      repl = all_repl_indices)
    A = A[-na_index, ]
    A
  })
  stack_mini = inla.stack(
    data = list(y = y_mini),
    A = list(spatial = A_inla_mini, 1),
    effects = effects_mini)

  formula_mini = y ~ -1 +
    f(spatial, model = spde_model_mini) +
    f(idx, model = a_model_mini)

  # Perform inference
  single_site_fits[[i]] = tryCatch({
    inla(
      formula = formula_mini,
      data = inla.stack.data(stack_mini),
      control.predictor = list(A = inla.stack.A(stack_mini)),
      control.mode = list(theta = est$par[c(6, 4:5, 3, 1:2)], restart = TRUE),
      only.hyperparam = TRUE,
      control.inla = list(control.vb = list(enable = FALSE)),
      #verbose = TRUE,
      silent = 2L,
      num.threads = 4,
      inla.mode = "experimental")
  }, error = function(e) NULL)
  single_site_fits[[i]]$model.matrix = NULL
  single_site_fits[[i]]$offset.linear.predictor = NULL

  pb$tick()
}
pb$terminate()

# ==============================================================================
# Estimate rho_b using an independence likelihood where
# Z_b is only Gaussian white noise
# ==============================================================================
iid_b_model = iid_generic_model_with_b_func(
  init = c(1, 1.5),
  priors = list(sigma = c(1, 5, .05), rho_b = c(1, 1.5, 3)),
  dist_to_s0 = dist_to_s0_inla)

# Fit the iid model
iid_fit = inla(
  formula = y ~ -1 + f(idx, model = a_model) + f(idx2, model = iid_b_model),
  data = data.frame(y = y_inla, idx = seq_along(y_inla), idx2 = seq_along(y_inla)),
  only.hyperparam = TRUE,
  control.inla = list(control.vb = list(enable = FALSE), int.strategy = "eb"),
  verbose = TRUE,
  num.threads = 4,
  inla.mode = "experimental")

# Extract the estimate for ρ_b, then delete the iid_fit to decrease memory requirements
log_rho_b = iid_fit$mode$theta[5]
rm(iid_fit)

# ==============================================================================
# Fit the global model using the fixed value of ρ_b
# ==============================================================================

# First, compute the maximum likelihood estimator
est2 = list(par = c(4, -.4, 4, .5, 4), convergence = 1)
while (est2$convergence != 0) {
  est2 = optim(
    par = est2$par,
    fn = loglik,
    y = data$y,
    y0 = data$y0,
    dist_to_s0 = data$dist_to_s0,
    dist_to_s0_from_mesh = data$dist_to_s0_from_mesh,
    A = lapply(data$obs_index, function(x) inla.spde.make.A(mesh, coords[x, ])),
    rho_b = exp(log_rho_b),
    n_cores = n_cores,
    control = list(fnscale = -1, maxit = 300, trace = 6, REPORT = 1))
}

# Then, redefine the spde_model for Z_b, so ρ_b is fixed
spde_priors$rho_b = c(0, log_rho_b)
spde_model = spde_generic_model_with_b_func(
  spde = spde,
  n = data$n,
  init = est2$par[3:4],
  priors = spde_priors,
  dist_to_s0 = do.call(rbind, data$dist_to_s0_from_mesh))

# Fit the model to data
fit_fixed = inla(
  formula = formula,
  data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack)),
  control.mode = list(theta = est2$par[c(5, 3:4, 1:2)], restart = TRUE),
  only.hyperparam = TRUE,
  control.inla = list(control.vb = list(enable = FALSE)),
  verbose = TRUE,
  num.threads = 4,
  inla.mode = "experimental")
fit_fixed$model.matrix = NULL
fit_fixed$offset.linear.predictor = NULL

# =================================================================
# Create a list of all the model fits, together with some additional
# information, and save the list into the file `filename`
# =================================================================
fits = list(
  fit = c(list(fit), single_site_fits, list(fit_fixed)),
  global = c(TRUE, rep(FALSE, length(single_site_fits)), TRUE),
  fixed = c(rep(FALSE, length(single_site_fits) + 1), TRUE),
  threshold = threshold,
  log_rho_b = log_rho_b,
  single_site_s0s = single_site_s0s)
saveRDS(fits, filename)

# =============================================================
# Load the model fits and prepare to evaluate them
# using the out-of-sample log-score
# =============================================================
fits = readRDS(filename)
log_rho_b = fits$log_rho_b

# Create the test data, based on all conditioning sites
# that were not used for performing inference with the global models
eval_data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = seq_len(n_loc)[-s0_index],
  threshold = threshold,
  n = 1,
  r = Inf)

# Compute distances to s0 from the mesh
eval_data$dist_to_s0_from_mesh = list()
for (i in seq_along(eval_data$s0)) {
  eval_data$dist_to_s0_from_mesh[[i]] = dist_euclid(mesh$loc[, 1:2], eval_data$s0[[i]])
}

# ================================================================
# Add some necessary extra information to the list of model fits
# ================================================================

# Which way should the parameters of the model fits be reordered to
# give the correct input to the log-likelihood functions?
# reordering1 is for the fits with estimated ρ_b while
# reordering2 is for the fits with fixed ρ_b
theta_reordering1 = c(5, 6, 4, 2, 3, 1)
theta_reordering2 = c(4, 5, 2, 3, 1)

# What are the names of the parameters from each model fit?
# This is mainly used for plotting the results later
theta_names1 = c("log_lambda", "log_kappa", "log_rho_b", "log_rho", "log_sigma", "log_tau")
theta_tex_names1 = paste0("$\\", c("lambda", "kappa", "rho_b", "rho", "sigma", "tau"), "$")
theta_names2 = c("log_lambda", "log_kappa", "log_rho", "log_sigma", "log_tau")
theta_tex_names2 = paste0("$\\", c("lambda", "kappa", "rho", "sigma", "tau"), "$")

# Add all this information to the list of model fits
fits$reordering = fits$names = fits$tex_names = list()
for (i in seq_along(fits$fit)) {
  if (fits$fixed[i]) {
    fits$reordering[[i]] = theta_reordering2
    fits$names[[i]] = theta_names2
    fits$tex_names[[i]] = theta_tex_names2
  } else {
    fits$reordering[[i]] = theta_reordering1
    fits$names[[i]] = theta_names1
    fits$tex_names[[i]] = theta_tex_names1
  }
}

# ===============================================================
# Compute H, J and C
# ===============================================================

# Compute gradients of the log-likelihood using numerical derivation
loglik_grad = function(theta, ...) {
  numDeriv::jacobian(loglik, theta, ...)
}

# Create the necessary projection matrices to compute the
# gradients of the log-likelihoods in-sample
A = inla.spde.make.A(mesh, coords)
A_mats = lapply(data$obs_index, function(x) A[x, ])
A_mats_single_site = lapply(single_site_data$obs_index, function(x) A[x, ])

# Compute C for each model fit, using only the data that the models were
# fitted to, i.e. we have more available data for computing C for the
# global model fit than we have for the single site fits
fits$C = list()
pb = progress_bar(length(fits$fit))
for (i in seq_along(fits$fit)) {

  # Find out which data and projection matrices correspond to the
  # given model fit
  if (fits$global[i]) {
    mydata = data
    my_A = A_mats
  } else {
    mydata = lapply(single_site_data, `[`, i - 1)
    my_A = A_mats_single_site[i - 1]
  }

  # Compute H using the estimate provided by R-INLA
  H = solve(fits$fit[[i]]$misc$cov.intern)[fits$reordering[[i]], fits$reordering[[i]]]

  # Compute all the gradients
  grads = loglik_grad(
    theta = fits$fit[[i]]$mode$theta[fits$reordering[[i]]],
    y = mydata$y,
    y0 = mydata$y0,
    dist_to_s0 = mydata$dist_to_s0,
    rho_b = if (fits$fixed[i]) exp(log_rho_b) else NULL,
    A = my_A,
    dist_to_s0_from_mesh = mydata$dist_to_s0_from_mesh,
    n_cores = n_cores,
    sum_terms = FALSE)

  # Estimate J using a sliding window approach where threshold exceedances are
  # considered as correlated if they happen less than 5 hours apart.
  times = unlist(mydata$time_index)
  window_radius = 5
  J = 0
  for (j in seq_along(times)) {
    index = which(
      times >= times[j] - window_radius &
        times <= times[j] + window_radius)
    for (k in index) {
      J = J + grads[j, ] %*% grads[k, , drop = FALSE]
    }
  }

  # Compute C given the estimates of H and J
  fits$C[[i]] = get_C(H, J)

  pb$tick()
}

# Save the list of model fits, now that we have added more information to it
saveRDS(fits, filename)

# ===============================================================
# Plot the adjusted/unadjusted posteriors for θ
# ===============================================================
fits = readRDS(filename)

n_posterior_samples = 1e5
plot_df = list()
set.seed(123)
for (i in seq_along(fits$fit)) {

  # Sample from the posteriors
  theta_uncorrected = inla.hyperpar.sample(
    n = n_posterior_samples,
    result = fits$fit[[i]],
    intern = TRUE)
  theta_uncorrected = theta_uncorrected[, fits$reordering[[i]]]

  # Switch between adjusting and not adjusting the posterior samples
  for (k in 1:2) {
    if (k == 1) {
      theta_corrected = theta_uncorrected
    } else {
      theta_corrected = matrix(
        rep(fits$fit[[i]]$mode$theta[fits$reordering[[i]]], each = n_posterior_samples),
        nrow = n_posterior_samples)
      theta_corrected = theta_corrected +
        (theta_uncorrected - theta_corrected) %*% t(fits$C[[i]])
    }
    theta_corrected = exp(theta_corrected)

    # Create a data.frame with relevant information
    df = as.data.frame(theta_corrected)
    names(df)[seq_along(fits$tex_names[[i]])] = fits$tex_names[[i]]
    df = tidyr::pivot_longer(df, tidyselect::everything()) |>
      dplyr::mutate(k = k, fixed_rho_b = fits$fixed[i], i = i) |>
      dplyr::mutate(name = factor(name, levels = fits$tex_names[[!!i]]))

    # Add the data.frame to the plot_df list
    plot_df[[length(plot_df) + 1]] = df
  }
}

# Create the final plot_df data.frame
plot_df = do.call(rbind, plot_df) |>
  dplyr::mutate(k = factor(k), i = factor(i)) |>
  dplyr::mutate(
    tag = case_when(
      i == 1 & k == 1 ~ 1,
      i == 1 & k == 2 ~ 2,
      i == 12 & k == 1 ~ 3,
      i == 12 & k == 2 ~ 4,
      k == 1 ~ 5,
      k == 2 ~ 6),
    group = paste0("i = ", i, ", k = ", k),
    tag = factor(
      tag,
      levels = 1:6,
      labels = c("Unadjusted global", "Adjusted global", "Unadjusted global 2",
                 "Adjusted global 2", "Unadjusted single site", "Adjusted single site")))

# Create the actual plot. First, we filter the plot_df and remove data that
# makes the final plot too cluttered or confusing. Explanationa for
# each filtering is given below
plot = plot_df |>
  # Select a subset of all the models
  dplyr::filter(i %in% c(1, 3, 4, 8, 9, 10)) |>
  # The density of the unadjusted global model with fixed ρ_b is too sharp,
  # and makes it difficult to focus on any of the other densities.
  # In addition, the KLD minimiser θ* for this model is different than that
  # of the models where ρ_b is estimated. So including the adjusted version
  # in the plot just makes it more difficult to extract the relevant information.
  # Therefore, we also remove the adjusted version of the model with fixed ρ_b
  dplyr::filter(! as.numeric(tag) %in% 3:4) |>
  # Remove observations with really low densities to make the plots more "readable"
  dplyr::filter(
    !(name == "$\\kappa$" & value > .9),
    !(name == "$\\kappa$" & value < .3),
    !(name == "$\\lambda$" & value > 30),
    !(name == "$\\lambda$" & value < 5),
    !(name == "$\\rho$" & value > 20),
    !(name == "$\\rho$" & value < 5),
    !(name == "$\\sigma$" & value > 2.5),
    !(name == "$\\sigma$" & value < 1.5),
    !(name == "$\\rho_b$" & value > .7),
    !(name == "$\\tau$" & value > 100)) |>
  # Define the plotting order for the different groups, so we get the densities
  # we care the most about on top of the others
  dplyr::mutate(
    group = factor(
      group,
      levels = unique(plot_df$group)[c(seq(3, 22, by = 2), seq(4, 22, by = 2), 1:2)])
  ) |>
  # Now, do the actual plotting
  ggplot(aes(x = value, col = tag, linetype = tag, group = group, size = tag)) +
  geom_density(trim = TRUE) +
  facet_wrap(~name, scales = "free") +
  scale_color_manual(values = c("black", "black", "gray", "dimgray")) +
  scale_linetype_manual(values = c("dotted", "solid", "dotted", "solid")) +
  scale_size_manual(values = c(1.5, 1.5, 1.2, .7)) +
  labs(color = "", y = "Density", x = "Value", linetype = "", size = "") +
  #scale_y_sqrt() +
  theme_light() +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(size = 15, colour = "black"),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
    legend.position = "top")

# Export the plot to pdf format
tikz_plot(file.path(image_dir(), "case-study-densities.pdf"),
          plot, width = 10, height = 6)

# ===============================================================
# Plot the adjusted/unadjusted posteriors for α(d) and ζ(d)
# ===============================================================
dist = c(5, 10, 20, 40, 60) # All distances of interest for plotting
n_posterior_samples = 1e5
plot_df = list()
set.seed(1)
for (i in seq_along(fits$fit)) {

  # Sample from the posteriors
  theta_uncorrected = inla.hyperpar.sample(
    n = n_posterior_samples,
    result = fits$fit[[i]],
    intern = TRUE)
  theta_uncorrected = theta_uncorrected[, fits$reordering[[i]]]

  # Switch between adjusting and not adjusting the posterior samples
  for (k in 1:2) {
    if (k == 1) {
      theta_corrected = theta_uncorrected
    } else {
      theta_corrected = matrix(
        rep(fits$fit[[i]]$mode$theta[fits$reordering[[i]]], each = n_posterior_samples),
        nrow = n_posterior_samples)
      theta_corrected = theta_corrected +
        (theta_uncorrected - theta_corrected) %*% t(fits$C[[i]])
    }
    theta_corrected = exp(theta_corrected)

    # Compute posteriors for α(d) and ζ(d) for all distances of interest
    for (j in seq_along(dist)) {
      alpha = exp(-(dist[j] / theta_corrected[, 1])^theta_corrected[, 2])
      if (fits$fixed[i]) {
        my_rho = fits$log_rho_b
      } else {
        my_rho = theta_corrected[, 3]
      }
      zeta = theta_corrected[, ncol(theta_corrected) - 1] *
        sqrt(1 - exp(-2 * dist[j] / my_rho))

      # Create a data.frame of relevant information
      plot_df[[length(plot_df) + 1]] = data.frame(
        value = c(alpha, zeta),
        name = rep(c("$\\alpha$", "$\\zeta$"), c(length(alpha), length(zeta))),
        dist = dist[j],
        k = k,
        fixed_rho_b = fits$fixed[i],
        global = fits$global[i],
        i = i)
    }
  }
}

# Create the final plot_df data.frame
plot_df = do.call(rbind, plot_df) |>
  dplyr::mutate(k = factor(k), i = factor(i), name = factor(name)) |>
  dplyr::mutate(
    tag = case_when(
      i == 1 & k == 1 ~ 1,
      i == 1 & k == 2 ~ 2,
      i == 12 & k == 1 ~ 3,
      i == 12 & k == 2 ~ 4,
      k == 1 ~ 5,
      k == 2 ~ 6),
    group = paste0("i = ", i, ", k = ", k),
    tag = factor(
      tag,
      levels = 1:6,
      labels = c("Unadjusted global", "Adjusted global", "Unadjusted global 2",
                 "Adjusted global 2", "Unadjusted single site", "Adjusted single site")))

plot = plot_df |>
  # Select a subset of all the models
  dplyr::filter(i %in% c(1, 3, 4, 8, 9, 10)) |>
  # The density of the unadjusted global model with fixed ρ_b is too sharp,
  # and makes it difficult to focus on any of the other densities.
  # In addition, the KLD minimiser θ* for this model is different than that
  # of the models where ρ_b is estimated. So including the adjusted version
  # in the plot just makes it more difficult to extract the relevant information.
  # Therefore, we also remove the adjusted version of the model with fixed ρ_b
  dplyr::filter(! as.numeric(tag) %in% 3:4) |>
  # Remove observations with really low densities to make the plots more "readable"
  dplyr::filter(
    !(name == "$\\alpha$" & dist == 5 & value < .45),
    !(name == "$\\alpha$" & dist == 5 & value > .75),
    !(name == "$\\alpha$" & dist == 10 & value < .35),
    !(name == "$\\alpha$" & dist == 10 & value > .65),
    !(name == "$\\alpha$" & dist == 20 & value < .15),
    !(name == "$\\alpha$" & dist == 20 & value > .5),
    !(name == "$\\alpha$" & dist == 40 & value < .05),
    !(name == "$\\alpha$" & dist == 40 & value > .35),
    !(name == "$\\alpha$" & dist == 60 & value > .3),
    !(name == "$\\zeta$" & value > 2.5),
    !(name == "$\\zeta$" & value < 1.5)) |>
  # Create the nametags that will be used as faceting labels in the final plot
  dplyr::mutate(
    nametag = paste0(
      sub("\\$$", "", name), "(",
      ifelse(name == "$\\zeta$", "x", paste(dist, "\\text{km}")), ")",
      ifelse(name == "$\\zeta$", ",\\ x > 2 \\text{km}", ""), "$"),
    nametag = factor(nametag, levels = unique(nametag)[c(2, 1, 3:6)]),
    nametag_nr = as.numeric(nametag)) |>
  # Define the plotting order for the different groups, so we get the densities
  # we care the most about on top of the others
  dplyr::mutate(
    group = factor(
      group,
      levels = unique(plot_df$group)[c(seq(3, 22, by = 2), seq(4, 22, by = 2), 1:2)])
  ) |>
  # Now, do the actual plotting
  ggplot(aes(x = value, col = tag, linetype = tag, group = group, size = tag)) +
  geom_density(trim = TRUE) +
  facet_wrap(~nametag, scales = "free", nrow = 2) +
  scale_color_manual(values = c("black", "black", "gray", "dimgray")) +
  scale_linetype_manual(values = c("dotted", "solid", "dotted", "solid")) +
  scale_size_manual(values = c(1.5, 1.5, 1.2, .7)) +
  labs(color = "", y = "Density", x = "Value", linetype = "", size = "") +
  theme_light() +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(size = 15, colour = "black"),
    strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
    legend.position = "top")

# Export the plot to pdf format
tikz_plot(file.path(image_dir(), "case-study-densities2.pdf"),
          plot, width = 10, height = 6)

# ===============================================================
# Compute out-of-sample log-scores
# ===============================================================

# First, simulate θ from all the different posteriors, and save all the samples
n_posterior_samples = 500
theta_samples = list()
set.seed(1)
for (i in seq_along(fits$fit)) {

  # Sample from the posteriors
  theta_uncorrected = inla.hyperpar.sample(
    n = n_posterior_samples,
    result = fits$fit[[i]],
    intern = TRUE)
  theta_uncorrected = theta_uncorrected[, fits$reordering[[i]]]

  # Switch between adjusting and not adjusting the posterior samples
  for (k in 1:2) {
    if (k == 1) {
      theta_corrected = theta_uncorrected
    } else {
      theta_corrected = matrix(
        rep(fits$fit[[i]]$mode$theta[fits$reordering[[i]]], each = n_posterior_samples),
        nrow = n_posterior_samples)
      theta_corrected = theta_corrected +
        (theta_uncorrected - theta_corrected) %*% t(fits$C[[i]])
    }

    # Save all the samples in a list, with some additional information
    theta_samples[[length(theta_samples) + 1]] = list(
      theta = theta_corrected,
      rho_b = if (fits$fixed[i]) exp(fits$log_rho_b) else NULL,
      model = paste0("nr. ", i, ", ", ifelse(k == 1, "unadjusted", "adjusted")))
  }
}

# Compute all the projection matrices that are necessary for computing
# out-of-sample log-scores with the data in eval_data
A_mats = lapply(eval_data$obs_index, function(x) inla.spde.make.A(mesh, coords[x, ]))

# Create a filename for storing the log-scores as we compute them
log_score_filename = file.path(results_dir(), "final-modelling-log-scores.rds")

# Preallocate a matrix for containing the log-scores of each threshold exceedance
# for each model fit
log_scores = matrix(-Inf, nrow = sum(eval_data$n), ncol = length(theta_samples))
colnames(log_scores) = sapply(theta_samples, `[[`, "model")
attributes(log_scores)$count = 0
attributes(log_scores)$i = NULL
# Compute the actual log-scores. This takes a lot of time, since we have many
# models and much data to evaluate them on
pb = progress_bar(n_posterior_samples * length(theta_samples))
for (i in seq_len(n_posterior_samples)) {
  for (j in seq_along(theta_samples)) {
    # Compute all the terms of the log-likelihood for sample nr i of posterior nr. j
    all_loglik_terms = loglik(
      theta = theta_samples[[j]]$theta[i, ],
      y = eval_data$y,
      y0 = eval_data$y0,
      dist_to_s0 = eval_data$dist_to_s0,
      dist_to_s0_from_mesh = eval_data$dist_to_s0_from_mesh,
      sum_terms = FALSE,
      A = A_mats,
      n_cores = n_cores,
      rho_b = theta_samples[[j]]$rho_b)
    for (k in seq_along(all_loglik_terms)) {
      # Update log_scores[k, j] to be equal to the logarithm
      # of exp(log_scores[k, j]) + exp(all_loglik_terms[k])
      log_scores[k, j] = log_sum_exp(
        x = c(log_scores[k, j], all_loglik_terms[k]),
        na.rm = TRUE)
    }
    pb$tick()
  }
  # Save the temporary results
  attributes(log_scores)$count = attributes(log_scores)$count + 1
  attributes(log_scores)$i = c(attributes(log_scores)$i, i)
  saveRDS(log_scores, log_score_filename)
}
pb$terminate()

# ===============================================================
# Plot the results from the model comparison using log-scores
# ===============================================================

# Read the data
log_scores = readRDS(log_score_filename)

# Rename the columns of `log_scores`, to make them better for plotting
ll_names = colnames(log_scores)
ll_names = sub("nr\\. 1,", "Global,", ll_names)
ll_names = sub("nr\\. 12,", "Global, fixed $\\\\rho_b$,", ll_names)
ll_names = sub("nr\\.", "Single site", ll_names)
for (i in 2:100) ll_names = sub(paste("Single site", i), paste("Single site", i - 1), ll_names)
colnames(log_scores) = ll_names

# Create the data.frame plot_data, that contains the ranking of all model fits
# for each of the conditioning sites used in eval_data
plot_data = as.data.frame(log_scores) |>
  dplyr::mutate(index = rep(seq_along(eval_data$n), eval_data$n)) |>
  tidyr::pivot_longer(-index) |>
  dplyr::group_by(index, name) |>
  dplyr::summarise(value = sum(value)) |>
  tidyr::pivot_wider(names_from = name, values_from = value) |>
  dplyr::ungroup() |>
  dplyr::select(-index) |>
  as.matrix()
for (i in seq_len(nrow(plot_data))) {
  plot_data[i, ] = rank(-plot_data[i, ])
}
model_names = unique(sub(", [^,]+$", "", colnames(plot_data)))
plot_data = lapply(
  X = model_names,
  FUN = function(name) {
    adjusted = plot_data[, which(colnames(plot_data) == paste0(name, ", adjusted"))]
    unadjusted = plot_data[, which(colnames(plot_data) == paste0(name, ", unadjusted"))]
    data.frame(
      value = c(adjusted, unadjusted),
      adjusted = rep(c("Adjusted", "Unadjusted"), each = nrow(plot_data)),
      model = name)
  })
plot_data = do.call(rbind, plot_data)

plot = plot_data |>
  dplyr::mutate(model = factor(model, levels = model_names[c(1:3, 5:12, 4)])) |>
  ggplot() +
  geom_histogram(
    aes(x = value, y = ..density.., col = adjusted, fill = adjusted),
    boundary = .5, binwidth = 1, position = "identity") +
  facet_wrap(~model, nrow = 2, scales = "free_x") +
  expand_limits(x = c(1, length(ll_names))) +
  scale_x_continuous(
    breaks = seq(1, length(ll_names), by = 3), expand = c(0, 0),
    minor_breaks = seq(2, length(ll_names), by = 1)) +
  labs(x = "Ranking", y = "Density", col = "", fill = "") +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = c("black", "gray")) +
  theme_light() +
  theme(text = element_text(size = 15)) +
  theme(strip.text = element_text(size = 13, colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
        legend.position = "top")

# Export the plot to pdf format
tikz_plot(file.path(image_dir(), "log-score-rankings.pdf"),
          plot, width = 12, height = 6)
