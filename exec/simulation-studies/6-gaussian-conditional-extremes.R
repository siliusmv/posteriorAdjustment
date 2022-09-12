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

filename = file.path(results_dir(), "gaussian-conditional-extremes.rds")
if (!file.exists(filename)) saveRDS(list(), filename)

# Parameter names used for creating a latex-friendly table of the results
theta_names = c("log_lambda", "log_kappa", "log_rho", "log_sigma", "log_precision")
theta_tex_names = paste0("$\\", c("lambda", "kappa", "rho", "sigma", "tau"), "$")

# ==============================================================================
# Load the case study data to create a similar setting in the simulation
# as for the case study
# ==============================================================================
radar = readRDS(file.path(downloads_dir(), "radar.rds"))
coords = st_coordinates(radar$coords)
n_loc = nrow(coords)
rm(radar)

coords = expand.grid(x = 275:290, y = 7090:7105) |>
  as.matrix()
n_loc = nrow(coords)

# Create the mesh and the SPDE
mesh = inla.mesh.2d(
  loc = coords,
  boundary = list(inla.nonconvex.hull(coords, convex = -.1),
                  inla.nonconvex.hull(coords, convex = -1)),
  max.edge = c(5, 30))
spde = inla.spde2.pcmatern(mesh, prior.range = c(60, .95), prior.sigma = c(4, .05))

# Create the projection matrix for the mesh and SPDE
A = inla.spde.make.A(mesh, coords)

# Sample a ton of Gaussian realisations
Q = inla.spde2.precision(spde, c(log(8), log(1)))

transformed_data = local({
  res = list()
  while (TRUE) {
    samples = as.matrix(A %*% rnorm_spde(1e4, Q))
    samples = qlaplace(pnorm(samples))
    good_index = which(apply(samples, 2, function(x) any(x > threshold)))
    res[[length(res) + 1]] = samples[, good_index, drop = FALSE]
    if (sum(sapply(res, ncol)) > 5e4) break
  }
  res = do.call(cbind, res)
  t(res)
})

loglik_pairwise = function(theta,
                           y,
                           y0,
                           alpha = NULL,
                           zeta = NULL,
                           beta = NULL,
                           gamma = NULL,
                           lower_lims = c(-Inf, 0, -Inf, -Inf),
                           upper_lims = rep(Inf, 4)) {
  i = 1
  if (is.null(alpha)) {
    alpha = theta[i]
    i = i + 1
  }
  if (is.null(zeta)) {
    zeta = theta[i]
    i = i + 1
  }
  if (is.null(beta)) {
    beta = theta[i]
    i = i + 1
  }
  if (is.null(gamma)) {
    gamma = theta[i]
    i = i + 1
  }

  # Test if the parameters are inside the upper/lower limits
  if (any(lower_lims > c(alpha, zeta, beta, gamma)) ||
        any(upper_lims < c(alpha, zeta, beta, gamma))) {
    return(-1e99)
  }

  # Compute the log-likelihood
  sum(dnorm(y, mean = alpha * y0 + gamma, sd = zeta * y0^beta, log = TRUE), na.rm = TRUE)
}

# ==============================================================================
# Fit pairwise models
# ==============================================================================
lower_lims = c(-Inf, 0, -Inf, -Inf)
upper_lims = c(Inf, 10, Inf, Inf)

n_pairs = 5000
cl = parallel::makeForkCluster(n_cores, mc.set.seed = TRUE)
parallel::clusterSetRNGStream(cl, 1)
pars = pbapply::pblapply(
  X = 1:n_pairs,
  cl = cl,
  FUN = function(i) {
    # Draw s1 and s0 randomly, then find the values of y and y0
    loc_index = sample.int(n_loc, 2)
    y0 = as.numeric(transformed_data[, loc_index[1]])
    y = as.numeric(transformed_data[, loc_index[2]])

    # Find all times where y0 is big and y is not NA
    y0_big_index = which(y0 > threshold)
    na_index = which(is.na(y))
    good_index = setdiff(y0_big_index, na_index)
    if (length(good_index) == 0) return(NULL)

    # Fit the different model where we vary between fixing or estimating gamma and beta
    p1 = optim(
      par = c(.5, 1),
      fn = loglik_pairwise,
      control = list(fnscale = -1),
      y = y[good_index],
      gamma = 0,
      beta = 0,
      lower_lims = lower_lims,
      upper_lims = upper_lims,
      y0 = y0[good_index])$par
    p2 = optim(
      par = c(.5, 1, .5),
      fn = loglik_pairwise,
      control = list(fnscale = -1),
      gamma = 0,
      lower_lims = lower_lims,
      upper_lims = upper_lims,
      y = y[good_index],
      y0 = y0[good_index])$par
    p3 = optim(
      par = c(.5, 1, 0),
      fn = loglik_pairwise,
      control = list(fnscale = -1),
      beta = 0,
      lower_lims = lower_lims,
      upper_lims = upper_lims,
      y = y[good_index],
      y0 = y0[good_index])$par
    p4 = optim(
      par = c(.5, 1, .5, 0),
      lower_lims = lower_lims,
      upper_lims = upper_lims,
      fn = loglik_pairwise,
      control = list(fnscale = -1),
      y = y[good_index],
      y0 = y0[good_index])$par

    # Return the results in a data.frame format
    res = data.frame(
      value = c(p1, p2, p3, p4),
      dist = as.numeric(dist(coords[loc_index, 1:2])),
      iter = i,
      n_y0 = length(good_index),
      par = c("alpha", "zeta",
              "alpha", "zeta", "beta",
              "alpha", "zeta", "gamma",
              "alpha", "zeta", "beta", "gamma"),
      model = c(1, 1,
                2, 2, 2,
                3, 3, 3,
                4, 4, 4, 4))

    res
  })
parallel::stopCluster(cl)

# ==============================================================================
# Plot the results
# ==============================================================================
do.call(rbind, pars) |>
  dplyr::mutate(
    model = paste("Model", model),
    par = factor(
      par,
      levels = c("alpha", "zeta", "beta", "gamma"),
      labels = paste0("$\\", c("alpha", "zeta", "beta", "gamma"), "$"))) |>
  ggplot() +
  geom_point(aes(x = dist, y = value), alpha = .05) +
  facet_grid(par ~ model, scales = "free") +
  theme_light() +
  theme(text = element_text(size = 15)) +
  theme(strip.text.y = element_text(angle = 0, size = 15, colour = "black"),
        strip.text.x = element_text(size = 12, colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
  labs(x = "Distance [km]", y = "Value")

# Extract all observations from when a threshold exceedance is observed,
# using every single location as a possible conditioning site
data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = 1:n_loc,
  threshold = threshold,
  n = 1,
  r = Inf)

# Compute the distances to the conditioning sites from the mesh nodes,
# and add the information to `data`
data$dist_to_s0_from_mesh = list()
for (i in seq_along(data$s0)) {
  data$dist_to_s0_from_mesh[[i]] = dist_euclid(mesh$loc[, 1:2], data$s0[[i]])
}

# Create functions for computing a(y0, d) = α(d) * y0, and ζ(d)
get_a_func = function(theta) {
  λ = exp(theta[1])
  κ = exp(theta[2])
  function(y0, dist) {
    alpha = exp(- (dist / λ)^κ)
    matrix(rep(y0, each = length(alpha)) * rep(alpha, length(y0)),
           nrow = length(dist),
           ncol = length(y0))
  }
}
get_zeta_func = function(theta) {
  zeta0 = exp(theta[1])
  rho = exp(theta[2])
  function(y0, dist) {
    n = length(y0)
    zeta = sqrt(1 - exp(-2 * dist / rho)) * zeta0
    matrix(rep(zeta, n), nrow = length(dist), ncol = n)
  }
}

# Create a function for computing the log-likelihood of the conditional extremes
# model, when Z is iid Gaussian white noise, meaning that all observations are
# independent, and Z_b is nonstationary Gaussian white noise
loglik_iid = function(theta, data) {
  a_func = get_a_func(theta[1:2])
  zeta_func = get_zeta_func(theta[3:4])
  res = 0
  for (i in seq_along(data$y)) {
    a = a_func(data$y0[[i]], data$dist_to_s0[[i]])
    zeta = zeta_func(data$y0[[i]], data$dist_to_s0[[i]])
    res = res + sum(dnorm(data$y[[i]], a, zeta, log = TRUE), na.rm = TRUE)
  }
  res
}

# Compute the maximum likelihood estimators for the parameters of α(d) and ζ(d)
res = optim(
  par = c(2.2724136, 0.1629746, 0.4168857, 1.1959553),
  fn = loglik_iid,
  control = list(fnscale = -1, trace = 6),
  data = data)

# Compute the maximum likelihood estimators of a(y0, d) and ζ(d)
a_func = get_a_func(res$par[1:2])
zeta_func = get_zeta_func(res$par[3:4])

rho_b = exp(res$par[4])

# =====================================================================
# Pairwise estimates with fixed alpha, zeta (without upper/lower limits)
# =====================================================================
n_pairs = 5000
cl = parallel::makeForkCluster(n_cores, mc.set.seed = TRUE)
parallel::clusterSetRNGStream(cl, 1)
pars2 = pbapply::pblapply(
  X = 1:n_pairs,
  cl = cl,
  FUN = function(i) {
    # Draw s1 and s0 randomly, then find the values of y and y0
    loc_index = sample.int(n_loc, 2)
    y0 = as.numeric(transformed_data[, loc_index[1]])
    y = as.numeric(transformed_data[, loc_index[2]])

    # Find all times where y0 is big and y is not NA
    y0_big_index = which(y0 > threshold)
    na_index = which(is.na(y))
    good_index = setdiff(y0_big_index, na_index)
    if (length(good_index) == 0) return(NULL)

    # Compute the values of alpha and zeta
    dist = as.numeric(dist(coords[loc_index, 1:2]))
    alpha = as.numeric(a_func(1, dist))
    zeta = as.numeric(zeta_func(1, dist))

    # Fit the different model where we fix alpha and zeta
    p1 = optim(
      par = .5,
      fn = loglik_pairwise,
      control = list(fnscale = -1),
      method = "Brent", lower = -1, upper = 2,
      alpha = alpha,
      gamma = 0,
      zeta = zeta,
      y = y[y0_big_index],
      y0 = y0[y0_big_index])$par
    p2 = optim(
      par = 0,
      fn = loglik_pairwise,
      control = list(fnscale = -1),
      method = "Brent", lower = -5, upper = 5,
      alpha = alpha,
      zeta = zeta,
      beta = 0,
      y = y[y0_big_index],
      y0 = y0[y0_big_index])$par
    p3 = optim(
      par = c(.5, 0),
      fn = loglik_pairwise,
      control = list(fnscale = -1),
      alpha = alpha,
      zeta = zeta,
      y = y[y0_big_index],
      y0 = y0[y0_big_index])$par

    # Return the results in a data.frame format
    res = data.frame(
      value = c(p1, p2, p3),
      dist = dist,
      iter = i,
      n_y0 = length(good_index),
      par = c("beta", "gamma", "beta", "gamma"),
      model = c(1, 2, 3, 3))

    res
  })
parallel::stopCluster(cl)

# ==============================================================================
# Plot the results
# ==============================================================================

# Examine the pairwise estimates of β and γ with α(d) and ζ(d) fixed
do.call(rbind, pars2) |>
  dplyr::mutate(
    model = paste("Model", model),
    par = factor(
      par,
      levels = c("alpha", "zeta", "beta", "gamma"),
      labels = paste0("$\\", c("alpha", "zeta", "beta", "gamma"), "$"))) |>
  ggplot() +
  geom_point(aes(x = dist, y = value), alpha = .05) +
  facet_grid(par ~ model, scales = "free") +
  theme_light() +
  theme(text = element_text(size = 15)) +
  theme(strip.text.y = element_text(angle = 0, size = 15, colour = "black"),
        strip.text.x = element_text(size = 12, colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
  labs(x = "Distance [km]", y = "Value")


# Plot the estimates of α(d) and ζ(d) together with their pairwise estimates
xx = seq(0, max(sapply(pars, `[[`, "dist")), length.out = 1000)
do.call(rbind, pars) |>
  dplyr::filter(model == 1) |>
  dplyr::mutate(par = factor(
    par,
    levels = c("alpha", "zeta", "beta", "gamma"),
    labels = paste0("$\\", c("alpha", "zeta", "beta", "gamma"), "$"))) |>
  ggplot() +
  geom_point(aes(x = dist, y = value), alpha = .1) +
  facet_wrap(~par, scales = "free", nrow = 1) +
  geom_line(data = data.frame(x = xx, y = a_func(1, xx), par = "$\\alpha$"),
            aes(x = x, y = y), col = "blue", size = 2) +
  geom_line(data = data.frame(x = xx, y = zeta_func(1, xx), par = "$\\zeta$"),
            aes(x = x, y = y), col = "blue", size = 2) +
  theme_light() +
  theme(text = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15, colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
  labs(x = "Distance [km]", y = "Value")


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

thinning = c(1, 2, 4, 6, 8, 16, 32)
data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n = thinning,
  r = cumsum(4 * thinning))

data$dist_to_s0_from_mesh = list()
for (i in seq_along(data$s0)) {
  data$dist_to_s0_from_mesh[[i]] = dist_euclid(mesh$loc[, 1:2], data$s0[[i]])
}

# Compute the maximum likelihood estimator, which is approximately equal to θ*.
# Use θ as initial values
est = list(par = c(res$par[c(1, 2, 4)], log(20), res$par[3], 4), convergence = 1)
est = list(par = c(res$par[c(1, 2, 4)], res$par[3], 4), convergence = 1)
while (est$convergence != 0) {
  est = optim(
    par = est$par,
    fn = loglik,
    y = data$y,
    y0 = data$y0,
    rho_b = rho_b,
    dist_to_s0 = data$dist_to_s0,
    dist_to_s0_from_mesh = data$dist_to_s0_from_mesh,
    A = lapply(data$obs_index, function(x) A[x, ]),
    n_cores = n_cores,
    control = list(fnscale = -1, maxit = 300, trace = 6, REPORT = 1))
}
theta_star = est$par

# > theta_star
# [1]  2.2726617  0.2282242 -0.3371791  2.5447459  0.3760144  5.6419177

# ===============================================================================
# Sample smaller amounts of data and fit the conditional extremes model
# many times, in order to estimate coverage percentages
# ===============================================================================

single_site_index = s0_index[6]

# Run all the n_repl experiments
cl = parallel::makeForkCluster(n_cores, mc.set.seed = TRUE)
parallel::clusterSetRNGStream(cl, 1)
pbapply::pblapply(
  X = 1:n_repl,
  cl = cl,
  FUN = function(i) {

    # Store the state of the random seed before we start simulating random data
    seed = .Random.seed

    transformed_data = local({
      samples = rnorm_spde(5e2, Q)
      samples = as.matrix(A %*% samples)
      samples = qlaplace(pnorm(samples))
      t(samples)
    })

    # Extract observations in a nonregular grid around conditioning sites that
    # exceed the threshold
    thinning = c(1, 2, 4, 6, 8, 16, 32)
    data = extract_extreme_fields(
      data = transformed_data,
      coords = coords,
      s0_index = s0_index,
      threshold = threshold,
      n = thinning,
      r = cumsum(4 * thinning))

    single_site_data = extract_extreme_fields(
      data = transformed_data,
      coords = coords,
      s0_index = single_site_index,
      threshold = threshold,
      n = thinning,
      r = cumsum(4 * thinning))

    rm(transformed_data)

    # Compute distances from the mesh to the conditioning sites
    data$dist_to_s0_from_mesh = list()
    for (j in seq_along(data$s0)) {
      data$dist_to_s0_from_mesh[[j]] = dist_euclid(mesh$loc[, 1:2], data$s0[[j]])
    }

    single_site_data$dist_to_s0_from_mesh = list()
    for (j in seq_along(single_site_data$s0)) {
      single_site_data$dist_to_s0_from_mesh[[j]] = dist_euclid(mesh$loc[, 1:2], single_site_data$s0[[j]])
    }

    # Estimate the maximum likelihood estimator to get good initial values for the
    # model fitting. We only allow maxit = 120 because this gives initial values
    # that are good enough, whitout spending too much time locating the actual maximum
    est = optim(
      par = theta_star,
      fn = loglik,
      y = data$y,
      y0 = data$y0,
      rho_b = rho_b,
      dist_to_s0 = data$dist_to_s0,
      dist_to_s0_from_mesh = data$dist_to_s0_from_mesh,
      A = lapply(data$obs_index, function(x) A[x, ]),
      n_cores = 1,
      control = list(fnscale = -1, maxit = 120))

    single_site_est = optim(
      par = theta_star,
      fn = loglik,
      y = single_site_data$y,
      y0 = single_site_data$y0,
      rho_b = rho_b,
      dist_to_s0 = single_site_data$dist_to_s0,
      dist_to_s0_from_mesh = single_site_data$dist_to_s0_from_mesh,
      A = lapply(single_site_data$obs_index, function(x) A[x, ]),
      n_cores = 1,
      control = list(fnscale = -1, maxit = 120))

    # R-INLA requires the observations y to be on a vector format
    y_inla = unlist(data$y)
    y_single_site_inla = unlist(single_site_data$y)

    # Our implementation of the cgeneric model for a requires y0 and dist_to_s0 as input.
    # However, it requires one value of y0 and dist_to_s0 for each of the observations y.
    # We therefore replicate y0 and dist_to_s0 in the correct way so we get to vectors with
    # equal length to y_inla, where y_inla[i] was observed when y(s0) was equal to
    # y0_inla[i], and the distance between them was dist_to_s0_inla[i].
    dist_to_s0_inla = data$dist_to_s0[rep(seq_along(data$n), data$n)]
    y0_inla = rep(unlist(data$y0), sapply(dist_to_s0_inla, length))
    dist_to_s0_inla = unlist(dist_to_s0_inla)

    dist_to_s0_single_site_inla = single_site_data$dist_to_s0[rep(seq_along(single_site_data$n), single_site_data$n)]
    y0_single_site_inla = rep(unlist(single_site_data$y0), sapply(dist_to_s0_single_site_inla, length))
    dist_to_s0_single_site_inla = unlist(dist_to_s0_single_site_inla)

    # Define the cgeneric model for a
    a_priors = list(lambda = c(1, 3, 3), kappa = c(1, 0, 3))
    a_model = a_generic_model(
      y0 = y0_inla,
      dist_to_s0 = dist_to_s0_inla,
      init = est$par[1:2],
      priors = a_priors)

    a_single_site_model = a_generic_model(
      y0 = y0_single_site_inla,
      dist_to_s0 = dist_to_s0_single_site_inla,
      init = single_site_est$par[1:2],
      priors = a_priors)

    # Define the cgeneric model for Z_b
    spde_priors = list(
      rho = c(1, 40, .95),
      sigma = c(1, 3, .05),
      #rho_b = c(1, log(6), 2))
      rho_b = c(0, log(rho_b)))
    spde_model = spde_generic_model_with_b_func(
      spde = spde,
      n = data$n,
      #init = est$par[c(4, 5, 3)],
      init = est$par[c(4, 5)],
      priors = spde_priors,
      dist_to_s0 = do.call(rbind, data$dist_to_s0_from_mesh))
    spde_single_site_model = spde_generic_model_with_b_func(
      spde = spde,
      n = single_site_data$n,
      #init = single_site_est$par[c(4, 5, 3)],
      init = single_site_est$par[c(4, 5)],
      priors = spde_priors,
      dist_to_s0 = do.call(rbind, single_site_data$dist_to_s0_from_mesh))

    # Create the necessary objects for running R-INLA with the conditional extremes model
    formula = y ~ -1 +
      f(spatial, model = spde_model) +
      f(idx, model = a_model)
    single_site_formula = y ~ -1 +
      f(spatial, model = spde_single_site_model) +
      f(idx, model = a_single_site_model)

    effects = list(
      spatial = seq_len(sum(data$n) * mesh$n),
      idx = seq_along(y_inla))
    single_site_effects = list(
      spatial = seq_len(sum(single_site_data$n) * mesh$n),
      idx = seq_along(y_single_site_inla))

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
    A_single_site_inla = local({
      all_location_indices = rep(single_site_data$obs_index, single_site_data$n)
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
    single_site_stack = inla.stack(
      data = list(y = y_single_site_inla),
      A = list(spatial = A_single_site_inla, 1),
      effects = single_site_effects)
    rm(A_inla, A_single_site_inla, effects, single_site_effects)

    # Fit the model
    fit = tryCatch({
      inla(
        formula = formula,
        data = inla.stack.data(stack),
        control.predictor = list(A = inla.stack.A(stack)),
        only.hyperparam = TRUE,
        #control.mode = list(theta = est$par[c(6, 4:5, 3, 1:2)], restart = TRUE),
        control.mode = list(theta = est$par[c(5, 3:4, 1:2)], restart = TRUE),
        control.inla = list(control.vb = list(enable = FALSE)),
        #verbose = TRUE,
        num.threads = 1,
        inla.mode = "experimental")
    }, error = function(e) NULL)
    if (is.null(fit)) next

    # Fit the model
    single_site_fit = tryCatch({
      inla(
        formula = single_site_formula,
        data = inla.stack.data(single_site_stack),
        control.predictor = list(A = inla.stack.A(single_site_stack)),
        only.hyperparam = TRUE,
        #control.mode = list(theta = single_site_est$par[c(6, 4:5, 3, 1:2)], restart = TRUE),
        control.mode = list(theta = single_site_est$par[c(5, 3:4, 1:2)], restart = TRUE),
        control.inla = list(control.vb = list(enable = FALSE)),
        #verbose = TRUE,
        num.threads = 1,
        inla.mode = "experimental")
    }, error = function(e) NULL)
    if (is.null(fit)) next

    # Which way should the parameters of the model fits be reordered to
    # give the correct input to the log-likelihood functions?
    theta_reordering = c(5, 6, 4, 2, 3, 1)
    theta_reordering = c(4, 5, 2, 3, 1)

    # Estimate H
    H = tryCatch(solve(fit$misc$cov.intern), error = \(e) NULL)
    if (is.null(H)) {
      warning("Inverting cov.intern failed for i = ", i)
      next
    }
    H = H[theta_reordering, theta_reordering]

    # Estimate single_site_H
    single_site_H = tryCatch(solve(single_site_fit$misc$cov.intern), error = \(e) NULL)
    if (is.null(single_site_H)) {
      warning("Inverting cov.intern failed for i = ", i)
      next
    }
    single_site_H = single_site_H[theta_reordering, theta_reordering]

    # Compute all terms of the gradient of the log-likelihood, so we can estimate J
    grads = loglik_grad(
      theta = fit$mode$theta[theta_reordering],
      y = data$y,
      y0 = data$y0,
      rho_b = rho_b,
      dist_to_s0 = data$dist_to_s0,
      dist_to_s0_from_mesh = data$dist_to_s0_from_mesh,
      A = lapply(data$obs_index, function(x) A[x, ]),
      n_cores = 1,
      sum_terms = FALSE)

    # Compute all terms of the gradient of the log-likelihood, so we can estimate J
    single_site_grads = loglik_grad(
      theta = single_site_fit$mode$theta[theta_reordering],
      y = single_site_data$y,
      y0 = single_site_data$y0,
      rho_b = rho_b,
      dist_to_s0 = single_site_data$dist_to_s0,
      dist_to_s0_from_mesh = single_site_data$dist_to_s0_from_mesh,
      A = lapply(single_site_data$obs_index, function(x) A[x, ]),
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

    # Estimate J. We know that observations from two different times are completely independent
    times = unlist(single_site_data$time_index)
    single_site_J = 0
    for (k in seq_along(times)) {
      index = which(times == times[k])
      for (j in index) {
        single_site_J = single_site_J + single_site_grads[k, ] %*% single_site_grads[j, , drop = FALSE]
      }
    }

    # Compute the estimate for C
    C = get_C(H, J)
    single_site_C = get_C(single_site_H, single_site_J)

    # Keep only the relevant parts of the model fit, to reduce requirements on file storage
    fit = list(misc = fit$misc, internal.marginals.hyperpar = fit$internal.marginals.hyperpar,
               mode = list(theta = fit$mode$theta))
    fit$misc$reordering = NULL
    fit$misc$family = NULL
    fit$misc$linkfunctions = NULL
    class(fit) = "inla"

    single_site_fit = list(
      misc = single_site_fit$misc,
      internal.marginals.hyperpar = single_site_fit$internal.marginals.hyperpar,
      mode = list(theta = single_site_fit$mode$theta))
    single_site_fit$misc$reordering = NULL
    single_site_fit$misc$family = NULL
    single_site_fit$misc$linkfunctions = NULL
    class(single_site_fit) = "inla"

    # Create a list containing all the data of interest
    model_fits = list(
      fit = fit,
      single_site_fit = single_site_fit,
      C = C,
      single_site_C = single_site_C,
      seed = seed,
      rho_b = rho_b,
      theta_reordering = theta_reordering,
      time_index = data$time_index,
      single_site_time_index = single_site_data$time_index,
      n = sum(data$n),
      single_site_n = sum(single_site_data$n))

    attributes(model_fits)$theta_star = theta_star

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
theta_star = attributes(res[[1]])$theta_star
theta_reordering = res[[1]]$theta_reordering
n_theta = length(theta_reordering)

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
  dplyr::mutate(name = factor(
    x = name,
    levels = theta_names,
    labels = theta_tex_names)) |>
  tidyr::pivot_wider(names_from = name, values_from = coverage)
#tmp = tmp[rep(2 * seq_len(nrow(tmp) / 2), each = 2) - c(0, 1), c(1:2, 4, 3, 7, 6, 8, 5)]
tmp = tmp[rep(2 * seq_len(nrow(tmp) / 2), each = 2) - c(0, 1), c(1:2, 4, 3, 6, 7, 5)]
print(tmp)

# ==============================================
# Reformat the results into latex tabular format
# ==============================================
table = paste(paste(c("Aim", "Method", names(tmp)[-c(1, 2)]), collapse = " & "), "\\\\")
table[2] = "\\midrule"
j = 3
for (i in 1:nrow(tmp)) {
  table[j] = paste(
    c(tmp$label[i], paste0("$", round(100 * tmp[i, -c(1, 2)], digits = 0), "\\%$")),
    collapse = " & ")
  if (i %% 2 == 1) {
    table[j] = paste(
      c(paste0("$", round(100 * tmp$prob[i], digits = 0), "\\%$"), table[j]), collapse = " & ")
  } else {
    table[j] = paste(c("", table[j]), collapse = " & ")
  }
  table[j] = paste(table[j], "\\\\")
  j = j + 1
  if (i %% 2 == 0 && i != nrow(tmp)) {
    table[j] = "\\midrule"
    j = j + 1
  }
}
table = paste(paste(table, collapse = "\n"), "\n")

cat(table)
