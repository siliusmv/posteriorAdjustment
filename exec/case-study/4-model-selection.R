devtools::load_all()
library(pbapply)
library(parallel)
library(ggplot2)
library(dplyr)
library(patchwork)
library(sf)

# ==============================================================================
# Fix necessary model parameters:
# ==============================================================================
threshold = qlaplace(.9975) # The threshold t for defining the conditional extremes model
n_cores = 6 # Run code in parallel
r = 5 # Radius used for computing aggregated empirical distribution functions

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
# Create the pairwise log likelihood function
# ==============================================================================
# theta: vector of parameters
# y: vector of observations at s1
# y0: vector of observations at s0
# alpha: if alpha = NULL, we estimate alpha. Else we fix alpha = alpha
# zeta: if zeta = NULL, we estimate zeta. Else we fix zeta = zeta
# beta: if beta = NULL, we estimate beta. Else we fix beta = beta
# gamma: if gamma = NULL, we estimate gamma. Else we fix gamma = gamma
# lower_lims: a vector of length 4 with lower limits for all the parameters.
#   This defaults to -Inf for everythng except zeta > 0.
# upper_lims: a vector of length 4 with upper limits for all the parameters.
#   This defaults to Inf for all parameters
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
      y0 = y0[good_index])$par
    p2 = optim(
      par = c(.5, 1, .5),
      fn = loglik_pairwise,
      control = list(fnscale = -1),
      gamma = 0,
      y = y[good_index],
      y0 = y0[good_index])$par
    p3 = optim(
      par = c(.5, 1, 0),
      fn = loglik_pairwise,
      control = list(fnscale = -1),
      beta = 0,
      y = y[good_index],
      y0 = y0[good_index])$par
    p4 = optim(
      par = c(.5, 1, .5, 0),
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
plot = do.call(rbind, pars) |>
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

# Export the plot to pdf format
tikz_plot(file.path(image_dir(), "model-selection1.pdf"),
          plot, width = 7, height = 7)

# ==============================================================================
# Some of the estimated parameter values are so large/small that it is difficult
# to detect the main patterns in the estimates. Repeat the entire procedure with
# stricter upper and lower limits for the parameters.
# ==============================================================================

lower_lims = c(-1, 0, -2, -Inf)
upper_lims = c(2, 10, 3, Inf)

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

# Plot the results
plot = do.call(rbind, pars) |>
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

# Export the plot to pdf format
tikz_plot(file.path(image_dir(), "model-selection1.pdf"),
          plot, width = 7, height = 7)

# =====================================================================
# Estimate α(d) and ζ(d)
# =====================================================================

# Extract all observations from when a threshold exceedance is observed,
# using every single location as a possible conditioning site
data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = 1:n_loc,
  threshold = threshold,
  n = 1,
  r = Inf)

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
  par = c(log(50), log(1), log(1), log(20)),
  fn = loglik_iid,
  control = list(fnscale = -1, trace = 6),
  data = data)

# Compute the maximum likelihood estimators of a(y0, d) and ζ(d)
a_func = get_a_func(res$par[1:2])
zeta_func = get_zeta_func(res$par[3:4])

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
plot1 = do.call(rbind, pars) |>
  dplyr::filter(model == 1) |>
  dplyr::mutate(par = factor(
    par,
    levels = c("alpha", "zeta", "beta", "gamma"),
    labels = paste0("$\\", c("alpha", "zeta", "beta", "gamma"), "$"))) |>
  ggplot() +
  geom_point(aes(x = dist, y = value), alpha = .2, col = "gray") +
  facet_wrap(~par, scales = "free", nrow = 1) +
  geom_line(data = data.frame(x = xx, y = a_func(1, xx), par = "$\\alpha$"),
            aes(x = x, y = y), size = 2) +
  geom_line(data = data.frame(x = xx, y = zeta_func(1, xx), par = "$\\zeta$"),
            aes(x = x, y = y), size = 2) +
  theme_light() +
  theme(text = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15, colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
  labs(x = "Distance [km]", y = "Value")

# Plot the joint pairwise estimates for β and γ with α(d) and ζ(d) fixed
plot2 = do.call(rbind, pars2) |>
  dplyr::mutate(model == 3) |>
  dplyr::mutate(par = factor(
    par,
    levels = c("alpha", "zeta", "beta", "gamma"),
    labels = paste0("$\\", c("alpha", "zeta", "beta", "gamma"), "$"))) |>
  ggplot() +
  geom_point(aes(x = dist, y = value), alpha = .1) +
  facet_wrap(~par, scales = "free") +
  theme_light() +
  theme(text = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15, colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
  labs(x = "Distance [km]", y = "Value")

# Combine the two plots into one
plot = patchwork::wrap_plots(plot1, plot2, nrow = 1)

# Export the plot to pdf format
tikz_plot(file.path(image_dir(), "model-selection2.pdf"),
          plot, width = 10, height = 5)

# ==================================================================
# Compute residuals and examine summary statistics for the residuals
# ==================================================================

residuals = list()
for (i in seq_along(data$y)) {
  a = a_func(data$y0[[i]], data$dist_to_s0[[i]])
  zeta = zeta_func(data$y0[[i]], data$dist_to_s0[[i]])
  residuals[[i]] = (data$y[[i]] - a) / zeta
}

# Summary statistics for all the residuals combined
mean(unlist(residuals), na.rm = TRUE)
var(unlist(residuals), na.rm = TRUE)

# Sort the residuals by all the unique values of dist_to_s0
unique_dists = unique(unlist(data$dist_to_s0))
dist_residuals = lapply(
  X = seq_along(unique_dists),
  FUN = function(i) {
    res = list()
    for (j in seq_along(data$dist_to_s0)) {
      index = which(data$dist_to_s0[[j]] == unique_dists[i])
      if (length(index) > 0) {
        res[[j]] = as.numeric(residuals[[j]][index, ])
      }
    }
    unlist(res)
  })

# Summary statistics for the residuals at each unique value of dist_to_s0
mean_vals = sapply(dist_residuals, mean, na.rm = TRUE)
var_vals = sapply(dist_residuals, var, na.rm = TRUE)

# Plot the summary statistics
data.frame(mean = mean_vals, var = var_vals, dist = unique_dists) |>
  tidyr::pivot_longer(-dist) |>
  ggplot() +
  geom_point(aes(x = dist, y = value)) +
  facet_wrap(~name, scales = "free")

# # =======================================================================================
# # Convert the "heavy" pdf files created in this script to "lighter" compressed jpeg files
# # =======================================================================================
# library(magick)
# image = magick::image_read_pdf(file.path(image_dir(), "model-selection1.pdf"))
# magick::image_write(
#   image = image,
#   path = file.path(image_dir(), "model-selection1.jpg"),
#   format = "jpg",
#   quality = 50)
# image = magick::image_read_pdf(file.path(image_dir(), "model-selection2.pdf"))
# magick::image_write(
#   image = image,
#   path = file.path(image_dir(), "model-selection2.jpg"),
#   format = "jpg",
#   quality = 50)
