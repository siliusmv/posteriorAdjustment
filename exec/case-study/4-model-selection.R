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
threshold = qlaplace(.95) # The threshold t for defining the conditional extremes model
n_cores = 6 # Run code in parallel
r = 5 # Radius used for computing aggregated empirical distribution functions

# ==============================================================================
# Load the radar data
# ==============================================================================
radar = readRDS(file.path(downloads_dir(), "radar.rds"))
coords = st_coordinates(radar$coords)
data = radar$data
n_loc = nrow(coords)
rm(radar)

# ==============================================================================
# Transform the data to Laplace margins
# ==============================================================================
cl = parallel::makeForkCluster(n_cores)
transformed_data = pbapply::pblapply(
  X = 1:n_loc,
  cl = cl,
  FUN = function(i) {
    F = aggregated_ecdf(data, coords, coords[i, ], radius = r)
    u = F(data[, i])
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
loglik_pairwise = function(theta,
                           y,
                           y0,
                           alpha = NULL,
                           zeta = NULL,
                           beta = NULL,
                           gamma = NULL) {
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
  if (zeta < 0) return(-1e199)
  if (zeta > 20) return(-1e199)
  if (beta < -2) return(-1e199)
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

tikz_plot(file.path(image_dir(), "model-selection1.pdf"),
          plot, width = 7, height = 7)

# =====================================================================
# Estimate alpha and zeta
# =====================================================================

data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = 1:n_loc,
  threshold = threshold,
  n = 1,
  r = Inf)

get_a_func = function(theta) {
  warning("change this to alpha_func")
  λ = exp(theta[1]); κ = exp(theta[2])
  function(y, dist) {
    alpha = exp(- (dist / λ)^κ)
    matrix(rep(y, each = length(alpha)) * rep(alpha, length(y)),
           nrow = length(dist),
           ncol = length(y))
  }
}

get_zeta_func = function(theta) {
  zeta0 = exp(theta[1])
  rho = exp(theta[2])
  function(dist, n) {
    zeta = sqrt(1 - exp(-2 * dist / rho)) * zeta0
    matrix(rep(zeta, n), nrow = length(dist), ncol = n)
  }
}

loglik = function(theta, data) {
  a_func = get_a_func(theta[1:2])
  zeta_func = get_zeta_func(theta[3:4])
  res = 0
  for (i in seq_along(data$y)) {
    a = a_func(data$y0[[i]], data$dist_to_s0[[i]])
    zeta = zeta_func(data$dist_to_s0[[i]], length(data$y0[[i]]))
    res = res + sum(dnorm(data$y[[i]], a, zeta, log = TRUE), na.rm = TRUE)
  }
  res
}

res = optim(
  par = c(log(50), log(1), log(1), log(20)),
  fn = loglik,
  control = list(fnscale = -1, trace = 6),
  data = data)

a_func = get_a_func(res$par[1:2])
zeta_func = get_zeta_func(res$par[3:4])

# =====================================================================
# Pairwise estimates with fixed alpha, zeta
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
    zeta = as.numeric(zeta_func(dist, 1))

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


xx = seq(0, max(sapply(pars, `[[`, "dist")), length.out = 1000)
plot1 = do.call(rbind, pars) |>
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
  geom_line(data = data.frame(x = xx, y = zeta_func(xx, 1), par = "$\\zeta$"),
            aes(x = x, y = y), col = "blue", size = 2) +
  theme_light() +
  theme(text = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15, colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
  labs(x = "Distance [km]", y = "Value")

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

plot = patchwork::wrap_plots(plot1, plot2, nrow = 1)

tikz_plot(file.path(image_dir(), "model-selection2.pdf"),
          plot, width = 10, height = 5)


residuals = list()
for (i in seq_along(data$y)) {
  a = a_func(data$y0[[i]], data$dist_to_s0[[i]])
  zeta = zeta_func(data$dist_to_s0[[i]], data$n[i])
  residuals[[i]] = as.numeric((data$y[[i]] - a) / zeta)
}
residuals = unlist(residuals)
mean(residuals, na.rm = TRUE)
var(residuals, na.rm = TRUE)

#library(magick)
#image = magick::image_read_pdf(file.path(image_dir(), "model-selection1.pdf"))
#magick::image_write(image, file.path(image_dir(), "model-selection1.jpg"), format = "jpg", quality = 50)
#image = magick::image_read_pdf(file.path(image_dir(), "model-selection2.pdf"))
#magick::image_write(image, file.path(image_dir(), "model-selection2.jpg"), format = "jpg", quality = 50)
