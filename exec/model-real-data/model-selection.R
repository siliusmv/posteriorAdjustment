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
p_gpd = 1
threshold = qlaplace(.999)
# How many cores to use in parallel?
num_cores = 15
r = 5

# ==============================================================================
# Load the radar data
# ==============================================================================
radar = readRDS(file.path(downloads_dir(), "radar-trondheim.rds"))
coords = st_coordinates(radar$coords)
n_loc = nrow(coords)

# ==============================================================================
# Restrict the data, for computational reasons
# ==============================================================================
bbox = c(xmin = 265, ymin = 7080, xmax = 295, ymax = 7110)
loc_index = which(bbox[1] <= coords[, 1] & coords[, 1] <= bbox[3] &
                    bbox[2] <= coords[, 2] & coords[, 2] <= bbox[4])
radar$data = radar$data[, loc_index]
radar$coords = radar$coords[loc_index, ]
coords = coords[loc_index, ]
n_loc = nrow(coords)

# ==============================================================================
# Transform the data to Laplace margins
# ==============================================================================
source(file.path(here::here(), "exec", "model-real-data", "extract-extremes-funcs.R"))

cl = parallel::makeForkCluster(num_cores)
transformed_data = pbapply::pblapply(
  X = 1:n_loc,
  cl = cl,
  FUN = function(i) {
    F = suppressWarnings({
      sliding_window_marginal(radar$data, radar$coords, radar$coords[i, ], radius = r, p = p_gpd)
    })
    tmp = F(radar$data[, i], remove_zeros = TRUE)
    # Ensure that we don't get infinities
    if (any(tmp == 1, na.rm = TRUE)) {
      tmp[which(tmp == 1)] = (1 + max(tmp[which(tmp != 1)])) / 2
    }
    qlaplace(tmp)
  })
parallel::stopCluster(cl)
transformed_data = do.call(cbind, transformed_data)

# ==============================================================================
# Pairwise log likelihood
# ==============================================================================
matern_corr = function(dist, ρ, ν = 1.5) {
  κ = sqrt(8 * ν) / ρ
  res = 2 ^ (1 - ν) / gamma(ν) * (κ * dist) ^ ν * besselK(κ * dist, ν)
  res[dist == 0] = 1
  res
}

ll_pairwise = function(θ, y, y0, dist, α = NULL, σ = NULL, β = NULL, γ = NULL, ρ = NULL) {
  i = 1
  if (is.null(α)) {
    α = θ[i]
    i = i + 1
  }
  if (is.null(σ)) {
    σ = θ[i]
    i = i + 1
  }
  if (is.null(β)) {
    β = θ[i]
    i = i + 1
  }
  if (is.null(γ)) {
    γ = θ[i]
    i = i + 1
  }
  if (is.null(ρ)) {
    ρ = θ[i]
    i = i + 1
  }
  if (σ < 0) return(-1e199)
  if (σ > 20) return(-1e199)
  if (ρ < 0) return(-1e199)
  if (ρ > 200) return(-1e199)
  corr = matern_corr(dist, ρ)
  sum(dnorm(y, mean = α * y0 + γ, sd = sqrt(1 - corr^2) * σ * y0^β, log = TRUE), na.rm = TRUE)
}

# ==============================================================================
# Fit pairwise models
# ==============================================================================
n_pairs = 5000
cl = parallel::makeForkCluster(num_cores, mc.set.seed = TRUE)
parallel::clusterSetRNGStream(cl, 1)
pars1 = pbapply::pblapply(
  X = 1:n_pairs,
  cl = cl,
  FUN = function(i) {
    loc_index = sample.int(n_loc, 2)
    x = as.numeric(transformed_data[, loc_index[1]])
    y = as.numeric(transformed_data[, loc_index[2]])
    na_index = which(is.na(x) | is.na(y))
    x_big_index = which(x > threshold)
    x_big_index = setdiff(x_big_index, na_index)
    if (length(x_big_index) == 0) return(NULL)
    dist = as.numeric(dist(coords[loc_index, 1:2]))
    p1 = optim(
      par = c(.5, 1),
      fn = ll_pairwise,
      control = list(fnscale = -1),
      y = y[x_big_index],
      γ = 0,
      ρ = 0.0001,
      β = 0,
      dist = dist,
      y0 = x[x_big_index])$par
    p2 = optim(
      par = c(.5, 1, .5),
      fn = ll_pairwise,
      control = list(fnscale = -1),
      ρ = 0.0001,
      γ = 0,
      dist = dist,
      y = y[x_big_index],
      y0 = x[x_big_index])$par
    p3 = optim(
      par = c(.5, 1, 0),
      fn = ll_pairwise,
      control = list(fnscale = -1),
      β = 0,
      ρ = 0.0001,
      dist = dist,
      y = y[x_big_index],
      y0 = x[x_big_index])$par
    p4 = optim(
      par = c(.5, 1, .5, 0),
      fn = ll_pairwise,
      control = list(fnscale = -1),
      ρ = 0.0001,
      dist = dist,
      y = y[x_big_index],
      y0 = x[x_big_index])$par
    p5 = optim(
      par = c(.5, 1, .5, 30),
      fn = ll_pairwise,
      control = list(fnscale = -1),
      γ = 0,
      dist = dist,
      y = y[x_big_index],
      y0 = x[x_big_index])$par
    #p6 = optim(
    #  par = c(.5, 1, .5, 0, 0),
    #  fn = ll_pairwise,
    #  control = list(fnscale = -1),
    #  y = y[x_big_index],
    #  y0 = x[x_big_index])$par
    data.frame(
      value = c(p1, p2, p3, p4, p5),# p6),
      dist = dist,
      iter = i,
      n_y0 = length(x_big_index),
      s0_index = loc_index[1],
      par = c("α", "σ",
              "α", "σ", "β",
              "α", "σ", "γ",
              "α", "σ", "β", "γ",
              "α", "σ", "β", "ρ"),
              #"α", "σ", "β", "μ",
              #"α", "σ", "β", "γ", "μ"),
      model = c(1, 1,
                2, 2, 2,
                3, 3, 3,
                4, 4, 4, 4,
                5, 5, 5, 5),
                #6, 6, 6, 6, 6),
      s_index = loc_index[2])
  })
parallel::stopCluster(cl)

plots = list()

plots$pairwise1 = do.call(rbind, pars1) |>
  dplyr::mutate(par = factor(
    par,
    levels = c("α", "σ", "β", "γ", "ρ"),
    labels = c("alpha", "sigma", "beta", "gamma", "rho"))) |>
  ggplot() +
  geom_point(aes(x = dist, y = value), alpha = .1) +
  facet_grid(par ~ model, scales = "free") +
  geom_smooth(aes(x = dist, y = value))

plots$pairwise1_improved = do.call(rbind, pars1) |>
  dplyr::filter(model != 5) |>
  dplyr::mutate(par = factor(
    par,
    levels = c("α", "σ", "β", "γ", "ρ"),
    labels = c("alpha", "sigma", "beta", "gamma", "rho"))) |>
  ggplot() +
  geom_point(aes(x = dist, y = value), alpha = .1) +
  facet_grid(par ~ model, scales = "free") +
  geom_smooth(aes(x = dist, y = value))
# There is a clear trend in α, but it is fighting a bit with γ.
# β and σ are also fighting a bit. We should start by estimating α
# and see what happens
# We seem to be unable to estimate ρ using the pairwise data

# =====================================================================
# Estimate α and σ
# =====================================================================

x_coords = coords[, 1] |> unique() |> sort()
y_coords = coords[, 2] |> unique() |> sort()
Δs0 = 1
s0_locs = expand.grid(
  x = x_coords[seq(Δs0, length(x_coords), by = Δs0)],
  y = y_coords[seq(Δs0, length(y_coords), by = Δs0)]) |>
  as.matrix()
n_s0 = nrow(s0_locs)
s0_index = lapply(
  X = 1:n_s0,
  FUN = function(i) {
    which(coords[, 1] == s0_locs[i, 1] & coords[, 2] == s0_locs[i, 2])
  })
s0_index = s0_index[sapply(s0_index, length) > 0] |>
  unlist() |>
  unname()
#plot(coords)
#points(s0_locs, col = "red")

data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n = 1,
  r = Inf)

get_a_func = function(θ) {
  λ = exp(θ[1]); κ = exp(θ[2])
  function(y, dist) {
    α = exp(- (dist / λ)^κ)
    matrix(rep(y, each = length(α)) * rep(α, length(y)),
           nrow = length(dist),
           ncol = length(y))
  }
}

get_σ_func = function(θ) {
  σ0 = exp(θ[1])
  ρ = exp(θ[2])
  function(dist, n) {
    σ = sqrt(1 - exp(-2 * dist / ρ)) * σ0
    matrix(rep(σ, n), nrow = length(dist), ncol = n)
  }
}

ll = function(θ, data) {
  a_func = get_a_func(θ[1:2])
  σ_func = get_σ_func(θ[3:4])
  res = 0
  for (i in seq_along(data$y)) {
    a = a_func(data$y0[[i]], data$dist_to_s0[[i]])
    σ = σ_func(data$dist_to_s0[[i]], length(data$y0[[i]]))
    res = res + sum(dnorm(data$y[[i]], a, σ, log = TRUE), na.rm = TRUE)
  }
  res
}

res = optim(
  par = c(log(50), log(1), log(1), log(20)),
  fn = ll,
  control = list(fnscale = -1, trace = 6),
  data = data)
est_par = res$par

a_func = get_a_func(est_par[1:2])
σ_func = get_σ_func(est_par[3:4])

xx = seq(0, max(sapply(pars1, `[[`, "dist")), length.out = 1000)
plots$alpha_and_sigma = do.call(rbind, pars1) |>
  dplyr::filter(model == 1) |>
  dplyr::mutate(par = factor(
    par,
    levels = c("α", "σ", "β", "γ", "ρ"),
    labels = c("alpha", "sigma", "beta", "gamma", "rho"))) |>
  ggplot() +
  geom_point(aes(x = dist, y = value), alpha = .5) +
  facet_wrap(~par, scales = "free", ncol = 1) +
  geom_line(data = data.frame(x = xx, y = a_func(1, xx), par = "alpha"),
            aes(x = x, y = y), col = "blue") +
  geom_line(data = data.frame(x = xx, y = σ_func(xx, 1), par = "sigma"),
            aes(x = x, y = y), col = "blue")

plot_df1 = list()
for (i in seq_along(data$y)) {
  a = a_func(data$y0[[i]], data$dist_to_s0[[i]])
  σ = σ_func(data$dist_to_s0[[i]], data$n[i])
  residuals = (data$y[[i]] - a) / σ
  plot_df1[[i]] = data.frame(
    y = as.numeric(residuals),
    dist = rep(data$dist_to_s0[[i]], ncol(residuals)),
    xcoord = rep(coords[data$loc_index[[i]], 1], ncol(residuals)),
    ycoord = rep(coords[data$loc_index[[i]], 2], ncol(residuals)))
}
plot_df1 = do.call(rbind, plot_df1)
ss = sample.int(nrow(plot_df1), 1e6)
plot_df1 = plot_df1[ss, ]

num_groups = 10
ii = round(seq(0, nrow(plot_df1), length = num_groups + 1), digits = 0)
dist_lims = sort(plot_df1$dist)[ii]
ff = function(x, lims) {
  res = x
  lims = sort(lims)
  for (i in seq_along(lims)) {
    res[x <= lims[length(lims) + 1 - i]] = length(lims) + 1 - i
  }
  res
}
plots$residuals = plot_df1 |>
  dplyr::mutate(dist_group = ff(dist, dist_lims)) |>
  dplyr::group_by(dist_group) |>
  dplyr::summarise(
    var = var(y, na.rm = TRUE),
    skew = moments::skewness(y, na.rm = TRUE),
    mean = mean(y, na.rm = TRUE)) |>
  tidyr::pivot_longer(c(var, skew, mean)) |>
  ggplot() +
  geom_line(aes(x = dist_group, y = value)) +
  facet_wrap(~name, scales = "free")

# =====================================================================
# Pairwise estimates with fixed α, σ
# =====================================================================
n_pairs = 5000
cl = parallel::makeForkCluster(num_cores, mc.set.seed = TRUE)
parallel::clusterSetRNGStream(cl, 1)
pars2 = pbapply::pblapply(
  X = 1:n_pairs,
  cl = cl,
  FUN = function(i) {
    loc_index = sample.int(n_loc, 2)
    x = as.numeric(transformed_data[, loc_index[1]])
    y = as.numeric(transformed_data[, loc_index[2]])
    na_index = which(is.na(x) | is.na(y))
    x_big_index = which(x > threshold)
    x_big_index = setdiff(x_big_index, na_index)
    if (length(x_big_index) == 0) return(NULL)
    dist = as.numeric(dist(coords[loc_index, 1:2]))
    α = as.numeric(a_func(1, dist))
    σ = as.numeric(σ_func(dist, 1))
    p1 = optim(
      par = .5,
      fn = ll_pairwise,
      control = list(fnscale = -1),
      method = "Brent", lower = -1, upper = 2,
      α = α,
      γ = 0,
      ρ = .0001,
      σ = σ,
      dist = dist,
      y = y[x_big_index],
      y0 = x[x_big_index])$par
    p2 = optim(
      par = 0,
      fn = ll_pairwise,
      control = list(fnscale = -1),
      method = "Brent", lower = -5, upper = 5,
      α = α,
      ρ = .0001,
      σ = σ,
      β = 0,
      dist = dist,
      y = y[x_big_index],
      y0 = x[x_big_index])$par
    data.frame(
      value = c(p1, p2),
      dist = dist,
      iter = i,
      n_y0 = length(x_big_index),
      s0_index = loc_index[1],
      par = c("β", "γ"),
      s_index = loc_index[2])
  })
parallel::stopCluster(cl)

plots$pairwise2 = do.call(rbind, pars2) |>
  dplyr::mutate(par = factor(
    par,
    levels = c("α", "σ", "β", "γ", "ρ"),
    labels = c("alpha", "sigma", "beta", "gamma", "rho"))) |>
  ggplot() +
  geom_point(aes(x = dist, y = value), alpha = .1) +
  geom_smooth(aes(x = dist, y = value)) +
  facet_wrap(~par, scales = "free")

jpeg()
for (name in names(plots)) {
  print(plots[[name]] + labs(title = name))
}
dev.off()

# ========================================================
# Plots for the paper
# ========================================================

plot1 = do.call(rbind, pars1) |>
  dplyr::filter(model != 5) |>
  dplyr::mutate(
    value = ifelse(par == "σ", pmin(value, 8), value),
    value = ifelse(par == "α", pmin(value, 2), value),
    value = ifelse(par == "α", pmax(value, -2), value),
    value = ifelse(par == "γ", pmax(value, -10), value),
    value = ifelse(par == "γ", pmin(value, 10), value),
    model = paste("Model", model),
    par = factor(
      par,
      levels = c("α", "σ", "β", "γ", "ρ"),
      labels = paste0("$\\", c("alpha", "sigma", "beta", "gamma", "rho"), "$"))) |>
  ggplot() +
  geom_point(aes(x = dist, y = value), alpha = .05) +
  facet_grid(par ~ model, scales = "free") +
  theme_light() +
  theme(text = element_text(size = 15)) +
  theme(strip.text.y = element_text(angle = 0, size = 15, colour = "black"),
        strip.text.x = element_text(size = 12, colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
  labs(x = "Distance [km]", y = "Value")

xx = seq(0, max(sapply(pars1, `[[`, "dist")), length.out = 1000)
plot2_1 = do.call(rbind, pars1) |>
  dplyr::filter(model == 1) |>
  dplyr::mutate(par = factor(
    par,
    levels = c("α", "σ", "β", "γ", "ρ"),
    labels = paste0("$\\", c("alpha", "sigma", "beta", "gamma", "rho"), "$"))) |>
  ggplot() +
  geom_point(aes(x = dist, y = value), alpha = .1) +
  facet_wrap(~par, scales = "free", nrow = 1) +
  geom_line(data = data.frame(x = xx, y = a_func(1, xx), par = "$\\alpha$"),
            aes(x = x, y = y), col = "blue", size = 2) +
  geom_line(data = data.frame(x = xx, y = σ_func(xx, 1), par = "$\\sigma$"),
            aes(x = x, y = y), col = "blue", size = 2) +
  theme_light() +
  theme(text = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15, colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
  labs(x = "Distance [km]", y = "Value")

plot2_2 = do.call(rbind, pars2) |>
  dplyr::mutate(par = factor(
    par,
    levels = c("α", "σ", "β", "γ", "ρ"),
    labels = paste0("$\\", c("alpha", "sigma", "beta", "gamma", "rho"), "$"))) |>
  ggplot() +
  geom_point(aes(x = dist, y = value), alpha = .1) +
  facet_wrap(~par, scales = "free") +
  theme_light() +
  theme(text = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15, colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0")) +
  labs(x = "Distance [km]", y = "Value")

plot2 = patchwork::wrap_plots(plot2_1, plot2_2, nrow = 1)

tikz_plot(file.path(image_dir(), "model-selection1.pdf"),
          plot1, width = 7, height = 7)

tikz_plot(file.path(image_dir(), "model-selection2.pdf"),
          plot2, width = 10, height = 5)


library(magick)
image = magick::image_read_pdf(file.path(image_dir(), "model-selection1.pdf"))
magick::image_write(image, file.path(image_dir(), "model-selection1.jpg"), format = "jpg", quality = 50)
image = magick::image_read_pdf(file.path(image_dir(), "model-selection2.pdf"))
magick::image_write(image, file.path(image_dir(), "model-selection2.jpg"), format = "jpg", quality = 50)
