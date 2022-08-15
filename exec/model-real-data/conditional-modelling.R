devtools::load_all()
library(sf)
library(scoringRules)
library(ggplot2)
library(dplyr)
library(tidyr)
library(INLA)
library(txtplot)

INLA::inla.setOption(pardiso.license = "~/.R/licences/pardiso.lic")

filename = file.path(tmp_dir(), "conditional-modelling.rds")
log_score_filename = file.path(tmp_dir(), "conditional-modelling-log-score.rds")

# ==============================================================================
# Fix necessary model parameters:
# ==============================================================================
p_gpd = 1
r = 5
threshold = qlaplace(.999)

# How many cores to use in parallel?
num_cores = 15
verbose = TRUE

# ==============================================================================
# Load the data
# ==============================================================================
radar = readRDS(file.path(downloads_dir(), "radar-trondheim.rds"))
coords = st_coordinates(radar$coords)

# ==============================================================================
# Restrict the data, for computational reasons
# ==============================================================================
bbox = c(xmin = 265, ymin = 7080, xmax = 295, ymax = 7110)
loc_index = which(bbox[1] <= coords[, 1] & coords[, 1] <= bbox[3] &
                    bbox[2] <= coords[, 2] & coords[, 2] <= bbox[4])
data = radar$data[, loc_index]
coords = st_coordinates(radar$coords)[loc_index, ]
n_loc = nrow(coords)
rm(radar)

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
      sliding_window_marginal(data, coords, coords[i, ], radius = r, p = p_gpd)
    })
    tmp = F(data[, i], remove_zeros = TRUE)
    # Ensure that we don't get infinities
    if (any(tmp == 1, na.rm = TRUE)) {
      tmp[which(tmp == 1)] = (1 + max(tmp[which(tmp != 1)])) / 2
    }
    qlaplace(tmp)
  })
parallel::stopCluster(cl)
transformed_data = do.call(cbind, transformed_data)

# ==============================================================================
# Place s0
# ==============================================================================
x_coords = coords[, 1] |> unique() |> sort()
y_coords = coords[, 2] |> unique() |> sort()
Δs0 = 4
s0_locs = expand.grid(
  x = x_coords[seq(Δs0, length(x_coords), by = Δs0)],
  y = y_coords[seq(Δs0, length(y_coords), by = Δs0)]) |>
  as.matrix()
s0_index = lapply(
  X = 1:nrow(s0_locs),
  FUN = function(i) {
    which(coords[, 1] == s0_locs[i, 1] & coords[, 2] == s0_locs[i, 2])
  })
s0_index = s0_index[sapply(s0_index, length) > 0] |>
  unlist() |>
  unname()
#plot(coords)
#points(s0_locs, col = "red")

# ==============================================================================
# Extract the extreme data
# ==============================================================================
my_n = c(1, 2, 4, 6, 8, 16, 32)
data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n = my_n,
  r = cumsum(4 * my_n))
  #n = 1,
  #r = Inf)

sum(data$n)
data$loc_index |> sapply(length) |> summary()
data$time_index |> unlist() |> table() |> as.numeric() |> summary()
# There are definitely some cases of multiple exceedences at the same time
mean(is.na(unlist(data$y)))
# a lot of the data are NAs

plot = local({
  i = 12
  df = as.data.frame(coords[data$loc_index[[i]], ]) |>
    dplyr::mutate(tag = "Used")
  df2 = as.data.frame(coords[-c(data$loc_index[[i]], data$s0_index[[i]]), ]) |>
    dplyr::mutate(tag = "Not used")
  s0_df = as.data.frame(coords[data$s0_index[[i]], , drop = FALSE]) |>
    dplyr::mutate(tag = "$\\bm s_0$")
  rbind(df, df2, s0_df) |>
    dplyr::mutate(tag = factor(tag, levels = c("$\\bm s_0$", "Used", "Not used"))) |>
    ggplot() +
    geom_point(aes(x = X, y = Y, shape = tag, size = tag)) +
    #scale_color_manual(values = c("blue", "red", "black")) +
    scale_shape_manual(values = c(8, 19, 18)) +
    scale_size_manual(values = c(2, 2, .3)) +
    theme_light() +
    guides(shape = "none", size = "none") +
    labs(x = "Easting", y = "Northing") +
    theme(axis.ticks = element_blank(), axis.text = element_blank())
})
tikz_plot(file.path(image_dir(), "design-of-experiment.pdf"),
          plot, width = 4, height = 4)



# ==============================================================================
# Create the mesh and spde
# ==============================================================================
mesh = inla.mesh.2d(
  loc = coords,
  boundary = list(inla.nonconvex.hull(coords, convex = -.1),
                  inla.nonconvex.hull(coords, convex = -1)),
  max.edge = c(5, 30))
mesh$n
plot(mesh)
points(coords)

data$dist_to_s0_from_mesh = list()
for (i in seq_along(data$s0)) {
  data$dist_to_s0_from_mesh[[i]] = dist_euclid(mesh$loc[, 1:2], data$s0[[i]])
}

spde = inla.spde2.pcmatern(mesh, prior.range = c(60, .95), prior.sigma = c(5, .05))

plot = local({
  i = 12
  df = as.data.frame(coords[data$loc_index[[i]], ]) |>
    dplyr::mutate(tag = "Used")
  df2 = as.data.frame(coords[-c(data$loc_index[[i]], data$s0_index[[i]]), ]) |>
    dplyr::mutate(tag = "Not used")
  s0_df = as.data.frame(coords[data$s0_index[[i]], , drop = FALSE]) |>
    dplyr::mutate(tag = "$\\bm s_0$")
  rbind(df, df2, s0_df) |>
    dplyr::mutate(tag = factor(tag, levels = c("$\\bm s_0$", "Used", "Not used"))) |>
    ggplot() +
    gg(mesh, int.color = "grey", ext.color = "grey") +
    geom_point(aes(x = X, y = Y, shape = tag, size = tag, color = tag)) +
    scale_color_manual(values = c("red", "black", "black")) +
    scale_shape_manual(values = c(19, 19, 19)) +
    scale_size_manual(values = c(2.4, 2.4, .5)) +
    #scale_size_manual(values = c(1.5, 1.5, .7)) +
    theme_light() +
    guides(shape = "none", size = "none", color = "none") +
    labs(x = "Easting", y = "Northing") +
    theme(axis.ticks = element_blank(), axis.text = element_blank(),
          text = element_text(size = 23))
})
tikz_plot(file.path(image_dir(), "design-of-experiment2.pdf"),
          plot, width = 12, height = 12)

# ==============================================================================
# Compute the MLE and use it for initial values in R-INLA
# ==============================================================================
get_a_func = function(θ) {
  λ = exp(θ[1])
  κ = exp(θ[2])
  function(y, dist) {
    α = exp(- (dist / λ)^κ)
    matrix(rep(y, each = length(α)) * rep(α, length(y)),
           nrow = length(dist),
           ncol = length(y))
  }
}
get_b_func = function(θ) {
  ρ = exp(θ)
  function(y, dist) {
    tmp = dist / ρ
    tmp[tmp < 1e-9] = 1e-9
    b = sqrt(1 - exp(-2 * tmp))
    matrix(rep(b, length(y)), nrow = length(dist), ncol = length(y))
  }
}
ll = function(θ, sum_terms = TRUE, verbose = FALSE, ρ2 = NULL, ...) {
  if (!is.null(ρ2)) {
    b_func = get_b_func(log(ρ2))
    Σ_func = get_Σ_func(θ[3:4], spde)$value
    τ = exp(θ[5])
  } else {
    b_func = get_b_func(θ[3])
    Σ_func = get_Σ_func(θ[4:5], spde)$value
    τ = exp(θ[6])
  }
  res = ll_conditional(
    a_func = get_a_func(θ[1:2]),
    b_func = b_func,
    Σ_func = Σ_func,
    no_beta = TRUE,
    τ = τ,
    ...)
  if (sum_terms) res = sum(res)
  if (verbose) message("Done")
  res
}

#est = list(par = c(4, -.4, 1.8, 4, .5, 4), convergence = 1)
#est = list(par = c(3.4, -.5, -1.3, 2.6, .6, 2.8), convergence = 1)
est = list(par = c(3.39, -.56, -1.24, 2.56, .62, 3.08), convergence = 1)
while (est$convergence != 0) {
  est = optim(
    par = est$par,
    fn = ll,
    y = data$y,
    y0 = data$y0,
    dist = data$dist_to_s0,
    dist_from_mesh = data$dist_to_s0_from_mesh,
    A = lapply(data$loc_index, function(x) inla.spde.make.A(mesh, coords[x, ])),
    num_cores = num_cores,
    control = list(fnscale = -1, trace = 6, REPORT = 1))
}

# ==============================================================================
# Create all the model components for R-INLA
# ==============================================================================

y_inla = unlist(data$y)
na_index = which(is.na(y_inla))

dist_to_s0_inla = data$dist_to_s0[rep(seq_along(data$n), data$n)]
y0_inla = rep(unlist(data$y0), sapply(dist_to_s0_inla, length))
dist_to_s0_inla = unlist(dist_to_s0_inla)

y_inla = y_inla[-na_index]
y0_inla = y0_inla[-na_index]
dist_to_s0_inla = dist_to_s0_inla[-na_index]

alpha_priors = list(lambda = c(1, 3, 3), kappa = c(1, -.4, 3))
alpha_model = alpha_generic_model(
  y0 = y0_inla,
  dist_to_s0 = dist_to_s0_inla,
  init = c(4, -.4),
  priors = alpha_priors)

spde_priors = list(
  rho = c(1, 60, .95),
  sigma = c(1, 5, .05),
  rho_2 = c(1, log(6), 2))
spde_model = spde_generic_model_with_b_func(
  spde = spde,
  n = data$n,
  init = c(4, .5, 1.8),
  priors = spde_priors,
  dist_to_s0 = do.call(rbind, data$dist_to_s0_from_mesh))

effects = list(
  spatial = seq_len(sum(data$n) * mesh$n),
  idx = seq_along(y_inla))

A_inla = local({
  all_location_indices = rep(data$loc_index, data$n)
  all_repl_indices = rep(seq_along(all_location_indices), sapply(all_location_indices, length))
  all_location_indices = unlist(all_location_indices)
  A = inla.spde.make.A(
    mesh,
    loc = coords,
    index = all_location_indices,
    repl = all_repl_indices)
  A = A[-na_index, ]
  A
})

stack = inla.stack(
  data = list(y = y_inla),
  A = list(spatial = A_inla, 1),
  effects = effects)
rm(A_inla, effects)

formula = y ~ -1 +
  f(spatial, model = spde_model) +
  f(idx, model = alpha_model)

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
fit$model.matrix = NULL
fit$offset.linear.predictor = NULL

# ==============================================================================
# Fit a model using only one conditioning site
# ==============================================================================

set.seed(123)
single_site_s0s = c(data$s0_index[which.max(data$n)], sample.int(nrow(coords), 9))
# local({
#   coords = as.data.frame(coords)
#   coords$tag = FALSE
#   coords$tag[single_site_s0s] = TRUE
#   ggplot(coords) +
#     geom_point(aes(x = X, y = Y, col = tag))
# })

single_site_data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = single_site_s0s,
  threshold = threshold,
  n = my_n,
  r = cumsum(4 * my_n))
single_site_data$dist_to_s0_from_mesh = list()
for (i in seq_along(single_site_data$s0)) {
  single_site_data$dist_to_s0_from_mesh[[i]] = dist_euclid(mesh$loc[, 1:2], single_site_data$s0[[i]])
}

single_site_fits = vector("list", length(single_site_s0s))
pb = progress_bar(length(single_site_s0s))
for (i in seq_along(single_site_s0s)) {
  y_mini = as.numeric(single_site_data$y[[i]])
  na_index = which(is.na(y_mini))
  mini_n = single_site_data$n[i]
  dist_to_s0_mini = rep(single_site_data$dist_to_s0[[i]], mini_n)
  y0_mini = rep(single_site_data$y0[[i]], each = length(single_site_data$dist_to_s0[[i]]))
  y_mini = y_mini[-na_index]
  y0_mini = y0_mini[-na_index]
  dist_to_s0_mini = dist_to_s0_mini[-na_index]

  alpha_model_mini = alpha_generic_model(
    y0 = y0_mini,
    dist_to_s0 = dist_to_s0_mini,
    init = c(4, -.4),
    priors = alpha_priors)
  spde_model_mini = spde_generic_model_with_b_func(
    spde = spde,
    n = mini_n,
    init = c(4, .5, 1.8),
    priors = spde_priors,
    dist_to_s0 = matrix(single_site_data$dist_to_s0_from_mesh[[i]], nrow = 1))
  effects_mini = list(
    spatial = seq_len(sum(mini_n) * mesh$n),
    idx = seq_along(y_mini))
  A_inla_mini = local({
    all_location_indices = rep(single_site_data$loc_index[i], mini_n)
    all_repl_indices = rep(seq_along(all_location_indices), sapply(all_location_indices, length))
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
  rm(A_inla_mini, effects_mini)

  formula_mini = y ~ -1 +
    f(spatial, model = spde_model_mini) +
    f(idx, model = alpha_model_mini)

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
# Estimate ρ2 using an independence likelihood
# ==============================================================================
iid_b_model = iid_generic_model_with_b_func(
  init = c(1, 1.5),
  priors = list(sigma = c(1, 5, .05), rho_2 = c(1, 1.5, 3)),
  dist_to_s0 = dist_to_s0_inla)

iid_fit = inla(
  formula = y ~ -1 + f(idx, model = alpha_model) + f(idx2, model = iid_b_model),
  data = data.frame(y = y_inla, idx = seq_along(y_inla), idx2 = seq_along(y_inla)),
  only.hyperparam = TRUE,
  control.inla = list(control.vb = list(enable = FALSE), int.strategy = "eb"),
  verbose = TRUE,
  num.threads = 4,
  inla.mode = "experimental")

log_ρ2 = iid_fit$mode$theta[5]
rm(iid_fit)

est2 = list(par = c(4, -.4, 4, .5, 4), convergence = 1)
while (est2$convergence != 0) {
  est2 = optim(
    par = est2$par,
    fn = ll,
    y = data$y,
    y0 = data$y0,
    dist = data$dist_to_s0,
    dist_from_mesh = data$dist_to_s0_from_mesh,
    A = lapply(data$loc_index, function(x) inla.spde.make.A(mesh, coords[x, ])),
    ρ2 = exp(log_ρ2),
    num_cores = num_cores,
    control = list(fnscale = -1, trace = 6, REPORT = 1))
}

# ==============================================================================
# Fit the global and single-site models using fixed ρ2
# ==============================================================================

spde_priors$rho_2 = c(0, log_ρ2)
spde_model = spde_generic_model_with_b_func(
  spde = spde,
  n = data$n,
  init = c(4, .5),
  priors = spde_priors,
  dist_to_s0 = do.call(rbind, data$dist_to_s0_from_mesh))

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

fits = list(
  fit = c(list(fit), single_site_fits, fit_fixed),
  global = c(TRUE, rep(FALSE, length(single_site_fits)), TRUE),
  fixed = c(rep(FALSE, length(single_site_fits) + 1), TRUE),
  threshold = threshold,
  log_ρ2 = log_ρ2,
  single_site_s0s = single_site_s0s)
saveRDS(fits, filename)

# ===============================================
# Load the model fits and evaluate them
# ===============================================
fits = readRDS(filename)
log_ρ2 = fits$log_ρ2

eval_data = extract_extreme_fields(
  data = transformed_data,
  coords = coords,
  s0_index = seq_len(n_loc)[-s0_index],
  threshold = threshold,
  n = 1,
  r = Inf)

eval_data$dist_to_s0_from_mesh = list()
for (i in seq_along(eval_data$s0)) {
  eval_data$dist_to_s0_from_mesh[[i]] = dist_euclid(mesh$loc[, 1:2], eval_data$s0[[i]])
}

# ===============================================================
# Compute H, J and C
# ===============================================================

ll_grad = function(θ, ...) {
  numDeriv::jacobian(ll, θ, ...)
}

A = inla.spde.make.A(mesh, coords)
A_mats = lapply(data$loc_index, function(x) A[x, ])
A_mats_single_site = lapply(single_site_data$loc_index, function(x) A[x, ])

θ_reordering1 = c(5, 6, 4, 2, 3, 1)
θ_reordering2 = c(4, 5, 2, 3, 1)
θ_names1 = c("log_lambda", "log_kappa", "log_rho_2", "log_rho", "log_sigma", "log_tau")
θ_tex_names1 = paste0("$\\", c("lambda", "kappa", "rho_2", "rho", "sigma", "tau"), "$")
θ_names2 = c("log_lambda", "log_kappa", "log_rho", "log_sigma", "log_tau")
θ_tex_names2 = paste0("$\\", c("lambda", "kappa", "rho", "sigma", "tau"), "$")

fits$reordering = fits$names = fits$tex_names = list()
for (i in seq_along(fits$fit)) {
  if (fits$fixed[i]) {
    fits$reordering[[i]] = θ_reordering2
    fits$names[[i]] = θ_names2
    fits$tex_names[[i]] = θ_tex_names2
  } else {
    fits$reordering[[i]] = θ_reordering1
    fits$names[[i]] = θ_names1
    fits$tex_names[[i]] = θ_tex_names1
  }
}

fits$C = list()
pb = progress_bar(length(fits$fit))
for (i in seq_along(fits$fit)) {
  if (fits$global[i]) {
    mydata = data
    my_A = A_mats
  } else {
    mydata = lapply(single_site_data, `[`, i - 1)
    my_A = A_mats_single_site[i - 1]
  }

  H = solve(fits$fit[[i]]$misc$cov.intern)[fits$reordering[[i]], fits$reordering[[i]]]
  grads = ll_grad(
    θ = fits$fit[[i]]$mode$theta[fits$reordering[[i]]],
    y = mydata$y,
    y0 = mydata$y0,
    dist = mydata$dist_to_s0,
    ρ2 = if (fits$fixed[i]) exp(log_ρ2) else NULL,
    A = my_A,
    dist_from_mesh = mydata$dist_to_s0_from_mesh,
    num_cores = num_cores,
    sum_terms = FALSE)

  times = unlist(mydata$time_index)
  window_radius = 5
  J = 0
  for (j in seq_along(times)) {
    index = which(times >= times[j] - window_radius & times <= times[j] + window_radius)
    for (k in index) {
      J = J + grads[j, ] %*% grads[k, , drop = FALSE]
    }
  }

  fits$C[[i]] = get_C(H, J)
  pb$tick()
}

saveRDS(fits, filename)

# ===============================================================
# Plot the corrected hyperparameter posteriors
# ===============================================================

fits = readRDS(filename)
log_ρ2 = fits$log_ρ2

n = 1e5
plot_df = local({
  df = list()
  for (i in seq_along(fits$fit)) {
    θ_uncorrected = inla.hyperpar.sample(n, fits$fit[[i]], intern = TRUE)[, fits$reordering[[i]]]
    for (k in 1:2) {
      if (k == 1) {
        θ_corrected = θ_uncorrected
      #} else if (!fits$global[i]) {
      #  next
      } else {
        θ_corrected = matrix(
          rep(fits$fit[[i]]$mode$theta[fits$reordering[[i]]], each = n),
          nrow = n)
        θ_corrected = θ_corrected + (θ_uncorrected - θ_corrected) %*% t(fits$C[[i]])
      }
      θ_corrected = exp(θ_corrected)
      my_names = fits$tex_names[[i]]

      df[[length(df) + 1]] = as.data.frame(θ_corrected)
      names(df[[length(df)]])[seq_along(fits$reordering[[i]])] = my_names
      df[[length(df)]] = tidyr::pivot_longer(df[[length(df)]], everything()) |>
        dplyr::mutate(k = k, fixed_ρ2 = fits$fixed[i], i = i) |>
        dplyr::mutate(name = factor(name, levels = my_names))
    }
  }
  do.call(rbind, df) |>
    dplyr::mutate(k = factor(k), i = factor(i)) |>
    dplyr::mutate(
      name = sub("rho_2", "rho_b", name),
      tag = case_when(
        i == 1 & k == 1 ~ 1,
        i == 1 & k == 2 ~ 2,
        i == 12 & k == 1 ~ 3,
        i == 12 & k == 2 ~ 4,
        k == 1 ~ 5,
        k == 2 ~ 6),
        #TRUE ~ 5),
      group = paste0("i = ", i, ", k = ", k),
      tag = factor(
        tag,
        levels = 1:6,
        labels = c("Unadjusted global", "Adjusted global", "Unadjusted global 2",
                   "Adjusted global 2", "Unadjusted local", "Adjusted local")))
})

plot = plot_df |>
  dplyr::filter(i != 12, !(i > 1 & k == 2)) |>
  #dplyr::mutate(group = factor(group, levels = unique(plot_df$group)[c(3:22, 1, 2)])) |>
  dplyr::mutate(group = factor(group, levels = unique(plot_df$group)[c(3:12, 1, 2)])) |>
  #dplyr::mutate(group = factor(group, levels = unique(plot_df$group)[c(3:12, 1, 2, 13, 14)])) |>
  dplyr::filter(
    !(name == "$\\kappa$" & value > 1.4),
    !(name == "$\\rho_b$" & value > 1.4),
    !(name == "$\\tau$" & value > 80)) |>
  ggplot(mapping = aes(x = value, col = tag, linetype = tag, group = group)) +
  geom_density(size = 1) +
  facet_wrap(~name, scales = "free") +
  scale_color_manual(values = c("black", "black", "gray")) +
  scale_linetype_manual(values = c("longdash", "solid", "dotted")) +
  labs(color = "", y = "Density", x = "Value", linetype = "") +
  theme_light() +
  theme(text = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15, colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))
tikz_plot(file.path(image_dir(), "case-study-densities.pdf"),
          plot, width = 10, height = 6)

dist = seq(0, 50, length.out = 200)
dist = c(5, 10, 20, 40, 60)
n = 1e5
plot_df = local({
  df = list()
  for (i in seq_along(fits$fit)) {
    θ_uncorrected = inla.hyperpar.sample(n, fits$fit[[i]], intern = TRUE)[, fits$reordering[[i]]]
    for (k in 1:2) {
      if (k == 1) {
        θ_corrected = θ_uncorrected
      } else if (!fits$global[i]) {
        next
      } else {
        θ_corrected = matrix(
          rep(fits$fit[[i]]$mode$theta[fits$reordering[[i]]], each = n),
          nrow = n)
        θ_corrected = θ_corrected + (θ_uncorrected - θ_corrected) %*% t(fits$C[[i]])
      }
      θ_corrected = exp(θ_corrected)
      my_names = fits$tex_names[[i]]
      for (j in seq_along(dist)) {
        α = exp(-(dist[j] / θ_corrected[, 1])^θ_corrected[, 2])
        if (fits$fixed[i]) {
          my_ρ = fits$log_ρ2
        } else {
          my_ρ = θ_corrected[, 3]
        }
        b = θ_corrected[, ncol(θ_corrected) - 1] * sqrt(1 - exp(-2 * dist[j] / my_ρ))
        df[[length(df) + 1]] = data.frame(
          value = c(α, b),
          name = rep(c("$\\alpha$", "$\\zeta$"), c(length(α), length(b))),
          dist = dist[j],
          k = k,
          fixed_ρ2 = fits$fixed[i],
          global = fits$global[i],
          i = i)
      }
    }
  }
  do.call(rbind, df) |>
    dplyr::mutate(
      k = factor(k), i = factor(i), name = factor(name),
      tag = case_when(
        i == 1 & k == 1 ~ 1,
        i == 1 & k == 2 ~ 2,
        i != 1 ~ 3),
      group = paste0("i = ", i, ", k = ", k),
      tag = factor(tag, levels = 1:3, labels = c("Unadjusted global", "Adjusted global", "Single site")))
})

plot = plot_df |>
  dplyr::mutate(group = factor(group, levels = unique(plot_df$group)[c(3:12, 1, 2)])) |>
  dplyr::mutate(
    nametag = paste0(
      sub("\\$$", "", name), "(", ifelse(name == "$\\zeta$", "x", paste(dist, "\\text{km}")), ")",
      ifelse(name == "$\\zeta$", ",\\ x > 2 \\text{km}", ""), "$"),
    nametag = factor(nametag, levels = unique(nametag)[c(2, 1, 3:6)]),
    nametag_nr = as.numeric(nametag)) |>
  dplyr::filter(!(nametag_nr == 1 & value > 2.5), !(nametag_nr == 4 & value < .28)) |>
  ggplot(mapping = aes(x = value, col = tag, linetype = tag, group = group)) +
  geom_density(size = 1) +
  facet_wrap(~nametag, scales = "free", nrow = 2) +
  scale_color_manual(values = c("black", "black", "gray")) +
  scale_linetype_manual(values = c("longdash", "solid", "dotted")) +
  labs(color = "", y = "Density", x = "Value", linetype = "") +
  theme_light() +
  theme(text = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15, colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))
tikz_plot(file.path(image_dir(), "case-study-densities2.pdf"),
          plot, width = 10, height = 7)

# ===============================================================
# Compute log-scores
# ===============================================================

n_sample = 500
set.seed(1)
θ_samples = list()
for (j in seq_along(fits$fit)) {
  θ_uncorrected = inla.hyperpar.sample(n_sample, fits$fit[[j]], intern = TRUE)[, fits$reordering[[j]]]
  for (k in 1:2) {
    if (k == 1) {
      θ_corrected = θ_uncorrected
    #} else if (!fits$global[j]) {
    #  next
    } else {
      θ_corrected = matrix(
        rep(fits$fit[[j]]$mode$theta[fits$reordering[[j]]], each = n_sample),
        nrow = n_sample)

      θ_corrected = θ_corrected + (θ_uncorrected - θ_corrected) %*% t(fits$C[[j]])
    }
    θ_samples[[length(θ_samples) + 1]] = list(
      θ = θ_corrected,
      ρ2 = if (fits$fixed[j]) exp(fits$log_ρ2) else NULL,
      model = paste0("nr. ", j, ", ", ifelse(k == 1, "unadjusted", "adjusted")))
  }
}

A_mats = lapply(eval_data$loc_index, function(x) inla.spde.make.A(mesh, coords[x, ]))

pb = progress_bar(n_sample * length(θ_samples))
log_scores = matrix(-Inf, nrow = sum(eval_data$n), ncol = length(θ_samples))
colnames(log_scores) = sapply(θ_samples, `[[`, "model")
attributes(log_scores)$count = 0
attributes(log_scores)$i = NULL
for (i in seq_len(n_sample)) {
  for (j in seq_along(θ_samples)) {
    tmp = ll(
      θ = θ_samples[[j]]$θ[i, ],
      y = eval_data$y,
      y0 = eval_data$y0,
      dist = eval_data$dist_to_s0,
      dist_from_mesh = eval_data$dist_to_s0_from_mesh,
      sum_terms = FALSE,
      A = A_mats,
      num_cores = num_cores,
      ρ2 = θ_samples[[j]]$ρ2)
    for (k in seq_along(tmp)) {
      log_scores[k, j] = log_sum_exp(c(log_scores[k, j], tmp[k]), na.rm = TRUE)
    }
    pb$tick()
  }
  attributes(log_scores)$count = attributes(log_scores)$count + 1
  attributes(log_scores)$i = c(attributes(log_scores)$i, i)
  saveRDS(log_scores, log_score_filename)
}
pb$terminate()

log_scores = readRDS(log_score_filename)

ll_names = colnames(log_scores)
ll_names = sub("nr\\. 1,", "Global,", ll_names)
ll_names = sub("nr\\. 12,", "Global 3,", ll_names)
ll_names = sub("nr\\.", "Single site", ll_names)
#ll_names[-c(1, 13)] = sub(", unadjusted", "", ll_names[-c(1, 13)])
for (i in 2:100) ll_names = sub(i, i - 1, ll_names)
colnames(log_scores) = ll_names

tmp = as.data.frame(log_scores) |>
  dplyr::mutate(index = rep(seq_along(eval_data$n), eval_data$n)) |>
  tidyr::pivot_longer(-index) |>
  dplyr::group_by(index, name) |>
  dplyr::summarise(value = sum(value)) |>
  tidyr::pivot_wider(names_from = name, values_from = value) |>
  dplyr::ungroup() |>
  dplyr::select(-index) |>
  as.matrix()

for (i in seq_len(nrow(tmp))) {
  tmp[i, ] = rank(-tmp[i, ])
}

plot = tmp |>
  as.data.frame() |>
  tidyr::pivot_longer(tidyselect::everything()) |>
  #dplyr::mutate(name = factor(name, levels = ll_names[c(2, 1, 14, 13, 3:12)])) |>
  dplyr::mutate(name = factor(name, levels = ll_names[c(2, 1, 24, 23, 3:22)])) |>
  ggplot() +
  geom_histogram(aes(x = value, y = ..density..), boundary = .5, binwidth = 1) +
  facet_wrap(~name) +
  scale_x_continuous(
    breaks = seq(1, ncol(tmp), by = 2), expand = c(0, 0),
    minor_breaks = seq(2, ncol(tmp), by = 2)) +
  labs(x = "Ranking", y = "Density") +
  theme_light() +
  theme(text = element_text(size = 15)) +
  theme(strip.text = element_text(size = 15, colour = "black"),
        strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"))

tikz_plot(file.path(image_dir(), "log-score-rankings.pdf"),
          plot, width = 12, height = 6)
