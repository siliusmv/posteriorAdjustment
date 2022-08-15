devtools::load_all()
library(INLA)
library(ggplot2)
library(mvtnorm)
library(mvnfast)
library(dplyr)
library(Matrix)
library(sf)

INLA::inla.setOption(pardiso.license = "~/.R/licences/pardiso.lic")

filename = file.path(tmp_dir(), "conditional-adjustment.rds")
overwrite = FALSE
if (!file.exists(filename)) saveRDS(list(), filename)

if (FALSE) {
  tmp = readRDS(filename)
  for (name in names(tmp)) assign(name, tmp[[name]])
}

real_data_filename = file.path(tmp_dir(), "conditional-modelling.rds")
ρ2_filename = file.path(tmp_dir(), "model-selection.rds")
truth_filename = file.path(tmp_dir(), "conditional-simulation-truth.rds")

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

mesh = inla.mesh.2d(
  loc = coords,
  boundary = list(inla.nonconvex.hull(coords, convex = -.1),
                  inla.nonconvex.hull(coords, convex = -1)),
  max.edge = c(5, 30))
spde = inla.spde2.pcmatern(mesh, prior.range = c(60, .95), prior.sigma = c(4, .05))
A = inla.spde.make.A(mesh, coords)

tmp = readRDS(real_data_filename)
ss = tmp$fit[[1]]$summary.hyperpar
τ = ss[1, 1]
ρ = exp(ss[2, 1])
σ = exp(ss[3, 1])
ρ2 = exp(ss[4, 1])
λ = exp(ss[5, 1])
κ = exp(ss[6, 1])
rm(tmp)

tmp = readRDS(ρ2_filename)
ρ2 = exp(tmp$pars$value[which(tmp$pars$name == "log_ρ2")])
rm(tmp)

info = list(
  τ = c(4, 1, 5e-5),
  ρ = c(20, 60, .95),
  σ = c(1, 5, .05),
  ρ2 = c(log(6), log(6), 2),
  λ = c(4, 3, 3),
  κ = c(-.4, -.4, 3))

n_sim = 100
n_cores = 15
threshold = qlaplace(.999)

θ_a = log(c(λ, κ))
θ_b = log(c(ρ2))
θ_Σ = log(c(ρ, σ))
θ = c(θ_a, θ_b, θ_Σ, log(τ))

get_a_func = function(θ) {
  λ = exp(θ[1]); κ = exp(θ[2])
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

a = get_a_func(θ_a)
b = get_b_func(θ_b)
Σ = get_Σ_func(θ_Σ, spde)

A = inla.spde.make.A(mesh, coords)

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
ll_grad = function(θ, sum_terms = TRUE, ...) {
  res = numDeriv::jacobian(ll, θ, ..., sum_terms = FALSE)
  if (sum_terms) res = apply(res, 2, sum)
  res
}

A_mats = lapply(1:n_loc, \(i) A[-i, ])
dist_to_s0 = lapply(1:n_loc, \(i) dist_euclid(coords[-i, ], coords[i, , drop = FALSE]))
dist_to_s0_from_mesh = lapply(1:n_loc, \(i) dist_euclid(mesh$loc[, 1:2], coords[i, , drop = FALSE]))

obs = keef_sampling(
  n = 5e4,
  a_func = a,
  b_func = b,
  Q = Σ$Q,
  τ = τ,
  threshold = threshold,
  dist_to_s0 = dist_to_s0,
  dist_to_s0_from_mesh = dist_to_s0_from_mesh,
  A = A_mats,
  verbose = TRUE)

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

source(file.path(here::here(), "exec", "model-real-data", "extract-extremes-funcs.R"))

my_n = c(1, 2, 4, 6, 8, 16, 32)

data = extract_extreme_fields(
  data = t(obs),
  coords = coords,
  s0_index = s0_index,
  threshold = threshold,
  n = my_n,
  r = cumsum(4 * my_n))
data$dist_to_s0_from_mesh = list()
for (i in seq_along(data$s0)) {
  data$dist_to_s0_from_mesh[[i]] = dist_euclid(mesh$loc[, 1:2], data$s0[[i]])
}

est = list(par = θ, convergence = 1)
while (est$convergence != 0) {
  est = optim(
    par = est$par,
    fn = ll,
    y = data$y,
    y0 = data$y0,
    dist = data$dist_to_s0,
    #ρ2 = ρ2,
    dist_from_mesh = data$dist_to_s0_from_mesh,
    A = lapply(data$loc_index, function(x) A[x, ]),
    num_cores = n_cores,
    control = list(fnscale = -1, maxit = 300, trace = 6, REPORT = 1))
}

saveRDS(est, truth_filename)

num_samples = 300
pb = progress_bar(num_samples)
n1 = 25
n1 = 60
for (i in 1:num_samples) {

  seed = .Random.seed

  obs = keef_sampling(
    n = n1,
    a_func = a,
    b_func = b,
    Q = Σ$Q,
    τ = τ,
    threshold = threshold,
    dist_to_s0 = dist_to_s0,
    dist_to_s0_from_mesh = dist_to_s0_from_mesh,
    A = A_mats)

  data = extract_extreme_fields(
    data = t(obs),
    coords = coords,
    s0_index = s0_index,
    threshold = threshold,
    n = my_n,
    r = cumsum(4 * my_n))
  data$dist_to_s0_from_mesh = list()
  for (j in seq_along(data$s0)) {
    data$dist_to_s0_from_mesh[[j]] = dist_euclid(mesh$loc[, 1:2], data$s0[[j]])
  }

  y_inla = unlist(data$y)
  dist_to_s0_inla = data$dist_to_s0[rep(seq_along(data$n), data$n)]
  y0_inla = rep(unlist(data$y0), sapply(dist_to_s0_inla, length))
  dist_to_s0_inla = unlist(dist_to_s0_inla)

  alpha_priors = list(lambda = c(1, info$λ[-1]), kappa = c(1, info$κ[-1]))
  alpha_model = alpha_generic_model(
    y0 = y0_inla,
    dist_to_s0 = dist_to_s0_inla,
    init = c(info$λ[1], info$κ[1]),
    priors = alpha_priors)

  spde_priors = list(
    rho = c(1, info$ρ[-1]),
    sigma = c(1, info$σ[-1]),
    rho_2 = c(1, info$ρ2[-1]))
    #rho_2 = c(0, log(ρ2)))
  spde_model = spde_generic_model_with_b_func(
    spde = spde,
    n = data$n,
    init = c(log(info$ρ[1]), log(info$σ[1]), log(info$ρ2[1])),
    #init = c(log(info$ρ[1]), log(info$σ[1])),
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
    A
  })

  stack = inla.stack(
    data = list(y = y_inla),
    A = list(spatial = A_inla, 1),
    effects = effects)

  formula = y ~ -1 +
    f(spatial, model = spde_model) +
    f(idx, model = alpha_model)

  est = optim(
    par = θ,
    fn = ll,
    y = data$y,
    y0 = data$y0,
    dist = data$dist_to_s0,
    dist_from_mesh = data$dist_to_s0_from_mesh,
    A = lapply(data$loc_index, function(x) A[x, ]),
    num_cores = n_cores,
    control = list(fnscale = -1, maxit = 120))

  fit = tryCatch({
    inla(
    formula = formula,
    data = inla.stack.data(stack),
    control.predictor = list(A = inla.stack.A(stack)),
    only.hyperparam = TRUE,
    #control.mode = list(theta = est$par[c(5, 3:4, 1:2)], restart = TRUE),
    control.mode = list(theta = est$par[c(6, 4:5, 3, 1:2)], restart = TRUE),
    control.inla = list(control.vb = list(enable = FALSE)),
    control.family = list(hyper = list(prec = list(param = info$τ[2:3]))),
    #verbose = TRUE,
    num.threads = n_cores,
    inla.mode = "experimental")
  }, error = function(e) NULL)
  if (is.null(fit)) next

  θ_reordering = c(5, 6, 4, 2, 3, 1)
  #θ_reordering = c(4, 5, 2, 3, 1)

  H = tryCatch(solve(fit$misc$cov.intern), error = \(e) NULL)
  if (is.null(H)) {
    warning("Inverting cov.intern failed for i = ", i)
    next
  }
  H = H[θ_reordering, θ_reordering]

  grads = ll_grad(
    θ = fit$mode$theta[θ_reordering],
    y = data$y,
    y0 = data$y0,
    dist = data$dist_to_s0,
    dist_from_mesh = data$dist_to_s0_from_mesh,
    A = lapply(data$loc_index, function(x) A[x, ]),
    num_cores = n_cores,
    sum_terms = FALSE)

  times = unlist(data$time_index)
  window_radius = 0
  J = 0
  for (k in seq_along(times)) {
    index = which(times >= times[k] - window_radius & times <= times[k] + window_radius)
    for (j in index) {
      J = J + grads[k, ] %*% grads[j, , drop = FALSE]
    }
  }

  # C = get_C(H, J)
  # θ_uncorrected = inla.hyperpar.sample(1e4, fit, intern = TRUE)
  # θ_uncorrected = θ_uncorrected[, θ_reordering]
  # θ_corrected = matrix(rep(fit$mode$theta[θ_reordering], each = 1e4), ncol = n_θ)
  # θ_corrected = θ_corrected + (θ_uncorrected - θ_corrected) %*% t(C)
  # apply(unname(θ_uncorrected), 2, quantile, probs = c(.025, .975))
  # apply(θ_corrected, 2, quantile, probs = c(.025, .975))
  # as.numeric(truth)

  fit = list(misc = fit$misc, internal.marginals.hyperpar = fit$internal.marginals.hyperpar,
             mode = list(theta = fit$mode$theta))
  fit$misc$reordering = NULL
  fit$misc$family = NULL
  fit$misc$linkfunctions = NULL
  class(fit) = "inla"

  model_fits = list(
    fit = fit,
    H = H,
    J = J,
    seed = seed,
    θ_reordering = θ_reordering,
    time_index = data$time_index,
    grads = grads,
    n = sum(data$n))

  tmp = readRDS(filename)
  tmp[[length(tmp) + 1]] = model_fits
  saveRDS(tmp, filename)

  pb$tick()
}
pb$terminate()

res = readRDS(filename)
truth = readRDS(truth_filename)$par
θ_reordering = res[[1]]$θ_reordering
n_θ = length(θ_reordering)

bad_index = which(sapply(res, length) == 0)
if (any(bad_index)) res = res[-bad_index]

for (i in seq_along(res)) {
  res[[i]]$C = get_C(res[[i]]$H, res[[i]]$J)
}

n_θ_per_res = 1e5
threshold = 5
probs = c(.9, .95, .99)
n_θ = length(θ_reordering)
pb = progress_bar(length(res))
if (length(truth) == 5) {
  θ_names = c("log_lambda", "log_kappa", "log_rho", "log_sigma", "log_precision")
  θ_tex_names = paste0("$\\", c("lambda", "kappa", "rho", "sigma", "tau"), "$")
} else {
  θ_names = c("log_lambda", "log_kappa", "log_rho_2", "log_rho", "log_sigma", "log_precision")
  θ_tex_names = paste0("$\\", c("lambda", "kappa", "rho_b", "rho", "sigma", "tau"), "$")
}
intervals = list()
for (i in seq_along(res)) {
  θ_uncorrected = inla.hyperpar.sample(n_θ_per_res, res[[i]]$fit, intern = TRUE)
  θ_uncorrected = θ_uncorrected[, θ_reordering]

  intervals[[i]] = list()
  for (k in 1:2) {
    if (k == 1) {
      θ_corrected = θ_uncorrected
    } else {
      θ_corrected = matrix(rep(res[[i]]$fit$mode$theta[θ_reordering], each = n_θ_per_res), ncol = n_θ)
      θ_corrected = θ_corrected + (θ_uncorrected - θ_corrected) %*% t(res[[i]]$C)
    }
    intervals[[i]][[k]] = rbind(
      as.data.frame(apply(θ_corrected, 2, quantile, probs = (1 - probs) / 2)) |>
        dplyr::mutate(k = k, prob = probs, lim = "lower"),
      as.data.frame(apply(θ_corrected, 2, quantile, probs = 1 - (1 - probs) / 2)) |>
        dplyr::mutate(k = k, prob = probs, lim = "upper"))
    row.names(intervals[[i]][[k]]) = NULL
    names(intervals[[i]][[k]])[1:n_θ] = θ_names
  }
  intervals[[i]] = do.call(rbind, intervals[[i]]) |>
    dplyr::mutate(iter = i)
  pb$tick()
}
intervals = do.call(rbind, intervals)
pb$terminate()

is_inside_interval = local({
  tmp1 = dplyr::filter(intervals, lim == "lower")
  tmp2 = dplyr::filter(intervals, lim == "upper")
  for (i in seq_along(θ_names)) {
    tmp1[[i]] = ifelse(tmp1[[i]] < truth[i], TRUE, FALSE)
    tmp2[[i]] = ifelse(tmp2[[i]] > truth[i], TRUE, FALSE)
  }
  tmp = (tmp1[, 1:n_θ] & tmp2[, 1:n_θ]) |>
    as.data.frame() |>
    dplyr::mutate(k = tmp1$k, prob = tmp1$prob, iter = tmp1$iter)
  tmp
})

is_inside_interval |>
  tidyr::pivot_longer(all_of(θ_names)) |>
  dplyr::group_by(prob, k, name) |>
  dplyr::summarise(coverage = mean(value)) |>
  #dplyr::mutate(k = paste("k =", k)) |>
  tidyr::pivot_wider(names_from = name, values_from = coverage) |>
  as.matrix() |>
  format(digits = 2) |>
  as.data.frame() |>
  dplyr::mutate(k = as.integer(k)) |>
  print(max = 5000)

tmp = is_inside_interval |>
  tidyr::pivot_longer(all_of(θ_names)) |>
  dplyr::group_by(prob, k, name) |>
  dplyr::summarise(coverage = mean(value)) |>
  dplyr::mutate(k = factor(k, levels = 1:2, labels = c("Unadjusted", "Adjusted"))) |>
  dplyr::mutate(name = factor(name, levels = θ_names, labels = θ_tex_names)) |>
  tidyr::pivot_wider(names_from = name, values_from = coverage)
if (length(truth) == 5) {
  tmp = tmp[, c(1, 2, 4, 3, 6, 7, 5)]
} else {
  tmp = tmp[, c(1, 2, 4, 3, 5, 7, 8, 6)]
}

table = paste(paste(c("Aim", "Method", names(tmp)[-c(1, 2)]), collapse = " & "), "\\\\")
table[2] = "\\midrule"
j = 3
for (i in 1:nrow(tmp)) {
  table[j] = paste(
    c(levels(tmp$k)[as.numeric(tmp$k[i])], paste0("$", round(100 * tmp[i, -c(1, 2)], digits = 0), "\\%$")),
    collapse = " & ")
  if (i %% length(unique(tmp$k)) == 1) {
    table[j] = paste(
      c(paste0("$", round(100 * tmp$prob[i], digits = 0), "\\%$"), table[j]), collapse = " & ")
  } else {
    table[j] = paste(c("", table[j]), collapse = " & ")
  }
  table[j] = paste(table[j], "\\\\")
  j = j + 1
  if (i %% length(unique(tmp$k)) == 0 && i != nrow(tmp)) {
    table[j] = "\\midrule"
    j = j + 1
  }
}
table = paste(paste(table, collapse = "\n"), "\n")

cat(table)
