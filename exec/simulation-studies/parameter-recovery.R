devtools::load_all()
library(INLA)
library(ggplot2)
library(sf)
library(dplyr)
library(tidyr)
library(scoringRules)

INLA::inla.setOption(pardiso.license = "~/.R/licences/pardiso.lic")

filename = file.path(tmp_dir(), "parameter_recovery.rds")
real_data_filename = file.path(tmp_dir(), "conditional-modelling.rds")
ρ2_filename = file.path(tmp_dir(), "model-selection.rds")

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

a_func = get_a_func(θ_a)
b_func = get_b_func(θ_b)
Σ = get_Σ_func(θ_Σ, spde)
Q = Σ$Q

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

s0 = matrix(c(281, 7092), nrow = 1)
s0_index = which(coords[, 1] == s0[1] & coords[, 2] == s0[2])

dist_to_s0 = dist_euclid(coords, s0)
dist_to_s0_from_mesh = dist_euclid(mesh$loc[, 1:2], s0)

# -----------------------------------------------------------------------
# Test out the model y = a(y0, |s - s0|) + b(y0, |s - s0|) * Z + ε
# -----------------------------------------------------------------------

n_train = 100
cl = parallel::makeForkCluster(n_cores, mc.set.seed = TRUE)
parallel::clusterSetRNGStream(cl, 1)
params = pbapply::pblapply(
  X = 1:n_sim,
  cl = cl,
  FUN = function(i) {

    seed = .Random.seed

    # Sample y0
    y0 = threshold + rexp(n_train)

    # Create y
    z = rnorm_spde(n = n_train, Q = Q)
    y = as.matrix(A %*% z)
    for (j in seq_along(y0)) {
      y[, j] = y[, j] * b_func(1, dist_to_s0) + a_func(y0[j], dist_to_s0)
    }
    y = y + rnorm(length(y), sd = τ^-.5)

    est = optim(
      par = θ,
      fn = function(x, ...) tryCatch(ll(x, ...), error = function(e) -9999999999),
      y = list(y),
      y0 = list(y0),
      dist = list(dist_to_s0),
      dist_from_mesh = list(dist_to_s0_from_mesh),
      A = list(A),
      control = list(fnscale = -1, maxit = 300))

    dist_to_s0_inla = rep(dist_to_s0, n_train)
    y0_inla = rep(y0, each = n_loc)
    y_inla = as.numeric(y)

    # Create the SPDE models
    priors = list(
      rho = c(1, info$ρ[2:3]),
      sigma = c(1, info$σ[2:3]),
      rho_2 = c(1, info$ρ2[2:3]))
    init = c(log(info$ρ[1]), log(info$σ[1]), info$ρ2[1])
    spde = inla.spde2.pcmatern(
      mesh,
      prior.range = priors$rho[-1],
      prior.sigma = priors$sigma[-1])
    spde_model = spde_generic_model_with_b_func(
      spde = spde,
      n = n_train,
      init = init,
      priors = priors,
      dist_to_s0 = dist_to_s0_from_mesh)

    # Create the alpha_model
    priors = list(
      lambda = c(1, info$λ[2:3]),
      kappa = c(1, info$κ[2:3]))
    alpha_model = alpha_generic_model(
      y0 = y0_inla,
      dist_to_s0 = dist_to_s0_inla,
      init = c(info$λ[1], info$κ[1]),
      priors = priors)

    formula = y ~ -1 +
      f(idx, model = alpha_model) +
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

    fit = tryCatch({
      inla(
        formula,
        data = inla.stack.data(stack),
        control.predictor = list(A = inla.stack.A(stack)),
        control.compute = list(waic = TRUE, dic = TRUE),
        num.threads = 1,
        only.hyperparam = TRUE,
        control.family = list(hyper = list(prec = list(
          #initial = info$τ[1],
          param = info$τ[2:3]))),
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

    hyper = inla.hyperpar.sample(10000, fit) |>
      (\(x) {x[is.infinite(x)] = NA; x})() |>
      apply(2, function(x, ...) c(mean(x, ...), quantile(x, c(.025, .975), ...)), na.rm = TRUE)

    res = data.frame(
      value = hyper[1, ] |> (\(x) c(1 / x[1], exp(x[2:3]), exp(x[4:6])))(),
      lower = hyper[2, ] |> (\(x) c(1 / x[1], exp(x[2:3]), exp(x[4:6])))(),
      upper = hyper[3, ] |> (\(x) c(1 / x[1], exp(x[2:3]), exp(x[4:6])))(),
      name = c("nugget", "λ", "κ", "ρ", "σ", "ρ2"),
      truth = c(1 / τ, λ, κ, ρ, σ, ρ2),
      cpu = fit$cpu[4],
      mlik = fit$mlik[1],
      i = i)
    res[1, 2:3] = res[1, 3:2]
    attributes(res)$seed = seed
    row.names(res) = NULL

    res
  })
parallel::stopCluster(cl)

if (FALSE) {
  saveRDS(params, filename)
}

bad_runs = which(sapply(params, length) == 2)
if (any(bad_runs)) params = params[-bad_runs]

params |>
  do.call(what = rbind) |>
  dplyr::mutate(err = truth - value, rel_err = err / truth) |>
  dplyr::mutate(is_included = truth > lower & truth < upper) |>
  dplyr::group_by(name) |>
  dplyr::summarise(mean = mean(err), mean_rel = mean(rel_err),
                   sd = sd(err), max = max(abs(err)),
                   mse = mean(err^2), mae = mean(abs(err)),
                   coverage = mean(is_included), mlik = mean(mlik))

params |>
  do.call(what = rbind) |>
  dplyr::mutate(
    name = factor(
      name,
      levels = c("nugget", "ρ", "σ", "ρ2", "λ", "κ"),
      labels = c("nugget", "rho", "sigma", "rho_b", "lambda", "kappa"))) |>
  ggplot() +
  geom_density(aes(x = value)) +
  facet_wrap(~name, scales = "free") +
  geom_vline(aes(xintercept = truth))
