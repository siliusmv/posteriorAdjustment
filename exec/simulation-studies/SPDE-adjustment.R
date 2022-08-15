devtools::load_all()
library(INLA)
library(mvtnorm)
library(ggplot2)

INLA::inla.setOption(pardiso.license = "~/.R/licences/pardiso.lic")

domain_size = 25
n_train = 200
n_loc = 400
n_cores = 7
verbose = TRUE
n_repl = 300

filename = file.path(tmp_dir(), "SPDE-adjustment.rds")

σ = 1
ρ = 12
τ = 100

set.seed(1)
loc = matrix(runif(n_loc * 2) * domain_size, ncol = 2)
dist = as.matrix(dist(loc))

matern_corr = function(dist, ρ, ν = 1.5) {
  κ = sqrt(8 * ν) / ρ
  res = 2 ^ (1 - ν) / gamma(ν) * (κ * dist) ^ ν * besselK(κ * dist, ν)
  res[dist == 0] = 1
  res
}
Σ = matern_corr(dist, ρ) * σ^2


# Create the mesh
mesh = inla.mesh.2d(
  loc.domain = loc,
  boundary = list(inla.nonconvex.hull(loc, convex = -.3),
                  inla.nonconvex.hull(loc, convex = -1)),
  cutoff = 4,
  max.edge = c(5, 50))
# plot(mesh); points(loc, cex = .5)

loglik = function(θ,
                  y,
                  spde,
                  A,
                  sum = TRUE) {
  τ = exp(θ[1])
  ρ = exp(θ[2])
  σ = exp(θ[3])
  Q = INLA::inla.spde2.precision(spde, theta = c(log(ρ), log(σ)))
  Σ = A %*% Matrix::solve(Q, Matrix::t(A))
  Σ = as.matrix(Σ) + diag(1 / τ, nrow(Σ))
  res = mvtnorm::dmvnorm(t(y), sigma = Σ, log = TRUE)
  if (sum) res = sum(res)
  res
}

n_truth = 1e5
set.seed(1)
y = mvtnorm::rmvnorm(n_truth, sigma = Σ) + rnorm(n_truth * n_loc, sd = τ^-.5)
y = t(y)

spde = inla.spde2.pcmatern(
  mesh,
  prior.range = c(ρ, 0.5),
  prior.sigma = c(σ, 0.5))
A = inla.spde.make.A(mesh, loc)

est = list(par = c(log(τ), log(ρ), log(σ)), convergence = 1)
while (est$convergence != 0) {
  est = optim(
    par = est$par,
    fn = loglik,
    y = y,
    spde = spde,
    A = A,
    control = list(fnscale = -1, trace = 6))
}
truth = est$par

set.seed(1)
model_fits = list()
pb = progress_bar(n_repl)
for (i in 1:n_repl) {

  # Sample from the SPDE
  y = mvtnorm::rmvnorm(n_train, sigma = Σ) + rnorm(n_train * n_loc, sd = τ^-.5)
  y = t(y)

  est = optim(
    par = c(log(τ), log(ρ), log(σ)),
    fn = loglik,
    y = y,
    spde = spde,
    A = A,
    control = list(fnscale = -1))

  spde = inla.spde2.pcmatern(
    mesh,
    prior.range = c(ρ, 0.5),
    prior.sigma = c(σ, 0.5))

  effects = list(spatial = inla.spde.make.index("spatial", n.spde = mesh$n, n.repl = n_train))
  y_inla = as.numeric(y)
  A_inla = inla.spde.make.A(
    mesh,
    loc,
    index = rep(1:n_loc, n_train),
    repl = rep(1:n_train, each = n_loc))
  stack = inla.stack(
    data = list(y = y_inla),
    A = list(spatial = A_inla),
    effects = effects)

  formula = y ~ -1 + f(spatial, model = spde, replicate = spatial.repl, nrep = n_train)

  fit = inla(
    formula = formula,
    data = inla.stack.data(stack),
    family = "gaussian",
    num.threads = n_cores,
    only.hyperparam = TRUE,
    control.mode = list(theta = est$par, restart = TRUE),
    control.inla = list(control.vb = list(enable = FALSE)),
    inla.mode = "experimental",
    control.predictor = list(A = inla.stack.A(stack)))

  θ = fit$mode$theta

  # --------------------------------------------------------------------------
  # Estimate H and J
  # --------------------------------------------------------------------------

  A = inla.spde.make.A(mesh, loc)

  H = solve(fit$misc$cov.intern)

  grads = numDeriv::jacobian(
    function(x) loglik(x, y, spde, A, sum = FALSE),
    x = θ)

  J = 0
  for (j in 1:n_train) {
    J = J + grads[j, ] %*% grads[j, , drop = FALSE]
  }

  fit = list(
    misc = fit$misc,
    internal.marginals.hyperpar = fit$internal.marginals.hyperpar,
    summary.hyperpar = fit$summary.hyperpar,
    model = list(theta = fit$mode$theta))
  fit$misc$reordering = NULL
  fit$misc$family = NULL
  fit$misc$linkfunctions = NULL
  class(fit) = "inla"

  model_fits[[i]] = list(fit = fit, J = J, H = H)
  saveRDS(model_fits, filename)
  pb$tick()
}
pb$terminate()

attributes(model_fits)$truth = truth
saveRDS(model_fits, filename)

res = readRDS(filename)
truth = attributes(res)$truth
n_θ = length(truth)
θ_names = c("log_precision", "log_rho", "log_sigma")

for (i in seq_along(res)) {
  res[[i]]$C = get_C(res[[i]]$H, res[[i]]$J)
}

n_θ_per_res = 1e5
threshold = 5
probs = c(.9, .95, .99)
pb = progress_bar(length(res))
intervals = list()
for (i in seq_along(res)) {
  intervals[[i]] = list()
  θ_uncorrected = inla.hyperpar.sample(n_θ_per_res, res[[i]]$fit, intern = TRUE)
  for (k in 1:2) {
    if (k == 1) {
      θ_corrected = θ_uncorrected
      label = "Unadjusted"
    } else {
      θ_corrected = matrix(rep(res[[i]]$fit$mode$theta, each = n_θ_per_res), ncol = n_θ)
      θ_corrected = θ_corrected + (θ_uncorrected - θ_corrected) %*% t(res[[i]]$C)
      label = "Adjusted"
    }
    intervals[[i]][[length(intervals[[i]]) + 1]] = rbind(
      as.data.frame(apply(θ_corrected, 2, quantile, probs = (1 - probs) / 2)) |>
        dplyr::mutate(label = label, prob = probs, lim = "lower"),
      as.data.frame(apply(θ_corrected, 2, quantile, probs = 1 - (1 - probs) / 2)) |>
        dplyr::mutate(label = label, prob = probs, lim = "upper"))
    row.names(intervals[[i]][[length(intervals[[i]])]]) = NULL
    names(intervals[[i]][[length(intervals[[i]])]])[1:n_θ] = θ_names
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
    dplyr::mutate(label = tmp1$label, prob = tmp1$prob, iter = tmp1$iter)
  tmp
})

tmp = is_inside_interval |>
  tidyr::pivot_longer(all_of(θ_names)) |>
  dplyr::group_by(prob, label, name) |>
  dplyr::summarise(coverage = mean(value)) |>
  dplyr::mutate(name = factor(
    x = name,
    levels = θ_names,
    labels = paste0("$\\", c("tau", "rho", "sigma"), "$"))) |>
  tidyr::pivot_wider(names_from = name, values_from = coverage)
tmp = tmp[rep(2 * seq_len(nrow(tmp) / 2), each = 2) - c(0, 1), ]
print(tmp)

table = paste(paste(c("\\(1 - \\alpha\\)", "Method", names(tmp)[-(1:2)]), collapse = " & "), "\\\\")
table[2] = "\\midrule"
j = 3
for (i in 1:nrow(tmp)) {
  table[j] = paste(
    c(tmp$label[i], paste0("$", round(100 * tmp[i, -(1:2)], digits = 0), "\\%$")),
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
