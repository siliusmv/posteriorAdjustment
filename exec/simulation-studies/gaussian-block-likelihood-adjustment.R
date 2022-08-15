devtools::load_all()
library(INLA)
library(mvtnorm)
library(ggplot2)

INLA::inla.setOption(pardiso.license = "~/.R/licences/pardiso.lic")

domain_size = 25
block_size = 5
n_blocks = (domain_size / block_size)^2
n_train = 100
n_loc = 400
n_cores = 7
verbose = TRUE
n_repl = 300

filename = file.path(tmp_dir(), "block-adjustment.rds")

σ = 1
ρ = 12
τ = 100

set.seed(1)
loc = matrix(runif(n_loc * 2) * domain_size, ncol = 2)

mesh = inla.mesh.2d(
  loc.domain = loc,
  boundary = list(inla.nonconvex.hull(loc, convex = -.3),
                  inla.nonconvex.hull(loc, convex = -1)),
  cutoff = 6,
  max.edge = c(9, 50))
#plot(mesh); points(loc, cex = .5)

block_index = list()
for (j in 1:(domain_size / block_size)) {
  for (i in 1:(domain_size / block_size)) {
    y_lims = c(i - 1, i) * block_size
    x_lims = c(j - 1, j) * block_size
    block_index[[length(block_index) + 1]] = which(
      x_lims[1] < loc[, 1] & loc[, 1] <= x_lims[2] &
      y_lims[1] < loc[, 2] & loc[, 2] <= y_lims[2])
  }
}
full_index = list(1:n_loc)

block_loglik = function(θ,
                        y,
                        block_index,
                        spde,
                        A,
                        return_mat = FALSE,
                        return_vec = FALSE) {
  τ = exp(θ[1])
  ρ = exp(θ[2])
  σ = exp(θ[3])
  Q = INLA::inla.spde2.precision(spde, theta = c(log(ρ), log(σ)))
  Σ = A %*% Matrix::solve(Q, Matrix::t(A))
  Σ = as.matrix(Σ) + diag(1 / τ, nrow(Σ))
  res = lapply(
    X = block_index,
    FUN = function(x) matrix(NA_real_, nrow = length(x), ncol = ncol(y)))
  for (i in seq_along(block_index)) {
    if (any(block_index[[i]])) {
      res[[i]] = mvtnorm::dmvnorm(
        t(y[block_index[[i]], ]),
        sigma = Σ[block_index[[i]], block_index[[i]]],
        log = TRUE)
    }
  }
  res = do.call(rbind, res)
  if (!return_mat) {
    res = apply(res, 2, sum)
    if (!return_vec) res = sum(res)
  }
  res
}

set.seed(1)
model_fits = list()
attributes(model_fits)$block_index = block_index
attributes(model_fits)$truth = log(c(τ, ρ, σ))
pb = progress_bar(n_repl)
for (i in 1:n_repl) {

  seed = .Random.seed

  # Sample from the SPDE
  spde_tmp = inla.spde2.pcmatern(
    mesh,
    prior.range = c(ρ, 0.5),
    prior.sigma = c(σ, 0.5))
  Q = inla.spde2.precision(spde_tmp, theta = c(log(ρ), log(σ)))
  A_minimal = inla.spde.make.A(mesh, loc)

  z = rnorm_spde(n_train, Q)
  y = as.matrix(A_minimal %*% z)
  y = y + rnorm(length(y), sd = τ^-.5)

  spde = inla.spde2.pcmatern(
    mesh,
    prior.range = c(ρ, 0.5),
    prior.sigma = c(σ, 0.5))

  effects_blocks = list(
    spatial = inla.spde.make.index("spatial", n.spde = mesh$n, n.repl = n_train * n_blocks))
  effects_full = list(
    spatial = inla.spde.make.index("spatial", n.spde = mesh$n, n.repl = n_train))

  y_inla_blocks = as.numeric(y[unlist(block_index), ])
  y_inla_full = as.numeric(y[unlist(full_index), ])

  A_blocks = list(spatial = inla.spde.make.A(
    mesh = mesh,
    loc = loc,
    index = rep(unlist(block_index), n_train),
    repl = rep(1:(n_blocks * n_train), rep(sapply(block_index, length), n_train))))
  A_full = list(spatial = inla.spde.make.A(
    mesh = mesh,
    loc = loc,
    index = rep(unlist(full_index), n_train),
    repl = rep(1:n_train, rep(sapply(full_index, length), n_train))))

  stack_blocks = inla.stack(data = list(y = y_inla_blocks), A = A_blocks, effects = effects_blocks)
  stack_full = inla.stack(data = list(y = y_inla_full), A = A_full, effects = effects_full)

  formula_blocks = y ~ -1 + f(spatial, model = spde, replicate = spatial.repl, nrep = n_train * n_blocks)
  formula_full = y ~ -1 + f(spatial, model = spde, replicate = spatial.repl, nrep = n_train)

  block_fit = inla(
    formula = formula_blocks,
    data = inla.stack.data(stack_blocks),
    family = "gaussian",
    num.threads = n_cores,
    only.hyperparam = TRUE,
    control.inla = list(control.vb = list(enable = FALSE)),
    inla.mode = "experimental",
    control.predictor = list(A = inla.stack.A(stack_blocks)))

  full_fit = inla(
    formula = formula_full,
    data = inla.stack.data(stack_full),
    family = "gaussian",
    num.threads = n_cores,
    only.hyperparam = TRUE,
    control.inla = list(control.vb = list(enable = FALSE)),
    inla.mode = "experimental",
    control.predictor = list(A = inla.stack.A(stack_full)))

  θ_block = block_fit$mode$theta
  θ_full = full_fit$mode$theta

  # --------------------------------------------------------------------------
  # Estimate H and J
  # --------------------------------------------------------------------------

  A = inla.spde.make.A(mesh, loc)

  H_block = solve(block_fit$misc$cov.intern)
  H_full = solve(full_fit$misc$cov.intern)

  grads_block = numDeriv::jacobian(
    function(x) block_loglik(x, y, block_index, spde, A, return_mat = TRUE),
    x = θ_block)
  grads_full = numDeriv::jacobian(
    function(x) block_loglik(x, y, full_index, spde, A, return_mat = TRUE),
    x = θ_full)

  J_block = 0
  J_full = 0
  for (j in 1:n_train) {
    for (k in 1:n_blocks) {
      for (l in 1:n_blocks) {
        J_block = J_block + grads_block[(j - 1) * n_blocks + k, ] %*%
          grads_block[(j - 1) * n_blocks + l, , drop = FALSE]
      }
    }
    J_full = J_full + grads_full[j, ] %*% grads_full[j, , drop = FALSE]
  }

  block_fit = list(
    misc = block_fit$misc,
    internal.marginals.hyperpar = block_fit$internal.marginals.hyperpar,
    model = list(theta = block_fit$mode$theta))
  block_fit$misc$reordering = NULL
  block_fit$misc$family = NULL
  block_fit$misc$linkfunctions = NULL
  class(block_fit) = "inla"

  full_fit = list(
    misc = full_fit$misc,
    internal.marginals.hyperpar = full_fit$internal.marginals.hyperpar,
    model = list(theta = full_fit$mode$theta))
  full_fit$misc$reordering = NULL
  full_fit$misc$family = NULL
  full_fit$misc$linkfunctions = NULL
  class(full_fit) = "inla"

  model_fits[[i]] = list(
    seed = seed,
    block_fit = block_fit,
    full_fit = full_fit,
    J_block = J_block,
    H_block = H_block,
    J_full = J_full,
    H_full = H_full)
  saveRDS(model_fits, filename)
  pb$tick()
}
pb$terminate()

res = readRDS(filename)
block_index = attributes(res)$block_index
truth = attributes(res)$truth
n_θ = length(truth)
θ_names = c("log_precision", "log_rho", "log_sigma")

for (i in seq_along(res)) {
  res[[i]]$C_block = get_C(res[[i]]$H_block, res[[i]]$J_block)
  res[[i]]$C_full = get_C(res[[i]]$H_full, res[[i]]$J_full)
}

n_θ_per_res = 1e5
threshold = 5
probs = c(.9, .95, .99)
n_θ = length(truth)
pb = progress_bar(length(res))
intervals = list()
for (i in seq_along(res)) {
  intervals[[i]] = list()
  for (j in 1:2) {
    if (j == 1) {
      fit = res[[i]]$block_fit
      C = res[[i]]$C_block
      model = "Block"
    } else {
      fit = res[[i]]$full_fit
      C = res[[i]]$C_full
      model = "Full"
    }
    θ_uncorrected = inla.hyperpar.sample(n_θ_per_res, fit, intern = TRUE)
    for (k in 1:2) {
      if (k == 1) {
        θ_corrected = θ_uncorrected
        label = "Unadjusted"
      } else {
        θ_corrected = matrix(rep(fit$mode$theta, each = n_θ_per_res), ncol = n_θ)
        θ_corrected = θ_corrected + (θ_uncorrected - θ_corrected) %*% t(C)
        label = "Adjusted"
      }
      intervals[[i]][[length(intervals[[i]]) + 1]] = rbind(
        as.data.frame(apply(θ_corrected, 2, quantile, probs = (1 - probs) / 2)) |>
          dplyr::mutate(model = model, label = label, prob = probs, lim = "lower"),
        as.data.frame(apply(θ_corrected, 2, quantile, probs = 1 - (1 - probs) / 2)) |>
          dplyr::mutate(model = model, label = label, prob = probs, lim = "upper"))
      row.names(intervals[[i]][[length(intervals[[i]])]]) = NULL
      names(intervals[[i]][[length(intervals[[i]])]])[1:n_θ] = θ_names
    }
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
    dplyr::mutate(model = tmp1$model, label = tmp1$label,
                  prob = tmp1$prob, iter = tmp1$iter)
  tmp
})

tmp = is_inside_interval |>
  tidyr::pivot_longer(all_of(θ_names)) |>
  dplyr::group_by(prob, label, model, name) |>
  dplyr::summarise(coverage = mean(value)) |>
  dplyr::mutate(name = factor(
    x = name,
    levels = θ_names,
    labels = paste0("$\\", c("tau", "rho", "sigma"), "$"))) |>
  tidyr::pivot_wider(names_from = name, values_from = coverage)
tmp = tmp[1:nrow(tmp) + rep(c(2, 2, -2, -2), 3), ]
print(tmp)

table = paste(paste(
  c("", "", "\\multicolumn{3}{c}{Block likelihood}", "\\multicolumn{3}{c}{Full likelihood}"),
  collapse = " & "), "\\\\")
table[2] = "\\cmidrule(r){3-5} \\cmidrule(l){6-8}"
table[3] = paste(paste(c("\\(1 - \\alpha\\)", "Method", rep(names(tmp)[-(1:3)], 2)), collapse = " & "), "\\\\")
table[4] = "\\midrule"
j = 5
for (i in seq(1, nrow(tmp) - 1, by = 2)) {
  table[j] = paste(
    c(tmp$label[i],
      paste0("$", round(100 * tmp[i, -(1:3)], digits = 0), "\\%$"),
      paste0("$", round(100 * tmp[i + 1, -(1:3)], digits = 0), "\\%$")),
    collapse = " & ")
  if (i %% 4 == 1) {
    table[j] = paste(
      c(paste0("$", round(100 * tmp$prob[i], digits = 0), "\\%$"), table[j]), collapse = " & ")
  } else {
    table[j] = paste(c("", table[j]), collapse = " & ")
  }
  table[j] = paste(table[j], "\\\\")
  j = j + 1
  if (i %% 4 == 3 && i != (nrow(tmp) - 1)) {
    table[j] = "\\midrule"
    j = j + 1
  }
}
table = paste(paste(table, collapse = "\n"), "\n")

cat(table)
