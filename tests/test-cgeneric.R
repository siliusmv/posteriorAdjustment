library(INLA)
library(Matrix)
library(MASS)
devtools::load_all()
INLA::inla.setOption(pardiso.license = "~/.R/licences/pardiso.lic")
INLA::inla.pardiso.check()

compile_cgeneric("alpha")
compile_cgeneric("spde")

# -----------------------------------------------------------------------
# Test out the model y = α * y(s0) + ε
# -----------------------------------------------------------------------

# Choose a design for the simulation study
loc = as.matrix(expand.grid(x = seq(2, 10, by = 2), y = seq(2, 10, by = 2)))
num_cores = 6
n_train = 100
σ_ε = .1
λ = 4.3
κ = 1.8
threshold = 4
verbose = TRUE

# Place s0
s0 = matrix(c(4.5, 5.5), nrow = 1)

# Compute distances
dist = as.matrix(dist(loc))
dist_to_s0 = as.matrix(dist(rbind(s0, loc)))[-1, 1]

# Sample y0
set.seed(1)
y0 = threshold + rexp(n_train)

a = get_a_funcs(log(c(λ, κ)))$value

# Create y
y = a(y0, dist_to_s0)
y = y + rnorm(length(y), sd = σ_ε)

# Create the alpha_model using rgeneric and cgeneric
init = c(log(λ), log(κ))
log_lambda_prior = c(log(λ), 1)
log_kappa_prior = c(log(κ), 1)
alpha_r_model = alpha_generic_model(
  y0 = y0,
  dist_to_s0 = dist_to_s0,
  n = length(y),
  init = init,
  log_lambda_prior = log_lambda_prior,
  log_kappa_prior = log_kappa_prior,
  use_c = FALSE)
alpha_c_model = alpha_generic_model(
  y0 = y0,
  dist_to_s0 = dist_to_s0,
  n = length(y),
  init = init,
  log_lambda_prior = log_lambda_prior,
  log_kappa_prior = log_kappa_prior,
  use_c = TRUE)

# Create everything else needed for running INLA
formula_r = y ~ -1 + f(idx, model = alpha_r_model)
formula_c = y ~ -1 + f(idx, model = alpha_c_model)

# Run INLA
res_r = inla(formula_r,
             data = data.frame(y = as.numeric(y), idx = seq_along(y)),
             num.threads = num_cores,
             inla.mode = "experimental",
             only.hyperparam = TRUE,
             control.inla = list(control.vb = list(enable = FALSE)),
             verbose = verbose)

res_c = inla(formula_c,
             data = data.frame(y = as.numeric(y), idx = seq_along(y)),
             num.threads = num_cores,
             inla.mode = "experimental",
             only.hyperparam = TRUE,
             control.inla = list(control.vb = list(enable = FALSE)),
             verbose = verbose)

σ_ε; res_r$summary.hyperpar[1, 1]^-.5; res_c$summary.hyperpar[1, 1]^-.5
λ; exp(res_r$summary.hyperpar[2, 1]); exp(res_c$summary.hyperpar[2, 1])
κ; exp(res_r$summary.hyperpar[3, 1]); exp(res_c$summary.hyperpar[3, 1])
res_r$mlik[1]; res_c$mlik[1]
res_r$cpu; res_c$cpu

# For some reason, the precision estimates are quite off.
# This seems to get improved if we add a lot more data, though

# -----------------------------------------------------------------------
# Test out the model y = α * y(s0) + ε, when using multiple s0's
# -----------------------------------------------------------------------

# Choose a design for the simulation study
loc = as.matrix(expand.grid(x = seq(2, 10, by = 2), y = seq(2, 10, by = 2)))
num_cores = 6
n_s0 = 5
n_train_per_s0 = 20
σ_ε = .1
λ = 4.3
κ = 1.8
threshold = 4
verbose = TRUE

# Place s0
s0 = matrix(runif(n_s0 * 2, 3, 9), ncol = 2)

# Compute distances
dist = as.matrix(dist(loc))
dist_to_s0 = lapply(1:n_s0, function(i) as.matrix(dist(rbind(s0[i, ], loc)))[-1, 1])

# Sample y0
y0 = lapply(1:n_s0, function(i) threshold + rexp(n_train_per_s0))

a = get_a_funcs(log(c(λ, κ)))$value

# Create y
y = list()
for (i in seq_along(y0)) {
  y[[i]] = a(y0[[i]], dist_to_s0[[i]])
  y[[i]] = y[[i]] + rnorm(length(y[[i]]), sd = σ_ε)
}

# Create the alpha_model using rgeneric and cgeneric
init = c(log(λ), log(κ))
log_lambda_prior = c(log(λ), 1)
log_kappa_prior = c(log(κ), 1)
alpha_r_model = alpha_generic_model(
  y0 = y0,
  dist_to_s0 = dist_to_s0,
  n = length(unlist(y)),
  init = init,
  log_lambda_prior = log_lambda_prior,
  log_kappa_prior = log_kappa_prior,
  use_c = FALSE)
alpha_c_model = alpha_generic_model(
  y0 = y0,
  dist_to_s0 = dist_to_s0,
  n = length(unlist(y)),
  init = init,
  log_lambda_prior = log_lambda_prior,
  log_kappa_prior = log_kappa_prior,
  use_c = TRUE)

# Create everything else needed for running INLA
formula_r = y ~ -1 + f(idx, model = alpha_r_model)
formula_c = y ~ -1 + f(idx, model = alpha_c_model)

# Run INLA
res_r = inla(formula_r,
             data = data.frame(y = unlist(y), idx = seq_along(unlist(y))),
             num.threads = num_cores,
             only.hyperparam = TRUE,
             control.inla = list(control.vb = list(enable = FALSE)),
             inla.mode = "experimental",
             verbose = verbose)

res_c = inla(formula_c,
             data = data.frame(y = unlist(y), idx = seq_along(unlist(y))),
             num.threads = num_cores,
             only.hyperparam = TRUE,
             control.inla = list(control.vb = list(enable = FALSE)),
             inla.mode = "experimental",
             verbose = verbose)

σ_ε; res_r$summary.hyperpar[1, 1]^-.5; res_c$summary.hyperpar[1, 1]^-.5
λ; exp(res_r$summary.hyperpar[2, 1]); exp(res_c$summary.hyperpar[2, 1])
κ; exp(res_r$summary.hyperpar[3, 1]); exp(res_c$summary.hyperpar[3, 1])
res_r$mlik[1]; res_c$mlik[1]
res_r$cpu; res_c$cpu


# -----------------------------------------------------------------------
# Test out the model y = Z + ε with only one repl
# -----------------------------------------------------------------------

loc = as.matrix(expand.grid(x = seq(2, 10, by = .6), y = seq(2, 10, by = .6)))
σ = 2
ρ = 5
n_train = 1

# Place s0
s0 = matrix(c(4.5, 5.5), nrow = 1)

# Create the mesh
mesh = inla.mesh.2d(
  loc.domain = loc,
  min.angle = 30,
  boundary = list(inla.nonconvex.hull(loc, convex = -.3),
                  inla.nonconvex.hull(loc, convex = -1)),
  #cutoff = 1,
  cutoff = .3,
  #max.edge = c(2, 5))
  max.edge = c(.6, 3))
#plot(mesh); points(loc)

# extract the precision matrix
spde_tmp = inla.spde2.pcmatern(mesh,
                               prior.range = c(ρ, 0.5),
                               prior.sigma = c(σ, 0.5))
Q = inla.spde2.precision(spde_tmp, theta = c(log(ρ), log(σ)))

# Sample y = Z + ε,
# with ε Gaussian white noise, and Z a Matern-field from the SPDE
A_minimal = inla.spde.make.A(mesh, loc)
constr = list(A = as.matrix(inla.spde.make.A(mesh, s0)), e = 0)
z = inla.qsample(n = n_train, Q = Q, constr = constr)
y = as.matrix(A_minimal %*% z)
y = y + rnorm(length(y), sd = σ_ε)

# Create the SPDE models
prior_range = c(ρ, .5)
prior_sigma = c(σ, .5)
init = c(log(ρ), log(σ))
spde_builtin = inla.spde2.pcmatern(
  mesh,
  prior.range = prior_range,
  prior.sigma = prior_sigma)
spde_r_model = spde_generic_model(
  spde = spde_builtin,
  init = init,
  rho_prior = prior_range,
  sigma_prior = prior_sigma,
  use_c = FALSE)
spde_c_model = spde_generic_model(
  spde = spde_builtin,
  init = init,
  rho_prior = prior_range,
  sigma_prior = prior_sigma,
  use_c = TRUE)
spde_r_model_with_replicates = spde_generic_model_with_replicates(
  spde = spde_builtin,
  init = init,
  n_repl = n_train,
  rho_prior = prior_range,
  sigma_prior = prior_sigma,
  use_c = FALSE)
spde_c_model_with_replicates = spde_generic_model_with_replicates(
  spde = spde_builtin,
  init = init,
  n_repl = n_train,
  rho_prior = prior_range,
  sigma_prior = prior_sigma,
  use_c = TRUE)

# Create everything else needed for running INLA
formula_builtin = y ~ -1 + f(spatial, model = spde_builtin, replicate = spatial.repl, nrep = n_train, extraconstr = constr)
formula_r_model = y ~ -1 + f(spatial, model = spde_r_model, replicate = spatial.repl, nrep = n_train, extraconstr = constr)
formula_c_model = y ~ -1 + f(spatial, model = spde_c_model, replicate = spatial.repl, nrep = n_train, extraconstr = constr)
formula_r_model_with_replicates = y ~ -1 + f(spatial, model = spde_r_model_with_replicates, extraconstr = constr)
formula_c_model_with_replicates = y ~ -1 + f(spatial, model = spde_c_model_with_replicates, extraconstr = constr)

effects = list(
  spatial = inla.spde.make.index(
    "spatial",
    n.spde = spde_builtin$n.spde,
    n.repl = n_train))

effects_with_replicates = list(spatial = seq_len(n_train * mesh$n))

A = inla.spde.make.A(
  mesh = mesh,
  loc = loc,
  index = rep(1:nrow(loc), n_train),
  repl = rep(1:n_train, each = nrow(loc)))

stack = inla.stack(
  data = list(y = as.numeric(y)),
  A = list(spatial = A),
  effects = effects)

stack_with_replicates = inla.stack(
  data = list(y = as.numeric(y)),
  A = list(spatial = A),
  effects = effects_with_replicates)

# Run INLA
res_builtin = inla(formula_builtin,
                   data = inla.stack.data(stack),
                   control.predictor = list(A = inla.stack.A(stack)),
                   inla.mode = "experimental",
                   only.hyperparam = TRUE,
                   control.inla = list(control.vb = list(enable = FALSE)),
                   num.threads = num_cores,
                   verbose = verbose)

res_r_model = inla(formula_r_model,
                   data = inla.stack.data(stack),
                   control.predictor = list(A = inla.stack.A(stack)),
                   only.hyperparam = TRUE,
                   control.inla = list(control.vb = list(enable = FALSE)),
                   inla.mode = "experimental",
                   num.threads = num_cores,
                   verbose = verbose)

res_c_model = inla(formula_c_model,
                   data = inla.stack.data(stack),
                   control.predictor = list(A = inla.stack.A(stack)),
                   only.hyperparam = TRUE,
                   control.inla = list(control.vb = list(enable = FALSE)),
                   inla.mode = "experimental",
                   num.threads = num_cores,
                   verbose = verbose)

res_r_model_with_replicates = inla(formula_r_model_with_replicates,
                   data = inla.stack.data(stack_with_replicates),
                   control.predictor = list(A = inla.stack.A(stack_with_replicates)),
                   only.hyperparam = TRUE,
                   control.inla = list(control.vb = list(enable = FALSE)),
                   inla.mode = "experimental",
                   num.threads = num_cores,
                   verbose = verbose)

res_c_model_with_replicates = inla(formula_c_model_with_replicates,
                   data = inla.stack.data(stack_with_replicates),
                   control.predictor = list(A = inla.stack.A(stack_with_replicates)),
                   inla.mode = "experimental",
                   only.hyperparam = TRUE,
                   control.inla = list(control.vb = list(enable = FALSE)),
                   num.threads = num_cores,
                   verbose = verbose)


σ_ε; res_builtin$summary.hyperpar[1, 1]^-.5; res_r_model$summary.hyperpar[1, 1]^-.5; res_c_model$summary.hyperpar[1, 1]^-.5; res_r_model_with_replicates$summary.hyperpar[1, 1]^-.5; res_c_model_with_replicates$summary.hyperpar[1, 1]^-.5
ρ; res_builtin$summary.hyperpar[2, 1]; exp(res_r_model$summary.hyperpar[2, 1]); exp(res_c_model$summary.hyperpar[2, 1]); exp(res_r_model_with_replicates$summary.hyperpar[2, 1]); exp(res_c_model_with_replicates$summary.hyperpar[2, 1])
σ; res_builtin$summary.hyperpar[3, 1]; exp(res_r_model$summary.hyperpar[3, 1]); exp(res_c_model$summary.hyperpar[3, 1]); exp(res_r_model_with_replicates$summary.hyperpar[3, 1]); exp(res_c_model_with_replicates$summary.hyperpar[3, 1])

res_builtin$cpu; res_r_model$cpu; res_c_model$cpu; res_r_model_with_replicates$cpu; res_c_model_with_replicates$cpu
res_builtin$mlik[1]; res_r_model$mlik[1]; res_c_model$mlik[1]; res_r_model_with_replicates$mlik[1]; res_c_model_with_replicates$mlik[1]

# -----------------------------------------------------------------------
# Test out the model y = Z + ε with multiple repls
# -----------------------------------------------------------------------

loc = as.matrix(expand.grid(x = seq(2, 10, by = .9), y = seq(2, 10, by = .9)))
σ = 2
ρ = 5
n_train = 50

# Create the mesh
mesh = inla.mesh.2d(
  loc.domain = loc,
  min.angle = 30,
  boundary = list(inla.nonconvex.hull(loc, convex = -.3),
                  inla.nonconvex.hull(loc, convex = -1)),
  #cutoff = 1,
  cutoff = .7,
  #max.edge = c(2, 5))
  max.edge = c(2, 8))
#plot(mesh); points(loc)

# Place s0
s0 = matrix(c(4.5, 5.5), nrow = 1)

# extract the precision matrix
spde_tmp = inla.spde2.pcmatern(mesh,
                               prior.range = c(ρ, 0.5),
                               prior.sigma = c(σ, 0.5))
Q = inla.spde2.precision(spde_tmp, theta = c(log(ρ), log(σ)))

# Compute distances
dist_mesh_to_s0 = as.matrix(dist(rbind(s0, mesh$loc[, -3])))[-1, 1]

# Sample y = Z + ε,
# with ε Gaussian white noise, and Z a Matern-field from the SPDE
A_minimal = inla.spde.make.A(mesh, loc)
constr = list(A = as.matrix(inla.spde.make.A(mesh, s0)), e = 0)
z = inla.qsample(n = n_train, Q = Q, constr = constr)
y = as.matrix(A_minimal %*% z)
y = y + rnorm(length(y), sd = σ_ε)

# Create the SPDE models
prior_range = c(ρ, .5)
prior_sigma = c(σ, .5)
init = c(log(ρ), log(σ))
spde_builtin = inla.spde2.pcmatern(
  mesh,
  prior.range = prior_range,
  prior.sigma = prior_sigma)
spde_r_model = spde_generic_model(
  spde = spde_builtin,
  init = init,
  rho_prior = prior_range,
  sigma_prior = prior_sigma,
  use_c = FALSE)
spde_c_model = spde_generic_model(
  spde = spde_builtin,
  init = init,
  rho_prior = prior_range,
  sigma_prior = prior_sigma,
  use_c = TRUE)
spde_r_model_with_replicates = spde_generic_model_with_replicates(
  spde = spde_builtin,
  init = init,
  n_repl = n_train,
  rho_prior = prior_range,
  sigma_prior = prior_sigma,
  use_c = FALSE)
spde_c_model_with_replicates = spde_generic_model_with_replicates(
  spde = spde_builtin,
  init = init,
  n_repl = n_train,
  rho_prior = prior_range,
  sigma_prior = prior_sigma,
  use_c = TRUE)

# Create everything else needed for running INLA
constr2 = list(A = as.matrix(inla.spde.make.A(mesh, s0, index = rep(1, n_train), repl = 1:n_train)), e = rep(0, n_train))
formula_builtin = y ~ -1 + f(spatial, model = spde_builtin, replicate = spatial.repl, nrep = n_train, extraconstr = constr)
formula_r_model = y ~ -1 + f(spatial, model = spde_r_model, replicate = spatial.repl, nrep = n_train, extraconstr = constr)
formula_c_model = y ~ -1 + f(spatial, model = spde_c_model, replicate = spatial.repl, nrep = n_train, extraconstr = constr)
formula_r_model_with_replicates = y ~ -1 + f(spatial, model = spde_r_model_with_replicates, extraconstr = constr2)
formula_c_model_with_replicates = y ~ -1 + f(spatial, model = spde_c_model_with_replicates, extraconstr = constr2)

effects = list(
  spatial = inla.spde.make.index(
    "spatial",
    n.spde = spde_builtin$n.spde,
    n.repl = n_train))

effects_with_replicates = list(spatial = seq_len(n_train * mesh$n))

A = inla.spde.make.A(
  mesh = mesh,
  loc = loc,
  index = rep(1:nrow(loc), n_train),
  repl = rep(1:n_train, each = nrow(loc)))

stack = inla.stack(
  data = list(y = as.numeric(y)),
  A = list(spatial = A),
  effects = effects)

stack_with_replicates = inla.stack(
  data = list(y = as.numeric(y)),
  A = list(spatial = A),
  effects = effects_with_replicates)

# Run INLA
res_builtin = inla(formula_builtin,
                   data = inla.stack.data(stack),
                   control.predictor = list(A = inla.stack.A(stack)),
                   inla.mode = "experimental",
                   only.hyperparam = TRUE,
                   control.inla = list(control.vb = list(enable = FALSE)),
                   num.threads = num_cores,
                   verbose = verbose)

res_r_model = inla(formula_r_model,
                   data = inla.stack.data(stack),
                   control.predictor = list(A = inla.stack.A(stack)),
                   inla.mode = "experimental",
                   only.hyperparam = TRUE,
                   control.inla = list(control.vb = list(enable = FALSE)),
                   num.threads = num_cores,
                   verbose = verbose)

res_c_model = inla(formula_c_model,
                   data = inla.stack.data(stack),
                   control.predictor = list(A = inla.stack.A(stack)),
                   inla.mode = "experimental",
                   only.hyperparam = TRUE,
                   control.inla = list(control.vb = list(enable = FALSE)),
                   num.threads = num_cores,
                   verbose = verbose)

res_r_model_with_replicates = inla(formula_r_model_with_replicates,
                   data = inla.stack.data(stack_with_replicates),
                   control.predictor = list(A = inla.stack.A(stack_with_replicates)),
                   inla.mode = "experimental",
                   only.hyperparam = TRUE,
                   control.inla = list(control.vb = list(enable = FALSE)),
                   num.threads = num_cores,
                   verbose = verbose)

res_c_model_with_replicates = inla(formula_c_model_with_replicates,
                   data = inla.stack.data(stack_with_replicates),
                   control.predictor = list(A = inla.stack.A(stack_with_replicates)),
                   inla.mode = "experimental",
                   only.hyperparam = TRUE,
                   control.inla = list(control.vb = list(enable = FALSE)),
                   num.threads = num_cores,
                   verbose = verbose)


σ_ε; res_builtin$summary.hyperpar[1, 1]^-.5; res_r_model$summary.hyperpar[1, 1]^-.5; res_c_model$summary.hyperpar[1, 1]^-.5; res_r_model_with_replicates$summary.hyperpar[1, 1]^-.5; res_c_model_with_replicates$summary.hyperpar[1, 1]^-.5
ρ; res_builtin$summary.hyperpar[2, 1]; exp(res_r_model$summary.hyperpar[2, 1]); exp(res_c_model$summary.hyperpar[2, 1]); exp(res_r_model_with_replicates$summary.hyperpar[2, 1]); exp(res_c_model_with_replicates$summary.hyperpar[2, 1])
σ; res_builtin$summary.hyperpar[3, 1]; exp(res_r_model$summary.hyperpar[3, 1]); exp(res_c_model$summary.hyperpar[3, 1]); exp(res_r_model_with_replicates$summary.hyperpar[3, 1]); exp(res_c_model_with_replicates$summary.hyperpar[3, 1])

res_builtin$cpu; res_r_model$cpu; res_c_model$cpu; res_r_model_with_replicates$cpu; res_c_model_with_replicates$cpu
res_builtin$mlik[1]; res_r_model$mlik[1]; res_c_model$mlik[1]; res_r_model_with_replicates$mlik[1]; res_c_model_with_replicates$mlik[1]

# The `with_replicates` models are much slower, which is a shame, since we have to use them when including β


# -----------------------------------------------------------------------
# Test out the model y = α * y0 + Z + ε
# -----------------------------------------------------------------------

# Choose a design for the simulation study
loc = as.matrix(expand.grid(x = seq(2, 10, by = 2), y = seq(2, 10, by = 2)))
num_cores = 6
n_train = 200
n_train = 20
σ_ε = .1
λ = 4.3
κ = 1.8
σ = 2
ρ = 4.2
threshold = 5
verbose = TRUE

# Place s0
s0 = matrix(c(4.5, 5.5), nrow = 1)

# Compute distances
dist = as.matrix(dist(loc))
dist_to_s0 = as.matrix(dist(rbind(s0, loc)))[-1, 1]

# Create the mesh
mesh = inla.mesh.2d(
  loc.domain = loc,
  min.angle = 30,
  boundary = list(inla.nonconvex.hull(loc, convex = -.3),
                  inla.nonconvex.hull(loc, convex = -1)),
  cutoff = 1,
  max.edge = c(2, 5))
# plot(mesh); points(loc)

# extract the precision matrix
spde_tmp = inla.spde2.pcmatern(mesh,
                               prior.range = c(ρ, 0.5),
                               prior.sigma = c(σ, 0.5))
Q = inla.spde2.precision(spde_tmp, theta = c(log(ρ), log(σ)))

# Sample y0
y0 = threshold + rexp(n_train)

a = get_a_funcs(log(c(λ, κ)))$value

# Create y
A_minimal = inla.spde.make.A(mesh, loc)
constr = list(A = as.matrix(inla.spde.make.A(mesh, s0)), e = 0)
z = inla.qsample(n = n_train, Q = Q, constr = constr)
y = as.matrix(A_minimal %*% z)
y = y + rnorm(length(y), sd = σ_ε)
y = y + a(y0, dist_to_s0)

# Create the alpha_model
alpha_c_model = alpha_generic_model(
  y0 = y0,
  dist_to_s0 = dist_to_s0,
  n = length(y),
  init = c(log(λ), log(κ)),
  log_lambda_prior = c(log(λ), 1),
  log_kappa_prior = c(log(κ), 1),
  use_c = TRUE)

# Create the SPDE models
prior_range = c(ρ, .5)
prior_sigma = c(σ, .5)
init = c(log(ρ), log(σ))
spde_builtin = inla.spde2.pcmatern(
  mesh,
  prior.range = prior_range,
  prior.sigma = prior_sigma)
spde_c_model = spde_generic_model(
  spde = spde_builtin,
  init = init,
  rho_prior = prior_range,
  sigma_prior = prior_sigma)
spde_c_model_with_replicates = spde_generic_model_with_replicates(
  spde = spde_builtin,
  init = init,
  n_repl = n_train,
  rho_prior = prior_range,
  sigma_prior = prior_sigma)

# Create everything else needed for running INLA
constr2 = list(A = as.matrix(inla.spde.make.A(mesh, s0, index = rep(1, n_train), repl = 1:n_train)), e = rep(0, n_train))
formula_builtin = y ~ -1 + f(idx, model = alpha_c_model) + f(spatial, model = spde_builtin, replicate = spatial.repl, nrep = n_train, extraconstr = constr)
formula_c = y ~ -1 + f(idx, model = alpha_c_model) + f(spatial, model = spde_c_model, replicate = spatial.repl, nrep = n_train, extraconstr = constr)
formula_c_model_with_replicates = y ~ -1 + f(idx, model = alpha_c_model) + f(spatial, model = spde_c_model_with_replicates, extraconstr = constr2)

effects = list(
  spatial = inla.spde.make.index(
    "spatial",
    n.spde = mesh$n,
    n.repl = n_train),
  idx = seq_along(y))

effects_with_replicates = list(
  spatial = seq_len(n_train * mesh$n),
  idx = seq_along(y))

A = inla.spde.make.A(
  mesh = mesh,
  loc = loc,
  index = rep(1:nrow(loc), n_train),
  repl = rep(1:n_train, each = nrow(loc)))

stack = inla.stack(data = list(y = as.numeric(y)), A = list(spatial = A, 1), effects = effects)
stack_with_replicates = inla.stack(data = list(y = as.numeric(y)), A = list(spatial = A, 1), effects = effects_with_replicates)

# Run INLA
res_builtin = inla(
  formula_builtin,
  data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack)),
  num.threads = num_cores,
  inla.mode = "experimental",
  verbose = verbose)

res_c = inla(
  formula_c,
  data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack)),
  num.threads = num_cores,
  inla.mode = "experimental",
  verbose = verbose)

res_c_with_replicates = inla(
  formula_c_model_with_replicates,
  data = inla.stack.data(stack_with_replicates),
  control.predictor = list(A = inla.stack.A(stack_with_replicates)),
  num.threads = num_cores,
  inla.mode = "experimental",
  verbose = verbose)

σ_ε; res_builtin$summary.hyperpar[1, 1]^-.5; res_c$summary.hyperpar[1, 1]^-.5;  res_c_with_replicates$summary.hyperpar[1, 1]^-.5
log(λ); res_builtin$summary.hyperpar[2, 1]; res_c$summary.hyperpar[2, 1]; res_c_with_replicates$summary.hyperpar[2, 1]
log(κ); res_builtin$summary.hyperpar[3, 1]; res_c$summary.hyperpar[3, 1]; res_c_with_replicates$summary.hyperpar[3, 1]
ρ; res_builtin$summary.hyperpar[4, 1]; exp(res_c$summary.hyperpar[4, 1]); exp(res_c_with_replicates$summary.hyperpar[4, 1])
σ; res_builtin$summary.hyperpar[5, 1]; exp(res_c$summary.hyperpar[5, 1]); exp(res_c_with_replicates$summary.hyperpar[5, 1])

res_builtin$cpu; res_c$cpu; res_c_with_replicates$cpu


# -----------------------------------------------------------------------
# Test out the model y = b(y0, |s - s0|) * Z + ε
# -----------------------------------------------------------------------

# Choose a design for the simulation study
loc = as.matrix(expand.grid(x = seq(2, 10, by = 2), y = seq(2, 10, by = 2)))
num_cores = 6
n_train = 100
σ_ε = .1
σ = 2
ρ = 4.2
λβ = 3.6
κβ = .5
threshold = 5
verbose = TRUE

# Place s0
s0 = matrix(c(4.5, 5.5), nrow = 1)

# Create the mesh
mesh = inla.mesh.2d(
  loc.domain = loc,
  min.angle = 30,
  boundary = list(inla.nonconvex.hull(loc, convex = -.3),
                  inla.nonconvex.hull(loc, convex = -1)),
  cutoff = 1,
  max.edge = c(2, 5))
# plot(mesh); points(loc)

# Compute distances
dist = as.matrix(dist(loc))
dist_to_s0 = as.matrix(dist(rbind(s0, loc)))[-1, 1]
dist_to_s0_from_mesh = as.matrix(dist(rbind(s0, mesh$loc[, -3])))[-1, 1]

# extract the precision matrix
spde_tmp = inla.spde2.pcmatern(mesh,
                               prior.range = c(ρ, 0.5),
                               prior.sigma = c(σ, 0.5))
Q = inla.spde2.precision(spde_tmp, theta = c(log(ρ), log(σ)))

# Sample y0
y0 = threshold + rexp(n_train)

b = get_b_funcs(log(c(λβ, κβ)))$value

# Create y
A_minimal = inla.spde.make.A(mesh, loc)
constr = list(A = as.matrix(inla.spde.make.A(mesh, s0)), e = 0)
z = inla.qsample(n = n_train, Q = Q, constr = constr)
y = as.matrix(A_minimal %*% z)
for (i in seq_along(y0)) {
  y[, i] = y[, i] * b(y0[i], dist_to_s0)
}
y = y + rnorm(length(y), sd = σ_ε)

# Create the SPDE models
priors = list(
  rho = c(ρ, .5),
  sigma = c(σ, .5),
  lambda = c(log(λβ), 1),
  kappa = c(log(κβ), 1))
init = log(c(ρ, σ, λβ, κβ))
spde = inla.spde2.pcmatern(
  mesh,
  prior.range = priors$rho,
  prior.sigma = priors$sigma)
spde_r_model = spde_generic_model_with_b_func(
  spde = spde,
  y0 = y0,
  init = init,
  priors = priors,
  dist_to_s0 = dist_to_s0_from_mesh,
  use_c = FALSE)
spde_c_model = spde_generic_model_with_b_func(
  spde = spde,
  y0 = y0,
  init = init,
  priors = priors,
  dist_to_s0 = dist_to_s0_from_mesh)
spde_c_model_const_sigma = spde_generic_model_with_b_func_const_sigma(
  spde = spde,
  log_sigma = log(σ),
  y0 = y0,
  init = init[-2],
  priors = priors[-2],
  dist_to_s0 = dist_to_s0_from_mesh)

# Create everything else needed for running INLA
constr2 = list(A = as.matrix(inla.spde.make.A(mesh, s0, index = rep(1, n_train), repl = 1:n_train)), e = rep(0, n_train))
formula_r = y ~ -1 + f(spatial, model = spde_r_model, extraconstr = constr2)
formula_c = y ~ -1 + f(spatial, model = spde_c_model, extraconstr = constr2)
formula_c_const_sigma = y ~ -1 + f(spatial, model = spde_c_model_const_sigma, extraconstr = constr2)

effects = list(spatial = seq_len(spde$n.spde * n_train))

A = inla.spde.make.A(
  mesh = mesh,
  loc = loc,
  index = rep(1:nrow(loc), n_train),
  repl = rep(1:n_train, each = nrow(loc)))

stack = inla.stack(data = list(y = as.numeric(y)), A = list(spatial = A), effects = effects)

# Run INLA
res_r = inla(
  formula_r,
  data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack)),
  num.threads = num_cores,
  inla.mode = "experimental",
  verbose = verbose)

res_c = inla(
  formula_c,
  data = inla.stack.data(stack),
  control.predictor = list(A = inla.stack.A(stack)),
  num.threads = num_cores,
  inla.mode = "experimental",
  verbose = verbose)

σ_ε; res_r$summary.hyperpar[1, 1]^-.5; res_c$summary.hyperpar[1, 1]^-.5
ρ; exp(res_r$summary.hyperpar[2, 1]); exp(res_c$summary.hyperpar[2, 1])
σ; exp(res_r$summary.hyperpar[3, 1]); exp(res_c$summary.hyperpar[3, 1])
λβ; exp(res_r$summary.hyperpar[4, 1]); exp(res_c$summary.hyperpar[4, 1])
κβ; exp(res_r$summary.hyperpar[5, 1]); exp(res_c$summary.hyperpar[5, 1])

res_r$cpu; res_c$cpu
res_r$mlik[1]; res_c$mlik[1]

res_c$summary.hyperpar
c(log(σ_ε^-2), log(ρ), log(σ), log(λβ), log(κβ))
