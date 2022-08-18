devtools::load_all()
library(microbenchmark)
library(mvtnorm)
library(mvnfast)
library(bayesm)

set.seed(123)
sigma = bayesm::rwishart(10, diag(8))$IW
means = rnorm(8)
X = mvtnorm::rmvnorm(900000, means, sigma)
X_t = t(X)
nugget = 2

vals = list(
  dmvnorm = mvtnorm::dmvnorm(X, means, sigma, log = TRUE),
  dmvnorm_arma = dmvnorm_arma(X_t, means, sigma, TRUE),
  dmvn = mvnfast::dmvn(X, means, sigma, log = TRUE))
sapply(vals, function(x) summary(as.numeric(x - vals[[1]])))

m1 = microbenchmark::microbenchmark(
  dmvnorm = mvtnorm::dmvnorm(X, means, sigma, log = TRUE),
  dmvnorm_arma = dmvnorm_arma(X_t, means, sigma, TRUE),
  dmvn = mvnfast::dmvn(X, means, sigma, log = TRUE))
m1

# ==============================================================================
# Test dconditional
# ==============================================================================

library(parallel)
library(pbapply)
library(INLA)
library(Matrix)
library(sf)

threshold = qlaplace(.95) # The threshold t for defining the conditional extremes model
n_cores = 6 # Run code in parallel
r = 5 # Radius used for computing aggregated empirical distribution functions
radar = readRDS(file.path(downloads_dir(), "radar.rds"))
coords = st_coordinates(radar$coords)
n_loc = nrow(coords)

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

delta_s0 = 4
x_coords = coords[, 1] |> unique() |> sort()
y_coords = coords[, 2] |> unique() |> sort()
s0_locs = expand.grid(
  x = x_coords[seq(delta_s0, length(x_coords), by = delta_s0)],
  y = y_coords[seq(delta_s0, length(y_coords), by = delta_s0)]) |>
  as.matrix()

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

mesh = inla.mesh.2d(
  loc = coords,
  boundary = list(inla.nonconvex.hull(coords, convex = -.1),
                  inla.nonconvex.hull(coords, convex = -1)),
  max.edge = c(5, 30))
spde = inla.spde2.pcmatern(mesh, prior.range = c(60, .95), prior.sigma = c(5, .05))

data$dist_to_s0_from_mesh = list()
for (i in seq_along(data$s0)) {
  data$dist_to_s0_from_mesh[[i]] = dist_euclid(mesh$loc[, 1:2], data$s0[[i]])
}

get_a_func = function(theta) {
  lambda = exp(theta[1])
  kappa = exp(theta[2])
  function(y, dist) {
    alpha = exp(- (dist / lambda)^kappa)
    matrix(rep(y, each = length(alpha)) * rep(alpha, length(y)),
           nrow = length(dist),
           ncol = length(y))
  }
}
get_b_func = function(theta) {
  rho = exp(theta)
  function(y, dist) {
    tmp = dist / rho
    tmp[tmp < 1e-9] = 1e-9
    b = sqrt(1 - exp(-2 * tmp))
    matrix(rep(b, length(y)), nrow = length(dist), ncol = length(y))
  }
}

theta = c(4, -.4, 1, 4, .5, 4)
y = data$y
y0 = data$y0
dist_to_s0 = data$dist_to_s0
dist_to_s0_from_mesh = data$dist_to_s0_from_mesh
A_mats = lapply(data$obs_index, function(x) inla.spde.make.A(mesh, coords[x, ]))
n_cores = n_cores

b_func = get_b_func(theta[3])
cov_mat_func = get_Σ_func(theta[4:5], spde)$value
tau = exp(theta[6])
a_func = get_a_func(theta[1:2])

lengths = c(length(y0), length(y), length(A_mats), length(dist_to_s0))
stopifnot(all(lengths == lengths[1]))
if (!is.null(dist_to_s0_from_mesh)) stopifnot(length(dist_to_s0_from_mesh) == length(y))

i = 1
b = b_func(1, dist_to_s0_from_mesh[[i]])
x = y[[i]] - a_func(y0[[i]], dist_to_s0[[i]])
A = A_mats[[i]]
B = Matrix::Diagonal(length(b), as.numeric(b))
B2 = do.call(cbind, rep(list(b), ncol(x)))
sigma0 = cov_mat_func()
nugget = 1 / tau

args = list(
  x = x,
  A = A,
  #B = B,
  b = b,
  sigma0 = sigma0,
  nugget = nugget,
  logd = TRUE,
  na_rm = TRUE)
args2 = args
args2$b = NULL
args2$B = B2

d1 = local({
  sigma = as.matrix(A %*% B %*% sigma0 %*% B %*% t(A)) + diag(nugget, nrow(A))
  sapply(
    X = 1:ncol(x),
    FUN = function(i) {
      good_index = which(!is.na(x[, i]))
      mvtnorm::dmvnorm(x[good_index, i], sigma = sigma[good_index, good_index], log = TRUE)
    })
})
d2 = local({
  sigma = as.matrix(A %*% B %*% sigma0 %*% B %*% t(A)) + diag(nugget, nrow(A))
  sapply(
    X = 1:ncol(x),
    FUN = function(i) {
      good_index = which(!is.na(x[, i]))
      mvnfast::dmvn(
        x[good_index, i],
        mu = rep(0, length(good_index)),
        sigma = sigma[good_index, good_index],
        log = TRUE)
    })
})
d3 = do.call(dconditional_no_beta, args)
d4 = do.call(dconditional, args2)
summary(as.numeric(d1 - d2))
summary(as.numeric(d1 - d3))
summary(as.numeric(d1 - d4))

m2 = microbenchmark::microbenchmark(
  d1 = local({
    sigma = as.matrix(A %*% B %*% sigma0 %*% B %*% t(A)) + diag(nugget, nrow(A))
    sapply(
      X = 1:ncol(x),
      FUN = function(i) {
        good_index = which(!is.na(x[, i]))
        mvtnorm::dmvnorm(
          x[good_index, i],
          sigma = sigma[good_index, good_index],
          log = TRUE)
      })
  }),
  d2 = local({
    sigma = as.matrix(A %*% B %*% sigma0 %*% B %*% t(A)) + diag(nugget, nrow(A))
    sapply(
      X = 1:ncol(x),
      FUN = function(i) {
        good_index = which(!is.na(x[, i]))
        mvnfast::dmvn(
          x[good_index, i],
          mu = rep(0, length(good_index)),
          sigma = sigma[good_index, good_index],
          log = TRUE)
      })
  }),
  d3 = do.call(dconditional_no_beta, args),
  d4 = do.call(dconditional, args2))
m2
