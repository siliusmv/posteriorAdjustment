# library(mvtnorm)
# library(numDeriv)
# library(INLA)
# library(spatialEVT)
# 
# library(RhpcBLASctl)
# RhpcBLASctl::blas_set_num_threads(1)
# RhpcBLASctl::omp_set_num_threads(1)
# 
# loc = expand.grid(seq(2, 10, by = 2), seq(2, 10, by = 2)) |>
#   as.matrix()
# s0 = c(4.5, 5.3)
# n_loc = nrow(loc)
# dist = as.matrix(dist(loc))
# dist_to_s0 = as.matrix(dist(rbind(s0, loc)))[-1, 1]
# 
# mesh = INLA::inla.mesh.2d(
#   loc = loc,
#   #min.angle = 30,
#   boundary = list(INLA::inla.nonconvex.hull(loc, convex = -.3),
#                   INLA::inla.nonconvex.hull(loc, convex = -1)),
#   cutoff = 3.2,
#   max.edge = c(4, 5))
# #plot(mesh); points(loc)
# 
# spde = INLA::inla.spde2.pcmatern(mesh, prior.range = c(1, .5), prior.sigma = c(1, .5))
# 
# ρ = runif(1, 1, 10)
# σ = runif(1, 1, 5)
# λ = runif(1, 2, 5)
# κ = runif(1, .7, 1.3)
# λβ = runif(1, 2, 5)
# κβ = runif(1, .2, .8)
# σ_ε = runif(1, .01, .5)
# τ = σ_ε^-2
# 
# θ_a = log(c(λ, κ))
# θ_b = log(c(λβ, κβ))
# θ_Σ = log(c(ρ, σ))
# 
# A = INLA::inla.spde.make.A(mesh, loc)
# A_constr = as.matrix(INLA::inla.spde.make.A(mesh, matrix(s0, nrow = 1)))
# threshold = 4
# 
# n_repl = 20
# a_funcs = get_a_funcs(θ_a, jacobian = TRUE)
# b_funcs = get_b_funcs(θ_b, jacobian = TRUE)
# Σ_funcs = get_Σ_funcs(θ_Σ, spde, jacobian = TRUE)
# 
# y0 = rexp(n_repl) + threshold
# y = rconditional(
#   y0 = y0,
#   a_func = a_funcs$value,
#   b_func = b_funcs$value,
#   Q = Σ_funcs$Q_unconstrained,
#   τ = τ,
#   dist = dist_to_s0,
#   A = A,
#   A_constr = A_constr)
# 
# # Insert some NAs into y
# ss = sample.int(length(y), 7)
# y[ss] = NA
# 
# test_that("derivative of Σ is correct", {
#   jac1 = get_Σ_funcs(θ_Σ, spde, TRUE)$jacobian(A, A_constr) |>
#                                     matrix(ncol = length(θ_Σ))
#   jac2 = numDeriv::jacobian(\(θ) get_Σ_funcs(θ, spde)$value(A, A_constr), x = θ_Σ)
#   expect_equal(jac1, jac2, tolerance = 1e-5)
# })
# 
# 
# test_that("derivative of b is correct", {
#   jac1 = get_b_funcs(θ_b, TRUE)$jacobian(y0, n_loc) |>
#                               matrix(ncol = length(θ_b))
#   jac2 = numDeriv::jacobian(\(θ) get_b_funcs(θ)$value(y0, n_loc), x = θ_b)
#   expect_equal(jac1, jac2, tolerance = 1e-5)
# })
# 
# test_that("computation of conditional ll is correct", {
#   ll1 = ll_conditional(
#     y = y,
#     y0 = y0,
#     a_func = a_funcs$value,
#     b_func = b_funcs$value,
#     Σ_func = Σ_funcs$value,
#     dist = dist_to_s0,
#     A = A,
#     A_constr = A_constr,
#     τ = τ)
#   ll2 = rep(NA, n_repl)
#   for (i in 1:n_repl) {
#     good_index = which(!is.na(y[, i]))
#     tmp_b = y0[i] ^ exp(- (dist_to_s0[good_index] / λβ) ^ κβ)
#     ll2[i] = mvtnorm::dmvnorm(
#       x = y[good_index, i],
#       mean = y0[i] * α_func(λ, κ, dist_to_s0[good_index]),
#       sigma = rep(tmp_b, each = length(good_index)) * rep(tmp_b, length(good_index)) *
#         Σ_funcs$value(A[good_index, ], A_constr) + diag(1 / τ, length(good_index)),
#       log = TRUE)
#   }
#   ll3 = ll_conditional(
#     y = y,
#     y0 = y0,
#     a_func = a_funcs$value,
#     b_func = b_funcs$value,
#     Σ_func = Σ_funcs$value,
#     dist = dist_to_s0,
#     A = A,
#     A_constr = A_constr,
#     τ = τ,
#     na.rm = FALSE)
#   ll4 = rep(NA, n_repl)
#   for (i in 1:n_repl) {
#     tmp_b = y0[i] ^ exp(- (dist_to_s0 / λβ) ^ κβ)
#     ll4[i] = mvtnorm::dmvnorm(
#       x = y[, i],
#       mean = y0[i] * α_func(λ, κ, dist_to_s0),
#       sigma = rep(tmp_b, each = length(dist_to_s0)) * rep(tmp_b, length(dist_to_s0)) *
#         Σ_funcs$value(A, A_constr) + diag(1 / τ, n_loc),
#       log = TRUE)
#   }
#   expect_equal(ll1, ll2, tolerance = 1e-5)
#   expect_equal(ll3, ll4, tolerance = 1e-5)
#   na_index = which(is.na(ll3))
#   expect_equal(ll1[-na_index], ll3[-na_index], tolerance = 1e-5)
# })
# 
# test_that("conditional ll gradient is correct", {
#   ll_grad1 = ll_conditional_grad(
#     y = y,
#     y0 = y0,
#     a_funcs = get_a_funcs(θ_a, TRUE),
#     b_funcs = get_b_funcs(θ_b, TRUE),
#     Σ_funcs = get_Σ_funcs(θ_Σ, spde, TRUE),
#     τ = τ,
#     dist = dist_to_s0,
#     A = A,
#     A_constr = A_constr)
#   ll_grad2 = numDeriv::jacobian(
#     func = function(θ) {
#       ll_conditional(
#         y = y,
#         y0 = y0,
#         a_func = get_a_funcs(θ[1:2])$value,
#         b_func = get_b_funcs(θ[3:4])$value,
#         Σ_func = get_Σ_funcs(θ[5:6], spde)$value,
#         dist = dist_to_s0,
#         A = A,
#         A_constr = A_constr,
#         τ = exp(θ[7]))
#     },
#     x = c(θ_a, θ_b, θ_Σ, log(τ)))
#   ll_grad3 = ll_conditional_grad(
#     y = y,
#     y0 = y0,
#     a_funcs = get_a_funcs(θ_a, TRUE),
#     b_funcs = get_b_funcs(θ_b, TRUE),
#     Σ_funcs = get_Σ_funcs(θ_Σ, spde, TRUE),
#     τ = τ,
#     dist = dist_to_s0,
#     A = A,
#     A_constr = A_constr,
#     na.rm = FALSE)
#   ll_grad4 = numDeriv::jacobian(
#     func = function(θ) {
#       ll_conditional(
#         y = y,
#         y0 = y0,
#         a_func = get_a_funcs(θ[1:2])$value,
#         b_func = get_b_funcs(θ[3:4])$value,
#         Σ_func = get_Σ_funcs(θ[5:6], spde)$value,
#         dist = dist_to_s0,
#         A = A,
#         A_constr = A_constr,
#         τ = exp(θ[7]),
#         na.rm = FALSE)
#     },
#     x = c(θ_a, θ_b, θ_Σ, log(τ)))
#   expect_equal(ll_grad1, ll_grad2, tolerance = 1e-5)
#   expect_equal(ll_grad3, ll_grad4, tolerance = 1e-5)
#   na_index = which(is.na(ll_grad3))
#   expect_equal(ll_grad1[-na_index], ll_grad3[-na_index], tolerance = 1e-5)
# })
# 
# 
# # Try again with multiple s0-locations
# # ---------------------------------------------------------------------------
# n_s0 = 8
# s0_locations = sort(sample.int(n_loc, n_s0))
# A = lapply(s0_locations, function(i) INLA::inla.spde.make.A(mesh, loc[-i, ]))
# A_constr = lapply(s0_locations, function(i) as.matrix(INLA::inla.spde.make.A(mesh, loc[i, , drop = FALSE])))
# dist_to_s0 = lapply(s0_locations, function(i) dist[-i, i])
# 
# n_repl_per_s0 = 4
# 
# a_funcs = get_a_funcs(θ_a, TRUE)
# b_funcs = get_b_funcs(θ_b, TRUE)
# Σ_funcs = get_Σ_funcs(θ_Σ, spde, TRUE)
# 
# y0 = lapply(1:n_s0, \(i) rexp(n_repl_per_s0) + threshold)
# y = rconditional(
#   y0 = y0,
#   a_func = a_funcs$value,
#   b_func = b_funcs$value,
#   Q = Σ_funcs$Q,
#   τ = τ,
#   dist = dist_to_s0,
#   A = A,
#   A_constr = A_constr)
# 
# test_that("conditional ll with multiple locations works", {
#   ll1 = ll_conditional(
#     y = y,
#     y0 = y0,
#     a_func = a_funcs$value,
#     b_func = b_funcs$value,
#     Σ_func = Σ_funcs$value,
#     A = A,
#     A_constr = A_constr,
#     dist = dist_to_s0,
#     τ = τ)
#   ll2 = NULL
#   for (i in 1:n_s0) {
#     if (!is.null(ncol(y[[i]]))) {
#       for (j in 1:n_repl_per_s0) {
#         tmp_b = y0[[i]][j] ^ exp(- (dist_to_s0[[i]] / λβ) ^ κβ)
#         ll2[length(ll2) + 1] = mvtnorm::dmvnorm(
#           x = y[[i]][, j],
#           mean = y0[[i]][j] * α_func(λ, κ, dist_to_s0[[i]]),
#           sigma = rep(tmp_b, each = length(dist_to_s0[[i]])) * rep(tmp_b, length(dist_to_s0[[i]])) *
#             as.matrix(A[[i]] %*% Σ_funcs$value(A_constr = A_constr[[i]]) %*% Matrix::t(A[[i]])) +
#             diag(1 / τ, n_loc - 1),
#           log = TRUE)
#       }
#     }
#   }
#   expect_equal(ll1, ll2, tolerance = 1e-5)
# })
# 
# test_that("conditional ll_grad with multiple locations works", {
#   ll_grad1 = ll_conditional_grad(
#     y = y,
#     y0 = y0,
#     a_funcs = get_a_funcs(θ_a, TRUE),
#     b_funcs = get_b_funcs(θ_b, TRUE),
#     Σ_funcs = get_Σ_funcs(θ_Σ, spde, TRUE),
#     τ = τ,
#     dist = dist_to_s0,
#     A = A,
#     A_constr = A_constr)
#   ll_grad2 = numDeriv::jacobian(
#     func = function(θ) {
#       ll_conditional(
#         y = y,
#         y0 = y0,
#         a_func = get_a_funcs(θ[1:2])$value,
#         b_func = get_b_funcs(θ[3:4])$value,
#         Σ_func = get_Σ_funcs(θ[5:6], spde)$value,
#         A = A,
#         dist = dist_to_s0,
#         A_constr = A_constr,
#         τ = exp(θ[7]))
#     },
#     x = c(θ_a, θ_b, θ_Σ, log(τ)))
#   expect_equal(ll_grad1, ll_grad2, tolerance = 1e-5)
# })
