
#' @export
spde_generic_model = function(spde,
                              init,
                              rho_prior,
                              sigma_prior,
                              debug = FALSE) {
  stopifnot(length(init) == 2)
  stopifnot(length(rho_prior) == 2)
  stopifnot(length(sigma_prior) == 2)
  args = list(debug = debug)
  args$model = "inla_cgeneric_spde_model"
  args$shlib = file.path(cgeneric_dir(), "b.so")
  args$n = spde$n.spde
  args$B0 = spde$param.inla$B0
  args$B1 = spde$param.inla$B1
  args$B2 = spde$param.inla$B2
  args$M0 = spde$param.inla$M0
  args$M1 = spde$param.inla$M1
  args$M2 = spde$param.inla$M2
  func = INLA::inla.cgeneric.define
  args$init = init
  args$rho_prior = rho_prior
  args$sigma_prior = sigma_prior

  do.call(func, args)
}

#' @export
spde_generic_model_with_replicates = function(spde,
                                              init,
                                              rho_prior,
                                              sigma_prior,
                                              n_repl,
                                              debug = FALSE) {
  stopifnot(length(init) == 2)
  stopifnot(length(rho_prior) == 2)
  stopifnot(length(sigma_prior) == 2)
  args = list(debug = debug)
  args$model = "inla_cgeneric_spde_model_with_replicates"
  args$shlib = file.path(cgeneric_dir(), "b.so")
  args$n = spde$n.spde * n_repl
  args$B0 = spde$param.inla$B0
  args$B1 = spde$param.inla$B1
  args$B2 = spde$param.inla$B2
  args$M0 = spde$param.inla$M0
  args$M1 = spde$param.inla$M1
  args$M2 = spde$param.inla$M2
  func = INLA::inla.cgeneric.define
  args$n_repl = as.integer(n_repl)
  args$init = init
  args$rho_prior = rho_prior
  args$sigma_prior = sigma_prior

  do.call(func, args)
}

#' @export
iid_generic_model_with_b_func = function(init,
                                         priors,
                                         dist_to_s0,
                                         debug = FALSE) {
  stopifnot(all(c("sigma", "rho_2") %in% names(priors)))
  stopifnot(length(priors) == 2)
  n_fixed = sum(sapply(priors, `[`, 1) == 0)
  stopifnot(length(init) == 2 - n_fixed)
  for (i in seq_along(priors)) if (priors[[i]][1] != 0) stopifnot(length(priors[[i]]) == 3)
  args = list(debug = debug)
  args$model = "inla_cgeneric_iid_model_with_b_func"
  args$shlib = file.path(cgeneric_dir(), "b.so")
  args$n = length(dist_to_s0)
  func = INLA::inla.cgeneric.define
  args$dist_to_s0 = dist_to_s0
  args$init = init
  args$sigma_prior = priors$sigma
  args$rho_2_prior = priors$rho_2
  do.call(func, args)
}

#' @export
iid_generic_model_with_a_and_b = function(init,
                                          y0,
                                          priors,
                                          dist_to_s0,
                                          debug = FALSE) {
  stopifnot(all(c("lambda", "kappa", "sigma", "rho_2") %in% names(priors)))
  stopifnot(length(priors) == 4)
  n_fixed = sum(sapply(priors, `[`, 1) == 0)
  stopifnot(length(init) == 4 - n_fixed)
  for (i in seq_along(priors)) if (priors[[i]][1] != 0) stopifnot(length(priors[[i]]) == 3)
  args = list(debug = debug)
  args$model = "inla_cgeneric_iid_model_with_a_and_b"
  args$shlib = file.path(cgeneric_dir(), "a.so")
  args$n = length(dist_to_s0)
  func = INLA::inla.cgeneric.define
  args$y0 = y0
  args$dist_to_s0 = dist_to_s0
  args$init = init
  args$lambda_prior = priors$lambda
  args$kappa_prior = priors$kappa
  args$sigma_prior = priors$sigma
  args$rho_2_prior = priors$rho_2
  do.call(func, args)
}


#' @export
spde_generic_model_with_b_func = function(spde,
                                          n,
                                          init,
                                          priors,
                                          dist_to_s0,
                                          debug = FALSE) {
  if (!is.matrix(dist_to_s0)) dist_to_s0 = matrix(dist_to_s0, nrow = 1)
  stopifnot(length(n) == nrow(dist_to_s0))
  stopifnot(spde$n.spde == ncol(dist_to_s0))
  s0_index = rep(seq_along(n), n)
  stopifnot(all(c("rho", "sigma", "rho_2") %in% names(priors)))
  stopifnot(length(priors) == 3)
  n_fixed = sum(sapply(priors, `[`, 1) == 0)
  stopifnot(length(init) == 3 - n_fixed)
  for (i in seq_along(priors)) if (priors[[i]][1] != 0) stopifnot(length(priors[[i]]) == 3)
  args = list(debug = debug)
  args$model = "inla_cgeneric_spde_model_with_b_func"
  args$shlib = file.path(cgeneric_dir(), "b.so")
  args$n = spde$n.spde * sum(n)
  args$B0 = spde$param.inla$B0
  args$B1 = spde$param.inla$B1
  args$B2 = spde$param.inla$B2
  args$M0 = spde$param.inla$M0
  args$M1 = spde$param.inla$M1
  args$M2 = spde$param.inla$M2
  func = INLA::inla.cgeneric.define
  args$init = init
  args$rho_prior = priors$rho
  args$sigma_prior = priors$sigma
  args$rho_2_prior = priors$rho_2
  args$dist_to_s0 = dist_to_s0
  args$s0_index = s0_index

  do.call(func, args)
}


#' @export
alpha_generic_model = function(y0,
                               dist_to_s0,
                               init,
                               priors,
                               debug = FALSE) {
  stopifnot(length(y0) == length(dist_to_s0))
  stopifnot(all(c("lambda", "kappa") %in% names(priors)))
  args = list(debug = debug)
  args$model = "inla_cgeneric_alpha_model"
  args$shlib = file.path(cgeneric_dir(), "a.so")
  func = INLA::inla.cgeneric.define
  for (i in seq_along(priors)) if (priors[[i]][1] != 0) stopifnot(length(priors[[i]]) == 3)
  args$n = length(y0)
  args$y0 = y0
  args$dist_to_s0 = dist_to_s0
  args$init = init
  args$lambda_prior = priors$lambda
  args$kappa_prior = priors$kappa

  do.call(func, args)
}
