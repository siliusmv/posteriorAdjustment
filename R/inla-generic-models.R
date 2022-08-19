#' Function for defininc a cgeneric model for a.
#'
#' The input variables are:
#' y0: A vector of equal length to the data (y_inla) used for running R-INLA,
#'   such that y0[i] contains the value at the conditioning site for the same
#'   time as when y_inla[i] was observed.
#' dist_to_s0: A vector of equal length to y_inla, such that dist_to_s0[i] contains
#'   the distance from the location where y_inla[i] was observed, to the conditioning
#'   site where y0[i] was observed.
#' init: Initial values for the subset of (log(λ), log(κ)) that will be estimated.
#'   See the documentation for the prior variable for more information.
#' priors: Priors for the two parameters of a. This is a list of length 2,
#'   with names "lambda" and "kappa". The two elements are vectors of length 2 or 3.
#'   If the first element of a vector is 0, then we will not estimate that parameter,
#'   but we will fix it equal to the value given in the second element of the vector.
#'   If the first element is not 0, then the logarithm of the parameter is given a
#'   Gaussian prior with mean equal to the second element of the vector, and standard
#'   deviation equal to the third element of the vector. If e.g.
#'   priors = list(lambda = c(0, 4), kappa = c(1, 1, 3)), then we fix log(λ) = 4,
#'   while we estimate log(κ) and give it the prior log(κ) ~ N(1, 3).
#' debug: A boolean stating if R-INLA should print debug information or not.
#' @export
a_generic_model = function(y0,
                           dist_to_s0,
                           init,
                           priors,
                           debug = FALSE) {
  stopifnot(length(y0) == length(dist_to_s0))
  stopifnot(all(c("lambda", "kappa") %in% names(priors)))
  stopifnot(length(priors) == 2)
  args = list(debug = debug)

  # The name and location of the required c-function for defining the cgeneric model
  args$model = "inla_cgeneric_a_model"
  args$shlib = file.path(cgeneric_dir(), "a.so")

  # Check that the priors have the correct lengths
  for (i in seq_along(priors)) {
    if (priors[[i]][1] == 0) {
      stopifnot(length(priors[[i]]) == 2)
    } else {
      stopifnot(length(priors[[i]]) == 3)
    }
  }

  # Check that the init vector has the correct length
  n_fixed = sum(sapply(priors, `[`, 1) == 0)
  stopifnot(length(init) == 2 - n_fixed)

  # Put all the arguments into args in the order defined in the c-function
  args$n = length(y0)
  args$y0 = y0
  args$dist_to_s0 = dist_to_s0
  args$init = init
  args$lambda_prior = priors$lambda
  args$kappa_prior = priors$kappa

  # Define the model
  do.call(INLA::inla.cgeneric.define, args)
}

#' Function for defininc a cgeneric model for Z_b.
#'
#' The input variables are:
#' spde: An inla.spde2 object that contains the matrices M0, M1, M2, B0, B1 and B2,
#'   that are necessary for computing the precision matrix of the Gaussian Matérn
#'   field used for creating the SPDE approximation.
#' n: A vector of length equal to the number of unique conditioning sites used
#'   for inference. Element nr. i of n contains the number of threshold exceedances
#'   found at conditioning site nr. i.
#' init: Initial values for the subset of (log(ρ), log(σ), log(ρ_b)) that will be
#'   estimated. See the documentation for the prior variable for more information.
#' priors: Priors for the three parameters of Z_b. This is a list of length 3,
#'   with names "rho", "sigma" and "rho_b". The three list elements are vectors of
#'   length 2 or 3. If the first element of a vector is 0, then we will not estimate
#'   that parameter, but we will fix it equal to the value given in the second element
#'   of the vector. If the first element is not 0, then the vector has length 3, and
#'   the two remaining values are used for placing a prior on that parameter.
#'   ρ and σ are given a joint PC prior, see ?inla.spde2.pcmatern for more information.
#'   log(ρ_b) is given a Gaussian prior. If e.g.
#'   priors = list(rho = c(1, 60, .95), sigma = c(1, 5, .05), rho_b = c(1, 2, 3)), then
#'   ρ and σ are given a PC prior such that P(ρ < 60) = .95 and P(σ > 5) = .05,
#'   while log(ρ_b) is given a Gaussian prior with mean 2 and standard deviation 3.
#'   If, on the other hand, priors$rho_b = c(0, 2), then the value of log(ρ_b) is
#'   fixed equal to 2.
#' dist_to_s0: A matrix of dimension (n x m), where m is the number of nodes in the
#'   triangular mesh, and n is the number of unique conditioning sites. Row nr. i
#'   of the matrix contains the distances from all mesh nodes to the ith conditioning
#'   site. dist_to_s0 can also be a vector, if only one conditioning site is used.
#' debug: A boolean stating if R-INLA should print debug information or not.
#' @export
spde_generic_model_with_b_func = function(spde,
                                          n,
                                          init,
                                          priors,
                                          dist_to_s0,
                                          debug = FALSE) {
  stopifnot(all(c("rho", "sigma", "rho_b") %in% names(priors)))
  if (!is.matrix(dist_to_s0)) dist_to_s0 = matrix(dist_to_s0, nrow = 1)
  stopifnot(length(n) == nrow(dist_to_s0))
  stopifnot(spde$n.spde == ncol(dist_to_s0))
  stopifnot(length(priors) == 3)
  args = list(debug = debug)

  # The name and location of the required c-function for defining the cgeneric model
  args$model = "inla_cgeneric_spde_model_with_b_func"
  args$shlib = file.path(cgeneric_dir(), "b.so")

  # Check that the priors have the correct lengths
  for (i in seq_along(priors)) {
    if (priors[[i]][1] == 0) {
      stopifnot(length(priors[[i]]) == 2)
    } else {
      stopifnot(length(priors[[i]]) == 3)
    }
  }

  # Check that the init vector has the correct length
  n_fixed = sum(sapply(priors, `[`, 1) == 0)
  stopifnot(length(init) == 3 - n_fixed)

  # Put all the arguments into args in the order defined in the c-function
  args$n = spde$n.spde * sum(n)
  args$B0 = spde$param.inla$B0
  args$B1 = spde$param.inla$B1
  args$B2 = spde$param.inla$B2
  args$M0 = spde$param.inla$M0
  args$M1 = spde$param.inla$M1
  args$M2 = spde$param.inla$M2
  args$init = init
  args$rho_prior = priors$rho
  args$sigma_prior = priors$sigma
  args$rho_b_prior = priors$rho_b
  args$dist_to_s0 = dist_to_s0
  args$s0_index = rep(seq_along(n), n)

  # Define the model
  do.call(INLA::inla.cgeneric.define, args)
}


#' @export
iid_generic_model_with_b_func = function(init,
                                         priors,
                                         dist_to_s0,
                                         debug = FALSE) {
  stopifnot(all(c("sigma", "rho_b") %in% names(priors)))
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
  args$rho_b_prior = priors$rho_b
  do.call(func, args)
}
