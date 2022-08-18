
#' Comput the composite log-likelihood for the global conditional extremes model.
#'
#' The input variables are:
#' y: A list of matrices. Matrix nr. i has dimension (d_i x n_i) and contains n_i
#'   d_i-dimensional vectors of observations for conditioning site nr. i.
#' y0: A list of vectors. Vector nr. i has dimension n_i and contains the threshold
#'   exceedances from conditioning site nr. i.
#' a_func: A function a(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of a(y0[j]), dist_to_s0[i]).
#' b_func: A function b(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of b(y0[j], dist_to_s0[i]).
#'   b_func is different from a_func in that a_func requires the distance to s0
#'   for all locations with observations, while b_func requires the distance to s0
#'   for all locations on the triangular mesh used for defining the SPDE.
#' sigma: The covariance matrix of the weights used for building the SPDE approximation.
#' tau: The precision of the nugget effect.
#' dist_to_s0: List containing the distances to s0 for all observation locations,
#'   used for computing a. The list contains one vector of distances for each
#'   conditioning site. Vector nr. i has dimension d_i.
#' dist_to_s0_from_mesh: List containing the distances to s0 for all locations in
#'   the mesh, used for computing b. The list contains one vector of distances for each
#'   conditioning site. Each of the distance vectors must have the same length.
#' A: List of projection matrices used for building the SPDE approximation. The list
#'   contains one projection matrix for each conditioning site.
#'   Each projection matrix must contain the same number of columns, equal to the
#'   lengths of dist_to_s0_from_mesh. Projection matrix nr. i has d_i rows.
#' n_cores: Number of cores to use for computing the log-likelihood in parallel. We
#'   cannot use more cores than the number of conditioning sites.
#' na.rm: Should NA's be removed or not? If na.rm = FALSE, then
#'   any column of x that returns an NA value will result in an NA value in the output.
#'   If na.rm = TRUE, then we remove the NA values before computing the likelihood.
#'   So if e.g. a column has 3 NA variables, then we remove these, and then we compute
#'   the likelihood for a (d-3)-dimensional Gaussian random variable.
#' @export
loglik_conditional = function(y,
                              y0,
                              a_func,
                              b_func,
                              sigma,
                              tau,
                              dist_to_s0,
                              dist_to_s0_from_mesh,
                              A,
                              n_cores = 1,
                              na.rm = TRUE) {
  lengths = c(length(y0), length(y), length(A), length(dist_to_s0), length(dist_to_s0_from_mesh))
  stopifnot(all(lengths == lengths[1]))
  stopifnot(all(sapply(A, nrow) == sapply(dist_to_s0, length)))
  stopifnot(all(sapply(A, ncol) == sapply(dist_to_s0_from_mesh, length)))
  res = parallel::mclapply(
    X = seq_along(A),
    mc.cores = n_cores,
    FUN = function(i) {
      b = b_func(y0[[i]], dist_to_s0_from_mesh[[i]])
      const_b = all(apply(b, 1, function(x) length(unique(x)) == 1))
      args = list(
        x = y[[i]] - a_func(y0[[i]], dist_to_s0[[i]]),
        A = A[[i]],
        sigma0 = sigma,
        nugget = 1 / tau,
        logd = TRUE,
        na_rm = na.rm)
      if (const_b) {
        args$b = b[, 1]
        func = dconditional_no_beta
      } else {
        args$B = B
        func = dconditional
      }

      as.numeric(do.call(func, args))
    })
  unlist(res)
}
