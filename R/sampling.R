
#' This is a wrapper for the inla.qsample() function for sampling from a Gaussian
#' random field with sparse precision matrix Q. The wrapper allows us to get
#' reproduceable results if the seed outside the function scope is known.
#' It also suppresses warnings and removes some unnecessary column and row names
#' @export
rnorm_spde = function(n, Q, ...) {
  res = suppressWarnings(INLA::inla.qsample(n, Q, seed = round(runif(1, 0, 1e6)), ...))
  colnames(res) = NULL
  rownames(res) = NULL
  res
}

#' Sample from the global conditional extremes distribution of [x | max(x) > t],
#' using the simulation algorithm of Keef et al. (2013).
#' The arguments of the function are:
#' n: Number of samples.
#' a_func: Function a(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of a(y0[j], dist_to_s0[i]).
#' b_func: Function b(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of a(y0[j], dist_to_s0[i]).
#'   b_func is different from a_func in that a_func requires the distance to s0
#'   for all locations with observations, while b_func requires the distance to s0
#'   for all locations on the triangular mesh used for defining the SPDE.
#' Q: Precision matrix of the weights used for building the SPDE approximation.
#' tau: Precision of the nugget effect.
#' threshold: The threshold t used for defining the conditional extremes distribution.
#' dist_to_s0: List containing the distances to s0 for all observation locations,
#'   used for computing a. The list contains one vector of distances for each
#'   possible conditioning site in the domain of interest. Each of the distance
#'   vectors must have the same length.
#' dist_to_s0_from_mesh: List containing the distances to s0 for all locations in
#'   the mesh, used for computing b. The list contains one vector of distances for each
#'   possible conditioning site in the domain of interest. Each of the distance
#'   vectors must have the same length.
#' A: List of projection matrices used for building the SPDE approximation. The list
#'   contains one projection matrix for each possible conditioning site in the
#'   domain of interest. Each projection matrix must contain the same number of
#'   rows (equal to the lengths of dist_to_s0) and columns (equal to the lengths of
#'   dist_to_s0_from_mesh).
#' verbose: bool that states if the function should display progress or not.
#' @export
keef_sampling = function(n,
                         a_func,
                         b_func,
                         Q,
                         tau,
                         threshold,
                         dist_to_s0,
                         dist_to_s0_from_mesh,
                         A,
                         verbose = FALSE) {
  n_s0 = length(A)
  # Check that the input has the correct length
  stopifnot(length(dist_to_s0) == n_s0 && length(dist_to_s0_from_mesh) == n_s0)
  # Check that the A matrices have the same number of rows as the lengths of dist_to_s0
  stopifnot(all(sapply(dist_to_s0, length) == sapply(A, nrow)))
  stopifnot(length(unique(sapply(dist_to_s0, length))) == 1)
  # Check that all the dist_to_s0_from_mesh distances have length equal to the
  # number of nodes in the triangular mesh
  stopifnot(all(sapply(dist_to_s0_from_mesh, length) == nrow(Q)))

  # Preallocate the result matrix
  res = matrix(NA_real_, nrow = n_s0, ncol = n)

  # Choose which locations the different simulations should use as conditioning sites
  s0_index = sort(sample.int(n_s0, n, replace = TRUE))

  # Index over all samples where y(s0) is not the maximum observed value
  y0_not_max = seq_len(n)

  if (verbose) message("Start simulation")

  # Sample from the conditional model until y(s0) is the maximum value for all
  # the replications
  while (any(y0_not_max)) {

    if (verbose) {
      message("Samples remaining: ", length(y0_not_max), " / ", n)
    }

    # For each possible conditioning site, sample y(s0) for the replications
    # where y(s0) is not the maximum value
    n_y0_per_s0 = sapply(1:n_s0, function(i) sum(s0_index[y0_not_max] == i))
    y0 = lapply(n_y0, function(n) rexp(n) + threshold)

    # Sample from the conditional extremes model given the values of y(s0)
    samples = rconditional(
      y0 = y0,
      a_func = a_func,
      b_func = b_func,
      Q = Q,
      tau = tau,
      dist = dist_to_s0,
      dist_from_mesh = dist_to_s0_from_mesh,
      A = A)

    # Add the value of y0 to the correct location in the samples
    for (i in seq_along(samples[-1])) {
      samples[[i]] = cbind(
        samples[[i]][, seq_len(i - 1)],
        y0[[i]],
        samples[[i]][, i:ncol(samples[[i]])])
    }
    samples[[n_s0]] = cbind(samples[[n_s0]], y0[[n_s0]])

    # Update the rows of `res` where y(s0) is still not the maximum
    res[, y0_not_max] = do.call(cbind, samples)

    # Find the indices where y(s0) is still not the maximum value
    y0_not_max = which(sapply(
      X = seq_len(n),
      FUN = function(i) any(res[-s0_index[i], i] >= res[s0_index[i], i])))
  }

  res
}

#' Function for simulating from the conditional extremes model given the values of y(s0).
#' The arguments of the function are:
#' y0: A list of vectors containing the values of y(s0) for each of the possible
#'   conditioning sites s0. The vectors can have length 0.
#' a_func: Function a(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of a(y0[j], dist_to_s0[i]).
#' b_func: Function b(y0, dist_to_s0) that takes in two vectors and returns a
#'   matrix where element (i, j) is the value of a(y0[j], dist_to_s0[i]).
#'   b_func is different from a_func in that a_func requires the distance to s0
#'   for all locations with observations, while b_func requires the distance to s0
#'   for all locations on the triangular mesh used for defining the SPDE.
#' Q: Precision matrix of the weights used for building the SPDE approximation.
#' tau: Precision of the nugget effect.
#' threshold: The threshold t used for defining the conditional extremes distribution.
#' dist_to_s0: List containing the distances to s0 for all observation locations,
#'   used for computing a. The list contains one vector of distances for each
#'   possible conditioning site in the domain of interest. Each of the distance
#'   vectors may be of different size.
#' dist_to_s0_from_mesh: List containing the distances to s0 for all locations in
#'   the mesh, used for computing b. The list contains one vector of distances for each
#'   possible conditioning site in the domain of interest. Each of the distance
#'   vectors must be of the same size.
#' A: List of projection matrices used for building the SPDE approximation. The list
#'   contains one projection matrix for each possible conditioning site in the
#'   domain of interest. Each projection matrix may contain a different number of
#'   rows, but all must have the same number of columns, which is equal to the number
#'   nodes in the triangular mesh.
#' @export
rconditional = function(y0,
                        a_func,
                        b_func,
                        Q,
                        tau,
                        dist_to_s0,
                        dist_to_s0_from_mesh,
                        A) {
  n = sapply(y0, length)
  n_s0 = length(y0)
  # Check that the input has correct lengths
  lengths = c(length(A), length(dist_to_s0), length(dist_to_s0_from_mesh))
  stopifnot(all(lengths == n_s0))
  # Check that the A matrices have the same number of rows as the lengths of dist_to_s0
  stopifnot(all(sapply(dist_to_s0, length) == sapply(A, nrow)))
  # Check that all the dist_to_s0_from_mesh dist_to_s0ances have length equal to the
  # number of nodes in the triangular mesh
  stopifnot(all(sapply(dist_to_s0_from_mesh, length) == nrow(Q)))

  # Preallocate the result
  res = vector("list", n_s0)
  count = 0
  z = rnorm_spde(n = sum(n), Q = Q)
  for (i in 1:n_s0) {
    if (n[i] == 0) next
    index = count + seq_len(n[i])
    count = count + n[i]

    # Compute a and b, and check that their dimensions are correct
    a = a_func(y0[[i]], dist_to_s0[[i]])
    b = b_func(y0[[i]], dist_to_s0_from_mesh[[i]])
    stopifnot(all(dim(z[, index]) == dim(b)))
    stopifnot(all(c(nrow(A[[i]]), length(index)) == dim(a)))

    # Sample from the conditional extremes model, given the values of a and b
    z[, index] = b * z[, index]
    res[[i]] = a + as.matrix(A[[i]] %*% z[, index]) + rnorm(length(a), sd = tau^-.5)
  }

  res
}
