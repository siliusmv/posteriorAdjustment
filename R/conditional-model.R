
#' @export
ll_conditional = function(y,
                          y0,
                          a_func,
                          b_func,
                          Σ_func,
                          tau,
                          dist,
                          A,
                          dist_from_mesh = NULL,
                          num_cores = 1,
                          no_beta = NULL,
                          na.rm = TRUE) {
  lengths = c(length(y0), length(y), length(A), length(dist))
  stopifnot(all(lengths == lengths[1]))
  if (!is.null(dist_from_mesh)) stopifnot(length(dist_from_mesh) == length(y))
  if (is.null(no_beta)) no_beta = (length(unique(as.numeric(b_func(1:20, 3)))) == 1)
  res = parallel::mclapply(
    X = seq_along(A),
    mc.cores = num_cores,
    FUN = function(i) {
      if (no_beta) {
        b = b_func(1, dist_from_mesh[[i]])
        res = dconditional_arma(
          x = y[[i]] - a_func(y0[[i]], dist[[i]]),
          A = A[[i]],
          B = Matrix::Diagonal(length(b), as.numeric(b)),
          sigma0 = Σ_func(),
          nugget = 1 / tau,
          logd = TRUE,
          na_rm = na.rm)
        res = as.numeric(res)
      } else {
        res = sapply(
          X = seq_along(y0[[i]]),
          FUN = function(j) {
            b = b_func(y0[[i]][j], dist_from_mesh[[i]])
            Σ = Σ_func(A[[i]], diag(as.numeric(b), length(b)))
            dconditional(
              y = y[[i]][, j, drop = FALSE],
              a = a_func(y0[[i]][j], dist[[i]]),
              Σ = Σ,
              tau = tau,
              log = TRUE,
              na.rm = na.rm)
          })
      }
      res
    })
  unlist(res)
}


#' @export
dconditional = function(y, a, Σ, tau, log = TRUE, num_cores = 1, na.rm = TRUE) {
  d = length(y)
  stopifnot(length(a) == d)
  stopifnot(all(dim(Σ) == d))
  nugget = diag(1 / tau, d)
  y = y - a # Standardise y
  if (na.rm && any(is.na(y))) {
    good_index = which(!is.na(y))
    if (length(good_index) == 0) return(NA_real_)
  } else {
    good_index = seq_len(d)
  }
  sigma = Σ + nugget
  res = tryCatch({
    mvnfast::dmvn(
      X = y[good_index],
      mu = rep(0, length(good_index)),
      sigma = sigma[good_index, good_index],
      log = TRUE)
  },
  error = function(e) NULL)
  if (is.null(res)) {
    res = mvtnorm::dmvnorm(
      x = y[good_index],
      sigma = sigma[good_index, good_index],
      log = log)
  }
  res
}


#' @export
get_Σ_func = function(theta, spde) {
  ρ = exp(theta[1]); σ = exp(theta[2])
  res = list()
  res$Q = INLA::inla.spde2.precision(spde, theta)
  Σ = as.matrix(Matrix::solve(res$Q))
  res$value = function(A = NULL, B = NULL) {
    if (!is.null(B)) {
      stopifnot(all(dim(B) == dim(Σ)))
      Σ = B %*% Σ %*% B
    }
    if (is.null(A)) {
      res = Σ
    } else {
      res = as.matrix(A %*% Σ %*% Matrix::t(A))
    }
    res
  }
  res
}
