
#' @export
ll_conditional = function(y,
                          y0,
                          a_func,
                          b_func,
                          Σ_func,
                          τ,
                          dist,
                          A,
                          dist_from_mesh = NULL,
                          A_constr = NULL,
                          num_cores = 1,
                          no_beta = NULL,
                          na.rm = TRUE) {
  lengths = c(length(y0), length(y), length(A), length(dist))
  stopifnot(all(lengths == lengths[1]))
  if (!is.null(dist_from_mesh)) stopifnot(length(dist_from_mesh) == length(y))
  if (!is.null(A_constr)) stopifnot(length(A_constr) == length(y))
  if (is.null(no_beta)) no_beta = (length(unique(as.numeric(b_func(1:20, 3)))) == 1)
  res = parallel::mclapply(
    X = seq_along(A),
    mc.cores = num_cores,
    FUN = function(i) {
      if (no_beta) {
        b = b_func(1, dist_from_mesh[[i]])
        if (!is.null(A_constr)) {
          Σ = Σ_func(A[[i]], A_constr[[i]], diag(as.numeric(b), length(b)))
          res = dconditional_arma(
            x = y[[i]] - a_func(y0[[i]], dist[[i]]),
            sigma = Σ,
            nugget = 1 / τ,
            logd = TRUE,
            na_rm = na.rm)
          res = as.numeric(res)
        } else {
          res = dconditional_arma2(
            x = y[[i]] - a_func(y0[[i]], dist[[i]]),
            A = A[[i]],
            B = Matrix::Diagonal(length(b), as.numeric(b)),
            sigma0 = Σ_func(),
            nugget = 1 / τ,
            logd = TRUE,
            na_rm = na.rm)
          res = as.numeric(res)
        }
       } else {
        res = sapply(
          X = seq_along(y0[[i]]),
          FUN = function(j) {
            b = b_func(y0[[i]][j], dist_from_mesh[[i]])
            Σ = Σ_func(A[[i]], A_constr[[i]], diag(as.numeric(b), length(b)))
            dconditional(
              y = y[[i]][, j, drop = FALSE],
              a = a_func(y0[[i]][j], dist[[i]]),
              Σ = Σ,
              τ = τ,
              log = TRUE,
              na.rm = na.rm)
          })
      }
      res
    })
  unlist(res)
}


#' @export
dconditional = function(y, a, Σ, τ, log = TRUE, num_cores = 1, na.rm = TRUE) {
  d = length(y)
  stopifnot(length(a) == d)
  stopifnot(all(dim(Σ) == d))
  nugget = diag(1 / τ, d)
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
rconditional = function(y0,
                        a_func,
                        b_func,
                        Q,
                        τ,
                        dist,
                        dist_from_mesh,
                        A,
                        A_constr = NULL) {
  if (is.list(A)) {
    func = rconditional_multiple_A
  } else {
    func = rconditional_one_A
  }
  args = list(
    y0 = y0,
    a_func = a_func,
    b_func = b_func,
    Q = Q,
    τ = τ,
    dist = dist,
    dist_from_mesh = dist_from_mesh,
    A = A,
    A_constr = A_constr)
  do.call(func, args)
}

rconditional_one_A = function(y0,
                              a_func,
                              b_func,
                              Q,
                              τ,
                              dist,
                              dist_from_mesh,
                              A,
                              A_constr = NULL) {
  z = rnorm_spde(n = length(y0), Q = Q)
  rconditional_given_z(
    z = z,
    a = a_func(y0, dist),
    b = b_func(y0, dist_from_mesh),
    τ = τ,
    A = A,
    Q = Q,
    A_constr = A_constr)
}

rconditional_multiple_A = function(y0,
                                   a_func,
                                   b_func,
                                   Q,
                                   τ,
                                   dist,
                                   dist_from_mesh,
                                   A,
                                   A_constr = NULL) {
  n = sapply(y0, length)
  n_A = length(A)
  lengths = c(length(y0), length(A), length(dist), length(dist_from_mesh))
  stopifnot(all(lengths == n_A))
  stopifnot(is.null(A_constr) || length(A_constr) == n_A)

  res = vector("list", n_A)
  count = 0
  z = rnorm_spde(n = sum(n), Q = Q)
  for (i in 1:n_A) {
    if (n[i] == 0) next
    index = count + seq_len(n[i])
    count = count + n[i]
    res[[i]] = rconditional_given_z(
      z = z[, index, drop = FALSE],
      a = a_func(y0[[i]], dist[[i]]),
      b = b_func(y0[[i]], dist_from_mesh[[i]]),
      A = A[[i]],
      A_constr = A_constr[[i]],
      Q = Q,
      τ = τ)
  }
  res
}

rconditional_given_z = function(z, a, b, τ, A, Q = NULL, A_constr = NULL) {
  stopifnot(all(dim(z) == dim(b)))
  stopifnot(nrow(a) == nrow(A))
  stopifnot(ncol(a) == ncol(z))
  z = b * z
  if (!is.null(A_constr)) {
    stopifnot(!is.null(Q))
    z = constrain_spde_to_zero(z, Q, A_constr)
  }
  res = a + as.matrix(A %*% z) + rnorm(length(a), sd = τ^-.5)
  res
}

#' @export
rnorm_spde = function(n, Q, ...) {
  # This allows us to get reproduceable sampling if we set the seed outside the function scope
  res = suppressWarnings(INLA::inla.qsample(n, Q, seed = round(runif(1, 0, 1e6)), ...))
  colnames(res) = NULL
  rownames(res) = NULL
  res
}

constrain_spde_to_zero = function(x, Q, A) {
  x - Matrix::solve(Q, t(A)) %*% Matrix::solve(A %*% Matrix::solve(Q, t(A)), A %*% x)
}

#' @export
get_Σ_func = function(θ, spde) {
  ρ = exp(θ[1]); σ = exp(θ[2])
  res = list()
  res$Q = INLA::inla.spde2.precision(spde, θ)
  Σ = as.matrix(Matrix::solve(res$Q))
  res$value = function(A = NULL, A_constr = NULL, B = NULL) {
    if (!is.null(B)) {
      stopifnot(all(dim(B) == dim(Σ)))
      Σ = B %*% Σ %*% B
    }
    if (!is.null(A_constr)) Σ = Σ_spde_constrained_func(Σ, A_constr)
    if (is.null(A)) {
      res = Σ
    } else {
      res = as.matrix(A %*% Σ %*% Matrix::t(A))
    }
    res
  }
  res
}

Σ_spde_constrained_func = function(Σ, A_constr) {
  stopifnot(all(dim(Σ) == ncol(A_constr)))
  Σ - Σ %*% t(A_constr) %*% solve(A_constr %*% Σ %*% t(A_constr), A_constr) %*% Σ
}
