
#' @export
keef_sampling = function(n,
                         a_func,
                         b_func,
                         Q,
                         τ,
                         threshold,
                         dist_to_s0,
                         dist_to_s0_from_mesh,
                         A,
                         A_constr = NULL,
                         verbose = FALSE) {
  n_loc = length(A)
  stopifnot(length(dist_to_s0) == n_loc)
  stopifnot(is.null(A_constr) || length(A_constr) == n_loc)
  s0_index = sort(sample.int(n_loc, n, replace = TRUE))
  res = matrix(0, nrow = n_loc, ncol = n)
  indx = which(sapply(1:n, \(i) any(res[-s0_index[i], i] >= res[s0_index[i], i])))
  while (TRUE) {
    n_y0 = sapply(1:n_loc, \(i) sum(s0_index[indx] == i))
    y0 = lapply(n_y0, \(n) rexp(n) + threshold)
    samples = rconditional(
      y0 = y0,
      a_func = a_func,
      b_func = b_func,
      Q = Q,
      τ = τ,
      dist = dist_to_s0,
      dist_from_mesh = dist_to_s0_from_mesh,
      A = A,
      A_constr = A_constr)
    samples = add_y0_to_rconditional(samples, y0, seq_along(y0))
    res[, indx] = do.call(cbind, samples)
    indx = which(sapply(1:n, \(i) any(res[-s0_index[i], i] >= res[s0_index[i], i])))
    if (length(indx) == 0) break
    if (verbose) message("Progress: ", length(indx), " / ", n, " remaining")
  }
  res
}

#' @export
wadsworth_sampling = function(n1,
                              n2,
                              a_func,
                              b_func,
                              Q,
                              τ,
                              threshold,
                              dist_to_s0,
                              dist_to_s0_from_mesh,
                              A,
                              A_constr = NULL,
                              replace = FALSE) {
  stopifnot(n2 > n1)
  n_loc = length(A)
  stopifnot(length(dist_to_s0) == n_loc)
  stopifnot(is.null(A_constr) || length(A_constr) == n_loc)
  s0_index = sort(sample.int(n_loc, n2, replace = TRUE))
  n_y0 = sapply(1:n_loc, \(i) sum(s0_index == i))
  y0 = lapply(n_y0, \(n) rexp(n) + threshold)
  res = rconditional(
    y0 = y0,
    a_func = a_func,
    b_func = b_func,
    Q = Q,
    τ = τ,
    dist = dist_to_s0,
    dist_from_mesh = dist_to_s0_from_mesh,
    A = A,
    A_constr = A_constr)
  res = add_y0_to_rconditional(res, y0, seq_along(y0))
  res = do.call(cbind, res)
  importance_weights = apply(res, 2, function(x) 1 / sum(x > threshold))
  sampling_index = sample.int(n2, n1, prob = importance_weights, replace = replace)
  res[, sampling_index]
}

add_y0_to_rconditional = function(y, y0, s0_index = rep(1, length(y0))) {
  if (is.list(y)) {
    y = add_y0_to_rconditional_multiple_y(y, y0, s0_index)
  } else {
    y = add_y0_to_rconditional_one_y(y, y0, s0_index)
  }
  y
}

add_y0_to_rconditional_multiple_y = function(y, y0, s0_index = rep(1, length(y0))) {
  lengths = c(length(y), length(y0), length(s0_index))
  stopifnot(all(lengths == lengths[1]))
  for (i in seq_along(y)) {
    if (is.null(y[[i]])) next
    y[[i]] = add_y0_to_rconditional_one_y(y[[i]], y0[[i]], rep(s0_index[i], length(y0[[i]])))
  }
  y
}

add_y0_to_rconditional_one_y = function(y, y0, s0_index = rep(1, length(y0))) {
  stopifnot(length(y0) == length(s0_index))
  stopifnot(length(y0) == ncol(y))
  n = nrow(y)
  y = rbind(NA, y)
  for (j in 1:ncol(y)) {
    y[, j] = c(y[seq_len(s0_index[j])[-1], j], y0[j], y[s0_index[j] + seq_len(n + 1 - s0_index[j]), j])
  }
  y
}

