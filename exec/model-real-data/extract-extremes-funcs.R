
locate_nonregular_grids = function(loc, s0_index, n = 1, r = Inf) {
  stopifnot(1 <= min(s0_index) && max(s0_index) <= nrow(loc))
  res = list(
    s0_index = list(),
    s0 = list(),
    loc_index = list(),
    dist_to_s0 = list())
  n_s0 = length(s0_index)
  for (i in seq_len(n_s0)) {
    res$loc_index[[i]] = locate_every_nth_cell_in_circles(
      coords = loc,
      center = loc[s0_index[i], , drop = FALSE],
      n = n,
      r = r,
      index_only = TRUE)
    res$loc_index[[i]] = res$loc_index[[i]][res$loc_index[[i]] != s0_index[i]]
    res$s0_index[[i]] = s0_index[i]
    res$s0[[i]] = loc[s0_index[i], , drop = FALSE]
    res$dist_to_s0[[i]] = as.numeric(dist_euclid(
      x = loc[s0_index[i], , drop = FALSE],
      y = loc[res$loc_index[[i]], , drop = FALSE]))
  }
  res$s0_index = unlist(res$s0_index)
  res
}

sliding_window_marginal = function(data, coords, center, radius, p) {
  window_index = locate_every_nth_cell_in_circles(
    coords = coords,
    center = center,
    n = 1,
    r = radius,
    index_only = TRUE)
  x = as.numeric(data[, window_index])
  marginal_distribution(x, p)
}

extract_extreme_fields = function(data,
                                  coords,
                                  s0_index,
                                  threshold,
                                  n,
                                  r,
                                  rm.inf = TRUE) {
  res = locate_nonregular_grids(
    loc = coords,
    s0_index = s0_index,
    n = n,
    r = r)

  res$y0 = res$y = res$time_index = vector("list", length(s0_index))
  good_index = NULL
  for (i in seq_along(s0_index)) {
    y0 = data[, res$s0_index[i]]
    res$time_index[[i]] = which(y0 > threshold)
    if (length(res$time_index[[i]]) == 0) next
    res$y0[[i]] = y0[res$time_index[[i]]]
    good_index = c(good_index, i)
    res$y[[i]] = data[res$time_index[[i]], res$loc_index[[i]], drop = FALSE]
    res$y[[i]] = t(res$y[[i]]) # we want each replicate of y to have one column
    if (rm.inf && any(is.infinite(res$y[[i]]))) {
      res$y[[i]][is.infinite(res$y[[i]])] = NA
    }
  }

  res$s0 = res$s0[good_index]
  res$s0_index = res$s0_index[good_index]
  res$loc_index = res$loc_index[good_index]
  res$dist_to_s0 = res$dist_to_s0[good_index]
  res$y = res$y[good_index]
  res$y0 = res$y0[good_index]
  res$time_index = res$time_index[good_index]

  res$n = sapply(res$y0, length)

  res
}
