
#' Compute the aggregated empirical cumulative distribution function for all
#' data inside a circle with a given radius. Input variables are:
#' data: An (n x d)-dimensional matrix of observations.
#' coords: A (d x 2)-dimensional matrix with the coordinates of the data.
#' center: A 2-dimensional vector containing the coordinates of the center we use
#'   for computing the aggregated ECDF.
#' radius: The radius of the circle used for computing the aggregated ECDF.
#' @export
aggregated_ecdf = function(data, coords, center, radius) {
  stopifnot(ncol(data) == nrow(coords))
  window_index = extract_thinned_out_circles(
    coords = coords,
    center = center,
    n = 1,
    r = radius,
    index_only = TRUE)
  x = as.numeric(data[, window_index])
  marginal_distribution(x)
}


#' Given a matrix of coordinates on a regular grid, extract a subset of the coordinates
#' that consists of circles centred around some center, with varying degrees of
#' densities. The input variables are:
#' coords: A matrix with two columns and multiple rows, containing (x, y) coordinates
#'   in a regular grid and with a distance-preserving projection.
#' center: One (x, y) tuple of coordinates, explaining the center of the circles we
#'   wish to extract.
#' n: A vector of densities for the different circles we wish to extract. If the ith
#'   value of n is k, then we only keep every kth of the original coords, both in
#'   the x-direction and in the y-direction.
#' r: A vector of radii for the different circles we wish to extract. If e.g.
#'   n = (1, 3) and r = (2, 5), then we will extract every coordinate that is closer
#'   to the center than a radius of 2. Then, we will only keep every 3rd coordinate,
#'   both in x-direction and y-direction, for all coordinates that have a distance
#'   between 2 and 5 from the center. Finally, all coords that are further away from
#'   the center than 5 will be dropped.
#' index_only: A bool describing if the function should return the extracted
#'   coordinates, or only their indices inside the coords object.
#' @export
extract_thinned_out_circles = function(coords,
                                       center,
                                       n,
                                       r = Inf,
                                       index_only = FALSE) {
  stopifnot(ncol(coords) == 2)
  stopifnot(length(n) == length(r))
  stopifnot(is.matrix(coords))

  # Locate all unique x-values and y-values
  x_vals = sort(unique(coords[, 1]))
  y_vals = sort(unique(coords[, 2]))

  # Compute distances and find the coord that is closest to the center
  dist_to_center = as.numeric(dist_euclid(center, coords))
  closest_to_center = coords[which.min(dist_to_center), ]
  x_center_index = which(x_vals == closest_to_center[1])
  y_center_index = which(y_vals == closest_to_center[2])

  # This will be the list of all indices we wish to extract from coords
  chosen_indices = NULL

  # Locate all the indices we wish to extract
  r = c(0, r)
  for (i in seq_along(n)) {
    # Indicies for all coords that have the correct distance from center
    dist_index = which(r[i] <= dist_to_center & dist_to_center <= r[i + 1])

    # Indices for all coords on a regular grid of resolution n[i] x n[i]
    thinned_x_index = c(seq(x_center_index, 1, by = -n[i]),
                        seq(x_center_index, length(x_vals), by = n[i]))
    thinned_y_index = c(seq(y_center_index, 1, by = -n[i]),
                        seq(y_center_index, length(y_vals), by = n[i]))
    thinned_x_vals = x_vals[thinned_x_index]
    thinned_y_vals = y_vals[thinned_y_index]
    thinned_index = which(coords[, 1] %in% thinned_x_vals & coords[, 2] %in% thinned_y_vals)

    # The indices of interest are the intersection of the two index sets above
    chosen_indices = c(chosen_indices, intersect(dist_index, thinned_index))
  }
  # Some times, duplicates may appear. Remove these
  chosen_indices = unique(chosen_indices)

  if (index_only) {
    res = chosen_indices
  } else {
    res = coords[chosen_indices, ]
  }
  res
}

#' Search through data after threshold exceedances at a set of conditioning sites.
#' Extract relevant information and observations for each of the threshold exceedances.
#'
#'
#' The input variables are:
#' data: An (n x d)-dimensional matrix of observations.
#' coords: A (d x 2)-dimensional matrix with the coordinates of the data.
#' s0_index: a vector of indices describing which of the d coordinates are chosen
#'   as conditioning sites.
#' threshold: The threshold t for the conditional extremes model.
#' n, r: Used for only extracting a subset of the observations in a circle around
#'   the conditioning site. See the docs for extract_thinned_out_circles for more info.
#'
#'
#' The output of the function is a list where each element is a list or vector
#' of length equal, containing these variables:
#' s0_index: This is a subset of all the s0_indices given as input such that each
#'   conditioning site has at least one threshold exceedance.
#' s0: A list where element nr. i is a 2-dimensional vector with the coordinates of
#'   conditioning site nr. i.
#' obs_index: A list where element nr. i is a vector of indices for all the
#'   locations in coords we use for performing inference using conditioning site nr. i.
#' dist_to_s0: A list where element nr. i is a vector with all the distances between
#'   conditioning site nr. i and the locations used for inference with conditioning
#'   site nr. i.
#' y0: A list where element nr. i is a vector of all threshold exceedances for
#'   conditioning site nr. i.
#' y: A list where element nr. i is a matrix consisting of one column per threshold
#'   exceedance in y0[[i]]. Each column contains observations from the locations
#'   described in obs_index[[i]].
#' time_index: A list where element nr. i contains the times of all threshold
#'   exceedances from conditioning site nr. i.
#' n: A vector where element nr. i is the number of threshold exceedances at
#'   conditioning site nr. i.
#' @export
extract_extreme_fields = function(data,
                                  coords,
                                  s0_index,
                                  threshold,
                                  n,
                                  r) {
  # s0_indices must be between 1 and nrow(coords)
  stopifnot(1 <= min(s0_index) && max(s0_index) <= nrow(coords))
  n_s0 = length(s0_index)

  # Initialise the output
  res = rep(list(vector("list", n_s0)), 6)
  names(res) = c("s0", "obs_index", "dist_to_s0", "y0", "y", "time_index")
  res$s0_index = s0_index

  # Look through the conditioning sites
  for (i in seq_len(n_s0)) {

    # Extract all variables at s0 nr. i, and locate threshold exceedances
    y0 = data[, res$s0_index[i]]
    res$time_index[[i]] = which(y0 > threshold)
    # If there are no threshold exceedances, skip to next conditioning site
    if (length(res$time_index[[i]]) == 0) next

    # Extract all threshold exceedances
    res$y0[[i]] = y0[res$time_index[[i]]]

    # Locate the thinned out observations locations
    res$obs_index[[i]] = extract_thinned_out_circles(
      coords = coords,
      center = coords[s0_index[i], , drop = FALSE],
      n = n,
      r = r,
      index_only = TRUE)
    res$obs_index[[i]] = res$obs_index[[i]][res$obs_index[[i]] != s0_index[i]]

    # Extract all observations of interest for the threshold exceedances
    res$y[[i]] = t(data[res$time_index[[i]], res$obs_index[[i]], drop = FALSE])

    # Fill in information about s0 and the dist_to_s0
    res$s0[[i]] = coords[s0_index[i], , drop = FALSE]
    res$dist_to_s0[[i]] = as.numeric(dist_euclid(
      x = coords[s0_index[i], , drop = FALSE],
      y = coords[res$obs_index[[i]], , drop = FALSE]))
  }

  # Only keep information about the conditioning sites with at least one threshold exceedance
  good_index = which(sapply(res$time_index, length) > 0)
  res = lapply(res, `[`, good_index)

  # Find the number of threshold exceedances for each conditioning site
  res$n = sapply(res$y0, length)

  res
}
