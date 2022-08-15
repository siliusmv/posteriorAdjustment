
#' @export
locate_every_nth_cell_in_circles = function(coords,
                                            center,
                                            n,
                                            r = Inf,
                                            index_only = FALSE) {
  stopifnot(length(n) == length(r))
  r = c(0, r)
  if (is(coords, "sf") || is(coords, "sfc")) {
    dist_to_center = as.numeric(sf::st_distance(center, coords))
    coord_mat = sf::st_coordinates(coords)
    center = sf::st_coordinates(center)
  } else {
    dist_to_center = as.numeric(dist_euclid(center, coords))
    coord_mat = coords
  }
  x_vals = sort(unique(coord_mat[, 1]))
  y_vals = sort(unique(coord_mat[, 2]))

  closest_to_center = coord_mat[which.min(dist_to_center), ]
  x_center_index = which(x_vals == closest_to_center[1])
  y_center_index = which(y_vals == closest_to_center[2])

  chosen_indices = NULL

  for (i in seq_along(n)) {
    dist_index = which(r[i] <= dist_to_center & dist_to_center <= r[i + 1])
    thinned_x_index = c(seq(x_center_index, 1, by = -n[i]),
                        seq(x_center_index, length(x_vals), by = n[i]))
    thinned_y_index = c(seq(y_center_index, 1, by = -n[i]),
                        seq(y_center_index, length(y_vals), by = n[i]))
    thinned_x_vals = x_vals[thinned_x_index]
    thinned_y_vals = y_vals[thinned_y_index]
    thinned_index = which(coord_mat[, 1] %in% thinned_x_vals & coord_mat[, 2] %in% thinned_y_vals)
    chosen_indices = c(chosen_indices, intersect(dist_index, thinned_index))
  }
  chosen_indices = unique(chosen_indices)

  if (index_only) {
    res = chosen_indices
  } else {
    res = coords[chosen_indices, ]
  }
  res
}
