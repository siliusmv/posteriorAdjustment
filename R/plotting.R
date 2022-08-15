
#' Plot a mesh as a ggplot object. This is necessary if we want to plot
#' the mesh on top of another ggplot object, like a map of Norway
#' @export
plot_mesh = function(mesh, coords = NULL) {
  plot = ggplot() +
    inlabru::gg(mesh, int.color = "black") +
    theme_light()
  if (!is.null(coords)) {
    if (is(coords, "sf") || is(coords, "sfc")) {
      mesh_bbox = c(min(mesh$loc[, 1]), min(mesh$loc[, 2]), max(mesh$loc[, 1]), max(mesh$loc[, 2]))
      plot = plot +
        geom_sf(data = coords) +
        add_norway_map(sf::st_crs(coords), mesh_bbox, col = "blue")
    } else {
      coords = as.data.frame(coords)
      plot = plot +
        geom_point(data = coords, aes(x = x, y = y))
    }
  }
  plot
}

#' This function takes an sf data.frame and plots it on a map
#' @export
plot_on_map = function(data, response_name = NULL, breaks = NULL, use_tex = FALSE, ...) {
  if (!is.null(breaks)) {
    data[[response_name]] = get_binned_data(data[[response_name]], breaks, use_tex = use_tex)
  }
  if (is.null(response_name)) {
    gg = ggplot2::ggplot(data) + ggplot2::geom_sf()
  } else {
    gg = ggplot2::ggplot(data) +
      ggplot2::geom_sf(ggplot2::aes_string(col = response_name)) +
      ggplot2::scale_color_viridis_c()
  }
  gg = style_map_plot(gg, data, use_tex, ...)
  gg
}

#' Take a ggplot object and add a map of Norway in the background. Also, do some styling
#' of the ggplot object.
#' @export
style_map_plot = function(plot, data, use_tex = FALSE, axis_text = TRUE, ...) {
  plot = plot +
    add_norway_map(sf::st_crs(data), sf::st_bbox(data), ...) +
    ggplot2::theme_light() +
    ggplot2::labs(x = "", y = "")
  if (!axis_text) {
    plot = plot + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  }
  if (use_tex) plot = latex_friendly_map_plot(plot)
  plot
}

latex_friendly_map_plot = function(x) {
  info = ggplot2::ggplot_build(x)$layout$panel_params[[1]]$graticule
  east_ticks = info$degree[info$type == "E"]
  north_ticks = info$degree[info$type == "N"]
  x +
    ggplot2::scale_x_continuous(breaks = east_ticks, labels = paste0(east_ticks, "$^\\circ$E"), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = north_ticks, labels = paste0(north_ticks, "$^\\circ$N"), expand = c(0, 0))
}

get_binned_data = function(x, breaks = NULL, digits = 1, use_tex = FALSE) {
  if (all(is.null(breaks))) breaks = seq(min(x), max(x), length.out = 5)
  if (!is.infinite(breaks[1])) breaks = c(-Inf, breaks)
  if (!is.infinite(breaks[length(breaks)])) breaks = c(breaks, Inf)
  breaks = round(breaks, digits)
  if (use_tex) {
    labels = c(
      paste("$<$", breaks[2]),
      paste(breaks[2:(length(breaks) - 2)], "-", breaks[3:(length(breaks) - 1)]),
      paste("$>$", breaks[length(breaks) - 1]))
  } else {
    labels = c(
      paste("<", breaks[2]),
      paste(breaks[2:(length(breaks) - 2)], "-", breaks[3:(length(breaks) - 1)]),
      paste(">", breaks[length(breaks) - 1]))
  }
  x = cut(x, breaks, labels = labels)
  x
}

#' @export
add_norway_map = function(crs = NULL, bbox = NULL, ...) {
  map = rnaturalearth::ne_countries(scale = 50, country = "Norway", returnclass = "sf")
  c(ggplot2::layer_sf(
    geom = ggplot2::GeomSf, data = map, mapping = ggplot2::aes(),
    stat = "sf", position = "identity", show.legend = NA,
    inherit.aes = TRUE,
    params = list(na.rm = FALSE, fill = NA, ...)),
    ggplot2::coord_sf(default = TRUE, crs = crs,
                      xlim = bbox[c(1, 3)], ylim = bbox[c(2, 4)]))
}

# #' Load data from the radar-data.rds object, and plot rasters from specified
# #' times or indices.
# #' @export
# plot_radar_data = function(index = NULL, time = NULL, use_unique_legends = TRUE, s0 = NULL) {
#   radar = readRDS(file.path(data_dir(), "radar-data.rds"))
#   data = radar$get_data()
#   on.exit(close(data))
# 
#   if (is.null(index)) {
#     if (is.null(time)) stop("Either index or time must be provided as input")
#     index = which(radar$times %in% time)
#   }
# 
#   df = sf::st_coordinates(radar$coords) %>%
#     cbind(t(data[index, ])) %>%
#     as.data.frame()
# 
#   if (use_unique_legends) {
#     plots = lapply(
#       X = seq_along(index),
#       FUN = function(i) {
#         plot = ggplot(df[, c(1, 2, 2 + i)]) +
#           geom_raster(aes_string(x = "X", y = "Y", fill = paste0("V", i + 2))) +
#           scale_fill_viridis_c() +
#           scale_x_continuous(expand = c(0, 0)) +
#           scale_y_continuous(expand = c(0, 0)) +
#           theme_light() +
#           theme(panel.grid = element_blank()) +
#           labs(fill = "mm/h")
#         if (!is.null(s0)) {
#           plot = plot + geom_point(x = s0[1, 1], y = s0[1, 2], col = "red")
#         }
#         plot
#       }) %>%
#       patchwork::wrap_plots()
#   } else {
#     plots = df %>%
#       tidyr::pivot_longer(-c(X, Y)) %>%
#       ggplot() +
#       geom_raster(aes(x = X, y = Y, fill = value)) +
#       facet_wrap(~name) +
#       scale_fill_viridis_c() +
#       scale_x_continuous(expand = c(0, 0)) +
#       scale_y_continuous(expand = c(0, 0)) +
#       theme_light() +
#       theme(panel.grid = element_blank()) +
#       labs(fill = "mm/h")
#     if (!is.null(s0)) {
#       plots = plots + geom_point(x = s0[1, 1], y = s0[1, 2], col = "red")
#     }
#   }
# 
#   plots
# }


