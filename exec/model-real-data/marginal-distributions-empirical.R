devtools::load_all()
library(pbapply)
library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)
library(parallel)

num_cores = 15

# ==============================================================================
# Load the radar data
# ==============================================================================
radar = readRDS(file.path(downloads_dir(), "radar-trondheim.rds"))
n_loc = nrow(radar$coords)
n_time = length(radar$times)
coords = st_coordinates(radar$coords)

# ==============================================================================
# Examine how to choose radius
# ==============================================================================
s0_index = locate_every_nth_cell_in_circles(coords, coords[1, ], n = 2, r = Inf, index_only = TRUE)
r = c(1, 2, 3, 5, 7, 10)
probs = c(.8, .9, .95, .99, .995, .999, .9995)

cl = parallel::makeForkCluster(num_cores)
df = pbapply::pblapply(
  X = seq_along(s0_index),
  cl = cl,
  FUN = function(i) {
    res = list()
    for (j in seq_along(r)) {
      ii = locate_every_nth_cell_in_circles(
        coords = coords,
        center = coords[s0_index[i], ],
        n = 1,
        r = r[j],
        index_only = TRUE)
      x = as.numeric(radar$data[, ii])
      quantiles = quantile(x[x > 0], probs = probs, na.rm = TRUE) |>
        unname()
      res[[j]] = data.frame(
        quantile = quantiles,
        prob = probs,
        i = i,
        r = r[j])
    }
    res = do.call(rbind, res)
    res
  })
stopCluster(cl)

quantile_data = do.call(rbind, df) |>
  tidyr::pivot_wider(names_from = "prob", values_from = "quantile") |>
  dplyr::select(-i)

plots = lapply(
  X = seq_along(probs),
  FUN = function(j) {
    lapply(
      X = seq_along(r),
      FUN = function(i) {
        coords[loc_index, ] |>
          as.data.frame() |>
          dplyr::mutate(value = dplyr::filter(quantile_data, r == !!r[i])[[j + 1]]) |>
          ggplot() +
          geom_raster(aes(x = X, y = Y, fill = value)) +
          scale_fill_viridis_c() +
          labs(fill = paste("r =", r[i]), title = paste("p =", probs[j]))
      })
  })

pdf("Rplots.pdf", width = 15, height = 8)
for (i in seq_along(plots)) print(patchwork::wrap_plots(plots[[i]]))
dev.off()

# tikz_plots = lapply(
#   X = seq_len(ncol(param_estimators) - 1),
#   FUN = function(j) {
#     names(param_estimators)[-1] = paste0("$", c("\\sigma", "\\xi", "u"), "$")
#     coords[loc_index, ] |>
#       as.data.frame() |>
#       dplyr::mutate(value = dplyr::filter(param_estimators, r == 5)[[j + 1]]) |>
#       ggplot() +
#       geom_raster(aes(x = X, y = Y, fill = value)) +
#       scale_fill_viridis_c() +
#       labs(fill = names(param_estimators)[j + 1],
#            x = "Easting", y = "Northing") +
#       theme_light() +
#       theme(axis.title.y = element_text(angle = 0, vjust = .5)) +
#       scale_x_continuous(expand = c(0, 0)) +
#       scale_y_continuous(expand = c(0, 0))
#   })
# 
# tikz_plot(file.path(image_dir(), "marginal-distributions.pdf"),
#           patchwork::wrap_plots(tikz_plots), width = 15, height = 8)
