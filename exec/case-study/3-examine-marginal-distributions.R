devtools::load_all()
library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)

# ==============================================================================
# Load the radar data
# ==============================================================================
radar = readRDS(file.path(downloads_dir(), "radar.rds"))
n_loc = nrow(radar$coords)
n_time = length(radar$times)
coords = st_coordinates(radar$coords)

# ==============================================================================
# Examine how to choose radius for the sliding window by plotting
# a set of quantiles of the estimated marginal distributions in a map
# ==============================================================================

# Examine different radii
r = c(0, 1, 2, 3, 5, 7)
# Examine different quantiles
probs = c(.5, .6, .7, .8, .9, .95, .99, .995, .999)

# Compute quantiles of the estimated marginal distributions for
# all combinations of r and probs
pb = progress_bar(n_loc)
df = lapply(
  X = seq_len(n_loc),
  FUN = function(i) {
    res = list()
    for (j in seq_along(r)) {
      F = aggregated_ecdf(radar$data, coords, coords[i, ], r[j])
      quantiles = quantile(environment(F)$F, probs = probs) |>
        unname()
      res[[j]] = data.frame(
        quantile = quantiles,
        prob = probs,
        i = i,
        r = r[j])
    }
    res = do.call(rbind, res)
    pb$tick()
    res
  })
pb$terminate()

# Reformat the data for plotting with ggplot
quantile_data = do.call(rbind, df) |>
  tidyr::pivot_wider(names_from = "prob", values_from = "quantile") |>
  dplyr::select(-i)

# Create all the plots
plots = lapply(
  X = seq_along(probs),
  FUN = function(j) {
    lapply(
      X = seq_along(r),
      FUN = function(i) {
        coords |>
          as.data.frame() |>
          dplyr::mutate(value = dplyr::filter(quantile_data, r == !!r[i])[[j + 1]]) |>
          ggplot() +
          geom_raster(aes(x = X, y = Y, fill = value)) +
          scale_fill_viridis_c() +
          labs(fill = paste("r =", r[i]), title = paste("p =", probs[j]))
      })
  })

# Plot all the maps into the file Rplots.pdf
pdf(file.path(image_dir(), "marginal-distributions.pdf"), width = 15, height = 8)
for (i in seq_along(plots)) print(patchwork::wrap_plots(plots[[i]]))
dev.off()
