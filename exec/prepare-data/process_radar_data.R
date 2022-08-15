devtools::load_all()
library(dplyr)
library(ggplot2)
library(sf)
library(lubridate)
library(terra)
library(raster)

# Load radar data
radar_file = file.path(downloads_dir(), "radar-trondheim.nc")
data = raster::brick(radar_file)

# Coordinates of the data
coords = raster::rasterToPoints(data[[1]], spatial = TRUE) |>
  sf::st_as_sf() |>
  {\(x) cbind(x, sf::st_coordinates(x))}() |>
  dplyr::select(geometry)

# Transform from [m] to [km]
myproj = st_crs(coords)[["input"]] |>
  sub("units=m", "units=km", x = _) |>
  st_crs()
coords = st_transform(coords, myproj)

# Add covariates (something is weird with the tif-file. the raster package crashes horribly when used)
height_raster = terra::rast(file.path(downloads_dir(), "dem.tif"))
transformed_coords = st_coordinates(st_transform(coords, crs(height_raster)))
coords$height = unlist(terra::extract(height_raster, transformed_coords))
coords$height = ifelse(is.na(coords$height), 0, coords$height)

# Detect the times that the Rissa radar is inactive
time_series = data[1, 132] # This is the upper right corner of the xy-plane
inactive_time_index = which(is.na(as.numeric(time_series)))
active_time_index = seq_along(as.numeric(time_series))[-inactive_time_index]
times = colnames(time_series) |>
  sub("X", "", x = _) |>
  lubridate::ymd_hms()
active_times = times[active_time_index]

# Location of the Rissa radar
rissa = st_point(c(10.203845, 63.690527)) |>
  st_sfc(crs = 4326) |>
  st_transform(st_crs(coords))

# Plot a time series of the data
for (i in 30:70) {plot(data[[i]]); Sys.sleep(.5)}

# Remove the Rissa radar and its closest neighbours, because
# there is a lot of artificial noise in the data close to the radar
dist_to_rissa = as.numeric(st_distance(rissa, coords))
bad_radius = 5
bad_location_index = which(dist_to_rissa < bad_radius)
good_location_index = seq_along(dist_to_rissa)[-bad_location_index]
good_coords = coords[good_location_index, ]
good_coords$index = seq_len(nrow(good_coords))

month_nrs = 6:8
time_index = which(lubridate::month(times) %in% month_nrs) |>
  intersect(active_time_index)

data_list = pbapply::pblapply(
  X = seq_along(time_index),
  cl = 20,
  FUN = function(i) data[[time_index[i]]][good_location_index])

data_matrix = do.call(rbind, data_list)

mydata = list(
  data = data_matrix,
  coords = good_coords,
  times = times[time_index])
filename = file.path(downloads_dir(), "radar-trondheim.rds")

saveRDS(mydata, filename)


# ======================================================================
# Plots for the paper
# ======================================================================

library(ggplot2)
library(patchwork)
library(gtable)
library(grid)

bbox = st_bbox(coords) |>
  st_as_sfc() |>
  st_transform(crs(height_raster))
heights = as.data.frame(terra::crop(height_raster, bbox), xy = TRUE)
circle_center = st_coordinates(st_transform(rissa, st_crs(height_raster)))
circle_df = data.frame(θ = seq(0, 2 * pi, length.out = 1000)) |>
  dplyr::mutate(x = circle_center[1] + bad_radius * cos(θ),
                y = circle_center[2] + bad_radius * sin(θ))
bbox2 = st_bbox(c(xmin = 265, ymin = 7080, xmax = 295, ymax = 7110), crs = st_crs(coords)) |>
  st_as_sfc()
bbox3 = st_bbox(c(xmin = 210, ymin = 7035, xmax = 250, ymax = 7065), crs = st_crs(coords)) |>
  st_as_sfc()

plot = ggplot(heights) +
  geom_raster(aes(x = x, y = y, fill = dem)) +
  geom_sf(data = bbox2, fill = NA, size = 3, col = "black") +
  #geom_sf(data = bbox3, fill = NA, size = 3, col = "black") +
  scale_fill_viridis_c()
plot = style_map_plot(plot, bbox, use_tex = FALSE, col = NA) +
  geom_point(data = as.data.frame(circle_center), aes(x = X, y = Y), shape = 20, size = 5) +
  labs(fill = "Altitude [m]") +
  geom_path(data = circle_df, aes(x = x, y = y), size = 3) +
  theme(text = element_text(size = 18)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
#plot


bbox4 = st_bbox(c(xmin = 4, ymin = 57, xmax = 31, ymax = 72))
plot1 = ggplot() +
  geom_blank() +
  geom_sf(data = bbox, fill = NA, size = 1, col = "black") +
  add_norway_map(crs = 4326, bbox = bbox4) +
  theme_light() +
  theme(text = element_text(size = 18))

plot2 = plot +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())

plot3 = patchwork::wrap_plots(plot1, plot2)
p = patchwork::patchworkGrob(plot3)

theme_nothing = theme(
  axis.text = element_blank(),
  axis.title = element_blank(),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.ticks.length = unit(0, "cm"),
  axis.ticks.margin = unit(0.01, "cm"),
  panel.margin = unit(0, "lines"),
  plot.margin = unit(c(0, 0, -.5, -.5), "lines"),
  plot.background = element_blank(),
  legend.position = "none",
  complete = TRUE)

arrow1 = ggplot() +
  geom_segment(data = data.frame(x = .28, y = .52, xend = 1, yend = 1),
               aes(x = x, y = y, xend = xend, yend = yend), size = 1) +
  theme_nothing +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
arrow1 = ggplotGrob(arrow1)
arrow2 = ggplot() +
  geom_segment(data = data.frame(x = .285, y = .445, xend = 1, yend = 0.02),
               aes(x = x, y = y, xend = xend, yend = yend), size = 1) +
  theme_nothing +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
arrow2 = ggplotGrob(arrow2)

pp = gtable_add_grob(p, arrow1, t = 11, l = 9, r = 23, name = "a1")
pp = gtable_add_grob(pp, arrow2, t = 11, l = 9, r = 23, name = "a2")

# png("Rplots.png", width = 800, height = 400)
# grid.newpage()
# grid.draw(pp)
# dev.off()
# system("open Rplots.png")

png(file.path(image_dir(), "data.png"), width = 800, height = 400)
grid.newpage()
grid.draw(pp)
dev.off()
