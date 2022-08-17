devtools::load_all()
library(dplyr)
library(sf)
library(lubridate)
library(raster)

# Load radar data
radar_file = file.path(downloads_dir(), "radar.nc")
data = raster::brick(radar_file)

# Get the coordinates of the data
coords = rasterToPoints(data[[1]], spatial = TRUE) |>
  st_as_sf() |>
  {\(x) cbind(x, st_coordinates(x))}() |>
  dplyr::select(geometry)

# Transform coordinates from [m] to [km] to remove unneccesary zeros
proj = st_crs(coords)[["input"]] |>
  sub(pattern = "units=m", replacement = "units=km") |>
  st_crs()
coords = st_transform(coords, proj)

# Get the times of the data
times = colnames(data[1, 1]) |>
  sub(pattern = "X", replacement = "") |>
  lubridate::ymd_hms()

# Extract the data into a matrix
data_matrix = raster::as.matrix(data)
colnames(data_matrix) = NULL
data_matrix = t(data_matrix)

mydata = list(
  data = data_matrix,
  coords = coords,
  times = times)
saveRDS(mydata, file.path(downloads_dir(), "radar.rds"))
