devtools::load_all()

# ==============================================================================
# Download radar data
# ==============================================================================

# Find all valid urls containing radar data
dates = seq(as.Date("2010-01-01"), as.Date("2021-12-31"), by = "1 month")
dates = dates[lubridate::month(dates) %in% c(6, 7, 8)]
urls = unlist(lapply(dates, get_radar_url))

# Arguments to CDO for downloading the data
args = c(
  # This selects an area of size 31x31 close to the radar
  "-selindexbox,561,591,1015,1045",
  # This is the name of the variable that contains hourly mean precipitation estimates
  "-selname,lwe_precipitation_rate")

filename = file.path(downloads_dir(), "radar.nc")

download(urls, args, filename)
