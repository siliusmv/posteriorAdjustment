devtools::load_all()

# ==============================================================================
# Download radar data
# ==============================================================================

bbox = c(xmin = 9, ymin = 63.4, xmax = 11.5, ymax = 64.6)
dates = seq(as.Date("2010-01-01"), as.Date("2021-12-31"), by = "1 month")
urls = lapply(dates, get_radar_url)

filename = file.path(downloads_dir(), "radar-trondheim.nc")
variable = "lwe_precipitation_rate"
args = c(get_cdo_bbox_args(bbox), get_cdo_varname_args(variable))
download(urls, args, filename)

# ============================================================================================
# Download DEM
# ============================================================================================
filename = file.path(downloads_dir(), "dem.tif")
url = "https://hoydedata.no/LaserInnsyn/Home/DownloadFile/56"
zipfile = file.path(downloads_dir(), "dtm-50m-utm33.zip")
download.file(url, zipfile)
temp_dir = file.path(downloads_dir(), "dtm-50m-utm33")
unzip(zipfile, exdir = temp_dir)
file.remove(zipfile)

files = list.files(temp_dir, full.names = TRUE)
tif_files = grep("*.tif$", files, value = TRUE)

# rast = terra::rast(tif_files[1])
# old_crs = as.character(rast@crs)
# new_crs = sub("units=m", "units=km", old_crs)

command = "gdalwarp"
args = c(#"-s_srs", shQuote(old_crs), "-t_srs", shQuote(new_crs), # Change crs
         #"-tr 1 1", # Change the resolution to 1x1km
         #"-r bilinear", # ... using bilinear interpolation
         "-overwrite", tif_files, filename)

execute_shell_script(command, args)
unlink(temp_dir)
