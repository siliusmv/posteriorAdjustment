
# The functions in this file are used for downloading climate data
# from online repositories.
# The main function is the download() function.
# The other functions are mainly helper functions

#' Given a list of lists of urls from https://thredds.met.no,
#' download data using the arguments args,
#' then concatenate all the data and save it into the file filename
#' @export
download = function(urls, args, filename) {
  tempfiles = get_tempfiles(filename, n = length(urls))
  for (i in seq_along(urls)) {
    file_already_exists = file.exists(tempfiles$files[i])
    if (!file_already_exists) {
      download_nc_files(urls[[i]], tempfiles$files[i], args)
    }
  }
  concatenate_nc_files(input = tempfiles$files, output = filename)
}

#' @export
concatenate_nc_files = function(input, output, args = NULL, ...) {
  execute_cdo_command(input, output, c("-cat", args), ...)
}

#' @export
download_nc_files = function(input, output, args = NULL, ...) {
  if (length(input) > 1) {
    tempfiles = get_tempfiles(output, length(input))
    for (i in seq_along(input)) {
      download_nc_files(input[i], tempfiles$files[i], args, ...)
    }
    concatenate_nc_files(tempfiles$files, output, args, ...)
    unlink(tempfiles$dir, recursive = TRUE)
  } else {
    execute_cdo_command(input, output, args, ...)
  }
}

download_one_nc_file = function(input, output, args = NULL, ...) {
  execute_cdo_command(input, output, args, ...)
}

#' @export
get_cdo_bbox_args = function(bbox) {
  paste0("-sellonlatbox,", bbox[1], ",", bbox[3], ",", bbox[2], ",", bbox[4])
}

#' @export
get_cdo_varname_args = function(varname) paste0("-selname,", varname)

#' @export
get_radar_url = function(date) {
  catalog_url = paste0(get_metno_base_url(), "remotesensingradaraccr/", lubridate::year(date), "/",
                       two_digit_month(date), "/catalog.html")
  files = locate_nc_files_in_metno_page(catalog_url)
  files = sort(files)
  get_full_url_of_metno_file(catalog_url, files)
}

locate_nc_files_in_metno_page = function(url) {
  page = xml2::read_html(url)
  links = rvest::html_nodes(page, "a")
  links = rvest::html_text(links)
  links = grep(".nc$", links, value = TRUE)
  links = grep(".nc.nc$", links, invert = TRUE, value = TRUE) # These files are bad, and unwanted
  links
}

get_full_url_of_metno_file = function(path, filename) {
  if (length(filename) == 0 || length(path) == 0) {
    warning("invalid metno_url")
    return(NULL)
  }
  path = sub("catalog/", "dodsC/", path)
  path = sub("catalog.html", "", path)
  paste0(path, filename)
}

two_digit_month = function(date) sprintf("%02d", lubridate::month(date))

get_metno_base_url = function() "https://thredds.met.no/thredds/catalog/"

get_tempfiles = function(filename, n) {
  tempdir = create_temp_dir(filename)
  tempfiles = paste0(tempdir, "/", seq_len(n), ".nc")
  list(dir = tempdir, files = tempfiles)
}

create_temp_dir = function(filename) {
  temp_dir = sub(".[^.]+$", "", filename)
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  temp_dir
}

execute_cdo_command = function(input, output, args = NULL) {
  warn_about_spaces(c(input, output))
  args = c(args, input, output)
  execute_shell_script("cdo", args)
}

warn_about_spaces = function(x) {
  if (any(grepl(" ", x))) {
    warning("There is a space in your file name. I don't think the code handles that correctly...")
  }
}
