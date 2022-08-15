
#' @export
data_dir = function() file.path(here::here(), "inst", "extdata")
#' @export
downloads_dir = function() file.path(data_dir(), "downloads")
#' @export
image_dir = function() file.path(data_dir(), "images")
#' @export
tmp_dir = function() file.path(data_dir(), "tmp")
#' @export
results_dir = function() file.path(data_dir(), "results")
#' @export
cgeneric_dir = function() file.path(here::here(), "cgeneric")

#' @export
compile_cgeneric = function(name) {
  files = file.path(cgeneric_dir(), paste0(name, c(".c", ".o", ".so")))
  shell_safe_files = shQuote(files)
  if (Sys.info()["sysname"] == "Darwin") {
    command = "gcc-11"
  } else {
    command = "gcc"
  }
  execute_shell_script(
    command = command,
    args = c(paste0("-Wall -fopenmp -fpic -g -O -c -o"), shell_safe_files[2], shell_safe_files[1]))
  execute_shell_script(
    command = command,
    args = c("-fopenmp -shared -o", shell_safe_files[3], shell_safe_files[2]))
  file.remove(files[2])
}

#' @export
log_sum_exp = function(x, na.rm = FALSE) {
  m = max(x, na.rm = na.rm)
  if (is.infinite(m) && m < 0) {
    res = -Inf
  } else {
    res = m + log(sum(exp(x - m), na.rm = na.rm))
  }
  res
}

#' @export
log_mean_exp = function(x, na.rm = FALSE) {
  log_sum_exp(x, na.rm) - log(sum(!is.na(x)))
}


#' @export
dist_euclid = function(x, y) {
  stopifnot(is.matrix(y))
  if (is.matrix(x)) {
    res = dist_euclid_mat_mat(x, y)
  } else {
    res = dist_euclid_vec_mat(x, y)
  }
  res
}

dist_euclid_vec_mat = function(x, y) {
  stopifnot(length(x) == ncol(y))
  tmp = 0
  for (i in seq_along(x)) {
    tmp = tmp + (x[i] - y[, i])^2
  }
  sqrt(as.numeric(tmp))
}

dist_euclid_mat_mat = function(x, y) {
  sapply(
    X = seq_len(nrow(x)),
    FUN = function(i) dist_euclid_vec_mat(x[i, ], y))
}

#' @export
execute_shell_script = function(command, args, ...) {
  output = system2(command, args, ...)
  success = (output == 0)
  if (!success) {
    formatted_args = paste(args, collapse = " ")
    stop("shell script ", command, " failed with arguments ", formatted_args)
  }
  0
}

#' @export
progress_bar = function(n) {
  pb = progress::progress_bar$new(
    format = ":percent [:bar] time elapsed: :elapsedfull, eta: :eta",
    total = n, width = 70, clear = FALSE)
  res = list()
  res$terminate = pb$terminate
  # pb throws an error if we tick too many times. I don't want that to happen
  res$tick = function(...) tryCatch(pb$tick(...), error = function(e) NULL)
  res
}

#' @export
tikz_plot = function(file, plot = NULL, expression = NULL, ...) {
  # Ensure that you are on an operating system that you have tested
  operating_system = Sys.info()[["sysname"]]
  if (operating_system == "Windows") {
    proceed = readline(paste("This function was written on a Mac,",
                             "I have no idea if it will work on Windows.",
                             "Proceed? (y/n) "))
    if (proceed != "y") return()
  }

  # Create a temporary file for the tikz-output
  tmp = tempfile(tmpdir = getwd())
  # Clean up after yourself on early interrupt
  on.exit(suppressWarnings(file.remove(tmp)), add = TRUE)

  # Extract default tex usepackages and add the bm package for bold greek letters
  opt = options()
  on.exit(options(opt)) #Reset global options on exit
  tikzDevice::setTikzDefaults(overwrite = FALSE)
  tex_packages = options()$tikzLatexPackages
  if (!any(grepl("usepackage\\{bm\\}", tex_packages))) {
    tex_packages = c(tex_packages, "\\usepackage{bm}\n")
  }
  if (!any(grepl("usepackage\\{amsmath\\}", tex_packages))) {
    tex_packages = c(tex_packages, "\\usepackage{amsmath}\n")
  }

  # Open a device for creating a tex-file
  tikzDevice::tikz(tmp, standAlone = TRUE, packages = tex_packages, ...)
  # Call dev.off() on exit in case of interruptions
  current_device = dev.cur()
  on.exit(dev.off(current_device))

  # Plot something into the tex-file
  if (!is.null(plot)) {
    if (any(class(plot) %in% c("gg", "ggplot", "patchwork"))) {
      print(plot)
    } else {
      for (p in plot) print(p)
    }
  } else {
    eval(substitute(expression), envir = parent.frame())
  }

  # Finish the creation of the tex-file
  dev.off()

  # Compile to pdf using lualatex
  system2("lualatex", shQuote(tmp))

  # Copy pdf file to final destination
  file.copy(paste0(tmp, ".pdf"), file, overwrite = TRUE)

  # Clean up all temporary files
  tmp_filename = tail(strsplit(tmp, "/")[[1]], 1)
  files_to_clean = grep(tmp_filename, list.files(full.names = TRUE), value = TRUE)
  file.remove(files_to_clean)
}
