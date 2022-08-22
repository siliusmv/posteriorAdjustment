devtools::load_all()
library(INLA)
library(Matrix)

# Create a small grid of locations
loc = as.matrix(expand.grid(x = 1:3, y = 1:2))

# Create a mesh on the grid
mesh = inla.mesh.2d(
  loc = loc,
  boundary = list(inla.nonconvex.hull(loc, convex = -.3)),
  cutoff = 1,
  max.edge = 4)
plot(mesh)
points(loc)
mesh$n

# Create the SPDE on the mesh
spde = inla.spde2.pcmatern(
  mesh,
  prior.range = c(rho, .5),
  prior.sigma = c(sigma, .5))

# ====================================================================
# Save the B- and M-matrices from the SPDE object in a file format
# that is readable for our c code, so we can compare precision matrices
# made by our c functions with matrices made using inla.spde2.precision()
# ====================================================================

# Functions for saving the matrices
save_mat = function(m, path) {
  file.create(path)
  cat("mat", nrow(m), ncol(m), sep = "\n", file = path)
  cat(t(m), sep = "\n", file = path, append = TRUE)
}
save_smat = function(m, path) {
  file.create(path)
  n = length(m@x)
  cat("smat", n, nrow(m), ncol(m), sep = "\n", file = path)
  for (i in seq_len(n)) {
    cat(m@i[i], m@j[i], m@x[i], sep = ";", file = path, append = TRUE)
    cat("\n", file = path, append = TRUE)
  }
}

# Create the directory for saving the matrices
path = file.path(here::here(), "cgeneric-data")
dir.create(path, recursive = TRUE)

for (i in 0:2) {
  # Save all the B-mats
  filename = file.path(path, paste0("B", i, ".txt"))
  save_mat(spde$param.inla[[paste0("B", i)]], filename)
  # Save all the M-smats
  filename = file.path(path, paste0("M", i, ".txt"))
  save_smat(spde$param.inla[[paste0("M", i)]], filename)
}


# ============================================================
# Compute Q using c and using inla.spde2.precision(), and
# compare the results
# ============================================================

execute_c_script = function(name, ...) {
  current_path = getwd()
  on.exit(setwd(current_path))
  setwd(cgeneric_dir())
  execute_shell_script(name, ...)
}

compare_r_and_c = function(log_rho, log_sigma) {
  Q = inla.spde2.precision(spde, c(log_rho, log_sigma))
  execute_c_script("./test1.o", c(log_rho, log_sigma))
  round(as.matrix(Q), digits = 6)
}

test_c_further = function(log_rho, log_sigma) {
  execute_c_script("./test2.o", c(log_rho, log_sigma))
}

make_cgeneric("test")

log_rho = runif(1, 0, 2)
log_sigma = runif(1, 0, 2)
compare_r_and_c(log_rho, log_sigma)

test_c_further(log_rho, log_sigma)

unlink(path)
