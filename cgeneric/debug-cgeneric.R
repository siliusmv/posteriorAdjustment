devtools::load_all()
library(INLA)
library(Matrix)

#loc = as.matrix(expand.grid(x = 1:4, y = 1:3))
loc = as.matrix(expand.grid(x = 1:3, y = 1:3))
σ = 2
ρ = 5

# Create the mesh
mesh = inla.mesh.2d(
  loc.domain = loc,
  min.angle = 30,
  boundary = list(inla.nonconvex.hull(loc, convex = -.3),
                  inla.nonconvex.hull(loc, convex = -1)),
  cutoff = 1,
  #cutoff = .4,
  #max.edge = c(1.3))
  max.edge = c(2))
plot(mesh); points(loc)

# Create the SPDE
spde = inla.spde2.pcmatern(
  mesh,
  prior.range = c(ρ, .5),
  prior.sigma = c(σ, .5))

# Functions for saving data in a format that is readable for the c code
save_mat = function(m, path) {
  file.create(path)
  cat("mat", nrow(m), ncol(m), sep = "\n", file = path)
  cat(m, sep = "\n", file = path, append = TRUE)
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

# Save all the necessary data
path = file.path(here::here(), "cgeneric-data")
dir.create(path)

for (i in 0:2) {
  # Save all the B-mats
  filename = file.path(path, paste0("B", i, ".txt"))
  save_mat(spde$param.inla[[paste0("B", i)]], filename)
  # Save all the M-smats
  filename = file.path(path, paste0("M", i, ".txt"))
  save_smat(spde$param.inla[[paste0("M", i)]], filename)
}

# inla.spde2.theta2phi0(spde, c(log(5), log(2)))
# inla.spde2.theta2phi1(spde, c(log(5), log(2)))
# inla.spde2.theta2phi2(spde, c(log(5), log(2)))

Q = inla.spde2.precision(spde, c(log(5), log(2)))
Q@x
length(Q@x)

Q[1:5, 1:5]
Q[20:28, 20:28]
tail(Q@x)
