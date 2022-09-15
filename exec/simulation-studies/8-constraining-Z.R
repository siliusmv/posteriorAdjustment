devtools::load_all()
library(ggplot2)
library(mvtnorm)

# Matérn correlation function
matern_corr = function(dist, rho, nu = 1.5) {
  kappa = sqrt(8 * nu) / rho
  res = 2 ^ (1 - nu) / gamma(nu) * (kappa * dist) ^ nu * besselK(kappa * dist, nu)
  res[dist == 0] = 1
  res
}

# Create a regular spacing of locations, with s0 = 0 in the middle
rho = 1 # Range parameter of the Gaussian random field
n_grid = 1000 # Number of points in the grid
s0_index = n_grid / 2 # Location of s0
loc = c(seq(-3.1 * rho, 0, length.out = n_grid / 2),
        seq(0, 3.1 * rho, length.out = n_grid / 2)[-1])

# Compute the correlation matrix for the unconstrained Gaussian field Z(s)
dist = as.matrix(dist(loc))
corr = matern_corr(dist, rho)

# Sample many realisations of Z(s)
n = 1e4
set.seed(1)
z = mvtnorm::rmvnorm(n, sigma = corr)

# Constrain Z(s) using subtraction
for (i in 1:n) {
  z[i, ] = z[i, ] - z[i, s0_index]
}

# Estimate the correlation structure of the constrained samples
constrained_corr = suppressWarnings(cor(z))

# Plot the correlation structure
plot = data.frame(
  x = rep(loc, length(loc)),
  y = rep(loc, each = length(loc)),
  corr = as.numeric(constrained_corr)) |>
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = corr)) +
  scale_fill_viridis_c(limits = c(-1, 1), option = "B") +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = rho * (-3:3),
    labels = c(paste0("$", -3:-1, "\\rho$"), "0", paste0("$", 1:3, "\\rho$"))) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = rho * (-3:3),
    labels = c(paste0("$", -3:-1, "\\rho$"), "0", paste0("$", 1:3, "\\rho$"))) +
  coord_equal() +
  labs(fill = "Correlation", x = "$x_1$", y = "$x_2$") +
  theme_light() +
  theme(axis.title.y = element_text(angle = 0, vjust = .5),
        text = element_text(size = 20))

tikz_plot(
  file = file.path(image_dir(), "constrained-correlation.pdf"),
  plot = plot, width = 8, height = 6)
