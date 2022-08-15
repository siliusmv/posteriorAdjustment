
set.seed(123)
sigma = bayesm::rwishart(10, diag(8))$IW
means = rnorm(8)
X = mvtnorm::rmvnorm(900000, means, sigma)
X_t = t(X)
nugget = 2

microbenchmark::microbenchmark(
  dmvnorm = mvtnorm::dmvnorm(X, means, sigma + diag(nugget, nrow(sigma)), log = TRUE),
  dmvnorm_arma = dmvnorm_arma(X_t, means, sigma + diag(nugget, nrow(sigma)), TRUE),
  dmvn = mvnfast::dmvn(X, means, sigma + diag(nugget, nrow(sigma)), log = TRUE),
  dconditional = dconditional_arma(X_t - means, sigma, nugget, TRUE))

y1 = mvtnorm::dmvnorm(X, means, sigma + diag(nugget, nrow(sigma)), log = TRUE)
y2 = dmvnorm_arma(X_t, means, sigma + diag(nugget, nrow(sigma)), TRUE)
y3 = mvnfast::dmvn(X, means, sigma + diag(nugget, nrow(sigma)), log = TRUE)
y4 = dconditional_arma(X_t - means, sigma, nugget, TRUE)
#y5 = sapply(
#  1:ncol(X_t),
#  function(i) {
#    dconditional(X_t[, i], a = means, Σ = sigma, τ = Inf, log = TRUE)
#  })

summary(y1 - y2)
summary(y2 - y3)
summary(y3 - y4)
#summary(y4 - y5)

X_t[3, 1] = NA
dmvnorm_arma(X_t[, 1, drop = FALSE], means, sigma + diag(nugget, nrow(sigma)), TRUE)
dconditional_arma(X_t[, 1, drop = FALSE] - means, sigma, nugget, na_rm = FALSE)
dconditional_arma(X_t[, 1, drop = FALSE] - means, sigma, nugget, TRUE)
mvtnorm::dmvnorm(X[1, -3, drop = FALSE], means[-3], sigma[-3, -3] + diag(nugget, nrow(sigma) - 1), TRUE)
