devtools::load_all()
library(INLA)
library(parallel)
library(numDeriv)
library(purrr)

σ = 1
n_train = 100
n_repl = 1000
n_cores = 15
n_posterior_samples = 10000
probs = c(.9, .95, .99)

filename = file.path(tmp_dir(), "misspecified-model-adjustment.rds")

ll = function(μ, y, σ, return_vec = FALSE) {
  res = dnorm(y, mean = μ, sd = σ, log = TRUE)
  if (!return_vec) res = sum(res)
  res
}

set.seed(1, kind = "L'Ecuyer-CMRG")
tmp = parallel::mclapply(
  X = seq_len(n_repl),
  mc.cores = n_cores,
  FUN = function(i) {
    y = rt(n_train, 1)
    res = inla(
      formula = y ~ 1,
      data = data.frame(y = y),
      control.compute = list(config = TRUE),
      control.fixed = list(prec.intercept = 0.001),
      control.family = list(hyper = list(prec = list(initial = log(σ^-2), fixed = TRUE))),
      control.inla = list(int.strategy = "grid", dz = .05, diff.logdens = 10),
      num.threads = 1)
    mode = inla.mmarginal(res$marginals.fixed[[1]])

    samples = inla.posterior.sample(
      n = n_posterior_samples,
      result = res,
      selection = list("(Intercept)" = 1),
      seed = 1L)
    samples = sapply(samples, `[[`, "latent")

    H = - numDeriv::hessian(ll, x = mode, y = y, σ = σ)
    grads = numDeriv::jacobian(ll, x = mode, y = y, σ = σ, return_vec = TRUE)
    J = sum(grads^2)
    C = get_C(H, J)

    samples = list(
      unadjusted = samples,
      adjusted = mode + (samples - mode) * as.numeric(C))

    if (i %% 10 == 0) message("Progress: ", i, " / ", n_repl)
    quantiles = list()
    for (name in names(samples)) {
      quantiles[[name]] = quantile(
        samples[[name]],
        probs = c((1 - probs) / 2, 1 - (1 - probs) / 2))
    }

    quantiles
  })
quantiles = purrr::transpose(tmp)
saveRDS(quantiles, filename)

quantiles = readRDS(filename)
for (name in names(quantiles)) {
  quantiles[[name]] = do.call(rbind, quantiles[[name]])
}

coverage = list()
for (name in names(quantiles)) {
  zero_above_lower = quantiles[[name]][, seq_along(probs)] < 0
  zero_below_upper = quantiles[[name]][, length(probs) + seq_along(probs)] > 0
  stopifnot(all(dim(zero_above_lower) == dim(zero_below_upper)))
  inside_interval = zero_below_upper & zero_above_lower
  coverage[[name]] = apply(inside_interval, 2, mean)
  names(coverage[[name]]) = paste(100 * probs, "%")
}
coverage

tmp = coverage |>
  do.call(what = cbind) |>
  (\(x) cbind(probs, x))() |>
  (\(x) round(x * 100, digits = 0))()
colnames(tmp) = c("Aim", "Unadjusted", "Adjusted")
rownames(tmp) = NULL

table = paste(paste(colnames(tmp), collapse = " & "), "\\\\")
table[2] = "\\midrule"
for (i in 1:nrow(tmp)) {
  table[[length(table) + 1]] = paste(paste0("$", tmp[i, ], "\\%$", collapse = " & "), "\\\\")
}
table = paste(paste(table, collapse = "\n"), "\n")

cat(table)
