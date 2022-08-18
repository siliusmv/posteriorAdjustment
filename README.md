# posteriorAdjustment


Start by calling 
```r
renv::load()
```

## order of exec scripts (and their belonging R/ functions)

1. simulation-studies/
   1. univariate: only needs get\_C()
   2. SPDE: only needs get\_C()
   3. block: get\_C() and rnorm\_spde()
   4. conditional parameter recovering (not looked at yet)
   5. conditional-extremes adjustment (not looked at yet)
   6. self-inconsistency (not looked at yet)
2. case-study/
   1. download-data: needs the functions in R/download-data.R
   2. process-data: doesn't need any R functions
   3. examine-marginal-distributions: data-processing.R
   4. model-selection: (working on it)
   5. final-modelling: (working on it)

R scripts

0. utils.R
1. compute-c-matrix.R
2. sampling.R (I have cleaned, but not tested this)
3. download-data.R
4. data-processing.R (not cleaned yet)
5. marginal-distributions.R
6. plotting.R (not looked at yet)
7. conditional-model.R (not looked at yet)
8. inla-generic-models.R (not looked at yet)

Rcpp scripts
1. dmvnorm.cpp

test scripts
1. test-Rcpp.R: the code runs, but is really ugly...
