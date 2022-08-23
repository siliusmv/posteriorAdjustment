# posteriorAdjustment

## Installation

Start by calling 
```r
renv::load()
```
to get the correct libraries and ensure reproducible results. Then you need to call
```r
devtools::install()
```
to install the package and compile the Rcpp functions correctly. It is also possible to use
```r
devtools::load_all()
```
instead of properly installing the R package.

Finally, you need to call
```r
make_cgeneric("all")
```
in order to compile and link the necessary cgeneric models. Before calling this function, you might
need to change the C compiler `GCC` (by default equal to `gcc-12`) used in the makefile
`cgeneric/Makefile` to a compiler that is available on your system.

## order of exec scripts (and their belonging R/ functions)

1. simulation-studies/
   1. univariate
   2. SPDE
   3. block
   4. conditional parameter recovering (working on it)
   5. conditional-extremes adjustment (not looked at yet)
   6. self-inconsistency (not looked at yet)

R scripts

0. utils.R
1. compute-c-matrix.R
2. sampling.R (I have cleaned, but not tested this)
3. download-data.R
4. data-processing.R
5. marginal-distributions.R
6. plotting.R (not looked at yet)
7. conditional-likelihood.R
8. inla-generic-models.R

Rcpp scripts
1. dmvnorm.cpp

test scripts
1. test-Rcpp.R: (working on it)
2. test-cgeneric (working on it)

