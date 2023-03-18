# funtimes: Functions for Time Series Analysis

[![R-CMD-check](https://github.com/vlyubchich/funtimes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vlyubchich/funtimes/actions/workflows/R-CMD-check.yaml)

This is the version of the package under development. To use the most recent stable version from CRAN, run
```r
install.packages("funtimes")
library(funtimes)
```

To get this version from GitHub, run
```r
devtools::install_github("vlyubchich/funtimes")
library(funtimes)
```

## Future work

-   Combine functions for tail comparison into one, without dots.
-   Remove dependence on Jmisc.

## CRAN version 9.1

-   Replace the `BIC` argument in AR estimation with `ic`, to easily switch to AIC. 
The directly affected functions are `notrend_test()`, `sync_test()`, and `wavk_test()`.

## CRAN version 9.0

- Added functions for out-of-sample Granger causality testing (based on prediction errors), with an option to restrict (disregard) near-contemporaneous lags: `causality_pred()` and `causality_predVAR()`.

## CRAN version 8.2

- Removed dependency on the package `TDA` by omitting the TDA-based clustering function `TopoCBN()`.
- Added a vignette on clusters of time series.

## CRAN version 8.1

- `mcusum_test` now runs faster, has the option to use heuristic shortened bootstrapping, and an option to set `m <= length(k)` (used to be `m = length(k)`) thanks to the contributions from Alexander Brenning.

## CRAN version 8.0

-   Deprecated functions (with dots) marked as defunct. <!-- https://devguide.ropensci.org/evolution.html  -->
-   New functions: `AuePolyReg_test()`, `ccf_boot()`, `cumsumCPA_test()`, `GombayCPA_test()`, and `TopoCBN()`.
-   Fixed a bug in `beales()` thanks to an email from Dave Lorenz.

## CRAN version 7.0

-   Use underscores in function names instead of dots. Mark those functions as deprecated. <!-- https://mirai-solutions.ch/news/2017/12/05/roxygen2-deprecate/ https://devguide.ropensci.org/evolution.html -->
-   Add vignettes for Beale's estimator and trend tests.
-   Add a package overview section in the documentation.
-   Format R code according to R style diagnostics (such as whitespaces).
-   Improve computational efficiency of `purity()` function, inspired by Brian Simmons's suggestions from 2019-06-12 about handling zero cases when all matches have been found.

## CRAN version 6.1 and below

See GitHub commits before 2019-04-15
