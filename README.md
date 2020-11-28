# Future work

* Mark deprecated functions as defunct. <!-- https://devguide.ropensci.org/evolution.html  -->
* Combine functions for tail comparison into one, without dots.
* Replace BIC argument in AR estimation with penalty, to easily switch to AIC or another penalty.
* Vignette trend tests: do need eval(parse...) in `wavk_test`? Add about `synch_test`.
* Add vignette on clusters of time series.


# CRAN version 7.0

* Use underscores in function names instead of dots. Mark those functions as deprecated. <!-- https://mirai-solutions.ch/news/2017/12/05/roxygen2-deprecate/ https://devguide.ropensci.org/evolution.html -->
* Add vignettes for Beale's estimator and trend tests.
* Add a package overview section in the documentation.
* Format R code according to R style diagnostics (such as whitespaces).
* Improve computational efficiency of `purity` function, inspired by Brian Simmons's suggestions from 2019-06-12 about handling zero-cases when all matches have been found.

# CRAN version 6.1 and below

See GitHub commits before 2019-04-15
