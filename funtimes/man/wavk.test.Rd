\name{wavk.test}
\alias{wavk.test}
\title{WAVK Trend Test}
\description{
Non-parametric test to detect (non-)monotonic parametric trends in time series.
}

\usage{
wavk.test(formula, factor.length=c("user.defined", "adaptive.selection"), 
          Window=NULL, q=3/4, j=c(8:11), B=1000, method=c("boot", "asympt"), 
         ar.order=NULL, ar.method = "HVK", BIC=TRUE, out=FALSE)
}
\arguments{
  \item{formula}{an object of class "formula", specifying the form of the parametric time trend to be tested. Variable \eqn{t} should be used to specify the form, where \eqn{t} is specified within the function as a regular sequence on the interval (0,1]. See 'Examples'.}

  \item{factor.length}{method to define the length of local windows (factors). Default option \code{"user.defined"} allows to set only one value of the argument \code{Window}. The option \code{"adaptive.selection"} sets \code{method = "boot"} and employs heuristic \eqn{m}-out-of-\eqn{n} subsampling algorithm (Bickel and Sakov, 2008) to select an optimal window from the set of possible windows \code{length(x)*q^j} whose values are mapped to the largest previous integer and greater than 2. \code{x} is the time series tested.}

  \item{Window}{length of the local window (factor), default is \code{round(0.1*length(x))}, where \code{x} is the time series tested. This argument is ignored if \code{factor.length = "adaptive.selection"}.}

  \item{q}{scalar from 0 to 1 to define the set of possible windows when \code{factor.length =} \code{"adaptive.selection"}. Default is \eqn{3/4}. This argument is ignored if \code{factor.length =} \code{"user.defined"}.}

  \item{j}{numeric vector to define the set of possible windows when \code{factor.length =} \code{"adaptive.selection"}. Default is \code{c(8:11)}. This argument is ignored if \code{factor.length =} \code{"user.defined"}.}

  \item{B}{number of bootstrap simulations to obtain empirical critical values.}

  \item{method}{method of obtaining critical values: from asymptotical (\code{"asympt"}) or bootstrap (\code{"boot"}) distribution. If \code{factor.length =} \code{"adaptive.selection"} the option \code{"boot"} is used.}

  \item{ar.order}{order of autoregressive filter when \code{BIC = FALSE}, or the maximal order for BIC-based filtering. Default is \code{round(10*log10(length(x)))}, where \code{x} is the time series tested.}
  
  \item{ar.method}{method of estimating autoregression coefficients. Default \code{"HVK"} delivers robust difference-based estimates by Hall and Van Keilegom (2003). Alternatively, options of \command{ar} function can be used, such as \code{"burg"}, \code{"ols"}, \code{"mle"}, and \code{"yw"}.}

  \item{BIC}{logical value indicates whether the order of autoregressive filter should be selected by Bayesian information criterion (BIC). If \code{TRUE} (default), models of orders 0,1,...,\code{ar.order} or 0,1,...,\code{round(10*log10(length(x)))} are considered, depending on whether \code{ar.order} is defined or not. \code{x} is the time series tested.}
  
  \item{out}{logical value indicates whether full output should be shown. Default is \code{FALSE}.}
}


\value{
A list with class htest containing the following components:

  \item{method}{name of the method.}

  \item{data.name}{name of the data.}

  \item{statistic}{value of the test statistic.}

  \item{p.value}{\eqn{p}-value of the test.}

  \item{alternative}{alternative hypothesis.}

  \item{parameter}{window that was used.}

  \item{estimate}{list with the following elements: estimated trend coefficients; user-defined or BIC-selected AR order; estimated AR coefficients; and, if \code{factor.length =} \code{"adaptive.selection"}, test results for all considered windows.}
}

\references{
Bickel, P. J. and Sakov, A. (2008) On the choice of \eqn{m} in the \eqn{m} out of \eqn{n} bootstrap and confidence bounds for extrema. \emph{Statistica Sinica} 18, 967--985.

Hall, P. and Van Keilegom, I. (2003) Using difference-based methods for inference in nonparametric regression with time series errors. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 65: 443--456. DOI: \href{http://dx.doi.org/10.1111/1467-9868.00395}{10.1111/1467-9868.00395}

Lyubchich, V., Gel, Y. R., and El-Shaarawi, A. (2013) On detecting non-monotonic trends in environmental time series: a fusion of local regression and bootstrap. \emph{Environmetrics} 24, 209--226. DOI: \href{http://dx.doi.org/10.1002/env.2212}{10.1002/env.2212}

Wang, L., Akritas, M. G., and Van Keilegom, I. (2008) An ANOVA-type nonparametric diagnostic test for heteroscedastic regression models. \emph{Journal of Nonparametric Statistics} 20(5), 365--382. DOI: \href{http://dx.doi.org/10.1080/10485250802066112}{10.1080/10485250802066112}

Wang, L. and Van Keilegom, I. (2007) Nonparametric test for the form of parametric regression with time series errors. \emph{Statistica Sinica} 17, 369--386.
}

\author{Yulia R. Gel, Vyacheslav Lyubchich, Ethan Schaeffer}

\seealso{\code{\link[stats]{ar}},
         \code{\link{HVK}}, 
         \code{\link{WAVK}},
         \code{\link{sync.test}}
}

\examples{

# Fix seed for reproduceable simulations:
set.seed(1)

#Simulate autoregressive time series of length n with smooth quadratic trend:
n <- 100
U <- 1 + 2*(1:n/n) + 4*(1:n/n)^2 + arima.sim(n=n, list(order=c(2,0,0), ar=c(-0.7, -0.1)))
plot.ts(U)


#Test H0 of a linear trend, with m-out-of-n selection of the local window:
\dontrun{
wavk.test(U~t, factor.length="adaptive.selection")}
# Sample output:
##	Trend test by Wang, Akritas, and Van Keilegom (bootstrap p-values)
##
##data:  U 
##WAVK test statistic = 5.3964, adaptively selected window = 4, p-value < 2.2e-16
##alternative hypothesis: trend is not of the form U ~ t.


#Test H0 of a quadratic trend, with m-out-of-n selection of the local window 
#and output of all results:
\dontrun{
wavk.test(U~poly(t,2), factor.length="adaptive.selection", out=TRUE)}
# Sample output:
##	Trend test by Wang, Akritas, and Van Keilegom (bootstrap p-values)
##
##data:  U 
##WAVK test statistic = 0.40083, adaptively selected window = 4, p-value = 0.576
##alternative hypothesis: trend is not of the form U ~ poly(t, 2).
##sample estimates:
##$trend_coefficients
##(Intercept) poly(t, 2)1 poly(t, 2)2 
##   3.408530   17.681422    2.597213 
##
##$AR_order
##[1] 1
##
##$AR_coefficients
##         phi_1 
##[1] -0.7406163
##
##$all_considered_windows
## Window WAVK-statistic p-value
##      4     0.40083181   0.576
##      5     0.06098625   0.760
##      7    -0.57115451   0.738
##     10    -1.02982929   0.360


# Test H0 of no trend (constant trend) using asymptotic distribution of statistic.
wavk.test(U~1, method="asympt")
# Sample output:
##	Trend test by Wang, Akritas, and Van Keilegom (asymptotic p-values)
##
##data:  U 
##WAVK test statistic = 25.999, user-defined window = 10, p-value < 2.2e-16
##alternative hypothesis: trend is not of the form U ~ 1.
}
\keyword{ts}
\keyword{htest}
\keyword{trend}

