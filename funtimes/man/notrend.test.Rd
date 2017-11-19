\name{notrend.test}
\alias{notrend.test}
\title{Sieve Bootstrap Based Test for the Null Hypothesis of no Trend}
\description{
A combination of time series trend tests for testing the null hypothesis of no trend, versus the alternative hypothesis of a linear trend (Student's t-test), or monotonic trend (Mann-Kendall test), or possibly non-monotonic trend (WAVK test).
}

\usage{
notrend.test(x, B = 1000, test = c("t", "MK", "WAVK"), 
            ar.method = "HVK", ar.order = NULL, BIC = TRUE, 
            factor.length = c("user.defined", "adaptive.selection"), 
            Window = NULL, q = 3/4, j = c(8:11))
}
\arguments{
  \item{x}{univariate time series. Missing values are not allowed.}
  
  \item{B}{number of bootstrap simulations to obtain empirical critical values.}
  
  \item{test}{trend test to implement: Student's t-test (\code{"t"}, default), Mann-Kendall test (\code{"MK"}), or WAVK test (\code{"WAVK"}, see \code{\link{WAVK}}).}

\item{ar.method}{method of estimating autoregression coefficients. Default \code{"HVK"} delivers robust difference-based estimates by Hall and Van Keilegom (2003). Alternatively, options of \command{ar} function can be used, such as \code{"burg"}, \code{"ols"}, \code{"mle"}, and \code{"yw"}.}

  \item{ar.order}{order of autoregressive filter when \code{BIC = FALSE}, or the maximal order for BIC-based filtering. Default is \code{round(10*log10(length(x)))}.}
  
  \item{BIC}{logical value indicates whether the order of autoregressive filter should be selected by Bayesian information criterion (BIC). If \code{TRUE} (default), models of orders 0,1,...,\code{ar.order} or 0,1,...,\code{round(10*log10(length(x)))} are considered, depending on whether \code{ar.order} is defined or not.}


  \item{factor.length}{method to define the length of local windows (factors). Used only if \code{test = "WAVK"}. Default option \code{"user.defined"} allows to set only one value of the argument \code{Window}. The option \code{"adaptive.selection"} sets \code{method = "boot"} and employs heuristic \eqn{m}-out-of-\eqn{n} subsampling algorithm (Bickel and Sakov, 2008) to select an optimal window from the set of possible windows \code{length(x)*q^j} whose values are mapped to the largest previous integer and greater than 2. \code{x} is the time series tested.}

  \item{Window}{length of the local window (factor), default is \code{round(0.1*length(x))}. Used only if \code{test = "WAVK"}. This argument is ignored if \code{factor.length = "adaptive.selection"}.}

  \item{q}{scalar from 0 to 1 to define the set of possible windows when \code{factor.length =} \code{"adaptive.selection"}. Used only if \code{test = "WAVK"}. Default is \eqn{3/4}. This argument is ignored if \code{factor.length =} \code{"user.defined"}.}

  \item{j}{numeric vector to define the set of possible windows when \code{factor.length =} \code{"adaptive.selection"}. Used only if \code{test = "WAVK"}. Default is \code{c(8:11)}. This argument is ignored if \code{factor.length =} \code{"user.defined"}.}

}

\details{
The current function tests the null hypothesis of no trend. For more general alternative shapes of the trend, use \code{\link{wavk.test}}. Notice that \code{\link{wavk.test}} employs hybrid bootstrap, which is alternative to the sieve bootstrap employed by the current function.
}

\value{
A list with class htest containing the following components:

  \item{method}{name of the method.}

  \item{data.name}{name of the data.}

  \item{statistic}{value of the test statistic.}

  \item{p.value}{\eqn{p}-value of the test.}

  \item{alternative}{alternative hypothesis.}

  \item{estimate}{list with the following elements: employed AR order and estimated AR coefficients.}

  \item{parameter}{window that was used in WAVK test, included in the output only if \code{test = "WAVK"}.}


}

\references{
Bickel, P. J. and Sakov, A. (2008). On the choice of \eqn{m} in the \eqn{m} out of \eqn{n} bootstrap and confidence bounds for extrema. \emph{Statistica Sinica} 18: 967--985.

Hall, P. and Van Keilegom, I. (2003). Using difference-based methods for inference in nonparametric regression with time series errors. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 65: 443--456. DOI: \href{http://dx.doi.org/10.1111/1467-9868.00395}{10.1111/1467-9868.00395}

Lyubchich, V., Gel, Y. R., and El-Shaarawi, A. (2013). On detecting non-monotonic trends in environmental time series: a fusion of local regression and bootstrap. \emph{Environmetrics} 24: 209--226. DOI: \href{http://dx.doi.org/10.1002/env.2212}{10.1002/env.2212}

Wang, L., Akritas, M. G., and Van Keilegom, I. (2008). An ANOVA-type nonparametric diagnostic test for heteroscedastic regression models. \emph{Journal of Nonparametric Statistics} 20(5): 365--382. DOI: \href{http://dx.doi.org/10.1080/10485250802066112}{10.1080/10485250802066112}

Wang, L. and Van Keilegom, I. (2007). Nonparametric test for the form of parametric regression with time series errors. \emph{Statistica Sinica} 17: 369--386.
}

\author{Yulia R. Gel, Vyacheslav Lyubchich}

\seealso{\code{\link[stats]{ar}},
         \code{\link{HVK}}, 
         \code{\link{WAVK}},
         \code{\link{wavk.test}}
}

\examples{
\dontrun{
# Fix seed for reproducible simulations:
set.seed(1)

#Simulate autoregressive time series of length n with smooth linear trend:
n <- 200
tsTrend <- 1 + 2*(1:n/n)
tsNoise <- arima.sim(n = n, list(order = c(2, 0, 0), ar = c(0.5, -0.1)))
U <- tsTrend + tsNoise
plot.ts(U)


#Use t-test
notrend.test(U)

#Use Mann-Kendall test and Yule-Walker estimates of the AR parameters
notrend.test(U, test = "MK", ar.method = "yw")


#Use WAVK test for the H0 of no trend, with m-out-of-n selection of the local window:
notrend.test(U, test = "WAVK", factor.length = "adaptive.selection")
# Sample output:
##	Sieve-bootstrap WAVK trend test
##
##data:  U
##WAVK test statistic = 21.654, moving window = 15, p-value < 2.2e-16
##alternative hypothesis: (non-)monotonic trend.
##sample estimates:
##$AR_order
##[1] 1
##
##$AR_coefficients
##    phi_1 
##0.4041848 
}
}

\keyword{ts}
\keyword{htest}
\keyword{trend}

