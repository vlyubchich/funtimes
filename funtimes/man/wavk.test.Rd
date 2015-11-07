\name{wavk.test}
\alias{wavk.test}
\title{WAVK trend test}
\description{
Non-parametric test to detect possibly non-monotonic parametric trend in a time series.
}

\usage{
wavk.test(x, factor.length=c("user.defined", "adaptive.selection"), 
    Window=round(0.1*length(x)), q=3/4, j=c(8:11), B=1000, 
    H0=c("no trend","linear"), method=c("boot", "asympt"), 
    ar.order=NULL, BIC=TRUE, robust=TRUE, out=FALSE)
}
\arguments{
  \item{x}{univariate time series. Missing values are not allowed.}

  \item{factor.length}{method to define the length of local windows (factors). Default option \code{"user.defined"} allows to set only one value of the argument \code{Window}. The option \code{"adaptive.selection"} sets \code{method = "boot"} and employs heuristic \eqn{m}-out-of-\eqn{n} subsampling algorithm (Bickel and Sakov, 2008) to select an optimal window from the set of possible windows \code{length(x)*q^j} whose values are mapped to the largest previous integer and greater than 2.}

  \item{Window}{length of the local window (factor), default is \code{round(0.1*length(x))}. This argument is ignored if \code{factor.length = "adaptive.selection"}.}

  \item{q}{scalar from 0 to 1 to define the set of possible windows when \code{factor.length =} \code{"adaptive.selection"}. Default is \eqn{3/4}. This argument is ignored if \code{factor.length =} \code{"user.defined"}.}

  \item{j}{numeric vector to define the set of possible windows when \code{factor.length =} \code{"adaptive.selection"}. Default is \code{c(8:11)}. This argument is ignored if \code{factor.length =} \code{"user.defined"}.}

  \item{B}{number of bootstrap simulations to obtain empirical critical values. Default is 1000.}

  \item{H0}{null hypothesis: \code{"no trend"} (default) for testing the absence of trend (in other words, constant trend) vs. any, possibly non-monotonic, trend; \code{"linear"} for testing the presence of parametric linear trend vs. alternative nonlinear trend.}

  \item{method}{method of obtaining critical values: from asymptotical (\code{"asympt"}) or bootstrap (\code{"boot"}) distribution. If \code{factor.length =} \code{"adaptive.selection"} the option \code{"boot"} is used.}

  \item{ar.order}{order of autoregressive filter when \code{BIC = FALSE}, or the maximal order for BIC-based filtering. Default is \code{floor(10*log10(length(x)))}.}

  \item{BIC}{logical value indicates whether the order of autoregressive filter should be selected by Bayesian information criterion (BIC). If \code{TRUE} (default), models of orders 1,...,\code{ar.order} or 1,...,\code{floor(10*log10(length(x)))} are be considered, depending on whether \code{ar.order} is defined or not.}
  
  \item{robust}{logical value indicates whether to use robust estimates of autoregression coefficients using \code{\link{HVK}} function (default), or to use Yule-Walker estimates delivered by \command{ar} function.}

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

  \item{estimate}{list, containing the estimated coefficients of linear trend (if \code{H0="linear"}); estimated AR coefficients; test results for all considered windows (if \code{factor.length =} \code{"adaptive.selection"}).}
}

\references{
Bickel, P. J. and Sakov, A. (2008) On the choice of \eqn{m} in the \eqn{m} out of \eqn{n} bootstrap and confidence bounds for extrema. \emph{Statistica Sinica} 18, 967--985.

Lyubchich, V., Gel, Y. R. and El-Shaarawi, A. (2013) On detecting non-monotonic trends in environmental time series: a fusion of local regression and bootstrap. \emph{Environmetrics} 24, 209--226.

Wang, L., Akritas, M. G. and Van Keilegom, I. (2008) An ANOVA-type nonparametric diagnostic test for heteroscedastic regression models. \emph{Journal of Nonparametric Statistics} 20(5), 365--382.

Wang, L. and Van Keilegom, I. (2007) Nonparametric test for the form of parametric regression with time series errors. \emph{Statistica Sinica} 17, 369--386.
}

\author{Yulia R. Gel, Vyacheslav Lyubchich}

\seealso{\code{\link{HVK}}, 
         \code{\link{WAVK}},
         \code{\link{sync.test}}.
}

\examples{

# Fix seed for reproduceable simulations.
set.seed(123)

# Simulate autoregressive time series of length n with linear trend 1+2*t, 
# where t is a regular sequence on the interval (0,1].
n <- 100
t <- c(1:n)/n
U <- 1+2*t + arima.sim(n=n, list(order = c(2,0,0), ar = c(-0.7, -0.1)))

# Test for linear trend with output of all results.
\dontrun{
wavk.test(U, factor.length = "adaptive.selection", H0="linear", out=TRUE, B=1000)
}
# Sample output:
##
##  Trend test by Wang, Akritas and Van Keilegom
##
##data:  U
##WAVK test statistic = 0.8562, adaptively selected window = 4, p-value = 0.356
##alternative hypothesis: presence of a nonlinear trend
##sample estimates:
##$linear_trend_coefficients
##(Intercept)           t 
##  0.9917251   2.0224272 
##
##$AR_coefficients
##     phi_1      phi_2 
##-0.6814546 -0.2404422 
##
##$all_considered_windows
## Window WAVK-statistic p-value
##      4      0.8561654   0.356
##      5      0.8620023   0.320
##      7      0.8691870   0.288
##     10      0.6837790   0.306


# Test H0 of absence of a trend using asymptotic distribution of statistic.
wavk.test(U, method="asympt")

# Sample output:
##
##        Trend test by Wang, Akritas and Van Keilegom
##
##data:  U
##WAVK test statistic = 18.4712, user-defined window = 10, p-value < 2.2e-16
##alternative hypothesis: presence of a trend
 
}
\keyword{ts}
\keyword{htest}
\keyword{trend}

