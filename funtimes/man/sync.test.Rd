\name{sync.test}
\alias{sync.test}
\title{Time Series Trend Synchronism Test}
\description{
Non-parametric test for synchronism of parametric trends in multiple time series. The method tests whether \eqn{N} observed time series exhibit the same trend of some pre-specified smooth parametric form. 
}

\usage{
sync.test(formula, B = 1000, Window = NULL, q = NULL, j = NULL, 
          ar.order = NULL, ar.method = "HVK", BIC = TRUE)
}

\arguments{
  \item{formula}{an object of class "formula", specifying the form of the common parametric time trend to be tested in a \eqn{T\times N} matrix of time series (time series in columns). Variable \eqn{t} should be used to specify the form of the trend, where \eqn{t} is specified within the function as a regular sequence on the interval (0,1]. See `Examples'.}

  \item{B}{number of bootstrap simulations to obtain empirical critical values.} 
  
  \item{Window}{scalar or \eqn{N}-vector with lengths of the local windows (factors). If only one value is set, the same \code{Window} is applied to each time series. \eqn{N}-vector specifies a particular window for each time series. If no \code{Window} is specified, the automatic algorithm for optimal window selection is performed as a default option (see `Details').}

  \item{q}{scalar from 0 to 1 to define the set of possible windows \code{T*q^j} and to automatically select an optimal window for each time series. Default is \eqn{3/4}. This argument is ignored if \code{Window} is set by user.}

  \item{j}{numeric vector to define the set of possible windows \code{T*q^j} and to automatically select an optimal window for each time series. Default is \code{c(8:11)}. This argument is ignored if \code{Window} is set by user.}
  
  \item{ar.order}{order of autoregressive filter when \code{BIC = FALSE}, or the maximal order for BIC-based filtering. Default is \code{round(10*log10(T))}. The \code{ar.order} can be a scalar or \eqn{N}-vector. If scalar, the same \code{ar.order} is applied to each time series. \eqn{N}-vector specifies a separate \code{ar.order} for each time series.}
  
    \item{ar.method}{method of estimating autoregression coefficients. Default \code{"HVK"} delivers robust difference-based estimates by Hall and Van Keilegom (2003). Alternatively, options of \command{ar} function can be used, such as \code{"burg"}, \code{"ols"}, \code{"mle"}, and \code{"yw"}.}

  \item{BIC}{logical value indicates whether the order of autoregressive filter should be selected by Bayesian information criterion (BIC). If \code{TRUE} (default), models of orders 0,1,...,\code{ar.order} or 0,1,...,\code{round(10*log10(T))} are considered, depending on whether \code{ar.order} is defined or not.}
}

\details{

Arguments \code{Window}, \code{j}, and \code{q} are used to set windows for the local regression. Current version of the function assumes two options: (1) user specifies one fixed window for each time series using the argument \code{Window} (if \code{Window} is set, \code{j} and \code{q} are ignored), and (2) user specifies a set of windows by \code{j} and \code{q} to apply this set to each time series and to select an optimal window using a heuristic \eqn{m}-out-of-\eqn{n} subsampling algorithm (Bickel and Sakov, 2008). The option of selecting windows automatically for some of the time series, while for other time series the window is fixed, is not available yet. If none of these three arguments is set, default \code{j} and \code{q} are used. Values \code{T*q^j} are mapped to the largest previous integer, then only those greater than 2 are used.
}

\value{
A list of class htest containing the following components:

  \item{method}{name of the method.}

  \item{data.name}{name of the data.}

  \item{statistic}{value of the test statistic.}
  
  \item{p.value}{\eqn{p}-value of the test.}

  \item{alternative}{alternative hypothesis.}

  \item{estimate}{list with elements \code{common_trend_estimates}, \code{ar_order_used}, \code{Window_used}, and \code{all_considered_windows}. The latter is a table with bootstrap and asymptotic test results for all considered windows, i.e., without adaptive selection of the local window.}
}

\references{
Bickel, P. J. and Sakov, A. (2008). On the choice of \eqn{m} in the \eqn{m} out of \eqn{n} bootstrap and confidence bounds for extrema. \emph{Statistica Sinica} 18: 967--985.

Lyubchich, V. and Gel, Y. R. (2016). A local factor nonparametric test for trend synchronism in multiple time series. \emph{Journal of Multivariate Analysis} 150: 91--104. DOI: \href{http://dx.doi.org/10.1016/j.jmva.2016.05.004}{10.1016/j.jmva.2016.05.004}

Lyubchich, V., Gel, Y. R., and El-Shaarawi, A. (2013). On detecting non-monotonic trends in environmental time series: a fusion of local regression and bootstrap. \emph{Environmetrics} 24: 209--226. DOI: \href{http://dx.doi.org/10.1002/env.2212}{10.1002/env.2212}

Wang, L., Akritas, M. G., and Van Keilegom, I. (2008). An ANOVA-type nonparametric diagnostic test for heteroscedastic regression models. \emph{Journal of Nonparametric Statistics} 20(5): 365--382. DOI: \href{http://dx.doi.org/10.1080/10485250802066112}{10.1080/10485250802066112}

Wang, L. and Van Keilegom, I. (2007). Nonparametric test for the form of parametric regression with time series errors. \emph{Statistica Sinica} 17: 369--386.
}

\author{Yulia R. Gel, Vyacheslav Lyubchich, Ethan Schaeffer, Xingyu Wang}

\seealso{\code{\link[stats]{ar}},
        \code{\link{HVK}}, 
        \code{\link{WAVK}},
		    \code{\link{wavk.test}}
}

\examples{

# Fix seed for reproduceable simulations:
set.seed(1)

# Simulate two autoregressive time series of length n without trend 
#(i.e., with zero or constant trend) 
# and arrange the series into a matrix:
n <- 200
y1 <- arima.sim(n=n, list(order=c(1,0,0), ar=c(0.6)))
y2 <- arima.sim(n=n, list(order=c(1,0,0), ar=c(-0.2)))
Y <- cbind(y1, y2)
plot.ts(Y)


#Test H0 of a common linear trend:
\dontrun{
sync.test(Y ~ t, B=500)}
# Sample output:
##	Non-parametric test for synchronism of parametric trends
##
##data:  Y 
##Test statistic = -0.0028999, p-value = 0.7
##alternative hypothesis: common trend is not of the form Y ~ t.
##sample estimates:
##$common_trend_estimates
##               Estimate Std. Error    t value  Pr(>|t|)
##(Intercept) -0.02472566  0.1014069 -0.2438261 0.8076179
##t            0.04920529  0.1749859  0.2811958 0.7788539
##
##$ar.order_used
##         y1 y2
##ar.order  1  1
##
##$Window_used
##       y1 y2
##Window 15  8
##
##$all_considered_windows
## Window    Statistic p-value Asympt. p-value
##      8 -0.000384583   0.728       0.9967082
##     11 -0.024994408   0.860       0.7886005
##     15 -0.047030164   0.976       0.6138976
##     20 -0.015078579   0.668       0.8714980


# Add a time series y3 with a different linear trend and re-apply the test:
y3 <- 1 + 3*((1:n)/n) + arima.sim(n=n, list(order=c(1,0,0), ar=c(-0.2)))
Y2 <- cbind(Y, y3)
plot.ts(Y2)
\dontrun{
sync.test(Y2 ~ t, B=500)}
# Sample output:
##	Non-parametric test for synchronism of parametric trends
##
##data:  Y2 
##Test statistic = 0.48579, p-value < 2.2e-16
##alternative hypothesis: common trend is not of the form Y2 ~ t.
##sample estimates:
##$common_trend_estimates
##              Estimate Std. Error  t value     Pr(>|t|)
##(Intercept) -0.3632963 0.07932649 -4.57976 8.219360e-06
##t            0.7229777 0.13688429  5.28167 3.356552e-07
##
##$ar.order_used
##         Y.y1 Y.y2 y3
##ar.order    1    1  0
##
##$Window_used
##       Y.y1 Y.y2 y3
##Window    8   11  8
##
##$all_considered_windows
## Window Statistic p-value Asympt. p-value
##      8 0.4930069       0    1.207378e-05
##     11 0.5637067       0    5.620248e-07
##     15 0.6369703       0    1.566057e-08
##     20 0.7431621       0    4.201484e-11


#Other hypothesized trend forms can be specified, for example:
\dontrun{
sync.test(Y ~ 1) #constant trend
sync.test(Y ~ poly(t, 2)) #quadratic trend
sync.test(Y ~ poly(t, 3)) #cubic trend}
}

\keyword{ts}
\keyword{htest}
\keyword{trend}
