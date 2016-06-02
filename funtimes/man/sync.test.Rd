\name{sync.test}
\alias{sync.test}
\title{Time series trend synchronism test}
\description{
Non-parametric test for synchronism of parametric trends in multiple time series. The method tests whether \eqn{N} observed time series exhibit the same trend of some pre-specified smooth parametric form. Current version of the code assumes a linear trend \eqn{f(\theta,t) = \theta_0 + \theta_1t}, where \eqn{\theta}'s are trend parameters, \eqn{t} is a regular sequence on the interval (0,1]. Note that an assumption of `no trend' among all \eqn{N} observed time series is a subcase of a linear trend.
}

\usage{
sync.test(formula, B=1000, Window=NULL, q=NULL, j=NULL, ar_order=NULL, BIC=TRUE, robust=TRUE)
}

\arguments{
  \item{formula}{an object of class "formula", specifying the form of the common parametric time trend to be tested. Variable 't' should be used to specify the form. See 'Details' and 'Examples' for more information.}

  \item{B}{number of bootstrap resamples.} 
  
  \item{Window}{scalar or \eqn{N}-vector with lengths of the local windows (factors). If only one value is set, the same \code{Window} is applied to each time series. \eqn{N}-vector specifies a particular window for each time series. If no \code{Window} is specified, the automatic algorithm for optimal window selection is performed as a default option (see details).}

  \item{q}{scalar from 0 to 1 to define the set of possible windows \code{T*q^j} and to automatically select an optimal window for each time series. Default is \eqn{3/4}. This argument is ignored if \code{Window} is set by user.}

  \item{j}{numeric vector to define the set of possible windows \code{T*q^j} and to automatically select an optimal window for each time series. Default is \code{c(8:11)}. This argument is ignored if \code{Window} is set by user.}
  
  \item{ar_order}{order of autoregressive filter when \code{BIC = FALSE}, or the maximal order for BIC-based filtering. Default is \code{floor(10*log10(T))}. \code{ar_order} can be a scalar or \eqn{N}-vector. If scalar, the same \code{ar_order} is applied to each time series. \eqn{N}-vector specifies a separate \code{ar_order} for each time series.}

  \item{BIC}{logical value indicates whether the order of autoregressive filter should be selected by Bayesian information criterion (BIC). If \code{TRUE} (default), models of orders 1,...,\code{ar_order} or 1,...,\code{floor(10*log10(T))} are considered, depending on whether \code{ar_order} is defined or not.}
  
  \item{robust}{logical value indicates whether to use robust estimates of autoregression coefficients using \code{\link{HVK}} function (default), or to use Yule-Walker estimates delivered by \command{ar} function.}
}

\details{

Arguments \code{Window}, \code{j} and \code{q} are used to set windows for the local regression. Current version of the function assumes two options: (1) user specifies one fixed window for each time series using the argument \code{Window} (if \code{Window} is set, \code{j} and \code{q} are ignored), and (2) user specifies a set of windows by \code{j} and \code{q} to apply this set to each time series and to select an optimal window using a heuristic \eqn{m}-out-of-\eqn{n} subsampling algorithm (Bickel and Sakov, 2008). The option of selecting windows automatically for some of the time series, while for other time series the window is fixed, is not available yet. If none of these three arguments is set, default \code{j} and \code{q} are used. Values \code{T*q^j} are mapped to the largest previous integer, then only those greater than 2 are used.
}

\value{
A list of class htest containing the following components:

  \item{method}{name of the method.}

  \item{data.name}{name of the data.}

  \item{statistic}{value of the test statistic.}
  
  \item{p.value}{\eqn{p}-value of the test.}

  \item{alternative}{alternative hypothesis.}

  \item{estimate}{list with elements \code{common_trend_estimates}, \code{ar_order_used}, \code{Window_used} and \code{all_considered_windows}. The latter is a table with bootstrap and asymptotic test results for all considered windows, i.e., without adaptive selection of the local window.}
}

\references{
Bickel, P. J. and Sakov, A. (2008) On the choice of \eqn{m} in the \eqn{m} out of \eqn{n} bootstrap and confidence bounds for extrema. \emph{Statistica Sinica} 18, 967--985.

Lyubchich, V., Gel, Y. R. and El-Shaarawi, A. (2013) On detecting non-monotonic trends in environmental time series: a fusion of local regression and bootstrap. \emph{Environmetrics} 24, 209--226.

Wang, L., Akritas, M. G. and Van Keilegom, I. (2008) An ANOVA-type nonparametric diagnostic test for heteroscedastic regression models. \emph{Journal of Nonparametric Statistics} 20(5), 365--382.

Wang, L. and Van Keilegom, I. (2007) Nonparametric test for the form of parametric regression with time series errors. \emph{Statistica Sinica} 17, 369--386.
}

\author{Yulia R. Gel, Vyacheslav Lyubchich, Xingyu Wang, Ethan Schaeffer}

\seealso{\code{\link{HVK}}, 
         \code{\link{WAVK}},
		 \code{\link{wavk.test}}.
}

\examples{

# Fix seed for reproduceable simulations.
set.seed(123)

# Simulate two autoregressive time series of length n without trend (i.e., with zero trend) 
# and apply the synchronism test.
n <- 200
y1 <- arima.sim(n=n, list(order=c(1,0,0), ar=c(0.6)))
y2 <- arima.sim(n=n, list(order=c(1,0,0), ar=c(-0.2)))
X1 <- cbind(y1, y2)

\dontrun{
sync.test(X1~t, B=1000)
}
# Sample output:
##
##	Non-parametric test for synchronism of parametric linear trends
##
##data:  X1 
##Test statistic = -0.071224, p-value = 0.742
##alternative hypothesis: trends are not synchronized.
##sample estimates:
##$common_trend_estimates
##               Estimate Std. Error    t value  Pr(>|t|)
##(Intercept)  0.02944134 0.09871156  0.2982563 0.7658203
##t           -0.05858974 0.17033482 -0.3439681 0.7312353
##
##$ar_order_used
##         y1 y2
##ar_order  1  1
##
##$Window_used
##       y1 y2
##Window  8 15
##
##$all_considered_windows
## Window   Statistic p-value Asympt. p-value
##      8 -0.09419625   0.498       0.3423957
##     11 -0.08168139   0.600       0.4103374
##     15 -0.08831680   0.590       0.3733687
##     20 -0.09337623   0.564       0.3466142


# Add a time series y3 with a different linear trend and apply the synchronism test.
y3 <- 1 + 3*((1:n)/n) + arima.sim(n=n, list(order=c(1,0,0), ar=c(-0.2)))
X2 <- cbind(y1, y3)

\dontrun{
sync.test(X2~t, B=1000)
}
# Sample output:
##
##  Non-parametric test for synchronism of parametric linear trends
##
##data:  X2 
##Test statistic = 0.39071, p-value < 2.2e-16
##alternative hypothesis: trends are not synchronized.
##sample estimates:
##$common_trend_estimates
##              Estimate Std. Error   t value     Pr(>|t|)
##(Intercept) -0.5612405 0.08850773 -6.341147 1.513739e-09
##t            1.1168966 0.15272728  7.313013 6.331373e-12
##
##$ar_order_used
##         y1 y3
##ar_order  1  1
##
##$Window_used
##       y1 y3
##Window  8 15
##
##$all_considered_windows
## Window Statistic p-value Asympt. p-value
##      8 0.3161412       0    3.508377e-05
##     11 0.3664799       0    1.614506e-06
##     15 0.4064099       0    1.042697e-07
##     20 0.4357487       0    1.176418e-08
}
\keyword{ts}
\keyword{htest}
\keyword{trend}
