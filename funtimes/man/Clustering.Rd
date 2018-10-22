\name{Clustering}
\alias{Clustering}

\title{Trend Synchronism Clustering}


\description{This function clusters  \eqn{N} time series with a common trend using sync.test funtion.}

\usage{
Clustering(formula, rate = 1, alpha = 0.05, ...)
}

\argument{
\item{formula}{An object of class "formula", to specify the type of common trend in a \eqn{T} by \eqn{N} matrix of time series (time series in colyumns). Formula can be linear, quadratic or complex. }
\item{rate}{rate of removal of time series.}

\item{alpha}{significance level for testing hypothesis of having a common trend of some specific parametric form in the selected time series. }

\item{...}{contains other input parameters for sync.test. For example, number of bootstrap (B), Window (Window), etc.}

\item{B}{number of bootstrap simulations to obtain empirical critical values.} 

\item{Window}{scalar or \eqn{N}-vector with lengths of the local windows (factors). If only one value is set, the same \code{Window} is applied to each time series. \eqn{N}-vector specifies a particular window for each time series. If no \code{Window} is specified, the automatic algorithm for optimal window selection is performed as a default option (see `Details').}


\item{q}{scalar from 0 to 1 to define the set of possible windows \code{T*q^j} and to automatically select an optimal window for each time series. Default is \eqn{3/4}. This argument is ignored if \code{Window} is set by user.}

\item{j}{numeric vector to define the set of possible windows \code{T*q^j} and to automatically select an optimal window for each time series. Default is \code{c(8:11)}. This argument is ignored if \code{Window} is set by user.}
  
\item{ar.order}{order of autoregressive filter when \code{BIC = FALSE}, or the maximal order for BIC-based filtering. Default is \code{round(10*log10(T))}. The \code{ar.order} can be a scalar or \eqn{N}-vector. If scalar, the same \code{ar.order} is applied to each time series. \eqn{N}-vector specifies a separate \code{ar.order} for each time series.}
  
\item{ar.method}{method of estimating autoregression coefficients. Default \code{"HVK"} delivers robust difference-based estimates by Hall and Van Keilegom (2003). Alternatively, options of \command{ar} function can be used, such as \code{"burg"}, \code{"ols"}, \code{"mle"}, and \code{"yw"}.}

\item{BIC}{logical value indicates whether the order of autoregressive filter should be selected by Bayesian information criterion (BIC). If \code{TRUE} (default), models of orders 0,1,...,\code{ar.order} or 0,1,...,\code{round(10*log10(T))} are considered, depending on whether \code{ar.order} is defined or not.}
}


\references{
Lyubchich, V. and Gel, Y. R. (2016). A local factor nonparametric test for trend synchronism in multiple time series. \emph{Journal of Multivariate Analysis} 150: 91--104. DOI: \href{http://dx.doi.org/10.1016/j.jmva.2016.05.004}{10.1016/j.jmva.2016.05.004}

Lyubchich, V., Gel, Y. R., and El-Shaarawi, A. (2013). On detecting non-monotonic trends in environmental time series: a fusion of local regression and bootstrap. \emph{Environmetrics} 24: 209--226. DOI: \href{http://dx.doi.org/10.1002/env.2212}{10.1002/env.2212}

}

\examples{
set.seed(1)

Y <- matrix(rnorm(5000),nrow=50,ncol=100)

#  Testing linear trend
Clustering(Y~t)
# $`Clusters`
# Lfinal
#   1 
# 100

# Testing quadratic trend with increased bootstrap simulations and different window size
Clustering(Y~poly(t,2),B=5000,Window=10)
# $`Clusters`
# Lfinal
#   1 
# 100
}

\keyword{trend}
\keyword{synchrony}
\keyword{clustering}
