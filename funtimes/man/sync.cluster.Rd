\name{sync.cluster}
\alias{sync.cluster}

\title{Time Series Clustering based on Trend Synchronism}

\description{Cluster time series with a common trend using the \code{sync.test} function.}

\usage{
sync.cluster(formula, rate = 1, alpha = 0.05, ...)
}

\arguments{
\item{formula}{an object of class "formula", specifying the type of common trend for clustering the time series in a \eqn{T} by \eqn{N} matrix of time series (time series in columns). Variable \eqn{t} should be used to specify the form of the trend, where \eqn{t} is specified within the function as a regular sequence on the interval (0,1]. See `Examples'.}

\item{rate}{rate of removal of time series. Default is 1 (i.e., one time series is removed at a time to re-test the remaining time series). Values above 1 are treated as number of time series to be removed. Values below 1 are treated as a fraction of time series to be removed.}

\item{alpha}{significance level for testing hypothesis of a common trend of the parametric form specified in \code{formula}. }

\item{...}{arguments passed to \code{sync.test}, for example, number of bootstrap replications (\code{B}).}
}

\details{
SL: In 1-2 sentences, describe what the function does.

}

\value{
SL: What is the output?
}


\references{
Ghahari, A., Gel, Y. R., Lyubchich, V., Chun, Y., and Uribe, D. (2017). On employing multi-resolution weather data in crop insurance. In \emph{Proceeding of the Workshop on Mining Big Data in Climate and Environment (MBDCE 2017), 17th SIAM International Conference on Data Mining (SDM 2017)}. April 27--29, 2017, Houston, Texas, USA.

Lyubchich, V. and Gel, Y. R. (2016). A local factor nonparametric test for trend synchronism in multiple time series. \emph{Journal of Multivariate Analysis} 150: 91--104. DOI: \href{http://dx.doi.org/10.1016/j.jmva.2016.05.004}{10.1016/j.jmva.2016.05.004}

}

\examples{
SL: make examples more meaningful, for example, show/simulate a case with several clusters and detect them; At the same time, demonstrate a couple of different trend shapes.

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
\keyword{cluster}
