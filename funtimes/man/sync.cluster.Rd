\name{sync.cluster}
\alias{sync.cluster}

\title{Time Series Clustering based on Trend Synchronism}

\description{Cluster time series with a common trend using the \code{\link{sync.test}} function.}

\usage{
sync.cluster(formula, rate = 1, alpha = 0.05, ...)
}

\arguments{
\item{formula}{an object of class "formula", specifying the type of common trend for clustering the time series in a \eqn{T} by \eqn{N} matrix of time series (time series in columns). Variable \eqn{t} should be used to specify the form of the trend, where \eqn{t} is specified within the function automatically as a regular sequence of length \eqn{N} on the interval (0,1]. See `Examples'.}

\item{rate}{rate of removal of time series. Default is 1 (i.e., if hypothesis of synchronism is rejected one time series is removed at a time to re-test the remaining time series). Integer values above 1 are treated as number of time series to be removed. Values from 0 to 1 are treated as a fraction of time series to be removed.}

\item{alpha}{significance level for testing hypothesis of a common trend (using \code{\link{sync.test}}) of the parametric form specified in \code{formula}.}

\item{...}{arguments passed to \code{\link{sync.test}}, for example, number of bootstrap replications (\code{B}).}
}

\details{
The \code{sync.cluster} function recursively clusters time series having a pre-specified common parametric trend until there are no time series left. Starting with the given \eqn{N} time series, the \code{\link{sync.test}} function is used to test for a common trend. If null hypothesis of common trend is not rejected by \code{\link{sync.test}}, the time series are grouped together (i.e., assigned to a cluster). Othewise, the time series with the largest contribution to the test statistics are temporarily removed (the number of time series to remove depends on the \code{rate} of removal) and \code{\link{sync.test}} is applied again. The largest contribution to the test statistic is assessed by WAVK test statistic. Output of the function is the number of clusters obtained, assigned cluster labels, and parameters used for the trend synchronism testing.
}

\value{
A list with the elements:
\item{Clusters}{Number of clusters obtained}
\item{Column.Index}{Index of columns of time series (based on the main matrix) in each cluster}
\item{Pval}{\code{p}-value of the \code{sync.test} statistics}
\item{TestStatistics}{\code{sync.test} test statistics}
\item{Estimate}{Parametric trend estimates of a cluster}
\item{AROrder}{AR filter order}
\item{WindowUsed}{Window used for the test}
\item{allConsideredWindow}{Different windows used in the \code{sync.test}}
\item{WAVKobs}{WAVK test statistics for each cluster}
}


\references{
Ghahari, A., Gel, Y. R., Lyubchich, V., Chun, Y., and Uribe, D. (2017). On employing multi-resolution weather data in crop insurance. In \emph{Proceeding of the Workshop on Mining Big Data in Climate and Environment (MBDCE 2017), 17th SIAM International Conference on Data Mining (SDM 2017)}. April 27--29, 2017, Houston, Texas, USA.

Lyubchich, V. and Gel, Y. R. (2016). A local factor nonparametric test for trend synchronism in multiple time series. \emph{Journal of Multivariate Analysis} 150: 91--104. DOI: \href{http://dx.doi.org/10.1016/j.jmva.2016.05.004}{10.1016/j.jmva.2016.05.004}

}

\examples{
set.seed(1)
## Simulating four autoregressive time series 
## Three having linear trend and one without trend

n=50
nc = 4
Y = matrix(NA,nrow = n, ncol = nc)

# Matrix of 20 autoregressive time series of length 

for ( i in 1:nc){
    if (i < 2){
        Y[,i] <-  arima.sim(n = n, list(order = c(1, 0, 0), ar = c(0.6)))
    } else {
        Y[,i] <- -10 + 0.5 * (1:n) + arima.sim(n = n, list(order = c(1, 0, 0), ar = c(0.6)))
    }
    
}

plot.ts(Y)
# Finding clusters with default parameters
\dontrun{
LinTrend <- sync.cluster(Y~t) 

## Sample Output:
## $`Clusters`
## Lfinal
## 1 
## 3 

}

## simulating seven auroregressive time series 
## Three have linear trend added and four have no trend
n=50
nc = 7
Y = matrix(NA,nrow = n, ncol = nc)

for ( i in 1:nc){
    if (i < 5){
        Y[,i] <-  arima.sim(n = n, list(order = c(1, 0, 0), ar = c(0.6)))
    } else {
        Y[,i] <- -10 + 0.5 * (1:n) + arima.sim(n = n, list(order = c(1, 0, 0), ar = c(0.6)))
    }
    
}

plot.ts(Y)
## Clustering with rate of removal = 5 and window = 15
\dontrun{
LinTrendR5W15 <- sync.cluster(Y~t, rate = 5, Window = 15) 

## Sample output:
# $`Clusters`
# Lfinal
# 1 
# 2 
}

## simulating five autoregressive time series to test for quadratic trend
## One has no trend, while rest of the series have quadratic trend
n=30
nc = 5
Y = matrix(NA,nrow = n, ncol = nc)
p <- 0.5
q <- seq(1,n,1)

for ( i in 1:nc){
    if (i < 2){
        Y[,i] <- arima.sim(n = n, list(order = c(1, 0, 0), ar = c(0.6)))
    } else {
        Y[,i] <- -10 + p*(q+10)^2 + arima.sim(n = n, list(order = c(1, 0, 0), ar = c(0.6)))
    }
    
}
plot.ts(Y)
# Clustering with default rate of removal
\dontrun{
QuadTrend <- sync.cluster(Y~poly(t,2))

## Sample output:
## $`Clusters`
## Lfinal
## 0 1 
## 2 2 
}
}

\keyword{trend}
\keyword{synchrony}
\keyword{cluster}
