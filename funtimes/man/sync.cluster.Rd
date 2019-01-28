\name{sync.cluster}
\alias{sync.cluster}

\title{Time Series Clustering based on Trend Synchronism}

\description{Cluster time series with a common trend using the \code{sync.test} function.}

\usage{
sync.cluster(formula, rate = 1, alpha = 0.05, ...)
}

\arguments{
\item{formula}{an object of class "formula", specifying the type of common trend for clustering the time series in a \eqn{T} by \eqn{N} matrix of time series (time series in columns). Variable \eqn{t} should be used to specify the form of the trend, where \eqn{t} is specified within the function as a regular sequence on the interval (0,1]. See `Examples'.}

\item{rate}{rate of removal of time series. Default is 1 (i.e., one time series is removed at a time to re-test the remaining time series). Values above 1 are treated as number of time series to be removed. Values from 0 to 1 are treated as a fraction of time series to be removed.}

\item{alpha}{significance level for testing hypothesis of a common trend of the parametric form specified in \code{formula}. }

\item{...}{arguments passed to \code{sync.test}, for example, number of bootstrap replications (\code{B}).}
}

\details{
%SL: In 1-2 sentences, describe what the function does.
The \code{sync.cluster} function recursively clusters the time series having a pre-specified common parametric trend until there are no time series in the matrix to test. With a given \eqn{N} time series, the \code{sync.test} function is used to check the common trend. In the presence of a common trend, time series are grouped together. Othewise, the time series with largest contribution to the test statistics is removed based on the rate of removal. The largest contribution to the test statistics is assessed by WAVK test statistics. Output of the function is the number of clusters obtained.
}

\value{
%SL: What is the output?
A list with the elements:
\item{Clusters}{Number of Clusters Obtained}
\item{Column.Index}{Index of columns of time series (based on the main matrix) with each cluster}
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
%SL: make examples more meaningful, for example, show/simulate a case with several clusters and detect them; At the same time, demonstrate a couple of different trend shapes.

#######  Testing linear trend
set.seed(1)
n=60
nc = 30
Y = matrix(NA,nrow = n, ncol = nc)

# Matrix of 20 autoregressive time series of length n
for ( i in 1:nc){
   
   Y[,i]<-  arima.sim(n = n, list(order = c(1, 0, 0), ar = c(0.6)))
    
}
# Rate of removal is default
LinTrend <- sync.cluster(Y~t)

## Sample output:
## $`Clusters`
## Lfinal
## 0  1  2  3  4 
## 2 21  2  3  2 

## Zero denoting that one time series is not tied to any other series
## 1, 2, 3 and 4 denote the number of clusters obtained. 
## The second row shows the number of time series in each cluster.


##### Testing Quadaratic Trend
## Rate of removal is default
QuadTrend <- sync.cluster(Y~poly(t,2))


## $`Clusters`
## Lfinal
## 0  1  2  3 
## 1 20  6  2 


# Rate of removal is 2
set.seed(1)
n=60
nc = 50
Y1= matrix(NA,nrow = n, ncol = nc)

# Matrix of 20 autoregressive time series of length n
for ( i in 1:nc){
    
    Y1[,i]<-  arima.sim(n = n, list(order = c(1, 0, 0), ar = c(0.6)))
    
}
LinTrendR2 <- sync.cluster(Y1~t, rate = 2)

## $`Clusters`
## Lfinal
## 1  2  3  4  5  6  7 
## 30  6  2  2  2  2  2 

QuadTrendR2 <- sync.cluster(Y1~poly(t,2), rate = 2)

## $`Clusters`
## Lfinal
## 1  2  3  4  5 
## 28  4  8  4  2
}

\keyword{trend}
\keyword{synchrony}
\keyword{cluster}
