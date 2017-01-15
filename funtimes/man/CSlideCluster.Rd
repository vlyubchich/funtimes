\name{CSlideCluster}
\alias{CSlideCluster}

\title{
Slide-Level Time Series Clustering
}

\description{This function clusters time series at a slide level, based on the Algorithm 1 of Ciampi et al. (2010).}

\usage{
CSlideCluster(X, Alpha = NULL, Beta = NULL, Delta = NULL, Theta = 0.8)
}

\arguments{
  \item{X}{a matrix of time series observed within a slide (time series in columns).}
  \item{Alpha}{lower limit of the time series domain.
  
  Default is \code{quantile(X)[2] - 1.5*(quantile(X)[4] - quantile(X)[2])}.}
  \item{Beta}{upper limit of the time series domain. 
  
  Default is \code{quantile(X)[2] + 1.5*(quantile(X)[4] - quantile(X)[2])}.}
  \item{Delta}{closeness parameter, a real value in [0,1]. Default is \code{0.1*(Beta - Alpha)}.}
  \item{Theta}{connectivity parameter, a real value in [0,1]. Default is 0.8.}
}

\value{A vector of length \code{dim(X)[2]} with cluster labels.}

\references{
Ciampi, A., Appice, A., and Malerba, D. (2010). Discovering trend-based clusters in spatially distributed data streams. In \emph{International Workshop of Mining Ubiquitous and Social Environments}, pages 107--122.
}

\author{
Vyacheslav Lyubchich
}

\seealso{\code{\link{CNeighbor}}, 
         \code{\link{CHomogeneity}},
         \code{\link{CExpandSlideCluster}},
         \code{\link{CExpandWindowCluster}},
         \code{\link{CWindowCluster}}
}

\examples{
set.seed(123)
X <- matrix(rnorm(50), 10, 5)
CSlideCluster(X)
}

\keyword{ts}
\keyword{trend}
