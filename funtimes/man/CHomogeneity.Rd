\name{CHomogeneity}
\alias{CHomogeneity}

\title{
Time Series Cluster Homogeneity
}

\description{
This is an auxillary function to check homogeneity of time series cluster, based on Definition 4 by Ciampi et al. (2010).}

\usage{
CHomogeneity(Bu, Bv, Alpha, Beta, Delta)
}

\arguments{
  \item{Bu}{bucket of time series already included in the cluster.}
  \item{Bv}{bucket of time series (neighbors) for potential inclusion in the cluster.}
  \item{Alpha}{lower limit of the time series domain.}
  \item{Beta}{upper limit of the time series domain.}
  \item{Delta}{closeness parameter, a real value in [0,1].}
}

\value{
A logical value indicating whether time series in \code{Bu} and \code{Bv} form a homogeneous cluster.}

\references{
Ciampi, A., Appice, A., and Malerba, D. (2010). Discovering trend-based clusters in spatially distributed data streams. In \emph{International Workshop of Mining Ubiquitous and Social Environments}, pages 107--122.
}

\author{
Vyacheslav Lyubchich
}

\seealso{\code{\link{CNeighbor}}, 
         \code{\link{CExpandSlideCluster}},
         \code{\link{CSlideCluster}},
         \code{\link{CExpandWindowCluster}},
         \code{\link{CWindowCluster}}
}

\examples{
Bu <- rnorm(10)
Bv <- rnorm(10)
Alpha <- min(c(Bu,Bv))
Beta <- max(c(Bu,Bv))
CHomogeneity(Bu, Bv, Alpha, Beta, Delta=0.5)
}

\keyword{ts}
\keyword{trend}
