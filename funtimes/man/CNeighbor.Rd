\name{CNeighbor}
\alias{CNeighbor}

\title{
Neighborhood of Time Series
}

\description{
This is an auxillary function to identify which time series in \code{Bv} are \eqn{E^\theta_\delta}-neighbors of \code{Bu}, based on Definition 2 by Ciampi et al. (2010).
}

\usage{
CNeighbor(Bu, Bv, Alpha, Beta, Delta, Theta)
}

\arguments{
  \item{Bu}{a time series vector for which the neighborhood is investigated.}
  \item{Bv}{a time series vector (of the same length as \code{Bu}) or a matrix (time series in columns) containing potential neighbors.}
  \item{Alpha}{lower limit of the time series domain.}
  \item{Beta}{upper limit of the time series domain.}
  \item{Delta}{closeness parameter, a real value in [0,1].}
  \item{Theta}{connectivity parameter, a real value in [0,1].}
}

\value{
A vector of logical values indicating which time series in \code{Bv} are \eqn{E^\theta_\delta}-neighbors of \code{Bu}.}

\references{
Ciampi, A., Appice, A., and Malerba, D. (2010) Discovering trend-based clusters in spatially distributed data streams. In \emph{International Workshop of Mining Ubiquitous and Social Environments}, pages 107--122.
}

\author{
Vyacheslav Lyubchich
}

\seealso{\code{\link{CHomogeneity}}, 
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
CNeighbor(Bu, Bv, Alpha, Beta, Delta=0.5, Theta=0.8)
}

\keyword{ts}
\keyword{trend}
