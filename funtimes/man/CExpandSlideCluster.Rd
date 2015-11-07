\name{CExpandSlideCluster}
\alias{CExpandSlideCluster}

\title{
Slide-level time series cluster expansion
}

\description{
This is an auxillary function to expand a slide-level time series cluster, based on Ciampi et al. (2010).}

\usage{
CExpandSlideCluster(u, Xuncl, Alpha, Beta, Delta, Theta)
}

\arguments{
  \item{u}{a time series vector --- a seed to expand the cluster.}
  \item{Xuncl}{a time series vector (of the same length as \code{u}) or a matrix (time series in columns) containing unclustered time series.}
  \item{Alpha}{lower limit of the time series domain.}
  \item{Beta}{upper limit of the time series domain.}
  \item{Delta}{closeness parameter, a real value in [0,1].}
  \item{Theta}{connectivity parameter, a real value in [0,1].}
}

\value{A vector of logical values indicating which time series in \code{Xuncl} should be included in the slide-level cluster with \code{u}.}

\references{
Ciampi, A., Appice, A. and Malerba, D. (2010) Discovering trend-based clusters in spatially distributed data streams. In \emph{International Workshop of Mining Ubiquitous and Social Environments}, pages 107--122.
}

\author{
Vyacheslav Lyubchich
}

\seealso{\code{\link{CNeighbor}}, 
         \code{\link{CHomogeneity}},
         \code{\link{CSlideCluster}},
         \code{\link{CExpandWindowCluster}},
         \code{\link{CWindowCluster}}.
}

\examples{
set.seed(123)
u <- rnorm(10)
Xuncl <- matrix(rt(50, 5), 10, 5)
Alpha <- min(cbind(u,Xuncl))
Beta <- max(cbind(u,Xuncl))
CExpandSlideCluster(u, Xuncl, Alpha, Beta, Delta=0.15, Theta=0.8)
}

\keyword{ts}
\keyword{trend}

