\name{CExpandWindowCluster}
\alias{CExpandWindowCluster}

\title{
Window-Level Time Series Cluster Expansion
}

\description{
This is an auxillary function to expand a window-level time series cluster, based on Ciampi et al. (2010).}

\usage{
CExpandWindowCluster(e, Euncl)
}

\arguments{
  \item{e}{a vector of logical values identifying which time series among \code{Euncl} were clustered together with \code{e} over at least \code{w*Epsilon} slides within a window (see Definition 7 by Ciampi et al., 2010). This is a seed for window-level clustering.}
  \item{Euncl}{a square matrix identifying the binary window cluster relation for yet unclustered time series.}
}

\value{A vector of logical values indicating which time series in \code{Euncl} should be included in the window-level cluster with \code{e}.}

\references{
Ciampi, A., Appice, A., and Malerba, D. (2010) Discovering trend-based clusters in spatially distributed data streams. In \emph{International Workshop of Mining Ubiquitous and Social Environments}, pages 107--122.
}

\author{
Vyacheslav Lyubchich
}

\seealso{\code{\link{CNeighbor}}, 
         \code{\link{CHomogeneity}},
         \code{\link{CSlideCluster}},
         \code{\link{CExpandSlideCluster}},
         \code{\link{CWindowCluster}}
}

\examples{
set.seed(123)
e <- sample(c(TRUE, FALSE), 5, replace=TRUE)
Euncl <- matrix(sample(c(TRUE, FALSE), 5, replace=TRUE), 5, 5)
CExpandWindowCluster(e, Euncl)
}

\keyword{ts}
\keyword{trend}
