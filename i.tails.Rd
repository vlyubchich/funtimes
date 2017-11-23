\name{i.tails}
\alias{i.tails}
\title{Interval-Based Tails Comparison}
\description{
This function compares right tails of two sample distributions using an interval-based approach (IBA). 
}

\usage{
i.tails(x0, x1, d = NULL)
}

\arguments{
  \item{x0,x1}{vectors of the same length (preferably). Tail in \code{x1} is compared against the tail in \code{x0}. }
  
  \item{d}{a threshold defining the tail. The threshold is the same for both \code{x0} and \code{x1}. Default is \code{quantile(x0, probs = 0.99)}.}
}

\details{
Sturges' formula is used to calculate number of intervals (\eqn{k}) for \code{x0} \eqn{\ge} \code{d}, then interval width is derived. The tails, \code{x0} \eqn{\ge} \code{d} and \code{x1} \eqn{\ge} \code{d}, are divided into the intervals. Number of \code{x1}-values within each interval is compared with the number of \code{x0}-values within the same interval (this difference is reported as \code{Nk}).
}

\value{
A list with two elements:
\item{Nk}{vector that tells how many more \code{x1}-values compared with \code{x0}-values there are within each interval.}
\item{Ck}{vector of intervals' centers.}
}

\references{
Chu, C., Gel, Y. R., and Lyubchich, V. (2015). Climate change from an insurance perspective: a case study of Norway. In J. G. Dy et al. (Eds.) \emph{Proceedings of the 5th International Workshop on Climate Informatics: CI2015}.  September 24--25, 2015, Boulder, Colorado, USA.

Lyubchich, V. and Gel, Y. R. (2017). Can we weather proof our insurance?
\emph{Environmetrics} 28(2): e2433. DOI: \href{http://dx.doi.org/10.1002/env.2433}{10.1002/env.2433}
}

\author{Calvin Chu, Yulia R. Gel, Vyacheslav Lyubchich}

\seealso{\code{\link{q.tails}}
}


\examples{
x0 <- rnorm(1000)
x1 <- rt(1000, 5)
i.tails(x0, x1)
}

\keyword{ts}
