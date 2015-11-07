\name{i.tails}
\alias{i.tails}
\title{Interval-based tails comparison}
\description{
This function compares tails of two sample distributions using an interval-based approach. 
}

\usage{
i.tails(x, y, d=NULL)
}

\arguments{
  \item{x,y}{vectors of the same length (preferably). Tail in \eqn{y} is compared against the tail in \eqn{x}. \eqn{x \ge d} is used to obtain interval width.}
  \item{d}{a threshold defining the tail. The threshold is the same for both \eqn{x} and \eqn{y}. Default is \code{quantile(x, probs=0.99)}.}
}

\details{
Sturges' formula is used to calculate number of intervals and interval width for \eqn{x \ge d}. The same interval width is applied to \eqn{y \ge d}. The tails, \eqn{x \ge d} and \eqn{y \ge d}, are divided into the intervals, and number of \eqn{y}-values is compared with the number of \eqn{x}-values within each interval. 
}

\value{
A list with two elements:
\item{Nk}{vector that tells how many more \eqn{y}-values compared with \eqn{x}-values there are within each interval.}
\item{Ck}{vector of intervals' centers.}
}

\author{Calvin Chu, Yulia R. Gel, Vyacheslav Lyubchich}

\seealso{\code{\link{q.tails}}
}


\examples{
x <- rnorm(1000); y <- rt(1000, 5)
i.tails(x, y)
}

\keyword{ts}
