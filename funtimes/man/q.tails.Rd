\name{q.tails}
\alias{q.tails}
\title{Quantile-based tails comparison}
\description{
This function compares right tails of two sample distributions using a quantile-based approach. 
}

\usage{
q.tails(x, y, q=0.99)
}

\arguments{
  \item{x,y}{vectors of the same length (preferably). Tail in \eqn{y} is compared against the tail in \eqn{x}.}
  \item{q}{a threshold defining the tail for both \eqn{x} and \eqn{y}, set as a quantile. Default is 99th percentile.}
}

\details{
Sturges' formula is used to calculate number of intervals (\eqn{k}) to split the upper \eqn{100(1-q)}\% portion of \eqn{x} and \eqn{y} (the right tail). Then, the tail is divided into equally-filled intervals with a quantile step \eqn{d=(1-q)/k}. \eqn{Pk} reports the difference between corresponding intervals' centers obtained from \eqn{x} and \eqn{y}.
}

\value{
A list with two elements:
\item{d}{the quantile step.}
\item{Pk}{vector of differences of intervals' centers.}
}

\references{
Soliman, M., Lyubchich, V., Gel, Y. R., Naser, D. and Esterby, S. (2015)
Evaluating the impact of climate change on dynamics of house insurance claims.
In V. Lakshmanan et al. (eds.) \emph{Machine Learning and Data Mining Approaches to Climate Science}, pages 175--183. Springer International Publishing Switzerland.
\url{http://www.springer.com/us/book/9783319172194#}

Soliman, M., Naser, D., Lyubchich, V., Gel, Y. R. and Esterby, S. (2014)
Evaluating the impact of climate change on dynamics of house insurance claims.
In \emph{Proceedings of the 4th International Workshop on Climate Informatics}.
\url{https://www2.image.ucar.edu/sites/default/files/35-%20InsuranceAbstractRevisedJan19.pdf#overlay-context=event/ci2014/poster35
}
}

\author{Vyacheslav Lyubchich, Yulia R. Gel}

\seealso{\code{\link{i.tails}}
}


\examples{
x <- rnorm(1000); y <- rt(1000, 5)
q.tails(x, y)
}

\keyword{ts}
