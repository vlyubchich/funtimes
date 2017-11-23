\name{q.tails}
\alias{q.tails}
\title{Quantile-Based Tails Comparison}
\description{
This function compares right tails of two sample distributions using a quantile-based approach (QBA). 
}

\usage{
q.tails(x0, x1, q = 0.99)
}

\arguments{
  \item{x0,x1}{vectors of the same length (preferably). Tail in \code{x1} is compared against the tail in \code{x0}.}
  \item{q}{a quantile defining the right tail for both \code{x0} and \code{x1}. Values above the  thresholds \code{quantile(x0, probs = q)} and \code{quantile(x1, probs = q)} are considered as the respective right tails.}
}

\details{
Sturges' formula is used to calculate number of intervals (\eqn{k}) to split the upper \eqn{100(1 - q)}\% portion of \code{x0} and \code{x1} (the right tails). Then, each tail is divided into equally-filled intervals with a quantile step \eqn{d=(1 - q)/k}. \code{Pk} reports the difference between corresponding intervals' centers obtained from \code{x0} and \code{x1}.
}

\value{
A list with two elements:
\item{d}{the quantile step.}
\item{Pk}{vector of differences of the intervals' centers.}
}

\references{
Lyubchich, V. and Gel, Y. R. (2017). Can we weather proof our insurance?
\emph{Environmetrics} 28(2): e2433. DOI: \href{http://dx.doi.org/10.1002/env.2433}{10.1002/env.2433}

Soliman, M., Lyubchich, V., Gel, Y. R., Naser, D., and Esterby, S. (2015).
Evaluating the impact of climate change on dynamics of house insurance claims.
Ch. 16 in V. Lakshmanan et al. (Eds.) \emph{Machine Learning and Data Mining Approaches to Climate Science}, pp. 175--183. Springer International Publishing.
DOI: \href{http://dx.doi.org/10.1007/978-3-319-17220-0_16}{10.1007/978-3-319-17220-0_16}

Soliman, M., Naser, D., Lyubchich, V., Gel, Y. R., and Esterby, S. (2014).
Evaluating the impact of climate change on dynamics of house insurance claims.
In \emph{Proceedings of the 4th International Workshop on Climate Informatics: CI2014}. September 25--26, 2014, Boulder, Colorado, USA.
}

\author{Vyacheslav Lyubchich, Yulia R. Gel}

\seealso{\code{\link{i.tails}}
}


\examples{
x0 <- rnorm(1000)
x1 <- rt(1000, 5)
q.tails(x0, x1)
}

\keyword{ts}
