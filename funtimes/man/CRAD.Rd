\name{crad}
\alias{crad}
\title{CRAD Clustering Algorithm}

\description{
A density-based clustering algorithm, named CRAD, which is based on an adaptive neighbor searching function and a robust data depth as the dissimilarity measure.
}

\usage{
crad(X, Nbin, StepSize = 1, Rcov = FALSE)
}

\arguments{
\item{X}{a \eqn{n \times k} matrix where \eqn{n} is a number of observations and \eqn{k} is a number of features. }

\item{Nbin}{the number of bins in generating histogram \code{H}.}

\item{StepSize}{the size of neighbor buckets in histogram \code{H} to compare for each \eqn{h_{i},i = 1, 1-width, 1-2*width, . . . , width}. Default value is 1.}

\item{Rcov}{a robust covariance matrix estimator, Minimum Covariance Determinant (MCD), will be used when \code{Rcov = TRUE}, or the sample covariance matrix \code{cov(X)} (a default option).}
}
\value{a \eqn{n \times 1} numeric vector indicating class memembership for each observation from CRAD clustering.}

\note{Since the code is calling Python from R, the data format requirement is that \code{dimnames} of argument \code{X} are to be \code{NULL}.}


\references{
Xin Huang, Yulia R. Gel. CRAD: Clustering with Robust Autocuts and Depth. In \emph{Proc. 17th IEEE International
Conference on Data Mining (ICDM).} 2017.
}



\author{
Huang, Xin, Gel, Yulia R.
}
\examples{
\dontrun{
require(mlbench)
set.seed(43)
data <- mlbench.smiley()
pclass = as.numeric(levels(data$classes))[data$classes]
p = data$x
# to visualize clusters of p
plot(p,pch=19)

# to set the dimnames of p as NULL, according to requirement for the data format
dimnames(p) <- NULL

# perform CRAD clustering
# here we set StepSize as 1, Nbin as 147,
# and without using robust covariance matrix estimator
res_cl = crad(p, 147, 1, FALSE)

# visulize clustering result
plot(p, col=as.factor(res_cl))

}
}
