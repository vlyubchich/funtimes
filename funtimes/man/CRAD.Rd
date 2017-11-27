\name{CRAD}
\alias{CRAD}
\title{CRAD Clustering Algorithm}

\description{
A density-based clustering algorithm, named CRAD, which is based on an adaptive neighbor searching function and a robust data depth as the dissimilarity measure.
}

\usage{
crad(X, StepSize, Nbin, Rcov = FALSE)
}

\arguments{
\item{X}{a n \eqn{\times} k matrix where n is a number of observations and k is a number of features. }

\item{Nbin}{the number of bins in generating histogram H.}

\item{StepSize}{the size of neighbor buckets in histogram H to compare for each \eqn{h_{i},i = 1, 1-width, 1-2*width, . . . , width}. Default value is 1.}

\item{Rcov}{a robust covariance matrix estimator, Minimum Covariance Determinant (MCD), will be used. If False (a default option), a sample covariance matrix \code{cov(X)} will used as the estimator.}
}
\value{a n \eqn{\times}1 numeric vector indicating class memembership for each observation from CRAD clustering.}

\note{Since the code is calling Python from R, the data format requirement is that \code{dimnames} of argument \code{X} are to be NULL.}


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
