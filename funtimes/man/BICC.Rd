\name{BICC}
\alias{BICC}

\title{BIC-Based Spatio-Temporal Clustering}

\description{
Apply the algorithm of unsupervised spatio-temporal clustering TRUST (Ciampi et al., 2010), with automatic selection of its tuning parameters \code{Delta} and \code{Epsilon} based on Bayesian information criterion (BIC).
}

\usage{
BICC(X, Alpha=NULL, Beta=NULL, Theta=0.8, p, w, s)
}

\arguments{
  \item{X}{a matrix of time series to be clustered (time series in columns).}
  \item{Alpha}{lower limit of the time series domain, passed to \code{\link{CSlideCluster}}.}
  \item{Beta}{upper limit of the time series domain, passed to \code{\link{CSlideCluster}}.}
  \item{Theta}{connectivity parameter, passed to \code{\link{CSlideCluster}}.}
  \item{p}{number of layers (time series observations) in each slide.}
  \item{w}{number of slides in each window.}
  \item{s}{step to shift a window, calculated in number of slides. Recommended values are 1 (overlapping windows) or equal to \code{w} (non-overlapping windows).}
}

\details{
This is the upper-level function for time series clustering. It exploits the functions \code{\link{CWindowCluster}} and \code{\link{CSlideCluster}} to cluster time series based on closeness and homogeneity measures. Clustering is performed multiple times with a range of equidistant values for the parameters \code{Delta} and \code{Epsilon}, then optimal parameters \code{Delta} and \code{Epsilon} along with the corresponding clustering results are shown (see Schaeffer et al., 2016, for more details).

The total length of time series (number of levels, i.e., \code{nrow(X)}) should be divisible by \code{p}.
}

\value{A list with the following elements:
\item{Delta.optimal}{`optimal' value for the clustering parameter \code{Delta}.}
\item{Epsilon.optimal}{`optimal' value for the clustering parameter \code{Epsilon}.}
\item{Clusters}{vector of length \code{ncol(X)} with cluster labels.}
\item{IC}{values of the information criterion (BIC) for each considered combination of \code{Delta} (rows) and \code{Epsilon} (columns).}
\item{Delta.all}{vector of considered values for \code{Delta}.}
\item{Epsilon.all}{vector of considered values for \code{Epsilon}.}
}

\references{
Ciampi, A., Appice, A., and Malerba, D. (2010). Discovering trend-based clusters in spatially distributed data streams. In \emph{International Workshop of Mining Ubiquitous and Social Environments}, pages 107--122.

Schaeffer, E. D., Testa, J. M., Gel, Y. R., and Lyubchich, V. (2016). On information criteria for dynamic spatio-temporal clustering. In A. Banerjee et al. (Eds.) \emph{Proceedings of the 6th International Workshop on Climate Informatics: CI 2016}. NCAR Technical Note NCAR/TN-529+PROC, September 2016, pages 5--8. DOI: \href{http://dx.doi.org/10.5065/D6K072N6}{10.5065/D6K072N6}
}

\author{
Vyacheslav Lyubchich, Ethan Schaeffer
}

\seealso{\code{\link{CSlideCluster}},
         \code{\link{CWindowCluster}},
         \code{\link{Purity}}
}

\examples{
<<<<<<< HEAD
# Fix seed for reproducible simulations:
=======
# Fix seed for reproduceable simulations:
>>>>>>> 3840f31c57e8f0298ddf920c8770d5e9e0ba02e8
set.seed(1)

##### Example 1
# Similar to Schaeffer et al. (2016), simulate 3 years of monthly data 
#for 10 locations and apply clustering:
# 1.1 Simulation
T <- 36 #total months
N <- 10 #locations
phi <- c(0.5) #parameter of autoregression
burn <- 300 #burn-in period for simulations
X <- sapply(1:N, function(x) 
            arima.sim(n=T+burn, list(order=c(length(phi), 0, 0), ar=phi)))[(burn+1):(T+burn),]
colnames(X) <- paste("TS", c(1:dim(X)[2]), sep="")

# 1.2 Clustering
# Assume that information arrives in year-long slides or data chunks
p <- 12 #number of time layers (months) in a slide
<<<<<<< HEAD
# Let the upper level of clustering (window) be the whole period of 3 years, so
=======
# Let the upper level of clistering (window) be the whole period of 3 years, so
>>>>>>> 3840f31c57e8f0298ddf920c8770d5e9e0ba02e8
w <- 3 #number of slides in a window
s <- w #step to shift a window, but it does not matter much here as we have only one window of data
tmp <- BICC(X, p=p, w=w, s=s)

# 1.3 Evaluate clustering
# In these simulations, it is known that all time series belong to one class,
#since they were all simulated the same way:
classes <- rep(1,10)
# Use the information on the classes to calculate clustering purity:
Purity(classes, tmp$Clusters[1,])

##### Example 2
# 2.1 Modify time series and update classes accordingly:
# Add a mean shift of to a half of the time series:
X2 <- X
X2[,1:(N/2)] <- X2[,1:(N/2)] + 3
classes2 <- rep(c(1,2), each=N/2)

# 2.2 Re-apply clustering procedure and evaluate clustering purity:
tmp2 <- BICC(X2, p=p, w=w, s=s)
tmp2$Clusters
Purity(classes2, tmp2$Clusters[1,])
}

\keyword{cluster}
\keyword{ts}
\keyword{trend}
