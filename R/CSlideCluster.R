#' Slide-Level Time Series Clustering
#' 
#' Cluster time series at a slide level, 
#' based on Algorithm 1 of \insertCite{Ciampi_etal_2010;textual}{funtimes}.
#' 
#' 
#' @param X a matrix of time series observed within a slide (time series in columns).
#' @param Alpha lower limit of the time series domain. Default is 
#' \code{quantile(X)[2] -}\cr \code{1.5*(quantile(X)[4] - quantile(X)[2])}.
#' @param Beta upper limit of the time series domain. 
#' Default is \code{quantile(X)[2] +}\cr \code{1.5*(quantile(X)[4] - quantile(X)[2])}.
#' @param Delta closeness parameter, a real value in \eqn{[0,1]}. 
#' Default is \code{0.1*(Beta - Alpha)}.
#' @param Theta connectivity parameter, a real value in \eqn{[0,1]}. Default is 0.8.
#' 
#' 
#' @return A vector of length \code{ncol(X)} with cluster labels.
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{CSlideCluster}}, \code{\link{CWindowCluster}}, 
#' and \code{\link{BICC}}
#' 
#' @keywords cluster ts trend
#' 
#' @author Vyacheslav Lyubchich
#' 
#' @export
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(50), 10, 5)
#' CSlideCluster(X)
#' 
CSlideCluster <- function(X, Alpha = NULL, Beta = NULL, Delta = NULL, Theta = 0.8)
{
    #X is a matrix of time series, with one slide
    N <- dim(X)[2]
    if (is.null(Alpha) | is.null(Beta)) {
        Alpha <- quantile(X)[2] - 1.5*(quantile(X)[4] - quantile(X)[2])
        Beta <- quantile(X)[4] + 1.5*(quantile(X)[4] - quantile(X)[2])
    }
    if (is.null(Delta)) {Delta <- 0.1*(Beta - Alpha)}
    TSclusters <- rep(NA, N)
    currentTS <- currentCluster <- 1
    TSclusters[currentTS] <- currentCluster
    Unclassified <- is.na(TSclusters)
    while (any(Unclassified)) {
        Include <- CExpandSlideCluster(X[,currentTS], X[,Unclassified], 
                                       Alpha, Beta, Delta, Theta)
        #Assign cluster number to the time series that were just grouped
        TSclusters[Unclassified][Include] <- currentCluster
        #Next cluster number will be:
        currentCluster <- currentCluster + 1
        #Still unclassified time series
        Unclassified <- is.na(TSclusters)
        #Select the next time series to start the cluster with, and immediately
        #assign it to the new cluster
        currentTS <- which(Unclassified)[1]
        TSclusters[currentTS] <- currentCluster
        Unclassified <- is.na(TSclusters)
    }
    TSclusters
}
