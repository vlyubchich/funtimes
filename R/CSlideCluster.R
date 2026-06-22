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
    N <- ncol(X)
    
    if (is.null(Alpha) || is.null(Beta)) {
        q <- quantile(X)
        IQR_X <- q[4] - q[2]
        Alpha <- q[2] - 1.5 * IQR_X
        Beta <- q[4] + 1.5 * IQR_X
    }
    
    if (is.null(Delta)) {
        Delta <- 0.1 * (Beta - Alpha)
    }
    
    ts_clusters <- rep(NA, N)
    cluster_label <- 1
    
    unclassified_indices <- 1:N
    
    while (length(unclassified_indices) > 0) {
        # Start a new cluster with the first unclassified time series
        seed_ts_idx <- unclassified_indices[1]
        ts_clusters[seed_ts_idx] <- cluster_label
        
        # Find other unclassified series to add to this cluster
        unclassified_subset_indices <- unclassified_indices[-1]
        
        if (length(unclassified_subset_indices) > 0) {
            series_to_include <- CExpandSlideCluster(X[, seed_ts_idx], 
                                                     X[, unclassified_subset_indices, drop = FALSE],
                                                     Alpha, Beta, Delta, Theta)
            
            # Assign the current cluster label to the included series
            ts_clusters[unclassified_subset_indices[series_to_include]] <- cluster_label
        }
        
        # Prepare for the next iteration
        cluster_label <- cluster_label + 1
        unclassified_indices <- which(is.na(ts_clusters))
    }
    
    ts_clusters
}
