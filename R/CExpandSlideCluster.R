#' Slide-Level Time Series Cluster Expansion
#' 
#' This is an auxiliary function to expand a slide-level time series cluster, 
#' based on \insertCite{Ciampi_etal_2010;textual}{funtimes}.
#' 
#' 
#' @param u a time series vector -- a seed to expand the cluster.
#' @param Xuncl a time series vector (of the same length as \code{u}) or a matrix 
#' (time series in columns) containing unclustered time series.
#' @inheritParams CNeighbor
#' 
#' 
#' @return A vector of logical values indicating which time series in \code{Xuncl} 
#' should be included in the slide-level cluster with \code{u}.
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{CSlideCluster}}, \code{\link{CWindowCluster}}, 
#' and \code{\link{BICC}}
#' 
#' @keywords cluster, ts, trend
#' 
#' @author Vyacheslav Lyubchich
#' 
#' @export
#' @examples
#' set.seed(123)
#' u <- rnorm(10)
#' Xuncl <- matrix(rt(50, 5), 10, 5)
#' Alpha <- min(cbind(u, Xuncl))
#' Beta <- max(cbind(u, Xuncl))
#' CExpandSlideCluster(u, Xuncl, Alpha, Beta, Delta = 0.15, Theta = 0.8)
#' 
CExpandSlideCluster <- function(u, Xuncl, Alpha, Beta, Delta, Theta)
{
    if(is.null(dim(Xuncl)[2])) {Xuncl <- matrix(Xuncl, ncol = 1)}
    ClusterInclude <- rep(FALSE, dim(Xuncl)[2])
    if(length(ClusterInclude) == 0) {return(ClusterInclude)}
    #Define neighbors of the time series among the unclussified unincluded time series Xuncl
    neighbors <- CNeighbor(u, Xuncl, Alpha, Beta, Delta, Theta)
    if (sum(neighbors) > 0) { #Proceed with all these things only if some neighbors were found
        neighborsID <- which(neighbors)
        #Check the homegeneity of the time series with its all neighbors
        if (CHomogeneity(u, Xuncl[,neighbors], Alpha, Beta, Delta)){ #The seed time series and all its neighbors are homogeneous
            #Add all neighbors to the cluster
            ClusterInclude[neighbors] <- TRUE
            #If all previously unclassified values were included in the cluster, return the result (FINISH).
            if (all(ClusterInclude)) {
                return(ClusterInclude)
            } else {
                #Otherwise expand cluster using each of the neighbors
                for (i in 1:sum(neighbors)){
                    tmp <- CExpandSlideCluster(Xuncl[,neighborsID[i]], 
                                               Xuncl[,!ClusterInclude], 
                                               Alpha, Beta, Delta, Theta)
                    ClusterInclude[!ClusterInclude][tmp] <- TRUE
                }
            }
        } else { #The cluster is not homogeneous
            #Check each neighbor one by one
            for (i in 1:sum(neighbors)){
                if (CHomogeneity(u, Xuncl[,neighborsID[i]], Alpha, Beta, Delta)){
                    ClusterInclude[neighborsID[i]] <- TRUE
                    tmp <- CExpandSlideCluster(Xuncl[,neighborsID[i]], 
                                               Xuncl[,!ClusterInclude], 
                                               Alpha, Beta, Delta, Theta)
                    ClusterInclude[!ClusterInclude][tmp] <- TRUE
                }
            }  
        }
    } 
    return(ClusterInclude)
}