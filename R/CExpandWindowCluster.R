#' Window-Level Time Series Cluster Expansion
#' 
#' This is an auxiliary function to expand a window-level time series cluster, 
#' based on \insertCite{Ciampi_etal_2010;textual}{funtimes}.
#' 
#' 
#' @param e a vector of logical values identifying which time series among \code{Euncl} 
#' were clustered together with \code{e} over at least \code{w*Epsilon} slides 
#' within a window \insertCite{@see Definition 7 by Ciampi_etal_2010}{funtimes}. 
#' This is a seed for window-level clustering.
#' @param Euncl a square matrix identifying the binary window cluster relation 
#' for yet unclustered time series.
#' 
#' 
#' @return A vector of logical values indicating which time series in \code{Euncl} 
#' should be included in the window-level cluster with \code{e}.
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
#' e <- sample(c(TRUE, FALSE), 5, replace = TRUE)
#' Euncl <- matrix(sample(c(TRUE, FALSE), 5, replace = TRUE), 5, 5)
#' CExpandWindowCluster(e, Euncl)
#' 
CExpandWindowCluster <- function(e, Euncl)
{
    if(is.null(dim(Euncl)[2])) {Euncl <- matrix(Euncl, ncol = 1)}
    ClusterInclude <- rep(FALSE, dim(Euncl)[2])
    if(length(ClusterInclude) == 0) {return(ClusterInclude)}
    if(any(e)) { #Perform it only if we have neighbors
        ClusterInclude[e] <- TRUE
        NeighID <- which(e)
        for (i in 1:sum(e)) {
            if (all(ClusterInclude)) {return(ClusterInclude)}
            tmp <- CExpandWindowCluster(Euncl[!ClusterInclude,NeighID[i]], 
                                        Euncl[!ClusterInclude, !ClusterInclude])
            ClusterInclude[!ClusterInclude][tmp] <- TRUE
        }
    }  
    return(ClusterInclude)
}
