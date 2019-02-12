#' Time Series Cluster Homogeneity
#' 
#' This is an auxiliary function to check homogeneity of time series cluster, based on 
#' Definition 4 by \insertCite{Ciampi_etal_2010;textual}{funtimes}.
#' 
#' 
#' @param Bu bucket of time series already included in the cluster.
#' @param Bv bucket of time series (neighbors) for potential inclusion in the cluster.
#' @inheritParams CNeighbor
#' 
#' 
#' @return A logical value indicating whether time series in \code{Bu} and \code{Bv} 
#' form a homogeneous cluster.
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
#' Bu <- rnorm(10)
#' Bv <- rnorm(10)
#' Alpha <- min(c(Bu, Bv))
#' Beta <- max(c(Bu, Bv))
#' CHomogeneity(Bu, Bv, Alpha, Beta, Delta = 0.5)
#' 
CHomogeneity <- function(Bu, Bv, Alpha, Beta, Delta)
{
    M <- cbind(Bu, Bv)
    med <- apply(M, 2, median)
    MED <- matrix(med, length(med), length(med))
    tmp <- abs(MED - t(MED))/(Beta - Alpha) <= Delta
    all(tmp)
}
