#' Neighborhood of Time Series
#' 
#' This is an auxiliary function to identify which time series in \code{Bv} 
#' are \eqn{E^\theta_\delta}-neighbors of \code{Bu}, based on Definition 2 
#' by \insertCite{Ciampi_etal_2010;textual}{funtimes}.
#' 
#' 
#' @param Bu a time series vector for which the neighborhood is investigated.
#' @param Bv a time series vector (of the same length as \code{Bu}) or a matrix 
#' (time series in columns) containing potential neighbors.
#' @param Alpha lower limit of the time series domain.
#' @param Beta upper limit of the time series domain.
#' @param Delta closeness parameter, a real value in \eqn{[0,1]}.
#' @param Theta connectivity parameter, a real value in \eqn{[0,1]}.
#' 
#' 
#' @return A vector of logical values indicating which time series in \code{Bv} 
#' are \eqn{E^\theta_\delta}-neighbors of \code{Bu}.
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
#' Bu <- rnorm(10)
#' Bv <- rnorm(10)
#' Alpha <- min(c(Bu, Bv))
#' Beta <- max(c(Bu, Bv))
#' CNeighbor(Bu, Bv, Alpha, Beta, Delta = 0.5, Theta = 0.8)
#' 
CNeighbor <- function(Bu, Bv, Alpha, Beta, Delta, Theta){
  p <- length(Bu)
  if(is.null(dim(Bv)[2])) {Bv <- matrix(Bv, ncol = 1)}
  colSums(abs(Bu - Bv) / (Beta - Alpha) <= Delta) >= Theta * p
}
