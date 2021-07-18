#' Time Series Clustering based on Trend Synchronism
#' 
#' @name sync.cluster-defunct
#' 
#' @seealso \code{\link{funtimes-defunct}}
#' 
#' @keywords internal
NULL

#' @rdname funtimes-defunct
#' @section \code{sync.cluster}:
#' For \code{sync.cluster}, use \code{\link{sync_cluster}}.
#' 
#' @author Srishti Vishwakarma, Vyacheslav Lyubchich
#' 
#' @export
#' 
sync.cluster <- function(...) 
{
    .Defunct("sync_cluster", msg = "sync.cluster is defunct (removed). Use sync_cluster instead.")
}
