#' Time Series Trend Synchronism Test
#'
#' @name sync.test-defunct
#'
#' @seealso \code{\link{funtimes-defunct}}
#'
#' @keywords internal
NULL

#' @rdname funtimes-defunct
#' @section \code{sync.test}:
#' For \code{sync.test}, use \code{\link{sync_test}}.
#'
#' @author Yulia R. Gel, Vyacheslav Lyubchich, Ethan Schaeffer, Xingyu Wang
#'
#' @export
#'
sync.test <- function(...)
{
    .Defunct("sync_test", msg = "sync.test is defunct (removed). Use sync_test instead.")
}
