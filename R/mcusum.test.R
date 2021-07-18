#' Change Point Test for Regression
#' 
#' @name mcusum.test-defunct
#'
#' @seealso \code{\link{funtimes-defunct}}
#' 
#' @keywords internal
NULL

#' @rdname funtimes-defunct
#' @section \code{mcusum.test}:
#' For \code{mcusum.test}, use \code{\link{mcusum_test}}.
#' 
#' @author Vyacheslav Lyubchich
#' 
#' @export
#' 
mcusum.test <- function(...)
{
    .Defunct("mcusum_test", msg = "mcusum.test is defunct (removed). Use mcusum_test instead.")
}
