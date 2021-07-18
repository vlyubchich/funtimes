#' Sieve Bootstrap Based Test for the Null Hypothesis of no Trend
#' 
#' @name notrend.test-defunct
#' 
#' @seealso \code{\link{funtimes-defunct}}
#' 
#' @keywords internal
NULL

#' @rdname funtimes-defunct
#' @section \code{notrend.test}:
#' For \code{notrend.test}, use \code{\link{notrend_test}}.
#' 
#' @author Vyacheslav Lyubchich, Yulia R. Gel
#' 
#' @export
#' 
notrend.test <- function(...)
{
    .Defunct("notrend_test", msg = "notrend.test is defunct (removed). Use notrend_test instead.")
}
