#' Quantile-Based Tails Comparison
#'
#' @name q.tails-defunct
#'
#' @seealso \code{\link{funtimes-defunct}}
#'
#' @keywords internal
NULL

#' @rdname funtimes-defunct
#' @section \code{q.tails}:
#' For \code{q.tails}, use \code{\link{tails_q}}.
#'
#' @author Vyacheslav Lyubchich, Yulia R. Gel
#'
#' @export
#'
q.tails <- function(...)
{
    .Defunct("tails_q", msg = "q.tails is defunct (removed). Use tails_q instead.")
}
