#' Interval-Based Tails Comparison
#'
#' @name i.tails-defunct
#'
#' @seealso \code{\link{funtimes-defunct}}
#'
#' @keywords internal
NULL

#' @rdname funtimes-defunct
#' @section \code{i.tails}:
#' For \code{i.tails}, use \code{\link{tails_i}}.
#'
#' @author Calvin Chu, Yulia R. Gel, Vyacheslav Lyubchich
#'
#' @export
#'
i.tails <- function(...)
{
    .Defunct("tails_i", msg = "i.tails is defunct (removed). Use tails_i instead.")
}
