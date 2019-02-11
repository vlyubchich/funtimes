#' Interval-Based Tails Comparison
#' 
#' This function compares right tails of two sample distributions using 
#' an interval-based approach (IBA); 
#' see \insertCite{Chu_etal_2015_insurance;textual}{funtimes}
#' and \insertCite{Lyubchich_Gel_2017_insurance;textual}{funtimes}. 
#' 
#' @details Sturges' formula is used to calculate number of intervals 
#' (\eqn{k}) for \code{x0} \eqn{\ge} \code{d}, then interval width is derived. 
#' The tails, \code{x0} \eqn{\ge} \code{d} and \code{x1} \eqn{\ge} \code{d}, 
#' are divided into the intervals. Number of \code{x1}-values within each interval 
#' is compared with the number of \code{x0}-values within the same interval 
#' (this difference is reported as \code{Nk}).
#' 
#' 
#' @param x0,x1 vectors of the same length (preferably). 
#' Tail in \code{x1} is compared against the tail in \code{x0}.
#' @param d a threshold defining the tail. The threshold is the same for both 
#' \code{x0} and \code{x1}. Default is \code{quantile(x0, probs = 0.99)}.
#' 
#' 
#' @return A list with two elements:
#' \item{Nk}{vector that tells how many more \code{x1}-values compared with 
#' \code{x0}-values there are within each interval.}
#' \item{Ck}{vector of intervals' centers.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{q.tails}}
#' 
#' @keywords ts
#' 
#' @author Calvin Chu, Yulia R. Gel, Vyacheslav Lyubchich
#' 
#' @export
#' @examples
#' x0 <- rnorm(1000)
#' x1 <- rt(1000, 5)
#' i.tails(x0, x1)
#' 
i.tails <- function(x0, x1, d = NULL)
{
    if(is.null(d)) d <- quantile(x0, probs = 0.99)
    n <- sum(x0 >= d)
    #Sturges' formula
    k <- ceiling(log2(n) + 1)
    width <- (max(x0) - d) / k
    MAX <- max(c(x0, x1))
    cutpoints <- seq(from = d, to = MAX + width, by = width)
    nx0 <- summary(cut(x0[x0 >= d], cutpoints, include.lowest = TRUE))
    nx1 <- summary(cut(x1[x1 >= d], cutpoints, include.lowest = TRUE))
    Nk <- nx1 - nx0
    Ck <- cutpoints[-1] - width/2
    return(list(Nk = Nk, Ck = Ck))
}
