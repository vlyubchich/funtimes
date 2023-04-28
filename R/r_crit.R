#' Critical Value for Correlation Coefficient
#'
#' Calculate critical value for correlation coefficient,
#' for a given sample size and confidence level.
#' If absolute value of an observed correlation coefficient is higher
#' than the critical value, then the correlation is statistically significant
#' with the specified confidence.
#' This approach is
#' identical to using textbook tables of critical values and
#' alternative to calculating \eqn{p}-values.
#'
#' @details
#' Using Student's \eqn{t}-distribution, the critical value is
#' \deqn{r_{crit} = \frac{t}{\sqrt{n - 2 + t^2}},}
#' where \eqn{t} is a quantile of \eqn{t}-distribution with
#' \eqn{n - 2} degrees of freedom for probability \eqn{1 - (1 - conf.level)/2}.
#'
#' @param n sample size(s) used to calculate correlations.
#' Values of \eqn{n < 4} are omitted.
#' @param conf.level confidence level for calculating the critical value(s).
#' Default is 0.95 (i.e., confidence of 95%).
#' If length of the input is higher than 1, only the first element is used.
#' @param method Method for calculating the critical values.
#' Currently only the method based on \eqn{t}-distribution is used.
#'
#' @return A vector of critical values.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link[stats]{cor.test}}, \code{\link[stats]{cor}}
#'
#' @keywords correlation threshold
#'
#' @author Vyacheslav Lyubchich
#'
#' @importFrom stats qt
#' @export
#' @examples
#' r_crit(120)
#' r_crit(20:30, conf.level = 0.9)
#'
r_crit <- function(n,
                   conf.level = 0.95,
                   method = c("t")) {
    if (any(n < 4)) {
        warning("Values of n < 4 were omitted.")
        n <- n[n >= 4]
    }
    alpha <- 1 - conf.level[1]
    t2 <- qt(p = 1 - alpha/2, df = n - 2)^2
    sqrt(t2 / (n - 2 + t2))
}
