#' Change Point Test for Regression
#' 
#' Apply change point test by \insertCite{Horvath_etal_2017;textual}{funtimes} 
#' for detecting at most \eqn{m} changes in regression coefficients, where
#' test statistic is a modified cumulative sum (CUSUM), and 
#' critical values are obtained with sieve bootstrap.
#' 
#' @details The sieve bootstrap is applied by approximating regression residuals \code{e} 
#' with an AR(\eqn{p}) model using function \code{\link{ARest}}, 
#' where the autoregressive coefficients are estimated with \code{ar.method},
#' and order \eqn{p} is selected based on \code{ar.order} and \code{BIC} settings
#' (see \code{\link{ARest}}). At the next step, \code{B} autoregressive processes 
#' are simulated under the null hypothesis of no change points.
#' The distribution of test statistics \eqn{M_T} computed on each of those 
#' bootstrapped series is used to obtain bootstrap-based \eqn{p}-values for the test. 
#' The test statistic corresponds to the maximal value of the modified CUSUM over
#' all combinations of hypothesized change points specified in \code{k}. The change 
#' points that correspond to that maximum are reported in \code{estimate$khat},
#' and their number is reported as \code{parameter}.
#' 
#'
#' @param e vector of regression residuals (a stationary time series). 
#' @param k an integer vector or scalar with hypothesized change point location(s) to 
#' test. The length of this vector is treated as \eqn{m}, that is, the number of change 
#' points being confirmed as statistially significant (from those 
#' specified in \code{k}) would be \eqn{\le m}.
#' @inheritParams wavk.test
#' @param ... additional arguments passed to \code{\link{ARest}}
#' (for example, \code{ar.method}).
#'
#'
#' @return A list of class \code{"htest"} containing the following components:
#' \item{method}{name of the method.}
#' \item{data.name}{name of the data.}
#' \item{statistic}{obseved value of the test statistic.}
#' \item{parameter}{\code{mhat} is the final number of change points, 
#' from those specified in the input \code{k},
#' for which the test statistic is reported. 
#' See the corresponding locations, \code{khat}, in the \code{estimate}.}
#' \item{p.value}{bootstrapped \eqn{p}-value of the test.}
#' \item{alternative}{alternative hypothesis.}
#' \item{estimate}{list with elements: \code{AR_order} and  
#' \code{AR_coefficients} (the autoregressive order and estimated autoregressive
#' coefficients used in sieve bootstrap procedure), \code{khat} (final change points, 
#' from those specified in the input \code{k} for which the test statistic is reported), 
#' and \code{B} (the number of bootstrap replications).} 
#' 
#' @references
#' \insertAllCited{}
#' 
#' @keywords changepoint htest ts
#' 
#' @author Vyacheslav Lyubchich
#' 
#' @export
#' @examples
#' ##### Model 1 with normal errors, by Horvath et al. (2017)
#' T <- 100 #length of time series
#' X <- rnorm(T, mean = 1, sd = 1)
#' E <- rnorm(T, mean = 0, sd = 1)
#' SizeOfChange <- 1.5
#' TimeOfChange <- 50
#' Y <- c(1 * X[1:TimeOfChange] + E[1:TimeOfChange], 
#'       (1 + SizeOfChange)*X[(TimeOfChange+1):T] + E[(TimeOfChange+1):T])
#' ehat <- lm(Y ~ X)$resid
#' mcusum.test(ehat, k = c(30, 50, 70))
#' 
mcusum.test <- function(e, k, B = 1000, ...)
{
    DNAME <- deparse(substitute(e))
    T <- length(e)
    e <- e - mean(e)
    phi <- ARest(e, ...)
    MTobs <- MTfun(e, k = k)
    if(length(phi) > 0) {
        e <-  as.vector(embed(e, length(phi) + 1L) %*% c(1, -phi))
        e <- e - mean(e)
    }
    MTboot <- sapply(1:B, function(b) 
        MTfun(arima.sim(T, model = list(order = c(length(phi), 0, 0), ar = phi), 
                        innov = sample(e, size = T, replace = TRUE)), k = k)$MT
    )
    STATISTIC <- MTobs$MT
    names(STATISTIC) <- "M_T"
    PARAMETER <- length(MTobs$k)
    names(PARAMETER) <- "mhat"
    ESTIMATE <- list(length(phi), phi, MTobs$k, B)
    names(ESTIMATE) <- c("AR_order", "AR_coefficients", "khat", "B")
    structure(list(method = "Test for at-most-m changes in linear regression model", 
                   data.name = DNAME, 
                   statistic = STATISTIC, 
                   parameter = PARAMETER, 
                   p.value = mean(MTboot >= MTobs$MT),
                   alternative = paste("at-most-", length(k), " changes exist", sep = ""), 
                   estimate = ESTIMATE), 
              class = "htest") 
}
