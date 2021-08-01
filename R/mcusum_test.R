#' Change Point Test for Regression
#' 
#' Apply change point test by \insertCite{Horvath_etal_2017;textual}{funtimes} 
#' for detecting at-most-\eqn{m} changes in regression coefficients, where
#' test statistic is a modified cumulative sum (CUSUM), and 
#' critical values are obtained with sieve bootstrap \insertCite{Lyubchich_etal_2020_changepoints}{funtimes}.
#' 
#' @details The sieve bootstrap is applied by approximating regression residuals \code{e} 
#' with an AR(\eqn{p}) model using function \code{\link{ARest}}, 
#' where the autoregressive coefficients are estimated with \code{ar.method},
#' and order \eqn{p} is selected based on \code{ar.order} and \code{BIC} settings
#' (see \code{\link{ARest}}). At the next step, \code{B} autoregressive processes 
#' are simulated under the null hypothesis of no change points.
#' The distribution of test statistics \eqn{M_T} computed on each of those 
#' bootstrapped series is used to obtain bootstrap-based \eqn{p}-values for the test 
#' \insertCite{Lyubchich_etal_2020_changepoints}{funtimes}. 
#' 
#' The test statistic corresponds to the maximal value of the modified CUSUM over
#' all up to \code{m} combinations of hypothesized change points specified in \code{k}. The change 
#' points that correspond to that maximum are reported in \code{estimate$khat},
#' and their number is reported as \code{parameter}.
#' 
#'
#' @param e vector of regression residuals (a stationary time series). 
#' @param k an integer vector or scalar with hypothesized change point location(s) to 
#' test.
#' @param m an integer specifying the maximum number of change 
#' points being confirmed as statistically significant (from those 
#' specified in \code{k}) would be \eqn{\le m}. Thus, \code{m} must
#' be in 1,...,\code{k}.
#' @inheritParams wavk_test
#' @param ksm logical value indicating whether a kernel smoothing to innovations in sieve 
#' bootstrap shall be applied (default is \code{FALSE}, that is, the original estimated 
#' innovations are bootstrapped, without the smoothing).
#' @param ksm.arg used only if \code{ksm = TRUE}. A list of arguments for kernel smoothing
#' to be passed to \code{\link[stats]{density}} function. Default settings specify the 
#' use of Gaussian kernel and the \code{"sj"} rule to choose the bandwidth.
#' @param shortboot if \code{TRUE} (and \code{ksm} is \code{FALSE}), then a heuristic
#' is used to perform the test with a reduced number of bootstrap replicates.
#' Specifically, \code{B/4} replicates are used, which may reduce computing time by
#' up to 75 percent when the number of retained null hypotheses is large. 
#' A p-value of 999 is reported whenever a null hypothesis
#' is retained as a result of this mechanism.
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
#' SizeOfChange <- 1
#' TimeOfChange <- 50
#' Y <- c(1 * X[1:TimeOfChange] + E[1:TimeOfChange], 
#'       (1 + SizeOfChange)*X[(TimeOfChange + 1):T] + E[(TimeOfChange + 1):T])
#' ehat <- lm(Y ~ X)$resid
#' mcusum_test(ehat, k = c(30, 50, 70))
#' 
#' #Same, but with bootstrapped innovations obtained from a kernel smoothed distribution:
#' mcusum_test(ehat, k = c(30, 50, 70), ksm = TRUE)
#' 
mcusum_test <- function(e, k,
                        m = length(k),
                        B = 1000, 
                        shortboot = FALSE,
                        ksm = FALSE, 
                        ksm.arg = list(kernel = "gaussian", bw = "sj"), ...)
{
    DNAME <- deparse(substitute(e))
    T <- length(e)
    e <- e - .Internal(mean(e))
    # AB: was:
    # e <- e - mean(e)
    # AB: Sort k so that it doesn't have to be sorted by Mfun and MTfun:
    k <- sort(unique(k))
    m <- min(m, length(k))
    
    phi <- ARest(e, ...)
    MTobs <- MTfun(e, k = k, m = m)
    if (length(phi) > 0) {
        e <-  as.vector(embed(e, length(phi) + 1L) %*% c(1, -phi))
        # AB: use more efficient .Internal(mean()):
        e <- e - .Internal(mean(e))
    }
    if (ksm) { #use e from a smoothed distribution of e
        ksm.arg$x <- e #append x to the arguments of the density function
        bw <- do.call(stats::density, ksm.arg) #estimate bandwidth
        bw <- bw$bw
        MTboot <- sapply(1:B, function(b) 
            MTfun(arima.sim(T, model = list(order = c(length(phi), 0, 0), ar = phi), 
                            innov = rnorm(T, mean = sample(e, size = T, replace = TRUE), sd = bw)), 
                  k = k, m = m)$MT
        )
        pval <- (1 + sum(MTboot >= MTobs$MT)) / (B + 1)
        
    } else {#use bootstrapped e
        
        if(shortboot){
            # AB:
            # Use a heuristic to reduce number of bootstrap replicates
            # needed to retain null hypothesis:
            # E.g. when using B = 1000 and 100 out of the first 250
            # bootstrapped statistics are >= the one obtained from the
            # actual sample, then the p-value must be >= 0.1 even if
            # the remaining 750 bootstrapped values were smaller.

            B_part <- ceiling(B/4) # use a quarter of B, but could also be less/more
            sig <- ceiling(B/10) # portion of bootstraps that has to be bigger than original for alpha = 0.1
            thr <- sig/B_part # prior threshold to discard time series which won't reach significant threshold anymore
            
            MTboot <- sapply(1:B_part, function(b)
                MTfun(as.vector(arima.sim(T, model = list(order = c(length(phi), 0, 0), ar = phi),
                                           innov = sample(e, size = T, replace = TRUE))),
                       k = k, m = m)$MT
            )
            
            if(mean(MTboot >= MTobs$MT) < thr){
                MTboot2 <- sapply(1:(B-B_part), function(b)
                    MTfun(as.vector(arima.sim(T, model = list(order = c(length(phi), 0, 0), ar = phi),
                                               innov = sample(e, size = T, replace = TRUE))),
                           k = k, m = m)$MT
                )
                MTboot <- c(MTboot, MTboot2)
                # AB:
                # p-value formula from Davison and Hinkley (1997), 
                # Bootstrap Methods and their Application, p. 141:
                pval <- (1 + sum(MTboot >= MTobs$MT)) / (B + 1)
            } else {
                # In this situation, we only know that the test is not significant:
                pval <- 999.0
            }
            
        } else { # no shortening of bootstrap
            MTboot <- sapply(1:B, function(b)
                MTfun(as.vector(arima.sim(T, model = list(order = c(length(phi), 0, 0), ar = phi),
                                           innov = sample(e, size = T, replace = TRUE))),
                       k = k, m = m)$MT
            )
            # AB:
            # p-value formula from Davison and Hinkley (1997), 
            # Bootstrap Methods and their Application, p. 141:
            pval <- (1 + sum(MTboot >= MTobs$MT)) / (B + 1)
        }
    }
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
                   p.value = pval,
                   alternative = paste("at-most-", m, " changes exist", sep = ""), 
                   estimate = ESTIMATE), 
              class = "htest") 
}
