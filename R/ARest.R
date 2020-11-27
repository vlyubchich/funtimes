#' Estimation of Autoregressive (AR) Parameters
#' 
#' Estimate parameters \eqn{\phi} of autoregressive time series model
#' \deqn{X_t = \sum_{i=1}^p\phi_iX_{t-i} + e_t,} 
#' by default using robust difference-based estimator and Bayesian information 
#' criterion (BIC) to select the order \eqn{p}. This function is employed 
#' for time series filtering in functions \code{\link{sync_test}} 
#' and \code{\link{wavk_test}}.
#' 
#' @details The same formula for BIC is used consistently for all methods:
#' \deqn{BIC=n\ln(\hat{\sigma}^2) + k\ln(n),} 
#' where \eqn{n} = \code{length(x)}, \eqn{k=p+1}.
#' 
#' 
#' @param x a vector containing a univariate time series. Missing values are not allowed.
#' @param ar.order order of autoregressive model when \code{BIC = FALSE}, or 
#' the maximal order for BIC-based filtering. Default is 
#' \code{round(10*log10(length(x)))}, where \code{x} is the time series.
#' @param ar.method method of estimating autoregression coefficients. 
#' Default \code{"HVK"} delivers robust difference-based estimates by 
#' \insertCite{Hall_VanKeilegom_2003;textual}{funtimes}. Alternatively, 
#' options of \command{ar} function can be used, such as \code{"burg"}, 
#' \code{"ols"}, \code{"mle"}, and \code{"yw"}.
#' @param BIC logical value indicates whether the order of autoregressive 
#' filter should be selected by Bayesian information criterion (BIC). 
#' If \code{TRUE} (default), models of orders \eqn{p=} 0,1,...,\code{ar.order} 
#' or \eqn{p=} 0,1,...,\code{round(10*log10(length(x)))} are considered, 
#' depending on whether \code{ar.order} is defined or not 
#' (\code{x} is the time series).
#' 
#' 
#' @return A vector of estimated AR coefficients. Returns \code{numeric(0)} if 
#' the final \eqn{p=0}.
#'  
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link[stats]{ar}}, \code{\link{HVK}}, 
#' \code{\link{sync_test}}, \code{\link{wavk_test}}
#' 
#' @keywords ts
#' 
#' @author Vyacheslav Lyubchich
#' 
#' @export
#' @examples
#' # Fix seed for reproducible simulations:
#' set.seed(1)
#' 
#' #Simulate some time series, possibly with trend:
#' n <- 100
#' Y <- arima.sim(n = n, list(order = c(2, 0, 0), ar = c(-0.7, -0.1)))
#' plot.ts(Y)
#' 
#' #Estimate the coefficients:
#' ARest(Y) #HVK by default
#' ARest(Y, ar.method = "yw") #Yule--Walker
#' ARest(Y, ar.method = "burg") #Burg
#' 
ARest <- function(x, ar.order = NULL, ar.method = "HVK", BIC = TRUE)
{
    x <- as.vector(x)
    n <- length(x)
    if (is.null(ar.order)) {
        ar.order <- round(10*log10(n))
    }
    bic <- rep(NA, ar.order + 1)
    bic[1] <- n*log(var(x)) #no AR-filtering (ar.order==0)
    pheta <- numeric(0) #if no AR-filtering, otherwise will be redefined below  
    if (ar.order > 0) {
        if (!BIC) { #BIC==FALSE, use fixed ar.order>0
            if (ar.method == "HVK") {
                pheta <- HVK(x, ar.order = ar.order)
            } else {
                a <- ar(x, aic = FALSE, order.max = ar.order, 
                        demean = TRUE, method = ar.method)
                pheta <- a$ar
            }
        } else {#BIC-based filtering
            for (i in 2:length(bic)) {
                if (ar.method == "HVK") {
                    pheta0 <- HVK(x, ar.order = i - 1)
                } else {
                    a <- ar(x, aic = FALSE, order.max = i - 1, 
                            demean = TRUE, method = ar.method)
                    pheta0 <- a$ar
                }
                tmp <- filter(x, pheta0, sides = 1)
                et <- x[i:n] - tmp[(i - 1):(n - 1)]
                bic[i] <- n*log(var(et)) + i*log(n) #here use i = ARorder p + 1 (variance)
            }
            if (which.min(bic) > 1) {
                if (ar.method == "HVK") {
                    pheta <- HVK(x, ar.order = (which.min(bic) - 1))
                } else {
                    a <- ar(x, aic = FALSE, order.max = (which.min(bic) - 1), 
                            demean = TRUE, method = ar.method)
                    pheta <- a$ar
                }
            }
        }
    }
    return(pheta)
}
