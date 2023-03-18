#' Estimation of Autoregressive (AR) Parameters
#'
#' Estimate parameters \eqn{\phi} of autoregressive time series model
#' \deqn{X_t = \sum_{i=1}^p\phi_iX_{t-i} + e_t,}
#' by default using robust difference-based estimator and Bayesian information
#' criterion (BIC) to select the order \eqn{p}. This function is employed
#' for time series filtering in the functions \code{\link{notrend_test}}, \code{\link{sync_test}},
#' and \code{\link{wavk_test}}.
#'
#' @details The formula for information criteria used consistently for all methods:
#' \deqn{IC=n\ln(\hat{\sigma}^2) + (p + 1)k,}
#' where \eqn{n} = \code{length(x)},
#' \eqn{p} is the autoregressive order (\eqn{p + 1} is the number of model parameters),
#' and \eqn{k} is the penalty (\eqn{k = \ln(n)} in BIC, and \eqn{k = 2} in AIC).
#'
#'
#' @param x a vector containing a univariate time series. Missing values are not allowed.
#' @param ar.order order of the autoregressive model when \code{ic = "none"}, or
#' the maximal order for IC-based filtering. Default is
#' \code{round(10*log10(length(x)))}, where \code{x} is the time series.
#' @param ar.method method of estimating autoregression coefficients.
#' Default \code{"HVK"} delivers robust difference-based estimates by
#' \insertCite{Hall_VanKeilegom_2003;textual}{funtimes}. Alternatively,
#' options of \command{ar} function can be used, such as \code{"burg"},
#' \code{"ols"}, \code{"mle"}, and \code{"yw"}.
#' @param ic information criterion used to select the order of autoregressive filter (AIC of BIC),
#' considering models of orders \eqn{p=} 0,1,...,\code{ar.order}.
#' If \code{ic = "none"}, the AR(\eqn{p}) model with \eqn{p=} \code{ar.order} is used,
#' without order selection.
#'
#'
#' @return A vector of estimated AR coefficients. Returns \code{numeric(0)} if
#' the final \eqn{p=0}.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link[stats]{ar}}, \code{\link{HVK}},
#' \code{\link{notrend_test}}, \code{\link{sync_test}}, \code{\link{wavk_test}}
#'
#' @keywords ts
#'
#' @author Vyacheslav Lyubchich
#'
#' @export
#' @examples
#' # Simulate a time series Y:
#' Y <- arima.sim(n = 200, list(order = c(2, 0, 0), ar = c(-0.7, -0.1)))
#' plot.ts(Y)
#'
#' # Estimate the coefficients:
#' ARest(Y) # HVK, by default
#' ARest(Y, ar.method = "yw") # Yule--Walker
#' ARest(Y, ar.method = "burg") # Burg
#'
ARest <- function(x,
                  ar.order = NULL,
                  ar.method = "HVK",
                  ic = c("BIC", "AIC", "none"))
{
    x <- as.vector(x)
    n <- length(x)
    if (is.null(ar.order)) {
        ar.order <- round(10*log10(n))
    }
    ic <- match.arg(ic)
    if (ic == "BIC") {
        useIC <- TRUE
        k <- log(n)
    } else if (ic == "AIC") {
        useIC <- TRUE
        k <- 2
    } else {
        useIC <- FALSE
    }
    ics <- rep(NA, ar.order + 1)
    ics[1] <- n*log(var(x)) # no AR-filtering (ar.order == 0)
    pheta <- numeric(0) # if no AR-filtering, otherwise will be redefined below
    if (ar.order > 0) {
        if (!useIC) { # useIC == FALSE, use fixed ar.order > 0
            if (ar.method == "HVK") {
                pheta <- HVK(x, ar.order = ar.order)
            } else {
                a <- ar(x, aic = FALSE, order.max = ar.order,
                        demean = TRUE, method = ar.method)
                pheta <- a$ar
            }
        } else { # IC-based filtering
            PHETAS <- list()
            for (i in 2:length(ics)) {
                if (ar.method == "HVK") {
                    PHETAS[[i]] <- HVK(x, ar.order = i - 1)
                } else {
                    a <- ar(x, aic = FALSE, order.max = i - 1,
                            demean = TRUE, method = ar.method)
                    PHETAS[[i]] <- a$ar
                }
                tmp <- filter(x, PHETAS[[i]], sides = 1)
                et <- x[i:n] - tmp[(i - 1):(n - 1)]
                ics[i] <- n*log(var(et)) + i*k # here use i = ARorder p + 1 (variance)
            }
            min_ic_index <- which.min(ics)
            if (min_ic_index > 1) {
                pheta <- PHETAS[[min_ic_index]]
            }
        }
    }
    return(pheta)
}
