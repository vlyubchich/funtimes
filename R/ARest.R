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
                  ar.method = c("HVK", "burg", "ols", "mle", "yw"),
                  ic = c("BIC", "AIC", "none"))
{
    x <- as.vector(x)
    if (!is.numeric(x) || anyNA(x)) {
        stop("x must be a numeric vector without missing values.")
    }
    n <- length(x)
    if (n < 5) {
        stop("x must contain at least 5 observations.")
    }

    if (is.null(ar.order)) {
        ar.order <- round(10 * log10(n))
    }
    if (length(ar.order) != 1L || !is.finite(ar.order) || ar.order < 0) {
        stop("ar.order must be a single non-negative number.")
    }
    ar.order <- as.integer(ar.order)

    ar.method <- match.arg(ar.method)
    ic <- match.arg(ic)

    x.var <- var(x)
    if (!is.finite(x.var) || x.var <= 0) {
        stop("x must have positive finite variance.")
    }

    # Handle the case where no information criterion is used for model selection
    if (ic == "none") {
        if (ar.order == 0) {
            return(numeric(0))
        }
        if (ar.method == "HVK") {
            pheta <- HVK(x, ar.order = ar.order)
        } else {
            a <- ar(x, aic = FALSE, order.max = ar.order,
                    demean = TRUE, method = ar.method)
            pheta <- a$ar
        }
        return(pheta)
    }

    # IC-based order selection
    k <- if (ic == "BIC") log(n) else 2

    # IC for order 0
    ics <- n * log(x.var)
    phetas_list <- list(numeric(0)) # For order p = 0

    if (ar.order > 0) {
        # Calculate IC for orders 1 to ar.order
        for (p in 1:ar.order) {
            if (ar.method == "HVK") {
                pheta <- HVK(x, ar.order = p)
            } else {
                a <- ar(x, aic = FALSE, order.max = p,
                        demean = TRUE, method = ar.method)
                pheta <- a$ar
            }
            phetas_list[[p + 1]] <- pheta

            # Calculate residuals
            res <- filter(x, pheta, sides = 1)
            et <- x[(p + 1):n] - res[p:(n - 1)]

            # Update IC vector
            ics[p + 1] <- n * log(var(et)) + (p + 1) * k
        }
    }

    # Find the best model order
    best_p_idx <- which.min(ics)
    return(phetas_list[[best_p_idx]])
}
