#' Sieve Bootstrap Based Test for the Null Hypothesis of no Trend
#'
#' A combination of time series trend tests for testing the null hypothesis of no trend,
#' versus the alternative hypothesis of a linear trend (Student's t-test),
#' or monotonic trend (Mann--Kendall test), or more general,
#' possibly non-monotonic trend (WAVK test).
#'
#' @details This function tests the null hypothesis of no trend
#' versus different alternatives.
#' To set some other shape of trend as the null hypothesis, use \code{\link{wavk_test}}.
#' Note that \code{\link{wavk_test}} employs hybrid bootstrap, which is an alternative
#' to the sieve bootstrap employed by the current function.
#'
#'
#' @inheritParams ARest
#' @inheritParams wavk_test
#' @param test trend test to implement: Student's t-test (\code{"t"}, default),
#' Mann--Kendall test (\code{"MK"}), or
#' WAVK test (\code{"WAVK"}, see \code{\link{WAVK}}).
#' @param factor.length method to define the length of local windows (factors).
#' Used only if \code{test = "WAVK"}. Default option \code{"user.defined"} allows
#' to set only one value of the argument \code{Window}. The option
#' \code{"adaptive.selection"} sets \code{method = "boot"} and employs
#' heuristic \eqn{m}-out-of-\eqn{n} subsampling algorithm
#' \insertCite{Bickel_Sakov_2008}{funtimes} to select an optimal window from the set
#' of possible windows \code{length(x)*q^j} whose values are mapped to the largest
#' previous integer and greater than 2. Vector \code{x} is the time series tested.
#' @param Window length of the local window (factor), default is
#' \code{round(0.1*length(x))}. Used only if \code{test = "WAVK"}.
#' This argument is ignored if\cr \code{factor.length = "adaptive.selection"}.
#' @param q scalar from 0 to 1 to define the set of possible windows when
#' \code{factor.length =} \code{"adaptive.selection"}.
#' Used only if \code{test = "WAVK"}. Default is \eqn{3/4}.
#' This argument is ignored if \code{factor.length =} \code{"user.defined"}.
#' @param j numeric vector to define the set of possible windows when
#' \code{factor.length =} \code{"adaptive.selection"}.
#' Used only if \code{test = "WAVK"}. Default is \code{c(8:11)}.
#' This argument is ignored if \code{factor.length =} \code{"user.defined"}.
#'
#'
#' @return A list with class \code{"htest"} containing the following components:
#' \item{method}{name of the method.}
#' \item{data.name}{name of the data.}
#' \item{statistic}{value of the test statistic.}
#' \item{p.value}{\eqn{p}-value of the test.}
#' \item{alternative}{alternative hypothesis.}
#' \item{estimate}{list with the following elements: employed AR order and estimated
#' AR coefficients.}
#' \item{parameter}{window that was used in WAVK test, included in the output only
#' if \code{test = "WAVK"}.}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link[stats]{ar}}, \code{\link{HVK}}, \code{\link{WAVK}},
#' \code{\link{wavk_test}}, \code{vignette("trendtests", package = "funtimes")}
#'
#' @keywords htest ts trend
#'
#' @author Vyacheslav Lyubchich
#'
#' @export
#' @examples
#' \dontrun{
#' # Fix seed for reproducible simulations:
#' set.seed(1)
#'
#' #Simulate autoregressive time series of length n with smooth linear trend:
#' n <- 200
#' tsTrend <- 1 + 2*(1:n/n)
#' tsNoise <- arima.sim(n = n, list(order = c(2, 0, 0), ar = c(0.5, -0.1)))
#' U <- tsTrend + tsNoise
#' plot.ts(U)
#'
#' #Use t-test
#' notrend_test(U)
#'
#' #Use Mann--Kendall test and Yule-Walker estimates of the AR parameters
#' notrend_test(U, test = "MK", ar.method = "yw")
#'
#' #Use WAVK test for the H0 of no trend, with m-out-of-n selection of the local window:
#' notrend_test(U, test = "WAVK", factor.length = "adaptive.selection")
#' # Sample output:
#' ##	Sieve-bootstrap WAVK trend test
#' ##
#' ##data:  U
#' ##WAVK test statistic = 21.654, moving window = 15, p-value < 2.2e-16
#' ##alternative hypothesis: (non-)monotonic trend.
#' ##sample estimates:
#' ##$AR_order
#' ##[1] 1
#' ##
#' ##$AR_coefficients
#' ##    phi_1
#' ##0.4041848
#' }
#'
notrend_test <- function(x, B = 1000, test = c("t", "MK", "WAVK"),
                         ar.method = "HVK", ar.order = NULL, ic = "BIC",
                         factor.length = c("user.defined", "adaptive.selection"),
                         Window = NULL, q = 3/4, j = c(8:11))
{
    ### Perform various checks.
    DNAME <- deparse(substitute(x))
    if (NCOL(x) > 1 || !is.numeric(x)) {
        stop("x is not a vector or univariate time series.")
    }
    if (any(is.na(x))) {
        stop("x contains missing values.")
    }
    x <- as.vector(x)
    n <- length(x)
    if (n < 5)
        stop("x must contain at least 5 observations.")

    if (!is.numeric(B) || length(B) != 1L || is.na(B))
        stop("B must be a single non-missing numeric value.")
    B <- as.integer(B)
    if (B <= 0)
        stop("number of bootstrap samples B must be positive.")

    if (!is.null(ar.order)) {
        if (length(ar.order) != 1L || is.na(ar.order))
            stop("ar.order must be a single non-missing value.")
        ar.order <- as.integer(ar.order)
        if (ar.order < 0)
            stop("ar.order must be non-negative.")
    }

    test <- match.arg(test)
    if (identical(test, "WAVK")) { #checks only for WAVK
        factor.length <- match.arg(factor.length)
        if (is.null(Window)) {
            Window <- round(0.1*n)
        }
        if (length(Window) != 1L || is.na(Window))
            stop("Window must be a single non-missing value.")
        Window <- as.integer(Window)
        if (Window < 2 || Window >= n)
            stop("Window must be in [2, length(x)).")

        if (length(q) != 1L || is.na(q))
            stop("q must be a single non-missing numeric value.")
        if (q <= 0 || q >= 1)
            stop("q must be in (0, 1).")

        if (!is.vector(j) || !is.numeric(j))
            stop("j must be a numeric vector.")
        if (any(is.na(j)))
            stop("j must not contain missing values.")

        if (factor.length == "user.defined") {
            kn <- Window
        } else {
            kn <- n*q^j
        }
        kn <- unique(sort(floor(kn)))
        kn <- kn[kn > 2 & kn < n]
        if (length(kn) == 0)
            stop("set a proper window. Check Window parameter or adjust q/j values.")
        if (factor.length == "adaptive.selection" && length(kn) < 3)
            stop("number of possible windows is not enough for adaptive selection. Change parameters 'q' and/or 'j'.")
    }
    
    ### Sieve bootstrap procedure.
    pheta <- ARest(x, ar.order = ar.order, ar.method = ar.method, ic = ic)
    if (length(pheta) > 0) {
        names(pheta) <- paste(rep("phi_", length(pheta)), 1:length(pheta), sep = "")
        tmp <- stats::filter(x, pheta, sides = 1)
        Z <- x[(length(pheta) + 1):n] - tmp[length(pheta):(n - 1)]
    } else {
        pheta <- numeric(0)
        Z <- x
    }
    Z_centered <- na.omit(Z) - mean(na.omit(Z))
    
    # Generate bootstrap samples in a single block
    Y <- vapply(1:B, function(i) {
        e <- sample(Z_centered, size = n, replace = TRUE)
        arima.sim(list(order = c(length(pheta), 0, 0), ar = pheta), n = n, innov = e)
    }, numeric(n))

    ESTIMATE <- list(length(pheta), pheta)
    names(ESTIMATE) <- c("AR_order", "AR_coefficients")
    
    t_seq <- 1:n / n
    
    # Perform the chosen test
    if (test == "t") {
        METHOD <- "Sieve-bootstrap Student's t-test for a linear trend"
        ALTERNATIVE <- "linear trend."
        STATISTIC <- summary(lm(x ~ t_seq))$coefficients["t_seq", "t value"]
        names(STATISTIC) <- "Student's t value"
        boot.stat <- apply(Y, 2, function(y_col) summary(lm(y_col ~ t_seq))$coefficients["t_seq", "t value"])
    } else if (test == "MK") {
        METHOD <- "Sieve-bootstrap Mann--Kendall's trend test"
        ALTERNATIVE <- "monotonic trend."
        STATISTIC <- mann_kendall_tau(x)
        names(STATISTIC) <- "Mann--Kendall's tau"
        boot.stat <- apply(Y, 2, mann_kendall_tau)
    } else { # WAVK test
        METHOD <- "Sieve-bootstrap WAVK trend test"
        ALTERNATIVE <- "(non-)monotonic trend."
        if (length(kn) < 3) {
            kn_opt <- kn[1]
            boot.stat <- apply(Y, 2, function(y_col) WAVK(y_col, kn_opt)$Tns)
        } else {
            s <- vapply(kn, function(k_i) apply(Y, 2, function(y_col) WAVK(y_col, k_i)$Tns), numeric(B))
            s_sorted <- apply(s, 2, sort) # Sort each column (for each kn)
            distances <- vapply(1:(ncol(s_sorted) - 1), function(i) dist(t(s_sorted[, i:(i + 1)])), numeric(1))
            argmin <- which.min(distances)
            kn_opt <- kn[argmin]
            boot.stat <- sort(s[, argmin]) # Use the sorted stats for the optimal window
        }
        STATISTIC <- WAVK(x, kn_opt)$Tns
        names(STATISTIC) <- "WAVK test statistic"
        PARAMETER <- kn_opt
        names(PARAMETER) <- "moving window"
    }
    
    P.VALUE <- mean(abs(boot.stat) >= abs(STATISTIC))
    
    # Construct the htest object
    result <- list(method = METHOD,
                   data.name = DNAME,
                   statistic = STATISTIC,
                   p.value = P.VALUE,
                   alternative = ALTERNATIVE,
                   estimate = ESTIMATE)
    if (test == "WAVK") {
        result$parameter <- PARAMETER
    }
    
    structure(result, class = "htest")
}
