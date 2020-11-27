#' Sieve Bootstrap Based Test for the Null Hypothesis of no Trend
#' 
#' A combination of time series trend tests for testing the null hypothesis of no trend, 
#' versus the alternative hypothesis of a linear trend (Student's t-test), 
#' or monotonic trend (Mann--Kendall test), or more general, 
#' possibly non-monotonic trend (WAVK test).
#' 
#' @details This function tests the null hypothesis of no trend 
#' versus different alternatives.
#' To set some other shape of trend as the null hypothesis, use \code{\link{wavk.test}}. 
#' Note that \code{\link{wavk.test}} employs hybrid bootstrap, which is alternative 
#' to the sieve bootstrap employed by the current function.
#' 
#' 
#' @inheritParams ARest
#' @inheritParams wavk.test
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
#' \code{\link{wavk.test}}
#' 
#' @keywords internal
NULL

#' @rdname funtimes-deprecated
#' @section \code{notrend.test}:
#' For \code{notrend.test}, use \code{\link{notrend_test}}.
#' 
#' @author Yulia R. Gel, Vyacheslav Lyubchich
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
#' notrend.test(U)
#'     
#' #Use Mann--Kendall test and Yule-Walker estimates of the AR parameters
#' notrend.test(U, test = "MK", ar.method = "yw")
#'     
#' #Use WAVK test for the H0 of no trend, with m-out-of-n selection of the local window:
#' notrend.test(U, test = "WAVK", factor.length = "adaptive.selection")
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
notrend.test <- function(x, B = 1000, test = c("t", "MK", "WAVK"), 
                         ar.method = "HVK", ar.order = NULL, BIC = TRUE, 
                         factor.length = c("user.defined", "adaptive.selection"), 
                         Window = NULL, q = 3/4, j = c(8:11))
{
    .Deprecated("notrend_test", msg = "notrend.test is deprecated and will be removed. Use notrend_test instead.")
    ### Perform various checks.
    DNAME <- deparse(substitute(x))  
    if (NCOL(x) > 1 | !is.numeric(x)) {
        stop("x is not a vector or univariate time series.")
    }
    if (any(is.na(x))) {
        stop("x contains missing values.")
    }
    x <- as.vector(x)
    n <- length(x)
    test <- match.arg(test)
    factor.length <- match.arg(factor.length)  
    if (is.null(Window)) {
        Window = round(0.1*n)
    }
    if (NCOL(q) > 1 | !is.numeric(q) | NROW(q) > 1) {
        stop("q is not a scalar.")
    }
    if (q >= 1 | q <= 0) {
        stop("q is out of range from 0 to 1.")
    }
    if (!is.vector(j) | !is.numeric(j)) {
        stop("j is not a numeric vector.")
    }
    if (factor.length == "user.defined") {
        kn <- Window[1]
    } else {
        kn <- n*q^j
    }  
    kn <- unique(sort(floor(kn)))
    kn <- kn[kn > 2 & kn < n]
    if (length(kn) == 0) {
        stop("set a proper window.")
    }
    if (factor.length == "adaptive.selection" & length(kn) < 3) {
        stop("number of possible windows is not enough for adaptive selection. Change parameters 'q' and/or 'j'.")
    }
    B <- round(B)
    if (B <= 0) {
        stop("number of bootstrap samples B must be positive.")
    }
    if (!is.null(ar.order) & (NCOL(ar.order) > 1 | !is.numeric(ar.order) | NROW(ar.order) > 1)) {
        stop("ar.order is not a scalar.")
    }
    if (!is.null(ar.order) && ar.order < 0) {
        stop("ar.order must be non-negative.")
    }
    ### Function.
    Y <- array(data = NA, c(n, B))
    t <- c(1:n)/n
    pheta <- ARest(x, ar.order = ar.order, ar.method = ar.method, BIC = BIC)
    if (length(pheta) > 0) {
        names(pheta) <- paste(rep("phi_", length(pheta)), c(1:length(pheta)), sep = "")
        tmp <- filter(x, pheta, sides = 1)
        Z <- x[(length(pheta) + 1):n] - tmp[length(pheta):(n - 1)]
        for (i in 1:B) {
            e <- sample(Z, size = n, replace = TRUE)
            Y[ ,i] <- arima.sim(list(order = c(length(pheta), 0, 0), ar = pheta), n = n, innov = e)
        }
    } else {
        Z <- x
        for (i in 1:B) {
            Y[ ,i] <- sample(Z, size = n, replace = TRUE)
        }
    }
    Z <- na.omit(Z) - mean(Z)
    ESTIMATE <- list(length(pheta), pheta)
    names(ESTIMATE) <- c("AR_order", "AR_coefficients")
    #If Student's t-test is used
    if (test == "t") {
        METHOD <- "Sieve-bootstrap Student's t-test for a linear trend"
        ALTERNATIVE <- "linear trend."
        STATISTIC <- summary(lm(x ~ t))$coefficients["t", "t value"]
        names(STATISTIC) <- "Student's t value"
        boot.stat <- sapply(1:dim(Y)[2], function(i) summary(lm(Y[,i] ~ t))$coefficients["t", "t value"])
    }
    #If Mann--Kendall's test is used
    if (test == "MK") {
        METHOD <- "Sieve-bootstrap Mann--Kendall's trend test"
        ALTERNATIVE <- "monotonic trend."
        STATISTIC <- MannKendall(x)$tau
        names(STATISTIC) <- "Mann--Kendall's tau"
        boot.stat <- sapply(1:dim(Y)[2], function(i) MannKendall(Y[,i])$tau)
    }
    #If WAVK test is used
    if (test == "WAVK") {
        METHOD <- "Sieve-bootstrap WAVK trend test"
        ALTERNATIVE <- "(non-)monotonic trend."
        if (length(kn) < 3) {
            kn_opt <- kn[1]
            boot.stat <- sapply(1:dim(Y)[2], function(j) WAVK(Y[,j], kn_opt)$Tns)
        } else {
            s <- array(data = NA, c(length(kn), B))
            for (i in 1:length(kn)) {
                s[i,] <- sapply(1:dim(Y)[2], function(j) WAVK(Y[,j], kn[i])$Tns)
            }
            s <- t(apply(s, 1, sort))
            distance <- sapply(1:(length(kn) - 1), function(x) dist(s[x:(x + 1),]))
            argmin <- which.min(distance)
            kn_opt <- kn[argmin]
            boot.stat <- s[argmin,]
        }
        STATISTIC <- WAVK(x, kn_opt)$Tns
        names(STATISTIC) <- "WAVK test statistic"
        PARAMETER <- kn_opt
        names(PARAMETER) <- "moving window"
    }
    P.VALUE <- mean(abs(boot.stat) >= abs(STATISTIC))
    if (test == "WAVK") {
        structure(list(method = METHOD, data.name = DNAME, statistic = STATISTIC, p.value = P.VALUE, 
                       alternative = ALTERNATIVE, estimate = ESTIMATE, parameter = PARAMETER), class = "htest") 
    } else {
        structure(list(method = METHOD, data.name = DNAME, statistic = STATISTIC, p.value = P.VALUE, 
                       alternative = ALTERNATIVE, estimate = ESTIMATE), class = "htest") 
    }
}
