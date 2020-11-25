#' WAVK Trend Test
#' 
#' Non-parametric test to detect (non-)monotonic parametric trends in time series
#' \insertCite{@based on @Lyubchich_etal_2013_wavk}{funtimes}.
#' 
#' @details See more details in \insertCite{Lyubchich_Gel_2016_synchronism;textual}{funtimes} 
#' and \insertCite{Lyubchich_2016_trends;textual}{funtimes}.
#' 
#' 
#' @param formula an object of class "\code{\link[stats]{formula}}", specifying the 
#' form of the parametric time trend to be tested. Variable \eqn{t} should be used 
#' to specify the form, where \eqn{t} is specified within the function as a regular 
#' sequence on the interval (0,1]. See `Examples'.
#' @param factor.length method to define the length of local windows (factors). 
#' Default option\cr \code{"user.defined"} allows to set only one value of the argument 
#' \code{Window}. The option \code{"adaptive.selection"} sets \code{method = "boot"} 
#' and employs heuristic \eqn{m}-out-of-\eqn{n} subsampling algorithm 
#' \insertCite{Bickel_Sakov_2008}{funtimes} to select an optimal window from the set 
#' of possible windows \code{length(x)*q^j} whose values are mapped to the largest 
#' previous integer and greater than 2. Vector \code{x} is the time series tested.
#' @param Window length of the local window (factor), default is 
#' \code{round(0.1*length(x))}, where \code{x} is the time series tested. 
#' This argument is ignored if\cr \code{factor.length = "adaptive.selection"}.
#' @param q scalar from 0 to 1 to define the set of possible windows when 
#' \code{factor.length =} \code{"adaptive.selection"}. Default is \eqn{3/4}. 
#' This argument is ignored if\cr \code{factor.length =} \code{"user.defined"}.
#' @param j numeric vector to define the set of possible windows when 
#' \code{factor.length =} \code{"adaptive.selection"}. Default is \code{c(8:11)}. 
#' This argument is ignored if\cr \code{factor.length = "user.defined"}.
#' @param B number of bootstrap simulations to obtain empirical critical values. 
#' Default is 1000.
#' @param method method of obtaining critical values: from asymptotical (\code{"asympt"}) 
#' or bootstrap (\code{"boot"}) distribution. 
#' If \code{factor.length =} \code{"adaptive.selection"} the option \code{"boot"} is used.
#' @inheritParams ARest
#' @param out logical value indicates whether full output should be shown. 
#' Default is \code{FALSE}.
#' 
#' 
#' @return A list with class \code{"htest"} containing the following components:
#' \item{method}{name of the method.}
#' \item{data.name}{name of the data.}
#' \item{statistic}{value of the test statistic.}
#' \item{p.value}{\eqn{p}-value of the test.}
#' \item{alternative}{alternative hypothesis.}
#' \item{parameter}{window that was used.}
#' \item{estimate}{list with the following elements: estimated trend coefficients; 
#' user-defined or BIC-selected AR order; estimated AR coefficients; and, 
#' if \code{factor.length =} \code{"adaptive.selection"}, 
#' test results for all considered windows.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link[stats]{ar}}, \code{\link{HVK}}, \code{\link{WAVK}}, 
#' \code{\link{sync.test}}
#' 
#' @keywords htest ts trend
#' 
#' @author Yulia R. Gel, Vyacheslav Lyubchich, Ethan Schaeffer
#' 
#' @export
#' @examples
#' # Fix seed for reproducible simulations:
#' set.seed(1)
#' 
#' #Simulate autoregressive time series of length n with smooth quadratic trend:
#' n <- 100
#' tsTrend <- 1 + 2*(1:n/n) + 4*(1:n/n)^2
#' tsNoise <- arima.sim(n = n, list(order = c(2, 0, 0), ar = c(-0.7, -0.1)))
#' U <- tsTrend + tsNoise
#' plot.ts(U)
#' 
#' #Test H0 of a linear trend, with m-out-of-n selection of the local window:
#' \dontrun{
#'     wavk.test(U ~ t, factor.length = "adaptive.selection")}
#' # Sample output:
#' ##	Trend test by Wang, Akritas, and Van Keilegom (bootstrap p-values)
#' ##
#' ##data:  U 
#' ##WAVK test statistic = 5.3964, adaptively selected window = 4, p-value < 2.2e-16
#' ##alternative hypothesis: trend is not of the form U ~ t.
#' 
#' #Test H0 of a quadratic trend, with m-out-of-n selection of the local window 
#' #and output of all results:
#' \dontrun{
#'     wavk.test(U ~ poly(t, 2), factor.length = "adaptive.selection", out = TRUE)}
#' # Sample output:
#' ##	Trend test by Wang, Akritas, and Van Keilegom (bootstrap p-values)
#' ##
#' ##data:  U 
#' ##WAVK test statistic = 0.40083, adaptively selected window = 4, p-value = 0.576
#' ##alternative hypothesis: trend is not of the form U ~ poly(t, 2).
#' ##sample estimates:
#' ##$trend_coefficients
#' ##(Intercept) poly(t, 2)1 poly(t, 2)2 
#' ##   3.408530   17.681422    2.597213 
#' ##
#' ##$AR_order
#' ##[1] 1
#' ##
#' ##$AR_coefficients
#' ##         phi_1 
#' ##[1] -0.7406163
#' ##
#' ##$all_considered_windows
#' ## Window WAVK-statistic p-value
#' ##      4     0.40083181   0.576
#' ##      5     0.06098625   0.760
#' ##      7    -0.57115451   0.738
#' ##     10    -1.02982929   0.360
#' 
#' # Test H0 of no trend (constant trend) using asymptotic distribution of statistic.
#' wavk.test(U ~ 1, method = "asympt")
#' # Sample output:
#' ##	Trend test by Wang, Akritas, and Van Keilegom (asymptotic p-values)
#' ##
#' ##data:  U 
#' ##WAVK test statistic = 25.999, user-defined window = 10, p-value < 2.2e-16
#' ##alternative hypothesis: trend is not of the form U ~ 1.
#' 
wavk.test <- function(formula, factor.length = c("user.defined", "adaptive.selection"), 
                      Window = NULL, q = 3/4, j = c(8:11), B = 1000, method = c("boot", "asympt"), 
                      ar.order = NULL, ar.method = "HVK", BIC = TRUE, out = FALSE)
{
    ### Perform various checks.
    frml <- deparse(substitute(formula))
    splt <- strsplit(frml, "~")[[1]]
    DNAME <- splt[1]
    x <- lm(formula = as.formula(paste0(DNAME, "~ 1"), env = parent.frame(n = 4)), 
            method = "model.frame")[,1]
    if (is.null(Window)) {
        Window = round(0.1*length(x))
    }
    if (NCOL(x) > 1 | !is.numeric(x)) {
        stop("x is not a vector or univariate time series.")
    }
    n <- length(x)
    if (any(is.na(x))) {
        stop("x contains missing values.")
    }
    factor.length <- match.arg(factor.length)
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
        kn <- length(x)*q^j
    }
    kn <- unique(sort(floor(kn)))
    kn <- kn[kn > 2 & kn < n]
    if (length(kn) == 0) {
        stop("set a proper window.")
    }
    if (factor.length == "adaptive.selection" & length(kn) < 3) {
        stop("number of possible windows is not enough for adaptive selection. Change parameters 'q' and/or 'j'.")
    }
    if (factor.length == "adaptive.selection") {
        method <- "boot"
    }
    B <- round(B)
    if (B <= 0) {
        stop("number of bootstrap samples B must be positive.")
    }
    method <- match.arg(method)
    if (!is.null(ar.order) & (NCOL(ar.order) > 1 | !is.numeric(ar.order) | NROW(ar.order) > 1)) {
        stop("ar.order is not a scalar.")
    }
    if (!is.null(ar.order) && ar.order < 0) {
        stop("ar.order must be non-negative.")
    }
    ### Function.
    t <- c(1:n)/n
    result <- matrix(NA, length(kn), 2)
    res <- matrix(NA, 1, 2)
    if (is.null(ar.order)) {
        ar.order <- floor(10*log10(n))
    }
    mod <- lm(as.formula(as.formula(paste0("x ~ ", splt[2]))))
    TrendCoeff <- mod$coefficients
    TS <- as.vector(mod$residuals)
    ALTERNATIVE <- paste("trend is not of the form ", frml, ".", sep = "")
    pheta <- ARest(TS, ar.order = ar.order, ar.method = ar.method, BIC = BIC)
    if (length(pheta) > 0) {
        names(pheta) <- paste(rep("phi_", length(pheta)), c(1:length(pheta)), sep = "")
        tmp <- filter(x, pheta, sides = 1)
        tmp2 <- filter(mod$fitted.values, pheta, sides = 1)
        Z <- (x[(length(pheta) + 1):n] - tmp[length(pheta):(n - 1)]) - 
            (mod$fitted.values[(length(pheta) + 1):n] - tmp2[length(pheta):(n - 1)])
    } else {
        Z <- TS
    }
    ESTIMATE <- list(TrendCoeff, length(pheta), pheta)
    names(ESTIMATE) <- c("trend_coefficients", "AR_order", "AR_coefficients")
    Z <- Z - mean(Z)
    sigma <- sqrt(sum(diff(Z)^2)/(2*(length(Z) - 1)))
    if (method == "asympt") {
        METHOD <- "Trend test by Wang, Akritas, and Van Keilegom (asymptotic p-values)"
        for (i in 1:length(kn)) {
            tmp <- WAVK(Z, kn[i])
            result[i,] <- c(tmp$Tns, tmp$p.value)
        }
        STATISTIC <- result[1,1]
        P.VALUE <- result[1,2]
        PARAMETER <- kn[1]
        names(PARAMETER) <- "user-defined window"
    } else {
        METHOD <- "Trend test by Wang, Akritas, and Van Keilegom (bootstrap p-values)"
        boot <- array(data = rnorm(n*B), c(n,B)) * sigma
        s <- array(data = NA, c(length(kn), B))
        for (i in 1:length(kn)) {
            s[i,] <- apply(boot, 2, function(x) WAVK(x, kn[i])$Tns)
            result[i, 1] <- WAVK(Z, kn[i])$Tns
            crit <- sum(result[i,1] < s[i,])/B
            if (crit < 0.5) {
                result[i, 2] <- 2*crit
            } else {
                result[i, 2] <- 2*(1 - crit)
            }
        }
        if (length(kn) < 3) {
            STATISTIC <- result[1,1]
            P.VALUE <- result[1,2]
            PARAMETER <- kn[1]
            names(PARAMETER) <- "user-defined window"
        } else {
            s <- t(apply(s, 1, sort))
            distance <- sapply(1:(length(kn) - 1), function(x) dist(s[x:(x + 1),]))
            kn_opt <- kn[which.min(distance)]
            res[1,] <- result[which.min(distance),]
            STATISTIC <- res[1,1]
            P.VALUE <- res[1,2]
            PARAMETER <- kn_opt
            names(PARAMETER) <- "adaptively selected window"
            if (out) {
                tmp <- names(ESTIMATE)
                ESTIMATE <- c(ESTIMATE, NA)
                names(ESTIMATE) <- c(tmp, "all_considered_windows")
                ESTIMATE[[length(ESTIMATE)]] <- cbind(kn, result)
                dimnames(ESTIMATE[[length(ESTIMATE)]]) <- list(rep("", length(kn)), 
                                                               c("Window", "WAVK-statistic", "p-value"))
            }
        }
    }
    names(STATISTIC) <- "WAVK test statistic"
    if (out) {
        structure(list(method = METHOD, data.name = DNAME, statistic = STATISTIC, 
                       p.value = P.VALUE, alternative = ALTERNATIVE, 
                       parameter = PARAMETER, estimate = ESTIMATE), class = "htest") 
    } else {
        structure(list(method = METHOD, data.name = DNAME, statistic = STATISTIC, 
                       p.value = P.VALUE, alternative = ALTERNATIVE, 
                       parameter = PARAMETER), class = "htest") 
    }
}
