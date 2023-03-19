#' Cross-Correlation of Autocorrelated Time Series
#'
#' Account for possible autocorrelation of time series when assessing the statistical significance
#' of their cross-correlation. A sieve bootstrap approach is used to generate multiple copies
#' of the time series with the same autoregressive dependence, under the null hypothesis of the
#' two time series under investigation being uncorrelated. The significance of cross-correlation
#' coefficients is assessed based on the distribution of their bootstrapped counterparts.
#' Both Pearson and Spearman types of coefficients are obtained, but a plot is provided for
#' only one type, with significant correlations shown using filled circles (see Examples).
#'
#' @details
#' Note that the smoothing of confidence bands is implemented purely for the look.
#' This smoothing is different from the
#' smoothing methods that can be applied to adjust bootstrap performance
#' \insertCite{DeAngelis_Young_1992}{funtimes}.
#' For correlations close to the significance bounds, the setting of \code{smooth} might
#' affect the decision on the statistical significance.
#' In this case, it is recommended to keep \code{smooth = FALSE} and set a higher \code{B}.
#'
#' @param x,y univariate numeric time-series objects or numeric vectors for which to
#' compute cross-correlation. Different time attributes in \code{ts} objects are
#' acknowledged, see Example 2 below.
#' @param lag.max maximum lag at which to calculate the cross-correlation. Will be
#' automatically limited as in \code{\link[stats]{ccf}}.
#' @param plot choose whether to plot results for Pearson correlation (default, or use
#' \code{plot = "Pearson"}), Spearman correlation (use \code{plot = "Spearman"}), or
#' suppress plotting (use \code{plot = "none"}). Both Pearson's and Spearman's results are
#' given in the output, regardless of the \code{plot} setting.
#' @param level confidence level, from 0 to 1. Default is 0.95, that is, 95% confidence.
#' @param B number of bootstrap simulations to obtain empirical critical values.
#' Default is 1000.
#' @param smooth logical value indicating whether the bootstrap confidence bands
#' should be smoothed across lags.
#' Default is \code{FALSE} meaning no smoothing.
#' @inheritParams causality_pred
#' @param ... other parameters passed to the function \code{\link{ARest}} to control
#' how autoregressive dependencies are estimated. The same set of parameters is used
#' separately on \code{x} and \code{y}.
#'
#'
#' @return A data frame with the following columns:
#' \item{Lag}{lags for which the following values were obtained.}
#' \item{r_P}{observed Pearson correlations.}
#' \item{lower_P, upper_P}{lower and upper confidence bounds (for the confidence level set by \code{level}) for Pearson correlations.}
#' \item{r_S}{observed Spearman correlations.}
#' \item{lower_S, upper_S}{lower and upper confidence bounds (for the confidence level set by \code{level}) for Spearman correlations.}
#'
#'
#' @seealso \code{\link{ARest}}, \code{\link[stats]{ar}}, \code{\link[stats]{ccf}},
#' \code{\link{HVK}}
#'
#' @keywords ts
#'
#' @author Vyacheslav Lyubchich
#'
#' @importFrom grDevices adjustcolor
#' @importFrom graphics grid lines matplot points polygon
#' @importFrom stats ccf
#' @export
#' @examples
#' \dontrun{
#' # Fix seed for reproducible simulations:
#' set.seed(1)
#'
#' # Example 1
#' # Simulate independent normal time series of same lengths
#' x <- rnorm(100)
#' y <- rnorm(100)
#' # Default CCF with parametric confidence band
#' ccf(x, y)
#' # CCF with bootstrap
#' tmp <- ccf_boot(x, y)
#' # One can extract results for both Pearson and Spearman correlations
#' tmp$rP
#' tmp$rS
#'
#' # Example 2
#' # Simulated ts objects of different lengths and starts (incomplete overlap)
#' x <- arima.sim(list(order = c(1, 0, 0), ar = 0.5), n = 30)
#' x <- ts(x, start = 2001)
#' y <- arima.sim(list(order = c(2, 0, 0), ar = c(0.5, 0.2)), n = 40)
#' y <- ts(y, start = 2020)
#' # Show how x and y are aligned
#' ts.plot(x, y, col = 1:2, lty = 1:2)
#' # The usual CCF
#' ccf(x, y)
#' # CCF with bootstrap confidence intervals
#' ccf_boot(x, y, plot = "Spearman")
#' # Notice that only +-7 lags can be calculated in both cases because of the small
#' # overlap of the time series. If we save these time series as plain vectors, the time
#' # information would be lost, and the time series will be misaligned.
#' ccf(as.numeric(x), as.numeric(y))
#'
#' # Example 3
#' # Box & Jenkins time series of sales and a leading indicator, see ?BJsales
#' plot.ts(cbind(BJsales.lead, BJsales))
#' # Each of the BJ time series looks as having a stochastic linear trend, so apply differences
#' plot.ts(cbind(diff(BJsales.lead), diff(BJsales)))
#' # Get cross-correlation of the differenced series
#' ccf_boot(diff(BJsales.lead), diff(BJsales), plot = "Spearman")
#' # The leading indicator "stands out" with significant correlations at negative lags,
#' # showing it can be used to predict the sales 2-3 time steps ahead (that is,
#' # diff(BJsales.lead) at times t-2 and t-3 is strongly correlated with diff(BJsales) at
#' # current time t).
#' }
#'
ccf_boot <- function(x,
                     y,
                     lag.max = NULL,
                     plot = c("Pearson", "Spearman", "none"),
                     level = 0.95,
                     B = 1000,
                     smooth = FALSE,
                     cl = 1L,
                     ...)
{
    ### Perform various checks
    namex <- deparse(substitute(x))[1L]
    namey <- deparse(substitute(y))[1L]
    if (is.matrix(x) || is.matrix(y)) {
        stop("x and y should be univariate time series only.")
    }
    if (any(is.na(x)) || any(is.na(y))) {
        stop("data should not contain missing values.")
    }
    nx <- length(x)
    ny <- length(y)
    B <- as.integer(B)
    if (B <= 0) {
        stop("number of bootstrap resamples B must be positive.")
    }
    plt <- match.arg(plot)
    bootparallel <- FALSE
    if (is.list(cl)) { #some other cluster supplied; use it but do not stop it
        bootparallel <- TRUE
        clStop <- FALSE
    } else {
        if (is.null(cl)) {
            cores <- parallel::detectCores()
        } else {
            cores <- cl
        }
        if (cores > 1) { #specified or detected cores>1; start a cluster and later stop it
            bootparallel <- TRUE
            cl <- parallel::makeCluster(cores)
            clStop <- TRUE
        }
    }
    ### Function
    xrank <- rank(x)
    yrank <- rank(y)
    attributes(xrank) <- attributes(x)
    attributes(yrank) <- attributes(y)
    tmp <- ccf(x, y, lag.max = lag.max, plot = FALSE)
    lags <- tmp$lag[,1,1]
    rP <- tmp$acf[,1,1]
    rS <- ccf(xrank, yrank, lag.max = lag.max, plot = FALSE)$acf[,1,1]
    phetax <- ARest(x, ...)
    phetay <- ARest(y, ...)
    if (length(phetax) > 0) {
        names(phetax) <- paste0(rep("phi_", length(phetax)), c(1:length(phetax)))
        tmp <- stats::filter(x, phetax, sides = 1)
        Zx <- x[(length(phetax) + 1):nx] - tmp[length(phetax):(nx - 1)]
    } else {
        Zx <- x
    }
    Zx <- Zx - mean(Zx)
    if (length(phetay) > 0) {
        names(phetay) <- paste0(rep("phi_", length(phetay)), c(1:length(phetay)))
        tmp <- stats::filter(y, phetay, sides = 1)
        Zy <- y[(length(phetay) + 1):ny] - tmp[length(phetay):(ny - 1)]
    } else {
        Zy <- y
    }
    Zy <- Zy - mean(Zy)
    ### Bootstrap
    if (bootparallel) {
        CCFs <- parallel::parSapply(cl, 1:B, function(b) {
            xboot <- arima.sim(list(order = c(length(phetax), 0, 0), ar = phetax), n = nx,
                               innov = sample(Zx, size = nx, replace = TRUE))
            yboot <- arima.sim(list(order = c(length(phetay), 0, 0), ar = phetay), n = ny,
                               innov = sample(Zy, size = ny, replace = TRUE))
            xrankboot <- rank(xboot)
            yrankboot <- rank(yboot)
            attributes(xboot) <- attributes(xrankboot) <- attributes(x)
            attributes(yboot) <- attributes(yrankboot) <- attributes(y)
            rPboot <- ccf(xboot, yboot, lag.max = lag.max, plot = FALSE)$acf[,1,1]
            rSboot <- ccf(xrankboot, yrankboot, lag.max = lag.max, plot = FALSE)$acf[,1,1]
            cbind(rPboot, rSboot)
        }, simplify = "array")
        if (clStop) {
            parallel::stopCluster(cl)
        }
    } else {
        CCFs <- sapply(1:B, function(b) {
            xboot <- arima.sim(list(order = c(length(phetax), 0, 0), ar = phetax), n = nx,
                               innov = sample(Zx, size = nx, replace = TRUE))
            yboot <- arima.sim(list(order = c(length(phetay), 0, 0), ar = phetay), n = ny,
                               innov = sample(Zy, size = ny, replace = TRUE))
            xrankboot <- rank(xboot)
            yrankboot <- rank(yboot)
            attributes(xboot) <- attributes(xrankboot) <- attributes(x)
            attributes(yboot) <- attributes(yrankboot) <- attributes(y)
            rPboot <- ccf(xboot, yboot, lag.max = lag.max, plot = FALSE)$acf[,1,1]
            rSboot <- ccf(xrankboot, yrankboot, lag.max = lag.max, plot = FALSE)$acf[,1,1]
            cbind(rPboot, rSboot)
        }, simplify = "array")
    } # end sequential bootstrap
    #CCFs has dimensions of nlags * 2 (Pearson and Spearman) * B
    ### Confidence regions
    alpha <- 1 - level
    # Percentile
    # crP <- apply(CCFs[,1,], 1, quantile, probs = c(alpha/2, 1 - alpha / 2))
    # crS <- apply(CCFs[,2,], 1, quantile, probs = c(alpha/2, 1 - alpha / 2))
    # Normal
    crP <- apply(CCFs[,1,], 1, function(x) qnorm(1 - alpha / 2, sd = sd(x)))
    if (smooth) {
        crP <- loess(crP ~ lags, span = 0.25)$fitted
    }
    crP <- rbind(-crP, crP)
    crS <- apply(CCFs[,2,], 1, function(x) qnorm(1 - alpha / 2, sd = sd(x)))
    if (smooth) {
        crS <- loess(crS ~ lags, span = 0.25)$fitted
    }
    crS <- rbind(-crS, crS)
    ### p-values
    # pP <- sapply(1:dim(CCFs)[1L], function(i) mean(abs(CCFs[i,1,]) > abs(rP[i])))
    # pS <- sapply(1:dim(CCFs)[1L], function(i) mean(abs(CCFs[i,2,]) > abs(rS[i])))
    RESULT <- data.frame(Lag = lags,
                         r_P = rP, #p_P = pP,
                         lower_P = crP[1,], upper_P = crP[2,], #Pearson
                         r_S = rS, #p_S = pS,
                         lower_S = crS[1,], upper_S = crS[2,]) #Spearman
    ### Plotting
    if (plt == "Pearson") {
        TMP <- RESULT[,grepl("_P", names(RESULT))]
    }
    if (plt == "Spearman") {
        TMP <- RESULT[,grepl("_S", names(RESULT))]
    }
    if (plt == "Pearson" || plt == "Spearman") {
        matplot(lags, TMP, type = "n",
                xlab = "Lag", ylab = "CCF",
                main = paste0(plt, " correlation of ", namex, "(t + Lag)", " and ", namey, "(t)\n",
                              "with ", level*100, "% bootstrap confidence region"),
                las = 1)
        grid(nx = 2, ny = NULL, lty = 1)
        polygon(x = c(lags, rev(lags)),
                y = c(TMP[,2], rev(TMP[,3])),
                col =  adjustcolor("deepskyblue", alpha.f = 0.80),
                border = NA)
        lines(lags, TMP[,1], type = "h")
        isoutside <- (TMP[,1] < TMP[,2]) | (TMP[,3] < TMP[,1])
        points(lags, TMP[,1], pch = c(1, 16)[1 + isoutside])
        return(invisible(RESULT))
    } else {
        return(RESULT)
    }
}
