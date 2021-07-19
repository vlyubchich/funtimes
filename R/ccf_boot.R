#' Cross-Correlation Function of Time Series with Sieve Bootstrap p-values
#' 
#' Account for possible autocorrelation of time series when assessing statistical significance 
#' of their cross-correlation. A sieve bootstrap approach is used to generate multiple copies
#' of the time series with the same autoregressive dependence, under the null hypothesis of the 
#' two time series under investigation being uncorrelated. Significance of cross-correlation
#' coefficients is assessed based on the distribution of their bootstrapped counterparts. 
#' Both Pearson and Spearman types of coefficients are obtained, but plot is provided for 
#' only one type, with significant correlations shown using filled circles.
#' 
#' 
#' @param x,y univariate numeric time series objects or numeric vectors for which to 
#' compute cross-correlation. Different time attributes in \code{ts} objects are 
#' acknowledged, see Example 2 below.
#' @param lag.max maximum lag at which to calculate the cross-correlation. Will be 
#' automatically limited as in \code{\link[stats]{ccf}}.
#' @param plot choose whether to plot results for Pearson correlation (default, or use
#' \code{plot = "Pearson"}), Spearman correlation (use \code{plot = "Spearman"}), or 
#' suppress plotting (use \code{plot = "none"}). Both Pearson and Spearman results are 
#' given in the output, irregardless of the \code{plot} setting.
#' @param level confidence level, from 0 to 1. Default is 0.95, that is, 95% confidence.
#' @param B number of bootstrap simulations to obtain empirical critical values. 
#' Default is 1000.
#' @param ... other parameters passed to the function \code{\link{ARest}} to control 
#' how autoregressive dependencies are estimated. The same set of parameters is used
#' separately on \code{x} and \code{y}.
#' 
#' 
#' @return A data frame with the following columns:
#' \item{Lag}{lags for which the following values were obtained.}
#' \item{rP}{observed Pearson correlations.}
#' \item{pP}{bootstrap p-value for Pearson correlations.}
#' \item{lowerP, upperP}{lower and upper quantiles (for the confidence level set by \code{level}) of the bootstrapped Pearson correlations.}
#' \item{rS}{observed Spearman correlations.}
#' \item{pS}{bootstrap p-value for Spearman correlations.}
#' \item{lowerS, upperS}{lower and upper quantiles (for the confidence level set by \code{level}) of the bootstrapped Spearman correlations.}
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
#' ccf(x, y) # default CCF with parametric confidence band
#' tmp <- ccf_boot(x, y) # CCF with bootstrap
#' tmp$rP; tmp$rS # can always extract results for both Pearson and Spearman correlations
#' 
#' # Example 2
#' # Simulated ts objects of different lengths and starts (incomplete overlap)
#' x <- arima.sim(list(order = c(1, 0, 0), ar = 0.5), n = 30)
#' x <- ts(x, start = 2001)
#' y <- arima.sim(list(order = c(2, 0, 0), ar = c(0.5, 0.2)), n = 40)
#' y <- ts(y, start = 2020)
#' ts.plot(x, y, col = 1:2, lty = 1:2) # show how x and y are aligned
#' ccf(x, y)
#' ccf_boot(x, y, plot = "Spearman") # CCF with bootstrap
#' # Notice that only +-7 lags can be calculated in both cases because of the small 
#' # overlap of the time series. If save these time series as plain vectors, the time
#' # information would be lost, and time series will be misaligned. 
#' ccf(as.numeric(x), as.numeric(y))
#' 
#' # Example 3
#' # Box & Jenkins time series of sales and a leading indicator, see ?BJsales
#' plot.ts(cbind(BJsales.lead, BJsales))
#' # Each of the BJ time series looks as having a stochastic linear trend, so apply differences:
#' plot.ts(cbind(diff(BJsales.lead), diff(BJsales)))
#' # Get cross-correlation of the differenced series:
#' ccf_boot(diff(BJsales.lead), diff(BJsales), plot = "Spearman")
#' # The leading indicator "stands out" with significant correlations at negative lags, 
#' # showing it can be used to predict the sales 2-3 time steps ahead (that is,
#' # diff(BJsales.lead) at times t-2 and t-3 is strongly correlated with diff(BJsales) at 
#' # current time t).
#' }
#' 
ccf_boot <- function(x, y, lag.max = NULL, 
                     plot = c("Pearson", "Spearman", "none"),
                     level = 0.95, B = 1000, ...)
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
    CCFs <- sapply(1:B, function(b) {
        xboot <- arima.sim(list(order = c(length(phetax), 0, 0), ar = phetax), n = nx, 
                           innov = sample(Zx, size = nx, replace = TRUE))
        yboot <- arima.sim(list(order = c(length(phetay), 0, 0), ar = phetay), n = ny, 
                           innov = sample(Zy, size = ny, replace = TRUE))
        attributes(xboot) <- attributes(x)
        attributes(yboot) <- attributes(y)
        xrankboot <- rank(xboot)
        yrankboot <- rank(yboot)
        attributes(xrankboot) <- attributes(x)
        attributes(yrankboot) <- attributes(y)
        rPboot <- ccf(xboot, yboot, lag.max = lag.max, plot = FALSE)$acf[,1,1]
        rSboot <- ccf(xrankboot, yrankboot, lag.max = lag.max, plot = FALSE)$acf[,1,1]
        cbind(rPboot, rSboot)
    }, simplify = "array")
    #CCFs has dimensions of nlags * 2 (Pearson and Spearman) * B
    ### Confidence regions
    alpha <- 1 - level
    crP <- apply(CCFs[,1,], 1, quantile, probs = c(alpha/2, 1 - alpha / 2))
    crS <- apply(CCFs[,2,], 1, quantile, probs = c(alpha/2, 1 - alpha / 2))
    ### p-values
    pP <- sapply(1:dim(CCFs)[1L], function(i) mean(abs(CCFs[i,1,]) > abs(rP[i])))
    pS <- sapply(1:dim(CCFs)[1L], function(i) mean(abs(CCFs[i,2,]) > abs(rS[i])))
    RESULT <- data.frame(Lag = lags,
                         rP = rP, pP = pP, lowerP = crP[1,], upperP = crP[2,], #Pearson
                         rS = rS, pS = pS, lowerS = crS[1,], upperS = crS[2,]) #Spearman
    ### Plotting
    if (plt == "Pearson") {
        TMP <- RESULT[,grepl("P", names(RESULT))]
    }
    if (plt == "Spearman") {
        TMP <- RESULT[,grepl("S", names(RESULT))]
    }
    if (plt == "Pearson" || plt == "Spearman") {
        matplot(lags, TMP[,-2], type = "n",
                xlab = "Lag", ylab = "CCF", 
                main = paste0(plt, " correlation of ", namex, "(t + Lag)", " and ", namey, "(t)\n",
                              "with ", level*100, "% bootstrapped confidence region"),
                las = 1)
        grid(nx = 2, ny = NULL, lty = 1)
        polygon(x = c(lags, rev(lags)),
                y = c(TMP[,3], rev(TMP[,4])),
                col =  adjustcolor("deepskyblue", alpha.f = 0.80),
                border = NA)
        lines(lags, TMP[,1], type = "h")
        points(lags, TMP[,1], pch = c(1, 16)[1 + (TMP[,2] < alpha)])
        return(invisible(RESULT))
    } else {
        return(RESULT)
    }
}
