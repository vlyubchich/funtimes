#' Beale's Estimator and Sample Size
#' 
#' Beale's ratio estimator \insertCite{Beale_1962}{funtimes} 
#' for estimating population total and  
#' confidence intervals, with an option of calculating sample size for a required
#' relative error (\code{p}) or margin of error (\code{d}).
#' 
#'
#' @param x a numeric vector with quantities of interest, such as river discharge 
#' per month. Missing values (\code{NA}) are allowed.
#' @param y a numeric vector with quantities of interest for which the total shall 
#' be estimated, such as total nutrient loads per month. 
#' Missing values (\code{NA}) are allowed. 
#' Lengths of \code{x} and \code{y} mush be the same.
#' @param level confidence level, from 0 to 1. 
#' Default is \code{0.95}, that is, 95% confidence.
#' @param N population size for which the estimate of the total \code{y} required.
#' By default, \code{length(x)} is used.
#' @param p optional argument specifying the required relative error, from 0 to 1,
#' for computing corresponding sample size. For example, \code{p = 0.15} defines
#' a 15% relative error.
#' @param d optional argument specifying the required margin of error
#' for computing corresponding sample size. If both \code{p} and \code{d} are specified,
#' only \code{p} is used.
#' @param verbose logical value defining whether the output should be printed out
#' in words. Default is set to \code{TRUE} to give such output.
#'
#' @return A list with the following components:
#' \item{estimate}{Beale's estimate of the population total for the variable \code{y}.}
#' \item{se}{standard error of the estimate.}
#' \item{CI}{a vector of length 2 with a confidence interval (lower and upper value) 
#' for the estimate.}
#' \item{level}{confidence level for the interval.}
#' \item{N}{population size.}
#' \item{n}{the actual sample size.}
#' \item{p}{the relative error used for sample size calculations. 
#' Reported only if \code{p} was specified in the input.}
#' \item{d}{the margin of error used for sample size calculations.
#' Reported only if \code{d} was specified and \code{p} was not specified in the input.}
#' \item{nhat}{estimated sample size for the given \code{level} and error 
#' (\code{p} or \code{d}).}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{vignette("beales", package = "funtimes")}
#' 
#' @keywords power sample ts
#' 
#' @author Vyacheslav Lyubchich, thanks to Dave Lorenz for pointing out an error in version 7 and below of the package
#' 
#' @export
#' @examples
#' #Some hypothetical data for monthly river discharge 
#' #and corresponding nutrient loads:
#' discharge <- c(NA, 50, 90, 100, 80, 90, 100, 90, 80, 70, NA, NA)
#' loads <- c(33, 22, 44, 48, NA, 44, 49, NA, NA, 36, NA, NA)
#' 
#' #Example 1:
#' #Estimate total annual load (12 months), 
#' #with 90% confidence intervals
#' beales(discharge, loads, level = 0.9)
#' 
#' #Example 2:
#' #Calculate sample size required for 90% confidence intervals 
#' #with a margin of error 30 units
#' beales(discharge, loads, level = 0.9, d = 30)
#' 
beales <- function(x, y, level = 0.95, N = NULL, p = NULL, d = NULL, verbose = TRUE){
    if (length(x) != length(y)) stop("Vectors 'x' and 'y' must be of the same length.")
    #Population and sample sizes:
    if (is.null(N)) N <- length(x)
    n <- sum(!is.na(y))
    if (n >= N || length(x) > N) stop("Population size 'N' must be bigger than the sample.")
    
    #Components of Beale's estimator:
    xbarprime <- mean(x, na.rm = TRUE)
    x[is.na(x)] <- xbarprime #if some x-values are missing, they are replaced with mean
    xbar <- mean(x[!is.na(y)])
    ybar <- mean(y, na.rm = TRUE)
    s2x <- var(x[!is.na(y)])
    X <- N * xbarprime
    theta <- 1/n - 1/N
    sxy <- cov(x, y, use = "complete.obs")
    
    #Beale's estimate of Y:
    Yhat <- X * ybar * (1 + theta*sxy/(xbar*ybar)) / (xbar * (1 + theta*s2x/xbar^2)) #output
    
    #Cumulants:
    i <- is.na(y)
    x <- x[!i]
    y <- y[!i]
    C11 <- sxy / (xbar * ybar)
    C20 <- var(x) / (xbar^2)
    C02 <- var(y) / (ybar^2)
    C21 <- (sum((x - xbar)^2*(y - ybar)) / (n - 1)) / (xbar^2 * ybar)
    C12 <- (sum((x - xbar)*(y - ybar)^2) / (n - 1)) / (xbar * ybar^2)
    C30 <- (sum((x - xbar)^3) / (n - 1)) / (xbar^3)
    a_hat <- 2*C20^2 - 4*C20*C11 + C11^2 + C20*C02
    b_hat <- C20 + C02 - 2*C11 + 2*(C30 - 2*C21 + C12)/N
    VarYhat <- X^2*ybar^2 * (theta * b_hat + theta^2 * a_hat ) / xbar^2
    #Confidence interval:
    z <- qnorm((1 - level)/2)
    seYhat <- sqrt(VarYhat) #output
    CI <- Yhat + c(1, -1) * z * seYhat #output
    if (verbose) {
        print(paste("Beale's estimate of the total (for population size ", N, ") is ", round(Yhat, 3), 
                    " with ", level*100, "% confidence interval from ", 
                    round(CI[1], 3), " to ", round(CI[2], 3), ".", sep = ""))
    }
    #Check if additional arguments are set, calculate sample size:
    if (is.null(p) && is.null(d)) {
        result <- list(estimate = Yhat,
                       se = seYhat,
                       CI = CI,
                       level = level,
                       N = N,
                       n = n)
    } else if (!is.null(p)) {#the relative error "p" is set by user
        c_hatprime <- p^2 / z^2
        theta_2 <- (-b_hat + sqrt(b_hat^2 + 4*a_hat*c_hatprime)) / (2*a_hat)
        n2 <- ceiling(1 / (theta_2 + 1/N))
        result <- list(estimate = Yhat,
                       se = seYhat,
                       CI = CI,
                       level = level,
                       N = N,
                       n = n,
                       p = p,
                       nhat = n2)
        if (verbose) {
            print(paste("To obtain a ", level*100, "% confidence interval with a relative error of ", 
                        p*100, "%, a sample of size ", n2, " is required.", sep = ""))
        }
    } else {#it means the margin of error "d" is set by user
        c_hat <- d^2 * xbar^2 / (z^2 * X^2 * ybar^2)
        theta_2 <- (-b_hat + sqrt(b_hat^2 + 4*a_hat*c_hat)) / (2*a_hat)
        n2 <- ceiling(1 / (theta_2 + 1/N) )
        result <- list(estimate = Yhat,
                       se = seYhat,
                       CI = CI,
                       level = level,
                       N = N,
                       n = n,
                       d = d,
                       nhat = n2)
        if (verbose) {
            print(paste("To obtain a ", level*100, "% confidence interval with a margin of error being ", 
                        round(d, 3), ", a sample of size ", n2, " is required.", sep = ""))
        }
    }
    return(result)
}
