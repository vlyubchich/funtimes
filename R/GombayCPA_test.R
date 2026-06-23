#' Change Point Detection in Autoregressive Time Series
#'
#' The function detects change points in autoregressive (AR) models for time series. Changes
#' can be detected in any of \code{p + 2} (mean, var, phi) autoregressive parameters where \code{p}
#' is the order of the AR model. The test statistic is based on the efficient score vector \insertCite{Gombay_2008}{funtimes}.
#'
#' @details The function tests for a temporary change and a change in specific model parameters.
#' Critical values can be estimated via asymptotic distribution \code{"asymptotic"} (i.e., the
#' default option) or sieve bootstrap \code{"bootstrap"}. The function employs internal
#' function \code{change.point} and sieve bootstrap \code{change.point.sieve} function.
#'
#'
#' @param y a vector that contains univariate time-series observations. Missing values are not allowed.
#' @param a.order order of the autoregressive model which must be a non-negative integer number.
#' @param alternatives a string parameter that specifies a type of the test (i.e., "two-sided",
#' "greater", "lesser", and "temporary").  The option "temporary" examines the temporary change
#' in one of the parameters \insertCite{Gombay_2008}{funtimes}.
#' @param crit.type method of obtaining critical values: "asymptotic" (default) or "bootstrap".
#' @param num.bootstrap number of bootstrap replications if \code{crit.type = "bootstrap"}.
#' The default number is 1000.
#'
#'
#' @return A list with the following components:
#' \item{index}{points of change for each parameter. The value of the \code{"alternatives"}
#'  determines the return:
#' "temporary" -- returns max, min, and abs.max points;
#' "greater" -- returns max points;
#' "lesser" --  returns min points;
#' "two-sided" -- returns abs.max.}
#' \item{stats}{test statistic values for change points in mean, var, phi.}
#' \item{p.values}{\code{p-value} of the change point test.}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{mcusum_test}}  change point test for regression and
#' \code{\link[Ecdat]{terrorism}} dataset used in the Example 2
#'
#' @keywords changepoint ts
#'
#' @author Palina Niamkova, Dorcas Ofori-Boateng, Yulia R. Gel
#'
#' @export
#' @examples
#' \dontrun{
#' #Example 1:
#'
#' #Simulate some time series:
#' series_1 = arima.sim(n = 100, list(order = c(2,0,0), ar = c(-0.7, -0.1)))
#' series_2 = arima.sim(n = 200, list(order = c(2,0,0), ar = c(0.1, -0.6)))
#' main_series = c(series_1, series_2)
#'
#' result11 = GombayCPA_test(series_1, 2, "two-sided")
#' result11 #== No change point ===#
#'
#' result12 = GombayCPA_test(main_series, 2, "two-sided")
#' result12  #=== One change at phi values ===#
#'
#' result13 = GombayCPA_test(main_series, 2, "two-sided", "bootstrap")
#' result13  #=== One change at phi values ===#
#'
#'
#'
#' #Example 2:
#'
#' #From the package 'Ecdat' consider a time series with annual world number of victims of
#' #terrorism in the US from 1970 till 2016:
#' c.data = Ecdat::terrorism['nkill.us']
#' nkill.us.ts <- ts(c.data, start = 1970, end = 2016)
#'
#' #Now perform a change point detection with one sided tests:
#' GombayCPA_test(nkill.us.ts, 0, "lesser")
#' GombayCPA_test(nkill.us.ts, 0, "greater")
#' nkill.us.ts[32]
#' year = 1970 + 31
#' print(year)
#' plot(nkill.us.ts)
#'
#' #In both cases we find that the change point is located at the position 31 or 32. We can
#' # examine it further by checking the value of this position (using: nkill.us.ts[32]) as well as
#' # by plotting the graph (using: plot(nkill.us.ts)). The detected change point corresponds to
#' #the year of 2001, when the 9/11 attack happened.
#' }
#'
.change_point_test <- function(y, a.order, alternatives = c("two-sided", "greater", "lesser", "temporary"), res = NULL) {
  if (a.order < 0) stop("a.order must be greater than or equal to 0.")
  n <- length(y)
  
  #=== demeaning ===#
  y.dem <- y - mean(y, na.rm = TRUE)
  ar.p <- arima0(y.dem, order = c(a.order, 0, 0), include.mean = TRUE, method = "ML")
  
  mu <- mean(y.dem, na.rm = TRUE) #=== this is equal to zero ===#
  sigma2 <- if (is.null(res)) ar.p$sigma2 else var(res)
  
  phi <- if (a.order > 0) ar.p$coef else NULL
  
  if (a.order > 0) {
    gamma <- matrix(NA, ncol = a.order, nrow = a.order)
    cov <- acf(y, plot = FALSE, type = "covariance")$acf[1:a.order]
    rownum <- rep(c(1:a.order), a.order)
    colnum <- rep(c(1:a.order), each = a.order)
    abdiff <- abs(rownum - colnum) + 1
    gamma <- matrix(cov[abdiff], ncol = a.order, nrow = a.order)
  }
  
  muvec <- double(n)
  muvec[1] <- y.dem[1]
  for (i in 2:n) {
    if (a.order > 0) {
      len <- min(a.order, (i - 1))
      muvec[i] <- y.dem[i] - sum(phi[1:len] * y.dem[(i - 1):(i - len)])
    } else {
      muvec[i] <- y.dem[i]
    }
  }
  
  mustat <- if (a.order > 0) {
    ((1 - sum(phi)) / sigma2) * cumsum(muvec)
  } else {
    1 / sigma2 * cumsum(muvec)
  }
  
  sigmastat <- -c(1:n) / (2 * sigma2) + 1 / (2 * sigma2^2) * cumsum((muvec^2))
  
  allstat <- if (a.order > 0) {
    phistat <- matrix(NA, nrow = a.order, ncol = n)
    for (s in 1:a.order) {
      phivec.s <- muvec * c(rep(0, s), y.dem[1:(n - s)])
      phistat.s <- 1 / sigma2 * cumsum(phivec.s)
      phistat[s, ] <- phistat.s
    }
    rbind(mustat, sigmastat, phistat)
  } else {
    rbind(mustat, sigmastat)
  }
  
  info <- matrix(0, nrow = (a.order + 2), ncol = (a.order + 2))
  if (a.order > 0) {
    info[c(3:(a.order + 2)), c(3:(a.order + 2))] <- 1 / sigma2 * gamma
    info[1, 1] <- (1 / sigma2) * (1 - sum(phi))^2
  } else {
    info[1, 1] <- 1 / sigma2
  }
  info[2, 2] <- 1 / (2 * (sigma2^2))
  
  svd.info <- svd(info)
  info.n <- svd.info$u %*% diag(svd.info$d^(-1/2)) %*% t(svd.info$v)
  
  b <- n^(-1/2) * info.n %*% allstat
  
  param_names <- if (a.order > 0) c("mean", "var", paste("phi", 1:a.order, sep = "")) else c("mean", "var")
  
  if (alternatives == "two-sided") {
    maxstats <- apply(abs(b), 1, max)
    points <- apply(abs(b[, -1]), 1, which.max) + 1
    stat_names <- "abs.max"
  } else if (alternatives == "greater") {
    maxstats <- apply(b, 1, max)
    points <- apply(b[, -1], 1, which.max) + 1
    stat_names <- "max"
  } else if (alternatives == "lesser") {
    maxstats <- apply(b, 1, min)
    points <- apply(b[, -1], 1, which.min) + 1
    stat_names <- "min"
  } else { # temporary
    supstats <- apply(b, 1, max)
    infstats <- apply(b, 1, min)
    maxstats <- supstats - infstats
    
    maxpoints <- apply(b[, -1], 1, which.max) + 1
    minpoints <- apply(b[, -1], 1, which.min) + 1
    abspoints <- apply(abs(b[, -1]), 1, which.max) + 1
    
    points <- rbind(maxpoints, minpoints, abspoints)
    rownames(points) <- c("max", "min", "abs.max")
    colnames(points) <- param_names
    
    maxstats <- rbind(maxstats)
    rownames(maxstats) <- "diff"
    colnames(maxstats) <- param_names
    
    return(list(stats = maxstats, points = points))
  }
  
  maxstats <- rbind(maxstats)
  rownames(maxstats) <- stat_names
  colnames(maxstats) <- param_names
  
  points <- rbind(points)
  rownames(points) <- stat_names
  colnames(points) <- param_names
  
  return(list(stats = maxstats, points = points))
}

GombayCPA_test <- function(y, a.order, alternatives = c("two-sided", "greater", "lesser", "temporary"), crit.type = c("asymptotic", "bootstrap"), num.bootstrap = 1000) {
  alternatives <- match.arg(alternatives)
  crit.type <- match.arg(crit.type)
  n <- length(y)
  
  # Original series analysis
  orig <- .change_point_test(y, a.order, alternatives)
  origstats <- orig$stats
  origpoints <- orig$points
  
  param_names <- if (a.order > 0) c("mean", "var", paste("phi", 1:a.order, sep = "")) else c("mean", "var")
  
  if (crit.type == "bootstrap") {
    y.dem <- y - mean(na.omit(y))
    ar.p <- arima0(y.dem, order = c(a.order, 0, 0), include.mean = TRUE, method = "ML")
    phi <- if (a.order > 0) ar.p$coef[1:a.order] else 0
    
    residuals <- na.omit(ar.p$residuals) - mean(na.omit(ar.p$residuals))
    
    bootstats <- matrix(NA, nrow = (a.order + 2), ncol = num.bootstrap)
    
    for (i in 1:num.bootstrap) {
      e <- sample(residuals, size = n + 100, replace = TRUE)
      
      y.sim <- if (a.order > 0) {
        arima.sim(list(order = c(a.order, 0, 0), ar = phi), n = n + 100, innov = e)[101:(n + 100)]
      } else {
        e[1:n] # Assuming e should be same length as y
      }
      bootstats[, i] <- .change_point_test(y.sim, a.order, alternatives, res = e)$stats
    }
    
    bootpvals <- double((a.order + 2))
    for (i in 1:(a.order + 2)) {
      if (alternatives == "two-sided") {
        bootpvals[i] <- mean(abs(bootstats[i, ]) > abs(origstats[i]))
      } else if (alternatives == "greater") {
        bootpvals[i] <- mean(bootstats[i, ] > origstats[i])
      } else if (alternatives == "lesser") {
        bootpvals[i] <- mean(bootstats[i, ] < origstats[i])
      } else { # temporary
        bootpvals[i] <- mean(bootstats[i, ] > origstats[i])
      }
    }
    
    bootpvals <- rbind(bootpvals)
    rownames(bootpvals) <- "p.value"
    colnames(bootpvals) <- param_names
    return(list(index = origpoints, stats = origstats, p.values = bootpvals))
    
  } else { # asymptotic
    asympvals <- double((a.order + 2))
    for (i in 1:(a.order + 2)) {
      SQtest_stat <- origstats[i]^2
      
      if (alternatives == "greater" || alternatives == "lesser") {
        # P(sup|W(t)| > x) approx exp(-2x^2)
        eqna <- exp(-2 * SQtest_stat)
        asympvals[i] <- round(eqna, 3)
      } else if (alternatives == "two-sided") {
        # More complex series for two-sided test
        k <- 1:100 # Truncate series for efficiency
        eqnb <- sum(((-1)^(k + 1)) * exp(-2 * SQtest_stat * k^2))
        asympvals[i] <- round(4 * eqnb, 3) # approximation
      } else { # temporary
        # P(sup W(t) - inf W(t) > x)
        k <- 1:100 # Truncate series for efficiency
        eqnd <- sum(2 * ((4 * k^2 * SQtest_stat) - 1) * exp(-2 * k^2 * SQtest_stat))
        asympvals[i] <- round(1 - eqnd, 3)
      }
    }
    
    asympvals <- rbind(asympvals)
    rownames(asympvals) <- "p.value"
    colnames(asympvals) <- param_names
    return(list(index = origpoints, stats = origstats, p.values = asympvals))
  }
}
