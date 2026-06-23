#' HVK Estimator
#'
#' Estimate coefficients in nonparametric autoregression using the difference-based
#' approach by \insertCite{Hall_VanKeilegom_2003;textual}{funtimes}.
#'
#' @details First, autocovariances are estimated
#' using formula (2.6) by \insertCite{Hall_VanKeilegom_2003;textual}{funtimes}:
#' \deqn{\hat{\gamma}(0)=\frac{1}{m_2-m_1+1}\sum_{m=m_1}^{m_2}
#' \frac{1}{2(n-m)}\sum_{i=m+1}^{n}\{(D_mX)_i\}^2,}
#' \deqn{\hat{\gamma}(j)=\hat{\gamma}(0)-\frac{1}{2(n-j)}\sum_{i=j+1}^n\{(D_jX)_i\}^2,}
#' where \eqn{n} = \code{length(X)} is sample size, \eqn{D_j} is a difference operator
#' such that \eqn{(D_jX)_i=X_i-X_{i-j}}. Then, Yule--Walker method is used to
#' derive autoregression coefficients.
#'
#'
#' @param X univariate time series. Missing values are not allowed.
#' @param m1,m2 subsidiary smoothing parameters. Default
#' \code{m1 = round(length(X)^(0.1))}, \code{m2 = round(length(X)^(0.5))}.
#' @param ar.order order of the nonparametric autoregression (specified by user).
#'
#'
#' @return Vector of length \code{ar.order} with estimated autoregression coefficients.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link[stats]{ar}}, \code{\link{ARest}}
#'
#' @keywords ts
#'
#' @author Yulia R. Gel, Vyacheslav Lyubchich, Xingyu Wang
#'
#' @export
#' @examples
#' X <- arima.sim(n = 300, list(order = c(1, 0, 0), ar = c(0.6)))
#' HVK(as.vector(X), ar.order = 1)
#'
HVK <- function(X, m1 = NULL, m2 = NULL, ar.order = 1) {
    # 1. Input Validation
    if (!is.numeric(X) || !is.vector(X)) {
        stop("input object should be a numeric vector.")
    }
    if (anyNA(X)) {
        stop("input vector should not contain missing values.")
    }

    n <- length(X)
    if (n < 5) stop("input vector is too short for HVK estimation.")

    if (length(ar.order) != 1L || is.na(ar.order))
        stop("ar.order must be a single non-missing value.")
    ar.order <- as.integer(ar.order)
    if (ar.order < 1 || ar.order >= n)
        stop("ar.order must satisfy 1 <= ar.order < length(X).")

    if (is.null(m1) || is.null(m2)) {
        m1 <- as.integer(round(n^(0.1)))
        m2 <- as.integer(round(n^(0.5)))
    } else {
        m1 <- as.integer(m1)
        m2 <- as.integer(m2)
    }
    if (anyNA(c(m1, m2)) || m1 < 1 || m2 < m1 || m2 >= n)
        stop("m1 and m2 must satisfy 1 <= m1 <= m2 < length(X).")

    # 2. Autocovariance at Lag 0
    m_lags <- m1:m2
    sum_sq_diffs_m <- sapply(m_lags, function(m) {
        # High-speed vector subtraction without calling diff()
        d <- X[(m + 1):n] - X[1:(n - m)]
        sum(d^2)
    })
    autocovar0 <- sum(sum_sq_diffs_m / (2 * (n - m_lags))) / (m2 - m1 + 1)

    # 3. Autocovariances for Lags 1 to ar.order
    lags_j <- 1:ar.order

    # Micro-optimization for the common ar.order = 1 case
    if (ar.order == 1L) {
        sum_sq_diffs_j <- sum((X[2:n] - X[1:(n - 1)])^2)
    } else {
        sum_sq_diffs_j <- sapply(lags_j, function(j) {
            d <- X[(j + 1):n] - X[1:(n - j)]
            sum(d^2)
        })
    }

    autocovarj <- autocovar0 - sum_sq_diffs_j / (2 * (n - lags_j))
    autocovar <- c(autocovar0, autocovarj)

    # 4. Create Toeplitz matrix
    G <- stats::toeplitz(autocovar[1:ar.order])

    # 5. Linear system solving: solve(A, b) instead of solve(A) %*% b
    coeffs <- solve(G, autocovar[2:(ar.order + 1)])

    return(as.vector(coeffs))
}