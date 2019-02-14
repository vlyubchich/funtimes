#' HVK Estimator
#' 
#' Estimate coefficients in non-parametric autoregression using the difference-based 
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
#' @param ar.order order of the non-parametric autoregression (specified by user).
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
    if (!is.numeric(X) | !is.vector(X)) {
        stop("input object should be a vector.")
    }    
    if (any(is.na(X))) {
        stop("input vector should not contain missing values.")
    }
    n <- length(X)
    if(is.null(m1)|is.null(m2)){
        m1 <- round(n^(0.1))
        m2 <- round(n^(0.5))
    }
    m <- c(m1:m2)
    tmp <- sapply(m, function(x) diff(X, lag = x))
    tmp <- sapply(1:length(tmp), function(x) sum(tmp[[x]]^2))
    autocovar0 <- sum(tmp/(2*(n-m)))/(m2-m1+1)
    tmp <- sapply(c(1:ar.order), function(x) diff(X, lag = x))
    if(ar.order == 1){
        tmp <- sum(tmp^2)
    }else{
        tmp <- sapply(1:length(tmp), function(x) sum(tmp[[x]]^2))
    }
    autocovarj <- autocovar0 - tmp/(2*(n-c(1:ar.order)))
    autocovar <- c(autocovar0, autocovarj)
    G <- matrix(NA, ar.order, ar.order)
    for (i in 1:(ar.order)){
        for (j in 1:(ar.order)){
            G[i,j] <- autocovar[abs(i-j) + 1]
        }
    }
    coeffs <- solve(G)%*%autocovar[2:(ar.order + 1)]
    return(as.vector(coeffs))
}
