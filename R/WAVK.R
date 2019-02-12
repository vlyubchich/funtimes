#' WAVK Statistic
#' 
#' Computes statistic for testing the parametric form of a regression function, 
#' suggested by \insertCite{Wang_etal_2008;textual}{funtimes}.
#' 
#' 
#' @param z pre-filtered univariate time series 
#' \insertCite{@see formula (2.1) by @Wang_VanKeilegom_2007}{funtimes}:
#' \deqn{Z_i=\left(Y_{i+p}-\sum_{j=1}^p{\hat{\phi}_{j,n}Y_{i+p-j}} \right)-
#' \left( f(\hat{\theta},t_{i+p})-
#' \sum_{j=1}^p{\hat{\phi}_{j,n}f(\hat{\theta},t_{i+p-j})} \right),}
#' where \eqn{Y_i} is observed time series of length \eqn{n}, \eqn{\hat{\theta}} 
#' is an estimator of hypothesized parametric trend \eqn{f(\theta, t)}, 
#' and \eqn{\hat{\phi}_p=(\hat{\phi}_{1,n}, \ldots, \hat{\phi}_{p,n})'} 
#' are estimated coefficients of an autoregressive filter of order \eqn{p}. 
#' Missing values are not allowed.
#' 
#' @param kn length of the local window.
#' 
#' 
#' @return A list with following components:
#' \item{Tn}{test statistic based on artificial ANOVA and defined 
#' by \insertCite{Wang_VanKeilegom_2007;textual}{funtimes}
#' as a difference of mean square for treatments (MST) and mean square for errors (MSE):
#' \deqn{T_n= MST - MSE =\frac{k_{n}}{n-1} \sum_{t=1}^T 
#' \biggl(\overline{V}_{t.}-\overline{V}_{..}\biggr)^2 - 
#' \frac{1}{n(k_{n}-1)} \sum_{t=1}^n \sum_{j=1}^{k_{n}}\biggl(V_{tj}-\overline{V}_{t.}\biggr)^2,}
#' where \eqn{\{V_{t1}, \ldots, V_{tk_n}\}=\{Z_j: j\in W_{t}\}}, \eqn{W_t} is a 
#' local window, \eqn{\overline{V}_{t.}} and \eqn{\overline{V}_{..}} are the mean 
#' of the \eqn{t}th group and the grand mean, respectively.}
#' 
#' \item{Tns}{standardized version of \code{Tn} according to 
#' Theorem 3.1 by \insertCite{Wang_VanKeilegom_2007;textual}{funtimes}:
#' \deqn{T_{ns} = \left( \frac{n}{k_n} \right)^{\frac{1}{2}}T_n \bigg/  
#' \left(\frac{4}{3}\right)^{\frac{1}{2}} \sigma^2,}{Tns = Tn*(n/kn)^0.5 / (sigma^2 * (4/3)^0.5),}
#' where \eqn{n} is length and \eqn{\sigma^2}{sigma^2} is variance of the time series. 
#' Robust difference-based Rice's estimator \insertCite{Rice_1984}{funtimes} 
#' is used to estimate \eqn{\sigma^2}{sigma^2}.}
#' 
#' \item{p.value}{\eqn{p}-value for \code{Tns} based on its 
#' asymptotic \eqn{N(0,1)} distribution.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{wavk.test}}
#' 
#' @keywords ts trend
#' 
#' @author Yulia R. Gel, Vyacheslav Lyubchich
#' 
#' @export
#' @examples
#' z <- rnorm(300)
#' WAVK(z, kn = 7)
#' 
WAVK <- function(z, kn = NULL) 
{
    if (!is.numeric(z) | !is.vector(z)) {
        stop("input object should be a vector.")
    }    
    if (any(is.na(z))) {
        stop("input vector should not contain missing values.")
    }
    kn <- round(kn)
    T <- length(z)
    ave_group <- sapply(c(1:(T-kn+1)), function(x) mean(z[x:(x+kn-1)]))
    ave_all <- mean(ave_group)
    MST <- sum((ave_group - ave_all)^2) * kn / (T - 1)
    MSE <- sum(sapply(c(1:(T-kn+1)), function(x) 
        sum((z[x:(x+kn-1)] - ave_group[x])^2))) / (T*(kn-1))
    Tn <- MST - MSE
    sigma2 <- sum(diff(z)^2)/(2 * (T - 1))
    Tns <- sqrt(T/kn) * Tn / (sqrt(4/3) * sigma2)
    crit <- pnorm(Tns, mean = 0, sd = 1)
    if (crit < 0.5) {
        p.value <- crit * 2
    } else {
        p.value <- (1 - crit) * 2
    }
    list(Tn = Tn, Tns = Tns, p.value = p.value)
}
