#' Time Series Trend Synchronism Test
#' 
#' Non-parametric test for synchronism of parametric trends in multiple time series
#' \insertCite{Lyubchich_Gel_2016_synchronism}{funtimes}. 
#' The method tests whether \eqn{N} observed time series exhibit the same trend 
#' of some pre-specified smooth parametric form.
#' 
#' @details Arguments \code{Window}, \code{j}, and \code{q} are used to set windows 
#' for the local regression. Current version of the function assumes two options: 
#' (1) user specifies one fixed window for each time series using the argument 
#' \code{Window} (if \code{Window} is set, \code{j} and \code{q} are ignored), 
#' and (2) user specifies a set of windows by \code{j} and \code{q} to apply 
#' this set to each time series and to select an optimal window using a heuristic 
#' \eqn{m}-out-of-\eqn{n} subsampling algorithm \insertCite{Bickel_Sakov_2008}{funtimes}.
#' The option of selecting windows automatically for some of the time series, 
#' while for other time series the window is fixed, is not available yet. 
#' If none of these three arguments is set, default \code{j} and \code{q} are used. 
#' Values \code{T*q^j} are mapped to the largest previous integer, then only 
#' those greater than 2 are used.
#' 
#' See more details in \insertCite{Lyubchich_Gel_2016_synchronism;textual}{funtimes} 
#' and \insertCite{Lyubchich_2016_trends;textual}{funtimes}.
#' 
#' 
#' @param formula an object of class "\code{\link[stats]{formula}}", 
#' specifying the form of the common parametric time trend to be tested 
#' in a \eqn{T} by \eqn{N} matrix of time series 
#' (time series in columns). Variable \eqn{t} should be used to specify the form of 
#' the trend, where \eqn{t} is specified within the function as a regular sequence 
#' on the interval (0,1]. See `Examples'.
#' @param Window scalar or \eqn{N}-vector with lengths of the local windows (factors). 
#' If only one value is set, the same \code{Window} is applied to each time series. 
#' An \eqn{N}-vector gives a specific window for each time series. 
#' If \code{Window} is not specified, an automatic algorithm for optimal 
#' window selection is applied as a default option (see `Details').
#' @param q scalar from 0 to 1 to define the set of possible windows \code{T*q^j} 
#' and to automatically select an optimal window for each time series. 
#' Default is \eqn{3/4}. This argument is ignored if \code{Window} is set by user.
#' @param j numeric vector to define the set of possible windows \code{T*q^j} 
#' and to automatically select an optimal window for each time series. 
#' Default is \code{c(8:11)}. This argument is ignored if \code{Window} is set by user.
#' @param ar.order order of autoregressive filter when \code{BIC = FALSE}, 
#' or the maximal order for BIC-based filtering. Default is \code{round(10*log10(T))}. 
#' The \code{ar.order} can be a scalar or \eqn{N}-vector. If scalar, the same 
#' \code{ar.order} is applied to each time series. An \eqn{N}-vector specifies 
#' a separate \code{ar.order} for each time series.
#' @inheritParams wavk.test
#' @inheritParams ARest
#' 
#' 
#' @return A list of class \code{"htest"} containing the following components:
#' \item{method}{name of the method.}
#' \item{data.name}{name of the data.}
#' \item{statistic}{value of the test statistic.}
#' \item{p.value}{\eqn{p}-value of the test.}
#' \item{alternative}{alternative hypothesis.}
#' \item{estimate}{list with elements \code{common_trend_estimates}, 
#' \code{ar_order_used}, \code{Window_used}, \code{wavk_obs}, and 
#' \code{all_considered_windows}. The latter is a table with bootstrap and 
#' asymptotic test results for all considered windows, that is, without adaptive 
#' selection of the local window.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link[stats]{ar}}, \code{\link{HVK}}, \code{\link{WAVK}}, 
#' \code{\link{wavk.test}}
#' 
#' @keywords htest trend ts synchrony
#' 
#' @author Yulia R. Gel, Vyacheslav Lyubchich, Ethan Schaeffer, Xingyu Wang
#' 
#' @export
#' @examples
#' #Fix seed for reproducible simulations:
#' set.seed(1)
#' 
#' # Simulate two autoregressive time series of length n without trend 
#' #(i.e., with zero or constant trend) 
#' # and arrange the series into a matrix:
#' n <- 200
#' y1 <- arima.sim(n = n, list(order = c(1, 0, 0), ar = c(0.6)))
#' y2 <- arima.sim(n = n, list(order = c(1, 0, 0), ar = c(-0.2)))
#' Y <- cbind(y1, y2)
#' plot.ts(Y)
#' 
#' 
#' #Test H0 of a common linear trend:
#' \dontrun{
#'     sync.test(Y ~ t, B = 500)
#' }
#' # Sample output:
#' ##	Non-parametric test for synchronism of parametric trends
#' ##
#' ##data:  Y 
#' ##Test statistic = -0.0028999, p-value = 0.7
#' ##alternative hypothesis: common trend is not of the form Y ~ t.
#' ##sample estimates:
#' ##$common_trend_estimates
#' ##               Estimate Std. Error    t value  Pr(>|t|)
#' ##(Intercept) -0.02472566  0.1014069 -0.2438261 0.8076179
#' ##t            0.04920529  0.1749859  0.2811958 0.7788539
#' ##
#' ##$ar.order_used
#' ##         y1 y2
#' ##ar.order  1  1
#' ##
#' ##$Window_used
#' ##       y1 y2
#' ##Window 15  8
#' ##
#' ##$all_considered_windows
#' ## Window    Statistic p-value Asympt. p-value
#' ##      8 -0.000384583   0.728       0.9967082
#' ##     11 -0.024994408   0.860       0.7886005
#' ##     15 -0.047030164   0.976       0.6138976
#' ##     20 -0.015078579   0.668       0.8714980
#' ##
#' ##$wavk_obs
#' ##[1]  0.05827148 -0.06117136
#' 
#' # Add a time series y3 with a different linear trend and re-apply the test:
#' y3 <- 1 + 3*((1:n)/n) + arima.sim(n = n, list(order = c(1, 0, 0), ar = c(-0.2)))
#' Y2 <- cbind(Y, y3)
#' plot.ts(Y2)
#' \dontrun{
#'     sync.test(Y2 ~ t, B = 500)}
#' # Sample output:
#' ##	Non-parametric test for synchronism of parametric trends
#' ##
#' ##data:  Y2 
#' ##Test statistic = 0.48579, p-value < 2.2e-16
#' ##alternative hypothesis: common trend is not of the form Y2 ~ t.
#' ##sample estimates:
#' ##$common_trend_estimates
#' ##              Estimate Std. Error  t value     Pr(>|t|)
#' ##(Intercept) -0.3632963 0.07932649 -4.57976 8.219360e-06
#' ##t            0.7229777 0.13688429  5.28167 3.356552e-07
#' ##
#' ##$ar.order_used
#' ##         Y.y1 Y.y2 y3
#' ##ar.order    1    1  0
#' ##
#' ##$Window_used
#' ##       Y.y1 Y.y2 y3
#' ##Window    8   11  8
#' ##
#' ##$all_considered_windows
#' ## Window Statistic p-value Asympt. p-value
#' ##      8 0.4930069       0    1.207378e-05
#' ##     11 0.5637067       0    5.620248e-07
#' ##     15 0.6369703       0    1.566057e-08
#' ##     20 0.7431621       0    4.201484e-11
#' ##
#' ##$wavk_obs
#' ##[1]  0.08941797 -0.07985614  0.34672734
#' 
#' #Other hypothesized trend forms can be specified, for example:
#' \dontrun{
#'     sync.test(Y ~ 1) #constant trend
#'     sync.test(Y ~ poly(t, 2)) #quadratic trend
#'     sync.test(Y ~ poly(t, 3)) #cubic trend
#' }
#' 
sync.test <- function(formula, B = 1000, Window = NULL, q = NULL, j = NULL, 
                      ar.order = NULL, ar.method = "HVK", BIC = TRUE)
{
    frml <- deparse(substitute(formula))
    splt <- strsplit(frml, "~")[[1]]
    DNAME <- splt[1]
    sh <- splt[2]
    X <- eval(parse(text = DNAME), parent.frame())
    n <- nrow(X)
    K <- ncol(X)
    t <- c(1:n)/n
    if(!is.null(Window)){ #if user set Window
        UseOneWindowPerTS <- TRUE
        ONEwindow <- FALSE
        if(length(Window) == 1){
            ONEwindow <- TRUE
            Window <- rep(Window, K)
        }else if(length(Window) != K){
            stop("number of windows does not match number of time series.")
        }
        if(!is.null(q)){warning("The parameter q was not used.")}
        if(!is.null(j)){warning("The parameter j was not used.")}
    }else{
        UseOneWindowPerTS <- FALSE
        if(!is.null(q)){
            if (NCOL(q) > 1 | !is.numeric(q) | NROW(q) > 1){
                stop("q is not a scalar.")
            }
            if (q >= 1 | q <= 0){
                stop("q is out of range from 0 to 1.")
            }
        }else{
            q <- 3/4
        }
        if(!is.null(j)){
            if (!is.vector(j) | !is.numeric(j)) {
                stop("j is not a numeric vector.")
            }
        }else{
            j <- c(8:11)
        }
        kn <- n*q^j
        kn <- unique(sort(floor(kn)))
        kn <- kn[kn > 2 & kn < n]
        if (length(kn) == 0) {
            stop("set proper q and/or j.")
        }
    } 
    if(!is.null(ar.order)){ #if user set ar.order
        maxARorder <- ar.order
        if(length(ar.order) == 1){
            maxARorder <- rep(ar.order, K)
        }else if(length(ar.order) != K){
            stop("number of elements in ar.order does not match number of time series.")
        }
    }else{
        maxARorder <- rep(round(10*log10(n)), K)
    }
    #Preallocate space:
    if(!UseOneWindowPerTS){
        s <- array(NA, dim = c(length(kn), B, K))
        wavk_obs_all <- matrix(NA, length(kn), K)        
    }
    wavk_boot_opt <- array(NA, c(B, K))
    wavk_obs <- rep(NA, K)
    sigma <- rep(NA, K)
    OutputARorder <- matrix(NA, 1, K, dimnames = list("ar.order", dimnames(X)[[2]]))
    OutputWindow <- matrix(NA, 1, K, dimnames = list("Window", dimnames(X)[[2]]))
    #Function:
    X <- scale(X)
    AveragedProcess <- apply(X, 1, mean)
    mod <- lm(as.formula(paste("AveragedProcess", sh, sep = "~"))) #common trend
    TrendCoeff <- summary(mod)$coefficients
    U <- demean(X - mod$fitted) #detrended time series
    for (k in 1:K){
        pheta <- ARest(U[,k], ar.order = ar.order, ar.method = ar.method, BIC = BIC)
        OutputARorder[1, k] <- length(pheta)
        if (length(pheta) > 0) {
            tmp <- filter(X[,k], pheta, sides = 1)
            tmp2 <- filter(mod$fitted, pheta, sides = 1)
            Z <- (X[(length(pheta)+1):n, k] - tmp[length(pheta):(n - 1)]) - 
                (mod$fitted[(length(pheta)+1):n] - tmp2[length(pheta):(n - 1)])
        } else {
            Z <- U[,k]
        }
        Z <- Z - mean(Z)
        sigma[k] <- sqrt(sum(diff(Z)^2)/(2*(length(Z) - 1)))
        boot <- array(data = rnorm(n*B), c(n,B))*sigma[k]
        if(!UseOneWindowPerTS){
            for (i in 1:length(kn)){
                s[i,,k] <- apply(boot, 2, function(x) WAVK(x, kn = kn[i])$Tn/sqrt(kn[i]))
            }
            if (length(kn) > 2){
                s1 <- t(apply(s[,,k], 1, sort))
                distance <- sapply(1:(length(kn)-1), function(x) dist(s1[x:(x+1),]))
                OutputWindow[1,k] <- kn[which.min(distance)]
                wavk_boot_opt[,k] <- s[which.min(distance),,k]
            }else{
                OutputWindow[1,k] <- kn[1]
                wavk_boot_opt[,k] <- s[1,,k]
            }
            wavk_obs[k] <- WAVK(Z, kn = OutputWindow[1,k])$Tn/sqrt(OutputWindow[1, k])
            wavk_obs_all[,k] <- sapply(kn, function(x) WAVK(Z, kn = x)$Tn/sqrt(x))
        }else{
            wavk_boot_opt[,k] <- apply(boot, 2, function(x) WAVK(x, kn = Window[k])$Tn/sqrt(Window[k]))
            OutputWindow[1,k] <- Window[k]
            wavk_obs[k] <- WAVK(Z, kn = OutputWindow[1,k])$Tn/sqrt(OutputWindow[1, k])
        }
    } #k=K
    
    #p-value for bootstrap with optimal window selected
    STATISTIC <- sum(wavk_obs)
    crit <-  sum(STATISTIC > apply(wavk_boot_opt, 1, sum))/B
    if (crit < 0.5) {
        P.VALUE <- 2*crit
    } else {
        P.VALUE  <- 2*(1 - crit)
    }
    if(!UseOneWindowPerTS){
        ST <- apply(wavk_obs_all, 1, sum)
        tmp_all <- apply(s, c(1,2), sum)
        crit.boot <- sapply(c(1:length(kn)), function(x) sum(ST[x]<tmp_all[x,]))/B
        p.value.boot.all <- 2 * crit.boot
        p.value.boot.all[crit.boot>0.5] <- 2 * (1 - crit.boot[crit.boot>0.5])
        #Asymptotic results
        StAssStandard <- ST*sqrt(n)/sqrt(4*sum(sigma^4)/3) 
        crit.ass <- pnorm(StAssStandard, mean = 0, sd = 1)
        p.value.ass <- crit.ass * 2
        p.value.ass[crit.ass>0.5] <- (1 - crit.ass[crit.ass>0.5]) * 2
        #
        ESTIMATE <- list(TrendCoeff, OutputARorder, OutputWindow, 
                         cbind(kn, ST, p.value.boot.all, p.value.ass), wavk_obs)
        names(ESTIMATE) <- list("common_trend_estimates", "ar.order_used", 
                                "Window_used", "all_considered_windows", "wavk_obs")
        dimnames(ESTIMATE[[4]]) <- list(rep("", NROW(ESTIMATE[[4]])), 
                                        c("Window", "Statistic", "p-value", "Asympt. p-value"))
    }else{
        p.value.boot.all <- P.VALUE
        ST <- sum(wavk_obs)
        #Asymptotic results
        StAssStandard <- ST*sqrt(n)/sqrt(4*sum(sigma^4)/3) 
        crit.ass <- pnorm(StAssStandard, mean = 0, sd = 1)
        p.value.ass <- crit.ass * 2
        p.value.ass[crit.ass>0.5] <- (1 - crit.ass[crit.ass>0.5]) * 2
        #
        if(ONEwindow){
            ESTIMATE <- list(TrendCoeff, OutputARorder, OutputWindow, 
                             cbind(Window[1], ST, p.value.boot.all, p.value.ass), wavk_obs)
            names(ESTIMATE) <- list("common_trend_estimates", "ar.order_used", 
                                    "Window_used", "all_considered_windows", "wavk_obs")
            dimnames(ESTIMATE[[4]]) <- list(rep("", NROW(ESTIMATE[[4]])), 
                                            c("Window", "Statistic", "p-value", "Asympt. p-value"))
        }else{
            ESTIMATE <- list(TrendCoeff, OutputARorder, OutputWindow, 
                             cbind(ST, p.value.boot.all, p.value.ass), wavk_obs)
            names(ESTIMATE) <- list("common_trend_estimates", "ar.order_used", 
                                    "Window_used", "all_considered_windows", "wavk_obs")
            dimnames(ESTIMATE[[4]]) <- list(rep("", NROW(ESTIMATE[[4]])), 
                                            c("Statistic", "p-value", "Asympt. p-value"))
        }
    }
    
    METHOD <- "Non-parametric test for synchronism of parametric trends"   
    names(STATISTIC) <- "Test statistic"
    ALTERNATIVE <- paste("common trend is not of the form ", frml, ".", sep="")
    structure(list(method = METHOD, data.name = DNAME, statistic = STATISTIC,  
                   p.value = P.VALUE,  alternative = ALTERNATIVE, estimate = ESTIMATE), 
              class = "htest") 	
}
