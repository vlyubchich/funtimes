#' Out-of-sample Tests of Granger Causality using (Restricted) Vector Autoregression
#'
#' Test for Granger causality using out-of-sample prediction errors from a vector
#' autoregression (VAR), where the original VAR can be restricted (see Details).
#' The tests include the MSE-t approach \insertCite{McCracken_2007}{funtimes} and
#' MSE-correlation test as in Chapter 9.3 of \insertCite{Granger_Newbold_2016;textual}{funtimes}.
#' The bootstrap is used to empirically derive distributions of the statistics.
#'
#'
#' @details The arguments specified in \code{...} are passed to the \code{\link[vars]{VAR}} function.
#' Additionally, \code{lag.restrict} can be specified to remove short-term lags from
#' consideration (\code{lag.restrict} is not an option in the original package \code{vars}).
#' Note that if \code{p} is specified, \code{lag.restrict} must be smaller
#' than \code{p} otherwise the default \code{lag.restrict = 0} will be used.
#' If \code{lag.max} is specified instead of \code{p}, VAR orders
#' \code{lag.restrict} + 1, \dots, \code{lag.max} will be considered using the training data
#' and the order \eqn{p} will be automatically selected according to the information criterion
#' (by default, AIC).
#'
#' In the current implementation, the bootstrapped \eqn{p}-value is calculated using equation 4.10 of
#' \insertCite{Davison_Hinkley_1997;textual}{funtimes}: \code{p.value} = (1 + \eqn{n}) / (\code{B} + 1),
#' where \eqn{n} is the number of bootstrapped statistics smaller or equal to the observed statistic.
#' In the fast bootstrap, \eqn{n} is the number of bootstrapped statistics greater or equal to 0.
#'
#' This function uses symmetric VAR with the same orders \eqn{p} for modeling both \eqn{Y} to \eqn{X}.
#' To select these orders more independently, consider using the function \code{\link{causality_pred}}.
#'
#'
#' @param y data frame or \code{ts} object for estimating VAR(\eqn{p}).
#' @param p an integer specifying the order \eqn{p} in VAR.
#' By default (if \code{p} is not specified),
#' \eqn{p} is selected based on the information criterion
#' \code{ic} (see \code{...} arguments; default \code{ic} is AIC).
#' @inheritParams causality_pred
#' @param ... other arguments passed to the function for VAR estimation.
#' The arguments include \code{lag.restrict} that is used to remove the first lags
#' in the cause variable from consideration (use restricted VAR to avoid testing for short-term causality);
#' default \code{lag.restrict = 0L}, i.e., no restrictions.
#' Other possible arguments are as in the \code{\link[vars]{VAR}} function.
#' Also, see Details and Examples.
#'
#'
#' @return Two lists (one for the fast bootstrap,
#' another for the bootstrap under the null hypothesis) each containing the following elements:
#' \item{result}{a table with the observed values of the test statistics and \eqn{p}-values.}
#' \item{cause}{the cause variable.}
#' \item{p}{the AR order used.}
#'
#' @references
#' \insertAllCited{}
#'
#' @keywords causality htest ts
#'
#' @author Vyacheslav Lyubchich
#'
#' @seealso \code{\link{causality_pred}}
#'
#' @importFrom stats predict residuals
#' @importFrom mlVAR simulateVAR
#' @export
#' @examples
#' \dontrun{
#' # Example 1: Canada time series (ts object)
#' Canada <- vars::Canada
#' causality_predVAR(Canada[,1:2], cause = "e", lag.max = 5)
#' causality_predVAR(Canada[,1:2], cause = "e", lag.restrict = 3, lag.max = 15)
#'
#' # Example 2 (run in parallel, initiate the cluster manually):
#' # Box & Jenkins time series
#' # of sales and a leading indicator, see ?BJsales
#'
#' # Initiate a local cluster
#' cores <- parallel::detectCores()
#' cl <- parallel::makeCluster(cores)
#' parallel::clusterSetRNGStream(cl, 123) # to make parallel computations reproducible
#'
#' D <- cbind(BJsales.lead, BJsales)
#' causality_predVAR(D, cause = "BJsales.lead", lag.max = 5, B = 1000, cl = cl)
#' causality_predVAR(D, cause = "BJsales.lead", lag.restrict = 3, p = 5, B = 1000, cl = cl)
#' parallel::stopCluster(cl)
#' }
#'
causality_predVAR <- function(y, p = NULL,
                              cause = NULL,
                              B = 500,
                              test = 0.3,
                              cl = 1,
                              ...)
{
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
    varnames <- colnames(y)
    if (length(varnames) != 2) stop("y must have 2 columns")
    if (is.null(cause)) {
        cause <- varnames[1]
    }
    dep <- setdiff(varnames, cause)[1] # dependent variable
    K <- 2L #length(varnames)

    # Define samples
    n <- nrow(y) # sample size
    if (test < 1) { # percentage split
        n_train <- round(n*(1 - test))
    } else { # use the last "test" observations as the testing set
        n_train <- n - test
    }
    n_test <- n - n_train

    # Estimate model on the training data to get the coefficient structure and
    # estimate p if using an information criterion
    x <- VAR(y[1:n_train,], p = p, ...)
    if (is.null(p)) {
        ptrain <- x$p
    } else {
        ptrain <- p
    }
    R2inv <- restrictions(x, cause)

    # # Recreate matrix of restrictions and overlay new restrictions for causality testing
    # co.names <- vars::Bcoef(x)
    # k <- which(gsub("\\.l\\d+", "", colnames(co.names)) %in% cause) # select cause regressors
    # l <- which(rownames(co.names) %in% cause) # select cause regressand
    # R2inv <- matrix(1, ncol = ncol(co.names), nrow = nrow(co.names))
    # R2inv[-l, k] <- 0 # select coef to be tested
    # # If the model already has restriction, overlay with the new ones
    # if (!is.null(x$restrictions)) {
    #     xr <- x$restrictions
    #     # match positions of variables
    #     xr <- xr[rownames(co.names), colnames(co.names)]
    #     # overlay
    #     R2inv <- xr * R2inv
    # }

    # Get 1-step ahead forecasts from the models (recursive)
    FCST <- sapply(1:n_test - 1, function(i) { # i = 0
        # estimate full model, VAR or restricted VAR (depends on lag.restrict)
        x <- VAR(y[(1 + i):(n_train + i),], p = ptrain, ...)
        ff <- predict(x, n.ahead = 1)$fcst[[dep]][1] # VAR predictions
        xres <- vars::restrict(x, method = "man", resmat = R2inv)
        fr <- predict(xres, n.ahead = 1)$fcst[[dep]][1] # VAR predictions
        c(ff, fr)
    })
    # Forecast errors
    efull <- y[(n_train + 1):n, dep] - FCST[1,]
    eres <- y[(n_train + 1):n, dep] - FCST[2,]
    # Observed test statistics
    OBS <- caustests(efull, eres)

    # Fast bootstrap (only the out-of-sample errors)
    BOOT <- #parallel::parSapply(cl, X = 1:B, FUN = function(b) {
        sapply(1:B, FUN = function(b) {
            # bootstrap prediction errors
            ind <- sample(n_test, replace = TRUE)
            efullb <- efull[ind]
            # ind <- sample(n_test, replace = TRUE)
            eresb <- eres[ind]
            # get bootstrapped statistics
            caustests(efullb, eresb)
        })
    Fast <- list(result = data.frame(MSEt = c(OBS["MSEt"],
                                              (sum(BOOT["MSEt",] >= 0) + 1) / (B + 1)),
                                     MSEcor = c(OBS["MSEcor"],
                                                (sum(BOOT["MSEcor",] >= 0) + 1) / (B + 1)),
                                     row.names = c("stat_obs", "p_boot")),
                 p = ptrain)

    # Bootstrap restricted model estimated on the full sample
    x <- VAR(y, p = p, ...)
    if (is.null(p)) {
        pfull <- x$p
    } else {
        pfull <- p
    }
    R2inv <- restrictions(x, cause)
    xres <- vars::restrict(x, method = "man", resmat = R2inv)
    xres_coef <- vars::Bcoef(xres)
    # disregard intercept and other coefficients, simulate just VAR
    keep <- which(gsub("\\.l\\d+", "", colnames(xres_coef)) %in% varnames)
    # convert the matrix of coefficients into a list
    xres_coef <- lapply(1:(length(keep)/K), function(i) xres_coef[,keep[((i - 1)*K + 1):(i*K)]])
    xres_cov <- cov(residuals(xres))
    if (bootparallel) {
        BOOT0 <- parallel::parSapply(cl, X = 1:B, FUN = function(b) {
            # sapply(1:B, FUN = function(b) {
            yb <- mlVAR::simulateVAR(pars = xres_coef, lags = 1:pfull,
                                     ,residuals = xres_cov
                                     ,Nt = n)
            names(yb) <- varnames
            # xb <- VAR(yb[1:n_train,], p = pfull, ...)
            # if (is.null(p)) {
            #     ptrainb <- xb$p
            # } else {
            #     ptrainb <- p
            # }
            # R2inv <- restrictions(xb, cause)
            FCST <- sapply(1:n_test - 1, function(i) { # i = 0
                # estimate full model, VAR or restricted VAR (depends on lag.restrict)
                x <- VAR(yb[(1 + i):(n_train + i),], p = pfull, ...)
                ff <- predict(x, n.ahead = 1)$fcst[[dep]][1] # VAR predictions
                xres <- vars::restrict(x, method = "man", resmat = R2inv)
                fr <- predict(xres, n.ahead = 1)$fcst[[dep]][1] # VAR predictions
                c(ff, fr)
            })
            # Forecast errors
            efullb <- yb[(n_train + 1):n, dep] - FCST[1,]
            eresb <- yb[(n_train + 1):n, dep] - FCST[2,]
            # test statistics
            caustests(efullb, eresb)
        })
        if (clStop) {
            parallel::stopCluster(cl)
        }
    } else {
        BOOT0 <- #parallel::parSapply(cl, X = 1:B, FUN = function(b) {
            sapply(1:B, FUN = function(b) {
                yb <- mlVAR::simulateVAR(pars = xres_coef, lags = 1:pfull,
                                         ,residuals = xres_cov
                                         ,Nt = n)
                names(yb) <- varnames
                # xb <- VAR(yb[1:n_train,], p = pfull, ...)
                # if (is.null(p)) {
                #     ptrainb <- xb$p
                # } else {
                #     ptrainb <- p
                # }
                # R2inv <- restrictions(xb, cause)
                FCST <- sapply(1:n_test - 1, function(i) { # i = 0
                    # estimate full model, VAR or restricted VAR (depends on lag.restrict)
                    x <- VAR(yb[(1 + i):(n_train + i),], p = pfull, ...)
                    ff <- predict(x, n.ahead = 1)$fcst[[dep]][1] # VAR predictions
                    xres <- vars::restrict(x, method = "man", resmat = R2inv)
                    fr <- predict(xres, n.ahead = 1)$fcst[[dep]][1] # VAR predictions
                    c(ff, fr)
                })
                # Forecast errors
                efullb <- yb[(n_train + 1):n, dep] - FCST[1,]
                eresb <- yb[(n_train + 1):n, dep] - FCST[2,]
                # test statistics
                caustests(efullb, eresb)
            })
    }#end sequential bootstrap
    FullH0 <- list(result = data.frame(MSEt = c(OBS["MSEt"],
                                                (sum(BOOT0["MSEt",] <= OBS["MSEt"]) + 1) / (B + 1)),
                                       MSEcor = c(OBS["MSEcor"],
                                                  (sum(BOOT0["MSEcor",] <= OBS["MSEcor"]) + 1) / (B + 1)),
                                       row.names = c("stat_obs", "p_boot")),
                   p = pfull)
    return(list(Fast = Fast, FullH0 = FullH0))
}
