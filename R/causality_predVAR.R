#' Out-of-sample Tests of Granger Causality using (Restricted) Vector Autoregression
#'
#' Test for Granger causality using out-of-sample prediction errors from a vector
#' autoregression (VAR), where the original VAR can be a restricted VAR (see Details).
#' The tests include the MSE-t approach \insertCite{McCracken_2007}{funtimes},
#' MSE-correlation test as in Chapter 9.3 of \insertCite{Granger_Newbold_2016;textual}{funtimes},
#' and difference of squared VAR prediction errors (Md statistic) when not using the cause
#' variable and with the cause variable.
#' Bootstrap is used to empirically derive distributions of the statistics.
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
#' where \eqn{n} is number of bootstrapped statistics greater or equal to 0.
#'
#'
#' @param y data frame or \code{ts} object for estimating VAR(\eqn{p}).
#' @param p an integer specifying the order \eqn{p} in VAR.
#' By default (if \code{p} is not specified),
#' \eqn{p} is selected based on the information criterion
#' \code{ic} (see \code{...} arguments; default \code{ic} is AIC).
#' @param cause name of the cause variable. If not specified, the first variable in
#' \code{y} is treated as the cause, and second -- as the dependent variable.
#' @param B number of bootstrap replications. Default is 100.
#' @param test a numeric value specifying the size of the testing set. If \code{test} < 1,
#' the value is treated as proportion of the sample size to be used as the testing set.
#' Otherwise, \code{test} is rounded and \code{test} values are used as the testing set.
#' Default is 0.3, which means that 30% of the sample are used for calculating
#' out-of-sample errors. The testing set is always at the end of the time series.
#' @param cl parameter to specify computer cluster for bootstrapping, passed to
#' the package \code{parallel} (default is \code{1}, meaning no cluster is used).
#' Possible values are:
#' \itemize{
#'   \item cluster object (list) produced by \link[parallel]{makeCluster}.
#'   In this case, new cluster is not started nor stopped;
#'   \item \code{NULL}. In this case, the function will detect
#'   available cores (see \link[parallel]{detectCores}) and, if there are
#'   multiple cores (\eqn{>1}), a cluster will be started with
#'   \link[parallel]{makeCluster}. If started, the cluster will be stopped
#'   after the computations are finished;
#'   \item positive integer defining the number of cores to start a cluster.
#'   If \code{cl = 1}, no attempt to create a cluster will be made.
#'   If \code{cl > 1}, cluster will be started (using \link[parallel]{makeCluster})
#'   and stopped afterwards (using \link[parallel]{stopCluster}).
#' }
#' @param ... other arguments passed to the function for VAR estimation.
#' The arguments include \code{lag.restrict} that is used to remove a number of first lags
#' in the cause variable from consideration (use restricted VAR to avoid testing for short-term causality);
#' default \code{lag.restrict = 0L}, i.e., no restrictions.
#' Other possible arguments are as in the \code{\link[vars]{VAR}} function.
#' Also see Details and Examples.
#'
#'
#' @return A list containing the following elements:
#' \item{MSEt}{observed value of the MSEt statistic.}
#' \item{MSEt_p}{bootstrapped \eqn{p}-value of the MSEt test.}
#' \item{MSEt_p_asympt}{asymptotic \eqn{p}-value of the MSEt test,
#' based on the left tail of the t distribution.}
#' \item{MSEcor}{observed value of the MSEcor statistic.}
#' \item{MSEcor_p}{bootstrapped \eqn{p}-value of the MSEcor test.}
#' \item{MSEcor_p_asympt}{asymptotic \eqn{p}-value of the MSEcor test,
#' based on the left tail of the t distribution.}
#' \item{Md}{observed value of the Md statistic.}
#' \item{Md_p}{bootstrapped \eqn{p}-value of the Md test.}
#' \item{p}{order \eqn{p} of the VAR(\eqn{p}) model.}
#'
#' @references
#' \insertAllCited{}
#'
#' @keywords causality htest ts
#'
#' @author Vyacheslav Lyubchich
#'
#' @importFrom stats predict residuals
#' @importFrom mlVAR simulateVAR
#' @export
#' @examples
#' \dontrun{
#' cores <- parallel::detectCores()
#' cl <- parallel::makeCluster(cores)
#' loadedfunc <- lsf.str()
#' parallel::clusterExport(cl, varlist = as.list(loadedfunc), envir = environment())

#' # Example 1: Canada time series (ts object)
#' Canada <- vars::Canada
#' causality_predVAR(Canada[,1:2], cause = "e", lag.max = 5, cl = cl)
#' causality_predVAR(Canada[,1:2], cause = "e", lag.restrict = 3, cl = cl, lag.max = 15)
#'
#' # Example 2: Box & Jenkins time series of sales and a leading indicator, see ?BJsales
#' D <- cbind(BJsales.lead, BJsales)
#' causality_predVAR(D, cause = "BJsales.lead", lag.max = 5, B = 100)
#' causality_predVAR(D, cause = "BJsales.lead", lag.restrict = 3, p = 5, B = 100)
#' }
#'
causality_predVAR <- function(y, p = NULL,
                              cause = NULL,
                              B = 100,
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
    FAST <- list(result = data.frame(MSEt = c(OBS["MSEt"],
                                              (sum(BOOT["MSEt",] >= 0) + 1) / (B + 1),
                                              stats::pt(OBS["MSEt"], n_test - 1, lower.tail = TRUE)),
                                     MSEcor = c(OBS["MSEcor"],
                                                (sum(BOOT["MSEcor",] >= 0) + 1) / (B + 1),
                                                stats::pt(OBS["MSEcor"], n_test - 1, lower.tail = TRUE)),
                                     row.names = c("stat_obs", "p_boot", "p_asympt")),
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
                                                (sum(BOOT0["MSEt",] <= OBS["MSEt"]) + 1) / (B + 1),
                                                stats::pt(OBS["MSEt"], n_test - 1, lower.tail = TRUE)),
                                       MSEcor = c(OBS["MSEcor"],
                                                  (sum(BOOT0["MSEcor",] <= OBS["MSEcor"]) + 1) / (B + 1),
                                                  stats::pt(OBS["MSEcor"], n_test - 1, lower.tail = TRUE)),
                                       row.names = c("stat_obs", "p_boot", "p_asympt")),
                   p = pfull)
    return(list(FAST = FAST, FullH0 = FullH0))
}
