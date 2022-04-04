#' Out-of-sample Tests of Granger Causality
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
#' @export
#' @examples
#' \dontrun{
#' # Example 1: Canada time series (ts object)
#' Canada <- vars::Canada
#' causality_pred(Canada[,1:2], cause = "e", lag.max = 5)
#' causality_pred(Canada[,1:2], cause = "e", lag.restrict = 3, lag.max = 15)
#'
#' # Example 2 (run in parallel): Box & Jenkins time series
#' # of sales and a leading indicator, see ?BJsales
#'
#' # Initiate a local cluster
#' cores <- parallel::detectCores()
#' cl <- parallel::makeCluster(cores)
#' parallel::clusterSetRNGStream(cl, 123) # to make parallel computations reproducible
#'
#' D <- cbind(BJsales.lead, BJsales)
#' causality_pred(D, cause = "BJsales.lead", lag.max = 5, B = 1000, cl = cl)
#' causality_pred(D, cause = "BJsales.lead", lag.restrict = 3, p = 5, B = 1000, cl = cl)
#' }
#'
causality_pred <- function(y, p = NULL,
                           cause = NULL,
                           B = 100L,
                           test = 0.3,
                           cl = 1L,
                           lag.max = NULL,
                           k = 2,
                           lag.restrict = 0L)
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
    dep <- setdiff(varnames, cause)[1] # dependent variable name

    # Define samples
    n <- nrow(y) # sample size
    if (test < 1) { # percentage split
        n_train <- round(n*(1 - test))
    } else { # use the last "test" observations as the testing set
        n_train <- n - test
    }
    n_test <- n - n_train

    # Estimate p if using an information criterion
    if (is.null(p) && is.null(lag.max)) {
        stop("Please specify p or lag.max.")
    }
    if (!is.null(lag.max)) { # then select p
        if (lag.restrict >= lag.max) {
            warning("lag.restrict >= lag.max. Using lag.restrict = 0 instead.")
            lag.restrict <- 0
        }
        if (lag.restrict > 0) {
            lagX <- embed(y[,cause], lag.max + 1)[, -c(1:(lag.restrict + 1)), drop = FALSE]
        } else {
            lagX <- embed(y[,cause], lag.max + 1)[, -1, drop = FALSE]
        }
        lagY <- embed(y[,dep], lag.max + 1)
        # Information criterion for the model
        # (use rows 1:(n_train - lag.max) to use training set only)
        IC <- sapply(1:ncol(lagX), function(s) {
            fit <- stats::lm.fit(x = cbind(1,
                                           lagY[, 2:(s + lag.restrict + 1)],
                                           lagX[, 1:s, drop = FALSE]),
                                 y = lagY[, 1])
            # see stats:::extractAIC.lm; but omit the scale option
            nfit <- length(fit$residuals)
            edf <- nfit - fit$df.residual
            RSS <- sum(fit$residuals^2, na.rm = TRUE) # stats:::deviance.lm(fit)
            dev <- nfit * log(RSS/nfit)
            dev + k * edf
        })
        p <- which.min(IC) + lag.restrict
    } # finish selection of p

    # Estimate model and get predictions on the testing set
    if (lag.restrict >= p) {
        warning("lag.restrict >= p. Using lag.restrict = 0 instead.")
        lag.restrict <- 0
    }
    if (lag.restrict > 0) {
        lagX <- embed(y[,cause], p + 1)[, -c(1:(lag.restrict + 1)), drop = FALSE]
    } else {
        lagX <- embed(y[,cause], p + 1)[, -1, drop = FALSE]
    }
    lagY <- embed(y[,dep], p + 1)
    n_actual <- nrow(lagY)
    # m_yx <- stats::lm(X1 ~ ., data = Dyx[1:(n_train - p),])
    # m_yx <- stats::lm.fit(x = cbind(1, lagY[1:(n_train - p), -1], lagX[1:(n_train - p), ]),
    #                     y = lagY[1:(n_train - p), 1])

    # Get 1-step ahead forecasts from the models (recursive)
    FCST <- sapply(1:n_test - 1, function(i) { # i = 0
        # estimate full and restricted models
        m_yx <- stats::lm.fit(x = cbind(1,
                                        lagY[(1 + i):(n_train - p + i), -1, drop = FALSE],
                                        lagX[(1 + i):(n_train - p + i), , drop = FALSE]),
                              y = lagY[(1 + i):(n_train - p + i), 1, drop = FALSE])
        m_y <- stats::lm.fit(x = cbind(1,
                                       lagY[(1 + i):(n_train - p + i), -1, drop = FALSE]),
                             y = lagY[(1 + i):(n_train - p + i), 1, drop = FALSE])
        # get predictions
        c(m_yx$coefficients %*% c(1, lagY[(n_train - p + i + 1), -1], lagX[(n_train - p + i + 1), ]),
          m_y$coefficients %*% c(1, lagY[(n_train - p + i + 1), -1]))
    })
    # Forecast errors
    efull <- y[(n_train + 1):n, dep] - FCST[1,]
    eres <- y[(n_train + 1):n, dep] - FCST[2,]
    # Observed test statistics
    OBS <- caustests(efull, eres) # caustests <- funtimes:::caustests

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
                 cause = cause,
                 p = p)

    # Bootstrap restricted model estimated on the full sample
    m_y <- stats::lm.fit(x = lagY[,-1, drop = FALSE],
                         y = lagY[, 1, drop = FALSE])
    m_y_fit <- m_y$fitted.values
    m_y_res <- m_y$residuals
    if (bootparallel) {
        BOOT0 <- parallel::parSapply(cl, X = 1:B, FUN = function(b) {
            dy_boot <- m_y_fit + sample(m_y_res, replace = TRUE) * rnorm(n_actual)
            FCST <- sapply(1:n_test - 1, function(i) { # i = 0
                # estimate full and restricted models
                m_yx <- stats::lm.fit(x = cbind(1,
                                                lagY[(1 + i):(n_train - p + i), -1, drop = FALSE],
                                                lagX[(1 + i):(n_train - p + i), , drop = FALSE]),
                                      y = dy_boot[(1 + i):(n_train - p + i)])
                m_y <- stats::lm.fit(x = cbind(1,
                                               lagY[(1 + i):(n_train - p + i), -1, drop = FALSE]),
                                     y = dy_boot[(1 + i):(n_train - p + i)])
                # get predictions
                c(m_yx$coefficients %*% c(1, lagY[(n_train - p + i + 1), -1], lagX[(n_train - p + i + 1), ]),
                  m_y$coefficients %*% c(1, lagY[(n_train - p + i + 1), -1]))
            })
            # Forecast errors
            efullb <- dy_boot[(n_train + 1):n - p] - FCST[1,]
            eresb <- dy_boot[(n_train + 1):n - p] - FCST[2,]
            # test statistics
            caustests(efullb, eresb)
        })
        if (clStop) {
            parallel::stopCluster(cl)
        }
    } else {
        BOOT0 <- #parallel::parSapply(cl, X = 1:B, FUN = function(b) {
            sapply(1:B, FUN = function(b) {
                dy_boot <- m_y_fit + sample(m_y_res, replace = TRUE) * rnorm(n_actual)
                FCST <- sapply(1:n_test - 1, function(i) { # i = 0
                    # estimate full and restricted models
                    m_yx <- stats::lm.fit(x = cbind(1,
                                                    lagY[(1 + i):(n_train - p + i), -1, drop = FALSE],
                                                    lagX[(1 + i):(n_train - p + i), , drop = FALSE]),
                                          y = dy_boot[(1 + i):(n_train - p + i)])
                    m_y <- stats::lm.fit(x = cbind(1,
                                                   lagY[(1 + i):(n_train - p + i), -1, drop = FALSE]),
                                         y = dy_boot[(1 + i):(n_train - p + i)])
                    # get predictions
                    c(m_yx$coefficients %*% c(1, lagY[(n_train - p + i + 1), -1], lagX[(n_train - p + i + 1), ]),
                      m_y$coefficients %*% c(1, lagY[(n_train - p + i + 1), -1]))
                })
                # Forecast errors
                efullb <- dy_boot[(n_train + 1):n - p] - FCST[1,]
                eresb <- dy_boot[(n_train + 1):n - p] - FCST[2,]
                # test statistics
                caustests(efullb, eresb)
            })
    } # end sequential bootstrap
    FullH0 <- list(result = data.frame(MSEt = c(OBS["MSEt"],
                                                (sum(BOOT0["MSEt",] <= OBS["MSEt"]) + 1) / (B + 1),
                                                stats::pt(OBS["MSEt"], n_test - 1, lower.tail = TRUE)),
                                       MSEcor = c(OBS["MSEcor"],
                                                  (sum(BOOT0["MSEcor",] <= OBS["MSEcor"]) + 1) / (B + 1),
                                                  stats::pt(OBS["MSEcor"], n_test - 1, lower.tail = TRUE)),
                                       row.names = c("stat_obs", "p_boot", "p_asympt")),
                   cause = cause,
                   p = p)
    return(list(FAST = FAST, FullH0 = FullH0))
}
