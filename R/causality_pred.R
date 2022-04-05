#' Out-of-sample Tests of Granger Causality
#'
#' Test for Granger causality using out-of-sample prediction errors from an
#' autoregression (AR) model, where some of the near-contemporaneous lags can be removed:
#' \deqn{Y_t = \sum_{i=1}^p\alpha_iY_{t-i} + \sum_{i=lag.restrict+1}^p\beta_iX_{t-i} + e_t,}
#' where \eqn{Y_t} is the dependent variable,
#' \eqn{X_t} is the cause variable,
#' \eqn{p} is the AR order,
#' \eqn{lag.restrict} is the number of restricted first lags (see the argument \code{lag.restrict}).
#'
#' The tests include the MSE-t approach \insertCite{McCracken_2007}{funtimes} and
#' MSE-correlation test as in Chapter 9.3 of \insertCite{Granger_Newbold_2016;textual}{funtimes}.
#' The bootstrap is used to empirically derive distributions of the statistics.
#'
#'
#' @details
#' Two versions of the bootstrap are used to derive the empirical distributions of the test statistics.
#' In the original version, residuals of a restricted model under the null hypothesis of no Granger
#' causality are bootstrapped to generate new data under the null. Then the full and restricted
#' models are re-estimated on the bootstrapped data to obtain new (bootstrapped) forecast errors.
#' In the fast bootstrap version, just the paired residuals from the full and restricted models are
#' bootstrapped.
#'
#' In the current implementation, the bootstrapped \eqn{p}-value is calculated using equation 4.10 of
#' \insertCite{Davison_Hinkley_1997;textual}{funtimes}: \code{p.value} = (1 + \eqn{n}) / (\code{B} + 1),
#' where \eqn{n} is the number of bootstrapped statistics smaller or equal to the observed statistic.
#' In the fast bootstrap, \eqn{n} is the number of bootstrapped statistics greater or equal to 0.
#'
#' This function allows using different orders \eqn{p} when testing the Granger causation
#' of \eqn{X} to \eqn{Y} and from \eqn{Y} to \eqn{X}
#' (need to run the function twice, with different argument \code{cause}).
#' To use the symmetric vector autoregression (VAR), use the function \code{\link{causality_predVAR}}.
#'
#' @param y matrix, data frame, or \code{ts} object with two columns
#' (an explanatory and the dependent time-series variable). Missing values are not allowed.
#' @param cause name of the cause variable. If not specified, the first variable in
#' \code{y} is treated as the cause and the second is treated as the dependent variable.
#' @param p an integer specifying the order \eqn{p} of autoregressive dependence
#' (same \eqn{p} used for both explanatory and dependent variables).
#' By default (if \code{p} is not specified); the user must specify \code{p} or \code{lag.max}.
#' If \code{lag.max} is specified, the argument \code{p} is ignored.
#' @param lag.restrict integer for the number of short-term lags in the cause variable
#' to remove from consideration (default is zero, meaning no lags are removed).
#' This setting does not affect the dependent variable lags that are always present.
#' @param lag.max integer for the highest lag order to explore.
#' The order is then selected using the Akaike information criterion (AIC; default),
#' see the argument \code{k} to change the criterion.
#' @param k numeric specifying the weight of the equivalent degrees of freedom part
#' in the AIC formula. Default \code{k = 2} corresponds to the traditional AIC.
#' Use \code{k = log(n)} to use the Bayesian information criterion instead
#' (see \code{\link[stats]{extractAIC}}).
#' @param B number of bootstrap replications. Default is 500.
#' @param test a numeric value specifying the size of the testing set. If \code{test} < 1,
#' the value is treated as a proportion of the sample size to be used as the testing set.
#' Otherwise, \code{test} is rounded and \code{test} values are used as the testing set.
#' Default is 0.3, which means that 30% of the sample is used for calculating
#' out-of-sample errors. The testing set is always at the end of the time series.
#' @param cl parameter to specify computer cluster for bootstrapping passed to
#' the package \code{parallel} (default \code{cl = 1}, means no cluster is used).
#' Possible values are:
#' \itemize{
#'   \item cluster object (list) produced by \link[parallel]{makeCluster}.
#'   In this case, a new cluster is not started nor stopped;
#'   \item \code{NULL}. In this case, the function will detect
#'   available cores (see \link[parallel]{detectCores}) and, if there are
#'   multiple cores (\eqn{>1}), a cluster will be started with
#'   \link[parallel]{makeCluster}. If started, the cluster will be stopped
#'   after the computations are finished;
#'   \item positive integer defining the number of cores to start a cluster.
#'   If \code{cl = 1} (default), no attempt to create a cluster will be made.
#'   If \code{cl} > 1, a cluster will be started (using \link[parallel]{makeCluster})
#'   and stopped afterward (using \link[parallel]{stopCluster}).
#' }
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
#' @seealso \code{\link{causality_predVAR}}
#' @export
#' @examples
#' \dontrun{
#' # Example 1: Canada time series (ts object)
#' Canada <- vars::Canada
#' causality_pred(Canada[,1:2], cause = "e", lag.max = 5)
#' causality_pred(Canada[,1:2], cause = "e", lag.restrict = 3, lag.max = 15)
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
#' causality_pred(D, cause = "BJsales.lead", lag.max = 5, B = 1000, cl = cl)
#' causality_pred(D, cause = "BJsales.lead", lag.restrict = 3, p = 5, B = 1000, cl = cl)
#' parallel::stopCluster(cl)
#' }
#'
causality_pred <- function(y, cause = NULL,
                           p = NULL,
                           lag.restrict = 0L,
                           lag.max = NULL,
                           k = 2,
                           B = 500L,
                           test = 0.3,
                           cl = 1L)
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
    Fast <- list(result = data.frame(MSEt = c(OBS["MSEt"],
                                              (sum(BOOT["MSEt",] >= 0) + 1) / (B + 1)),
                                     MSEcor = c(OBS["MSEcor"],
                                                (sum(BOOT["MSEcor",] >= 0) + 1) / (B + 1)),
                                     row.names = c("stat_obs", "p_boot")),
                 cause = cause,
                 p = p)

    # Bootstrap restricted model estimated on the full sample
    m_y <- stats::lm.fit(x = lagY[,-1, drop = FALSE],
                         y = lagY[, 1, drop = FALSE])
    m_y_fit <- m_y$fitted.values
    m_y_res <- m_y$residuals
    if (bootparallel) {
        BOOT0 <- parallel::parSapply(cl, X = 1:B, FUN = function(b) {
            dy_boot <- m_y_fit + sample(m_y_res, replace = TRUE)
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
                dy_boot <- m_y_fit + sample(m_y_res, replace = TRUE)
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
                                                (sum(BOOT0["MSEt",] <= OBS["MSEt"]) + 1) / (B + 1)),
                                       MSEcor = c(OBS["MSEcor"],
                                                  (sum(BOOT0["MSEcor",] <= OBS["MSEcor"]) + 1) / (B + 1)),
                                       row.names = c("stat_obs", "p_boot")),
                   cause = cause,
                   p = p)
    return(list(Fast = Fast, FullH0 = FullH0))
}
