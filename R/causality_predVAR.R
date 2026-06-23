#' @importFrom stats predict residuals
#' @importFrom mlVAR simulateVAR

.get_recursive_VAR_fcsts <- function(y, p, n_train, n_test, dep, R2inv, ...) {
    sapply(1:n_test - 1, function(i) {
        # estimate full model, VAR or restricted VAR (depends on lag.restrict)
        x <- VAR(y[(1 + i):(n_train + i), ], p = p, ...)
        ff <- predict(x, n.ahead = 1)$fcst[[dep]][1] # VAR predictions
        xres <- vars::restrict(x, method = "man", resmat = R2inv)
        fr <- predict(xres, n.ahead = 1)$fcst[[dep]][1] # VAR predictions
        c(ff, fr)
    })
}

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
                              B = 500L,
                              test = 0.3,
                              cl = 1L,
                              ...)
{
    # Input validation
    y <- as.data.frame(y)
    if (ncol(y) != 2) stop("y must have 2 columns")
    if (!all(vapply(y, is.numeric, logical(1L)))) stop("Both columns in y must be numeric.")
    if (anyNA(y)) stop("Missing values are not allowed in y.")
    if (!is.null(p) && any(p < 1L)) stop("p must be positive integer(s).")
    if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1) stop("B must be a single positive integer.")
    B <- as.integer(B)
    if (!is.numeric(test) || length(test) != 1L || is.na(test) || test <= 0) stop("test must be a single positive numeric value.")

    # Setup parallel processing
    bootparallel <- FALSE
    if (is.list(cl)) {
        bootparallel <- TRUE
    } else {
        cores <- if (is.null(cl)) parallel::detectCores() else cl
        if (cores > 1) {
            bootparallel <- TRUE
            cl <- parallel::makeCluster(cores)
            on.exit(parallel::stopCluster(cl), add = TRUE)
        }
    }

    # Identify cause and dependent variables
    varnames <- colnames(y)
    if (is.null(cause)) {
        cause <- varnames[1]
    } else {
        cause <- as.character(cause)[1L]
        if (!(cause %in% varnames)) stop("cause must match one of the column names in y.")
    }
    dep <- setdiff(varnames, cause)[1]
    K <- 2L

    # Define samples
    n <- nrow(y)
    n_train <- if (test < 1) round(n * (1 - test)) else n - as.integer(test)
    n_test <- n - n_train
    if (n_train < 3) stop("Training sample is too short for model estimation.")
    if (n_test < 2) stop("Testing sample must contain at least 2 observations.")

    # Estimate model on the training data to get p and coefficient restrictions
    x_train <- VAR(y[1:n_train, ], p = p, ...)
    ptrain <- if (is.null(p)) x_train$p else p
    if (n_train <= ptrain) stop("Training sample is too short for selected lag order.")
    if (n <= ptrain) stop("Sample is too short for selected lag order.")
    R2inv <- restrictions(x_train, cause)

    # Get 1-step ahead forecasts from the models (recursive)
    FCST <- .get_recursive_VAR_fcsts(y, ptrain, n_train, n_test, dep, R2inv, ...)
    efull <- y[(n_train + 1):n, dep] - FCST[1, ]
    eres <- y[(n_train + 1):n, dep] - FCST[2, ]
    OBS <- caustests(efull, eres)

    # Fast bootstrap (on out-of-sample errors)
    BOOT <- sapply(1:B, FUN = function(b) {
        ind <- sample(n_test, replace = TRUE)
        caustests(efull[ind], eres[ind])
    })
    Fast <- list(result = data.frame(MSEt = c(OBS["MSEt"], (sum(BOOT["MSEt",] >= 0) + 1) / (B + 1)),
                                     MSEcor = c(OBS["MSEcor"], (sum(BOOT["MSEcor",] >= 0) + 1) / (B + 1)),
                                     row.names = c("stat_obs", "p_boot")),
                 p = ptrain)

    # Bootstrap under the null hypothesis
    x_full <- VAR(y, p = p, ...)
    pfull <- if (is.null(p)) x_full$p else p
    R2inv_full <- restrictions(x_full, cause)
    xres_full <- vars::restrict(x_full, method = "man", resmat = R2inv_full)
    xres_coef_mat <- vars::Bcoef(xres_full)
    keep <- which(gsub("\\.l\\d+", "", colnames(xres_coef_mat)) %in% varnames)
    xres_coef_list <- lapply(1:(length(keep) / K), function(i) xres_coef_mat[, keep[((i - 1) * K + 1):(i * K)]])
    xres_cov <- cov(residuals(xres_full))

    boot_worker <- function(b) {
        yb <- mlVAR::simulateVAR(pars = xres_coef_list, lags = 1:pfull,
                                 residuals = xres_cov, Nt = n)
        names(yb) <- varnames
        FCST_boot <- .get_recursive_VAR_fcsts(yb, pfull, n_train, n_test, dep, R2inv_full, ...)
        efullb <- yb[(n_train + 1):n, dep] - FCST_boot[1, ]
        eresb <- yb[(n_train + 1):n, dep] - FCST_boot[2, ]
        caustests(efullb, eresb)
    }

    BOOT0 <- if (bootparallel) {
        parallel::parSapply(cl, 1:B, FUN = boot_worker)
    } else {
        sapply(1:B, FUN = boot_worker)
    }

    FullH0 <- list(result = data.frame(MSEt = c(OBS["MSEt"], (sum(BOOT0["MSEt",] <= OBS["MSEt"]) + 1) / (B + 1)),
                                       MSEcor = c(OBS["MSEcor"], (sum(BOOT0["MSEcor",] <= OBS["MSEcor"]) + 1) / (B + 1)),
                                       row.names = c("stat_obs", "p_boot")),
                   p = pfull)

    return(list(Fast = Fast, FullH0 = FullH0))
}
