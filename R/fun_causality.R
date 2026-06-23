# Internal Helper Functions for Granger Causality Testing
#
# This file contains internal (non-exported) helper functions for Granger causality
# testing with restricted VAR models. These functions are used by causality_pred()
# and causality_predVAR().
#
# Functions:
#   - restrictions(): Creates/overlays restriction matrix for VAR models
#   - caustests(): Computes multiple causality test statistics
#
# @keywords internal
# @noRd

# Function to recreate matrix of restrictions and overlay new restrictions for causality testing
restrictions <- function(x, cause) {
    if (missing(cause) || length(cause) < 1L || any(is.na(cause)))
        stop("'cause' must contain at least one non-missing variable name.")
    cause <- as.character(cause)

    co.names <- vars::Bcoef(x)
    if (is.null(colnames(co.names)) || is.null(rownames(co.names)))
        stop("Could not determine coefficient names for restrictions().")

    k <- which(gsub("\\.l\\d+", "", colnames(co.names)) %in% cause) # select cause regressors
    l <- which(rownames(co.names) %in% cause) # select cause regressand
    if (length(k) == 0L)
        stop("No lagged regressors matched 'cause'.")
    if (length(l) == 0L)
        stop("No response variables matched 'cause'.")

    R2inv <- matrix(1, ncol = ncol(co.names), nrow = nrow(co.names))
    R2inv[-l, k] <- 0 # select coef to be tested
    dimnames(R2inv) <- dimnames(co.names)

    # If the model already has restriction, overlay with the new ones
    if (!is.null(x$restrictions)) {
        xr <- x$restrictions
        if (is.null(rownames(xr)) || is.null(colnames(xr)))
            stop("Existing restriction matrix must have row and column names.")
        if (!all(rownames(co.names) %in% rownames(xr)) ||
            !all(colnames(co.names) %in% colnames(xr))) {
            stop("Existing restriction matrix is incompatible with model coefficients.")
        }
        # match positions of variables
        xr <- xr[rownames(co.names), colnames(co.names)]
        if (!identical(dim(xr), dim(R2inv)))
            stop("Existing restriction matrix has incompatible dimensions.")
        # overlay
        R2inv <- xr * R2inv
    }
    R2inv
}

# Function to calculate causality statistics
caustests <- function(efull, eres) {
    if (!is.numeric(efull) || !is.numeric(eres))
        stop("efull and eres must be numeric vectors.")
    if (length(efull) == 0L || length(eres) == 0L)
        stop("efull and eres must be non-empty.")
    if (length(efull) != length(eres))
        stop("efull and eres must have the same length.")
    if (any(is.na(efull)) || any(is.na(eres)))
        stop("efull and eres must not contain missing values.")

    # MSEt test
    dt <- efull^2 - eres^2
    mlmt <- lm(dt ~ 1)
    # MSEcor test (or MSEreg)
    diff_res <- efull - eres
    sum_res <- efull + eres
    mlm <- lm(diff_res ~ sum_res - 1)
    # ENC-NEW from (3) in Clark and McCracken (2001)
    P <- length(efull)
    L2full <- sum(efull^2)
    if (!is.finite(L2full) || L2full <= 0)
        stop("Sum of squared efull values must be positive.")
    EN <- P * sum(eres^2 - eres * efull) / L2full
    # OOS-F from (3) in McCracken (2007)
    OOSF <- -P * sum(dt) / L2full
    # output
    stats <- c(
        summary(mlmt)$coefficients[1, "t value"],
        summary(mlm)$coefficients[1, "t value"],
        OOSF,
        EN
    )
    names(stats) <- c("MSEt", "MSEcor", "OOSF", "EN")
    stats
}

.select_lags_ar <- function(y, cause, dep, n_train, lag.max, lag.restrict, p.free, k) {
    if (length(lag.max) == 1) {
        lag.max <- rep(lag.max, 2)
    }
    maxl <- max(lag.max)
    if (lag.restrict >= lag.max[2]) {
        warning("lag.restrict >= lag.max or lag.max[2]. Using lag.restrict = 0 instead.")
        lag.restrict <- 0
    }

    # Embed both X and Y up to maxl to have the results of the same length
    lagX_train <- if (lag.restrict > 0) {
        embed(y[1:n_train, cause], maxl + 1)[, -c(1:(lag.restrict + 1)), drop = FALSE]
    } else {
        embed(y[1:n_train, cause], maxl + 1)[, -1, drop = FALSE]
    }
    lagY_train <- embed(y[1:n_train, dep], maxl + 1)

    # Information criterion for the model
    if (p.free) {
        best_ic <- Inf
        p <- c(NA, NA)
        for (p1 in 1:lag.max[1]) {
            for (p2 in (lag.restrict + 1):lag.max[2]) {
                fit <- stats::lm.fit(x = cbind(1,
                                               lagY_train[, 2:(p1 + 1), drop = FALSE],
                                               lagX_train[, 1:(p2 - lag.restrict), drop = FALSE]),
                                     y = lagY_train[, 1])
                # see stats:::extractAIC.lm() but omit the scale option
                nfit <- length(fit$residuals)
                edf <- nfit - fit$df.residual
                RSS <- sum(fit$residuals^2, na.rm = TRUE)
                dev <- nfit * log(RSS / nfit)
                fit_ic <- dev + k * edf
                if (fit_ic < best_ic) {
                    best_ic <- fit_ic
                    p <- c(p1, p2)
                }
            }
        }
    } else {
        IC <- sapply((lag.restrict + 1):lag.max[2], function(s) {
            fit <- stats::lm.fit(x = cbind(1,
                                           lagY_train[, 2:(s + 1), drop = FALSE],
                                           lagX_train[, 1:(s - lag.restrict), drop = FALSE]),
                                 y = lagY_train[, 1])
            # see stats:::extractAIC.lm; but omit the scale option
            nfit <- length(fit$residuals)
            edf <- nfit - fit$df.residual
            RSS <- sum(fit$residuals^2, na.rm = TRUE)
            dev <- nfit * log(RSS/nfit)
            dev + k * edf
        })
        p_val <- which.min(IC) + lag.restrict
        p <- rep(p_val, 2)
    }
    return(p)
}
