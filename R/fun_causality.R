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
    tmp1 <- efull - eres
    tmp2 <- efull + eres
    mlm <- lm(tmp1 ~ tmp2 - 1)
    # ENC-NEW from (3) in Clark and McCracken (2001)
    P <- length(efull)
    L2full <- sum(efull^2)
    if (!is.finite(L2full) || L2full <= 0)
        stop("Sum of squared efull values must be positive.")
    EN <- P * sum(eres^2 - eres * efull) / L2full
    # OOS-F from (3) in McCracken (2007)
    OOSF <- -P * sum(dt) / L2full
    # output
        c(MSEt = summary(mlmt)$coefficients[1, "t value"],
            MSEcor = summary(mlm)$coefficients[1, "t value"],
      OOSF = OOSF,
      EN = EN)
}
