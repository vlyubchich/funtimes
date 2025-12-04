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
    co.names <- vars::Bcoef(x)
    k <- which(gsub("\\.l\\d+", "", colnames(co.names)) %in% cause) # select cause regressors
    l <- which(rownames(co.names) %in% cause) # select cause regressand
    R2inv <- matrix(1, ncol = ncol(co.names), nrow = nrow(co.names))
    R2inv[-l, k] <- 0 # select coef to be tested
    # If the model already has restriction, overlay with the new ones
    if (!is.null(x$restrictions)) {
        xr <- x$restrictions
        # match positions of variables
        xr <- xr[rownames(co.names), colnames(co.names)]
        # overlay
        R2inv <- xr * R2inv
    }
    R2inv
}

# Function to calculate causality statistics
caustests <- function(efull, eres) {
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
    EN <- P * sum(eres^2 - eres * efull) / L2full
    # OOS-F from (3) in McCracken (2007)
    OOSF <- -P * sum(dt) / L2full
    # output
    c(MSEt = summary(mlmt)$coefficients[3],
      MSEcor = summary(mlm)$coefficients[3],
      OOSF = OOSF,
      EN = EN)
}
