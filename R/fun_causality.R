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
    MSEt <- summary(mlmt)$coefficients[3]
    # MSEcor test (or MSEreg)
    tmp1 <- efull - eres
    tmp2 <- efull + eres
    mlm <- lm(tmp1 ~ tmp2 - 1)
    MSEcor <- summary(mlm)$coefficients[3]
    # output
    c(MSEt = MSEt, MSEcor = MSEcor)
}
