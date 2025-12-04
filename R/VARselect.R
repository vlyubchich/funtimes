# Internal VAR Order Selection Function (Modified from vars package)
#
# This is a modified copy of the VARselect() function from the vars package by Bernhard Pfaff.
# The modifications were submitted via pull request on 2021-08-22 but were not merged:
# https://github.com/bpfaff/vars/pull/10/files
#
# The vars package v.1.5-6 (2021-09-17) and later versions do not include these changes.
# This modified version is maintained internally to support restricted VAR estimation
# with lag restrictions for Granger causality testing.
#
# Key Modifications (marked with "#VL" comments):
# - Added lag.restrict parameter for excluding near-contemporaneous lags
# - Validation checks for lag.restrict parameter
# - Integration with modified VAR() function
#
# @keywords internal
# @noRd

VARselect <- function(y,
                      lag.max = 10,
                      lag.restrict = 0L, #VL
                      type = c("const", "trend", "both", "none"),
                      season = NULL,
                      exogen = NULL)
{
    y <- as.matrix(y)
    if (any(is.na(y)))
        stop("\nNAs in y.\n")
    colnames(y) <- make.names(colnames(y))
    K <- ncol(y)
    lag.max <- abs(as.integer(lag.max))
    #VL: check that lag.restrict is an N0 number below lag.max [0,1,...,lag.max)
    lag.restrict <- abs(as.integer(lag.restrict))
    if (lag.restrict >= lag.max) { #VL
        warning("lag.restrict >= lag.max. Using lag.restrict = 0 instead.")
        lag.restrict <- 0
    }
    type <- match.arg(type)
    lag <- abs(as.integer(lag.max + 1))
    ylagged <- embed(y, lag)[, -c(1:K)]
    yendog <- y[-c(1:lag.max), ]
    sample <- nrow(ylagged)
    rhs <- switch(type, const = rep(1, sample), trend = seq(lag.max + 1,
                                                            length = sample), both = cbind(rep(1, sample), seq(lag.max + 1, length = sample)), none = NULL)
    if (!(is.null(season))) {
        season <- abs(as.integer(season))
        dum <- (diag(season) - 1/season)[, -season]
        dums <- dum
        while (nrow(dums) < sample) {
            dums <- rbind(dums, dum)
        }
        dums <- dums[1:sample, ]
        rhs <- cbind(rhs, dums)
    }
    if (!(is.null(exogen))) {
        exogen <- as.matrix(exogen)
        if (!identical(nrow(exogen), nrow(y))) {
            stop("\nDifferent row size of y and exogen.\n")
        }
        if (is.null(colnames(exogen))) {
            colnames(exogen) <- paste("exo", 1:ncol(exogen),
                                      sep = "")
            warning(paste("No column names supplied in exogen, using:",
                          paste(colnames(exogen), collapse = ", "), ", instead.\n"))
        }
        colnames(exogen) <- make.names(colnames(exogen))
        rhs <- cbind(rhs, exogen[-c(1:lag.max), ])
    }
    idx <- seq(K, K * lag.max, K)
    if (!is.null(rhs)) {
        detint <- ncol(as.matrix(rhs))
    } else {
        detint <- 0
    }
    criteria <- matrix(NA, nrow = 4, ncol = lag.max - lag.restrict) #VL
    rownames(criteria) <- c("AIC(n)", "HQ(n)", "SC(n)", "FPE(n)")
    colnames(criteria) <- (lag.restrict + 1):lag.max #VL
    ii <- 1 #VL
    for (i in (lag.restrict + 1):lag.max) { #VL
        if (lag.restrict == 0) { #VL
            ys.lagged <- cbind(ylagged[, c(1:idx[i])], rhs)
            nstar <- ncol(ys.lagged)
            resids <- stats::lm.fit(x = ys.lagged, y = yendog)$residuals
        } else { #VL
            resids <- numeric()
            for (k in 1:K) {
                #VL: restrict to 0 coeffs 1,...,lag.restrict for all variables except k-th
                tmp <- rep(TRUE, K)
                tmp[k] <- FALSE #do not restrict for the current variable
                irestr <- which(rep(tmp, lag.restrict))
                ys.lagged <- cbind(ylagged[, c(1:idx[i])[-irestr]], rhs)
                resids <- cbind(resids, stats::lm.fit(x = ys.lagged, y = yendog[,k])$residuals)
            }
            nstar <- ncol(ys.lagged)
        }
        sigma.det <- det(crossprod(resids)/sample)
        criteria[1, ii] <- log(sigma.det) + (2/sample) * (i * K^2 + K * detint)
        criteria[2, ii] <- log(sigma.det) + (2 * log(log(sample))/sample) * (i * K^2 + K * detint)
        criteria[3, ii] <- log(sigma.det) + (log(sample)/sample) * (i * K^2 + K * detint)
        criteria[4, ii] <- ((sample + nstar)/(sample - nstar))^K * sigma.det
        ii <- ii + 1
    }
    order <- apply(criteria, 1, which.min) + lag.restrict #VL
    list(selection = order, criteria = criteria)
}
