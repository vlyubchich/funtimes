# Internal VAR Estimation Function (Modified from vars package)
#
# This is a modified copy of the VAR() function from the vars package by Bernhard Pfaff.
# The modifications were submitted via pull request on 2021-08-22 but were not merged:
# https://github.com/bpfaff/vars/pull/10/files
#
# The vars package v.1.5-6 (2021-09-17) and later versions do not include these changes.
# This modified version is maintained internally to support restricted VAR estimation
# and Granger causality testing with lag restrictions.
#
# Key Modifications (marked with "#VL" comments):
# - Added lag.restrict parameter to exclude near-contemporaneous lags
# - Enhanced integration with causality_pred() and causality_predVAR() functions
# - Support for restricted coefficient matrices in VAR models
#
# @keywords internal
# @noRd

VAR <- function(y, p = 1, type = c("const", "trend", "both", "none"),
    season = NULL, exogen = NULL, lag.max = NULL,
    lag.restrict = 0L, #VL
    ic = c("AIC", "HQ", "SC", "FPE"))
{
  y <- as.matrix(y)
    if (any(is.na(y)))
        stop("\nNAs in y.\n")
    if (ncol(y) < 2)
        stop("The matrix 'y' should contain at least two variables. For univariate analysis consider ar() and arima() in package stats.\n")
    if (is.null(colnames(y))) {
        colnames(y) <- paste("y", 1:ncol(y), sep = "")
        warning(paste("No column names supplied in y, using:",
            paste(colnames(y), collapse = ", "), ", instead.\n"))
    }
    colnames(y) <- make.names(colnames(y))
    y.orig <- y
    type <- match.arg(type)
    obs <- dim(y)[1]
    K <- dim(y)[2]
    if (!is.null(lag.max) && is.null(p)) { #VL
      lag.max <- abs(as.integer(lag.max))
      ic <- paste(match.arg(ic), "(n)", sep = "")
      p <- VARselect(y, lag.max = lag.max, lag.restrict = lag.restrict, #VL
                     type = type, season = season, exogen = exogen)$selection[ic]
    }
    if (p <= lag.restrict) { #VL
      warning("lag.restrict >= p. Using lag.restrict = 0 instead.")
      lag.restrict <- 0
    }
    sample <- obs - p
    ylags <- embed(y, dimension = p + 1)[, -(1:K)]
    temp1 <- NULL
    for (i in 1:p) {
        temp <- paste(colnames(y), ".l", i, sep = "")
        temp1 <- c(temp1, temp)
    }
    colnames(ylags) <- temp1
    yend <- y[-c(1:p), ]
    if (type == "const") {
        rhs <- cbind(ylags, rep(1, sample))
        colnames(rhs) <- c(colnames(ylags), "const")
    }
    else if (type == "trend") {
        rhs <- cbind(ylags, seq(p + 1, length = sample))
        colnames(rhs) <- c(colnames(ylags), "trend")
    }
    else if (type == "both") {
        rhs <- cbind(ylags, rep(1, sample), seq(p + 1, length = sample))
        colnames(rhs) <- c(colnames(ylags), "const", "trend")
    }
    else if (type == "none") {
        rhs <- ylags
        colnames(rhs) <- colnames(ylags)
    }
    if (!(is.null(season))) {
        season <- abs(as.integer(season))
        dum <- (diag(season) - 1/season)[, -season]
        dums <- dum
        while (nrow(dums) < obs) {
            dums <- rbind(dums, dum)
        }
        dums <- dums[1:obs, ]
        colnames(dums) <- paste("sd", 1:ncol(dums), sep = "")
        rhs <- cbind(rhs, dums[-c(1:p), ])
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
        tmp <- colnames(rhs)
        rhs <- cbind(rhs, exogen[-c(1:p), ])
        colnames(rhs) <- c(tmp, colnames(exogen))
    }
    datamat <- as.data.frame(rhs)
    colnames(datamat) <- colnames(rhs)
    equation <- list()
    if (lag.restrict == 0) { #VL
      restrictions = NULL
    } else { #VL
      restrictions <- matrix(1, K, ncol(rhs))
      dimnames(restrictions) <- list(colnames(y.orig), colnames(rhs))
    }
    for (i in 1:K) {
      y <- yend[, i]
      if (lag.restrict == 0) { #VL
        equation[[colnames(yend)[i]]] <- lm(y ~ -1 + ., data = datamat)
      } else { #VL
        tmp <- rep(TRUE, K)
        tmp[i] <- FALSE #do not restrict for the current variable
        irestr <- which(rep(tmp, lag.restrict))
        restrictions[i, irestr] <- 0
        equation[[colnames(yend)[i]]] <- lm(y ~ -1 + ., data = datamat[,-irestr])
      }
      if (any(c("const", "both") %in% type)) {
        attr(equation[[colnames(yend)[i]]]$terms, "intercept") <- 1
      }
    }
  call <- match.call()
  if ("season" %in% names(call)) call$season <- eval(season)
    result <- list(varresult = equation, datamat = data.frame(cbind(yend,
        rhs)), y = y.orig, type = type, p = p, K = K, obs = sample,
        totobs = sample + p,
        restrictions = restrictions, #VL
        call = call)
    class(result) <- "varest"
    return(result)
}
