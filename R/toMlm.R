# Internal VAR to MLM Conversion Function (Modified from vars package)
#
# This is a modified copy of toMlm() methods from the vars package by Bernhard Pfaff.
# The modifications were submitted via pull request on 2021-08-22 but were not merged:
# https://github.com/bpfaff/vars/pull/10/files
#
# The vars package v.1.5-6 (2021-09-17) and later versions do not include these changes.
# This modified version is maintained internally to support conversion of restricted
# VAR models to multivariate linear models (MLM) for coefficient estimation and testing.
#
# Key Modifications (marked with "#VL" comments):
# - Support for x$restrictions matrix in varest objects
# - Proper handling of restricted coefficients in conversion to MLM
# - Individual equation fitting respecting coefficient restrictions
#
# @keywords internal
# @noRd

toMlm <- function(x, ...){
    UseMethod("toMlm")
}

#' @exportS3Method toMlm default
toMlm.default <- function(x, ...){
    lm(x$model)
}

#' @exportS3Method toMlm varest
toMlm.varest <- function(x, ...){
    ix <- 1:x$K
    X <- x$datamat
    int <- "-1" #automatic intercept is removed by default
    xr <- x$restrictions #VL
    if (x$type %in% c("const", "both")) { #VL
        #if x contains a constant, remove it from data X but allow to be added by lm automatically
        #remove constant in datamat
        X <- X[, -grep("const", colnames(X))]
        int <- ""
        if (!is.null(xr)) {
            xr <- xr[, -grep("const", colnames(xr))]
        }
    }
    #construct formula
    left <- paste(names(X)[ix], collapse = ",")
    fo <- as.formula(paste("cbind(", left, ") ~ .", int))
    #apply lm
    res <- eval(substitute(lm(fo, X), list(fo = fo))) #code suggested by Gabor Groothendick
    if (!is.null(xr)) { #VL
        for (i in ix) {
            foi <- as.formula(paste(names(X)[i], "~",
                                    paste(colnames(X)[-ix][xr[i,] == 1], collapse = "+"),
                                    int))
            mi <- eval(substitute(lm(fo, X), list(fo = foi)))
            #replace elements in the res object
            res$residuals[,i] <- mi$residuals
            res$effects[,i] <- mi$effects
            res$fitted.values[,i] <- mi$fitted.values
            #coefficients vectors are of different lengths, match them first
            #the code is bulky to avoid using other packages such as plyr::join
            rc <- data.frame(name = rownames(res$coefficients), rc = res$coefficients[,i])
            rc$id <- 1:nrow(rc)
            mc <- data.frame(name = names(mi$coefficients), mc = mi$coefficients)
            L <- merge(rc, mc, all.x = TRUE)
            L <- L[order(L$id),]
            L[is.na(L)] <- 0L
            res$coefficients[,i] <- as.numeric(L$mc)
        }
    }
    return(res)
}

coeftest.varest <- function(x, ...){
    lmtest::coeftest(toMlm.varest(x), ...)
}

bread.varest <- function(x, ...){
    sandwich::bread(toMlm.varest(x), ...)
}

#' @exportS3Method vcov varest
vcov.varest <- function(object, ...){
    stats::vcov(toMlm.varest(object), ...)
}

vcovHC.varest <- function(x, ...){
    sandwich::vcovHC(toMlm.varest(x), ...)
}

estfun.varest <- function(x, ...){
    sandwich::estfun(toMlm.varest(x), ...)
}
