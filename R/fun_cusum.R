#Auxiliary functions related to the modified CUSUM statistic 
#by Horvath et al. (2017) for testing changepoints.


##### The function for M(k1,...,km) on p. 553:
Mfun <- function(e, k) {
    T <- length(e)
    # AB: Disabled, since sorting is redundant and time-consuming -
    #     k already comes pre-sorted from the calling function.
    # k <- unique(sort(k))
    m <- length(k)
    sume <- sum(e)
    m1 <- (sum(e[1:k[1]]) - k[1]*sume/T) / sqrt(k[1])
    mm <- (sum(e[(k[m] + 1):T]) - sume*(T - k[m])/T) / sqrt(T - k[m])
    mi <- numeric()
    if (m > 1) {
        sT <- sqrt(T)
        for (i in 2:m) {
            mi[i - 1] <- (sum(e[(k[i - 1] + 1):k[i]]) - sume*(k[i] - k[i - 1])/T) / sT
        }
    }
    sum(abs(c(m1, mi, mm)))
}


# AB: same as Mfun, but only for at most one change point:
M1fun <- function(x, e) { 
    T <- length(e)
    me <- .Internal(mean(e)) # 'naked' mean function
    #me <- mean(e)
    m1 <- (sum(e[1:x]) - x*me) / sqrt(x)
    mm <- (sum(e[(x + 1):T]) - me*(T - x)) / sqrt(T - x)
    sum(abs(c(m1, mm)))
}


##### The function for MT on p. 554, or as the first equation on p. 568 (if k = NULL):
# AB: ignore the argument x, which is a dummy variable needed for sapply
MTfun <- function(e, m, k = NULL, x = NULL) {
    if (is.null(k)) { #explore all combinations of at-most-m change points
        T <- length(e)
        M <- K <- as.list(rep(NA, m))
        for (i in 1:m) {
            K[[i]] <- combn(T - 1, i)
            M[[i]] <- sapply(1:ncol(K[[i]]), function(x) Mfun(e, K[[i]][,x]))
        }
    } else {#explore only the pre-defined k's
        # AB: Sorting is redundant because it is already sorted:
        # k <- unique(sort(k))
        # AB: m can be independent of length(k):
        # m <- length(k)
        M <- K <- as.list(rep(NA, m))
        # AB: m can now be <length(k); in the case where both are 1,
        #     we can optimize the execution by avoiding matrix and sapply:
        if (m == 1) {
            if (length(k) == 1) {
                K[[1]] <- k
                M[[1]] <- M1fun(x = k, e = e)
            } else {
                K[[1]] <- matrix(k, nrow = 1)
                M[[1]] <- sapply(k, M1fun, e = e)
            }
        } else {
            K[[1]] <- matrix(k, nrow = 1)
            M[[1]] <- sapply(k, M1fun, e = e)
            for (i in 2:m) {
                K[[i]] <- combn(k, i)
                M[[i]] <- apply(K[[i]], 2, function(x) Mfun(e, x))
            }
        }
    }
    mhat <- which.max(sapply(M, max))
    tmp <- which.max(M[[mhat]])
    if (m == 1) {
        khat <- k[tmp]
        # AB: was:
        # khat <- k
    } else {
        khat <- K[[mhat]][,tmp]
    }
    MT <- M[[mhat]][tmp]
    list(MT = MT, m = mhat, k = khat)
}