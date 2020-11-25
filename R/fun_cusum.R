#Auxiliary functions related to the modified CUSUM statistic 
#by Horvath et al. (2017) for testing changepoints.


##### The function for M(k1,...,km) on p. 553:
Mfun <- function(e, k) {
    T <- length(e)
    k <- unique(sort(k))
    m <- length(k)
    sume <- sum(e)
    m1 <- (sum(e[1:k[1]]) - k[1]*sume/T) / sqrt(k[1])
    mm <- (sum(e[(k[m] + 1):T]) - sume*(T - k[m])/T) / sqrt(T - k[m])
    mi <- numeric()
    if (m > 1) {
        for (i in 2:m) {
            mi[i - 1] <- (sum(e[(k[i - 1] + 1):k[i]]) - sume*(k[i] - k[i - 1])/T) / sqrt(T)
        }
    }
    sum(abs(c(m1, mi, mm)))
}

##### The function for MT on p. 554, or as the first equation on p. 568 (if k = NULL):
MTfun <- function(e, m, k = NULL) {
    if (is.null(k)) { #explore all combinations of at-most-m change points
        T <- length(e)
        M <- K <- as.list(rep(NA, m))
        for (i in 1:m) {
            K[[i]] <- combn(T - 1, i)
            M[[i]] <- sapply(1:ncol(K[[i]]), function(x) Mfun(e, K[[i]][,x]))
        }
    } else {#explore only the pre-defined k's
        k <- unique(sort(k))
        m <- length(k)
        M <- K <- as.list(rep(NA, m))
        if (m == 1) {
            K[[1]] <- k
            M[[1]] <- Mfun(e, k)
        } else {
            for (i in 1:m) {
                K[[i]] <- combn(k, i)
                M[[i]] <- sapply(1:ncol(K[[i]]), function(x) Mfun(e, K[[i]][,x]))
            }
        }
    }
    mhat <- which.max(sapply(M, max))
    tmp <- which.max(M[[mhat]])
    if (m == 1) {
        khat <- k
    } else {
        khat <- K[[mhat]][,tmp]
    }
    MT <- M[[mhat]][tmp]
    return(list(MT = MT, m = mhat, k = khat))
}