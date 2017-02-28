notrend.test <- function(x, B = 1000, test = c("t", "MK", "WAVK"), 
                          ar.method = "HVK", ar.order = NULL, BIC = TRUE, 
                          factor.length = c("user.defined", "adaptive.selection"), 
                          Window = NULL, q = 3/4, j = c(8:11))
{
  ### Perform various checks.
  DNAME <- deparse(substitute(x))  
  if (NCOL(x) > 1 | !is.numeric(x)) {
    stop("x is not a vector or univariate time series.")
  }
  if (any(is.na(x))) {
    stop("x contains missing values.")
  }
  x <- as.vector(x)
  n <- length(x)
  test <- match.arg(test)
  factor.length <- match.arg(factor.length)  
  if (is.null(Window)) {
    Window = round(0.1*n)
  }
  if (NCOL(q) > 1 | !is.numeric(q) | NROW(q) > 1) {
    stop("q is not a scalar.")
  }
  if (q >= 1 | q <= 0) {
    stop("q is out of range from 0 to 1.")
  }
  if (!is.vector(j) | !is.numeric(j)) {
    stop("j is not a numeric vector.")
  }
  if (factor.length == "user.defined") {
    kn <- Window[1]
  } else {
    kn <- n*q^j
  }  
  kn <- unique(sort(floor(kn)))
  kn <- kn[kn > 2 & kn < n]
  if (length(kn) == 0) {
    stop("set a proper window.")
  }
  if (factor.length == "adaptive.selection" & length(kn) < 3) {
    stop("number of possible windows is not enough for adaptive selection. Change parameters 'q' and/or 'j'.")
  }
  B <- round(B)
  if (B <= 0) {
    stop("number of bootstrap samples B must be positive.")
  }
  if (!is.null(ar.order) & (NCOL(ar.order) > 1 | !is.numeric(ar.order) | NROW(ar.order) > 1)) {
    stop("ar.order is not a scalar.")
  }
  if (!is.null(ar.order) && ar.order < 0) {
    stop("ar.order must be non-negative.")
  }
  
  ### Function.
  Y <- array(data = NA, c(n, B))
  t <- c(1:n)/n
  pheta <- ARest(x, ar.order = ar.order, ar.method = ar.method, BIC = BIC)
  if (length(pheta) > 0) {
    names(pheta) <- paste(rep("phi_", length(pheta)), c(1:length(pheta)), sep = "")
    tmp <- filter(x, pheta, sides = 1)
    Z <- x[(length(pheta)+1):n] - tmp[length(pheta):(n-1)]
    for (i in 1:B){
      e <- sample(Z, size = n, replace = TRUE)
      Y[ ,i] <- arima.sim(list(order = c(length(pheta), 0, 0), ar = pheta), n = n, innov = e)
    }
  } else {
    Z <- x
    for (i in 1:B){
      Y[ ,i] <- sample(Z, size = n, replace = TRUE)
    }
  }
  Z <- na.omit(Z) - mean(Z)
  ESTIMATE <- list(length(pheta), pheta)
  names(ESTIMATE) <- c("AR_order", "AR_coefficients")
  
  #If Student's t-test is used
  if(test == "t"){
    METHOD <- "Sieve-bootstrap Student's t-test for a linear trend"
    ALTERNATIVE <- "linear trend."
    STATISTIC <- summary(lm(x ~ t))$coefficients["t", "t value"]
    names(STATISTIC) <- "Student's t value"
    boot.stat <- sapply(1:dim(Y)[2], function(i) summary(lm(Y[,i] ~ t))$coefficients["t", "t value"])
  }
  #If Mann-Kendall's test is used
  if(test == "MK"){
    METHOD <- "Sieve-bootstrap Mann-Kendall's trend test"
    ALTERNATIVE <- "monotonic trend."
    STATISTIC <- MannKendall(x)$tau
    names(STATISTIC) <- "Mann-Kendall's tau"
    boot.stat <- sapply(1:dim(Y)[2], function(i) MannKendall(Y[,i])$tau)
  }
  #If WAVK test is used
  if(test == "WAVK"){
    METHOD <- "Sieve-bootstrap WAVK trend test"
    ALTERNATIVE <- "(non-)monotonic trend."
    if (length(kn) < 3) {
      kn_opt <- kn[1]
      boot.stat <- sapply(1:dim(Y)[2], function(j) WAVK(Y[,j], kn_opt)$Tns)
    } else {
      s <- array(data=NA, c(length(kn), B))
      for (i in 1:length(kn)){
        s[i,] <- sapply(1:dim(Y)[2], function(j) WAVK(Y[,j], kn[i])$Tns)
      }
      s <- t(apply(s, 1, sort))
      distance <- sapply(1:(length(kn)-1), function(x) dist(s[x:(x+1),]))
      argmin <- which.min(distance)
      kn_opt <- kn[argmin]
      boot.stat <- s[argmin,]
    }
    STATISTIC <- WAVK(x, kn_opt)$Tns
    names(STATISTIC) <- "WAVK test statistic"
    PARAMETER <- kn_opt
    names(PARAMETER) <- "moving window"
  }
  
  P.VALUE <- mean(abs(boot.stat) >= abs(STATISTIC))
  if(test == "WAVK"){
    structure(list(method = METHOD, data.name = DNAME, statistic = STATISTIC, p.value = P.VALUE, 
                   alternative = ALTERNATIVE, estimate = ESTIMATE, parameter = PARAMETER), class = "htest") 
  } else {
    structure(list(method = METHOD, data.name = DNAME, statistic = STATISTIC, p.value = P.VALUE, 
                   alternative = ALTERNATIVE, estimate = ESTIMATE), class = "htest") 
  }
}


