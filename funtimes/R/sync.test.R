sync.test <- function(formula, B = 1000, Window = NULL, q = NULL, j = NULL, 
                      ar.order = NULL, ar.method = "HVK", BIC = TRUE)
{
  frml <- deparse(substitute(formula))
  splt <- strsplit(frml, "~")[[1]]
  DNAME <- splt[1]
  sh <- splt[2]
  X <- eval(parse(text = DNAME))
  n <- nrow(X)
  K <- ncol(X)
  t <- c(1:n)/n
  if(!is.null(Window)){ #if user set Window
    UseOneWindowPerTS <- TRUE
    ONEwindow <- FALSE
    if(length(Window)==1){
      ONEwindow <- TRUE
      Window <- rep(Window, K)
    }else if(length(Window) != K){
      stop("number of windows does not match number of time series.")
    }
    if(!is.null(q)){warning("The parameter q was not used.")}
    if(!is.null(j)){warning("The parameter j was not used.")}
  }else{
    UseOneWindowPerTS <- FALSE
    if(!is.null(q)){
      if (NCOL(q) > 1 | !is.numeric(q) | NROW(q) > 1){
        stop("q is not a scalar.")
      }
      if (q >= 1 | q <= 0){
        stop("q is out of range from 0 to 1.")
      }
    }else{
      q <- 3/4
    }
    if(!is.null(j)){
      if (!is.vector(j) | !is.numeric(j)) {
        stop("j is not a numeric vector.")
      }
    }else{
      j <- c(8:11)
    }
    kn <- n*q^j
    kn <- unique(sort(floor(kn)))
    kn <- kn[kn > 2 & kn < n]
    if (length(kn) == 0) {
      stop("set proper q and/or j.")
    }
  } 
  if(!is.null(ar.order)){ #if user set ar.order
    maxARorder <- ar.order
    if(length(ar.order)==1){
      maxARorder <- rep(ar.order, K)
    }else if(length(ar.order)!=K){
      stop("number of elements in ar.order does not match number of time series.")
    }
  }else{
    maxARorder <- rep(round(10*log10(n)), K)
  }
  #Preallocate space:
  if(!UseOneWindowPerTS){
    s <- array(NA, dim=c(length(kn), B, K))
    wavk_obs_all <- matrix(NA, length(kn), K)        
  }
  wavk_boot_opt <- array(NA, c(B, K))
  wavk_obs <- rep(NA, K)
  sigma <- rep(NA, K)
  OutputARorder <- matrix(NA, 1, K, dimnames=list("ar.order", dimnames(X)[[2]]))
  OutputWindow <- matrix(NA, 1, K, dimnames=list("Window", dimnames(X)[[2]]))
  #Function:
  X <- scale(X)
  AveragedProcess <- apply(X, 1, mean)
  mod <- lm(as.formula(paste("AveragedProcess", sh, sep = "~"))) #common trend
  TrendCoeff <- summary(mod)$coefficients
  U <- demean(X - mod$fitted) #detrended time series
  for (k in 1:K){
    pheta <- ARest(U[,k], ar.order=ar.order, ar.method=ar.method, BIC=BIC)
    OutputARorder[1,k] <- length(pheta)
    if (length(pheta)>0) {
      tmp <- filter(X[,k], pheta, sides=1)
      tmp2 <- filter(mod$fitted, pheta, sides=1)
      Z <- (X[(length(pheta)+1):n,k] - tmp[length(pheta):(n - 1)]) - (mod$fitted[(length(pheta)+1):n] - tmp2[length(pheta):(n - 1)])
    } else {
      Z <- U[,k]
    }
    Z <- Z - mean(Z)
    sigma[k] <- sqrt(sum(diff(Z)^2)/(2*(length(Z)-1)))
    boot <- array(data = rnorm(n*B), c(n,B))*sigma[k]
    if(!UseOneWindowPerTS){
      for (i in 1:length(kn)){
        s[i,,k] <- apply(boot, 2, function(x) WAVK(x, kn=kn[i])$Tn/sqrt(kn[i]))
      }
      if (length(kn)>2){
        s1 <- t(apply(s[,,k], 1, sort))
        distance <- sapply(1:(length(kn)-1), function(x) dist(s1[x:(x+1),]))
        OutputWindow[1,k] <- kn[which.min(distance)]
        wavk_boot_opt[,k] <- s[which.min(distance),,k]
      }else{
        OutputWindow[1,k] <- kn[1]
        wavk_boot_opt[,k] <- s[1,,k]
      }
      wavk_obs[k] <- WAVK(Z, kn=OutputWindow[1,k])$Tn/sqrt(OutputWindow[1,k])
      wavk_obs_all[,k] <- sapply(kn, function(x) WAVK(Z, kn=x)$Tn/sqrt(x))
    }else{
      wavk_boot_opt[,k] <- apply(boot, 2, function(x) WAVK(x, kn=Window[k])$Tn/sqrt(Window[k]))
      OutputWindow[1,k] <- Window[k]
      wavk_obs[k] <- WAVK(Z, kn=OutputWindow[1,k])$Tn/sqrt(OutputWindow[1,k])
    }
  } #k=K
  
  #p-value for bootstrap with optimal window selected
  STATISTIC <- sum(wavk_obs)
  crit <-  sum(STATISTIC > apply(wavk_boot_opt, 1, sum))/B
  if (crit < 0.5) {
    P.VALUE <- 2*crit
  } else {
    P.VALUE  <- 2*(1 - crit)
  }
  
  if(!UseOneWindowPerTS){
    ST <- apply(wavk_obs_all, 1, sum)
    tmp_all <- apply(s, c(1,2), sum)
    crit.boot <- sapply(c(1:length(kn)), function(x) sum(ST[x]<tmp_all[x,]))/B
    p.value.boot.all <- 2 * crit.boot
    p.value.boot.all[crit.boot>0.5] <- 2 * (1 - crit.boot[crit.boot>0.5])
    #Asymptotic results
    StAssStandard <- ST*sqrt(n)/sqrt(4*sum(sigma^4)/3) 
    crit.ass <- pnorm(StAssStandard, mean = 0, sd = 1)
    p.value.ass <- crit.ass * 2
    p.value.ass[crit.ass>0.5] <- (1 - crit.ass[crit.ass>0.5]) * 2
    #
    ESTIMATE <- list(TrendCoeff, OutputARorder, OutputWindow, cbind(kn, ST, p.value.boot.all, p.value.ass))
    names(ESTIMATE) <- list("common_trend_estimates", "ar.order_used", "Window_used", "all_considered_windows")
    dimnames(ESTIMATE[[4]]) <- list(rep("", NROW(ESTIMATE[[4]])), c("Window", "Statistic", "p-value", "Asympt. p-value"))
  }else{
    p.value.boot.all <- P.VALUE
    ST <- sum(wavk_obs)
    #Asymptotic results
    StAssStandard <- ST*sqrt(n)/sqrt(4*sum(sigma^4)/3) 
    crit.ass <- pnorm(StAssStandard, mean = 0, sd = 1)
    p.value.ass <- crit.ass * 2
    p.value.ass[crit.ass>0.5] <- (1 - crit.ass[crit.ass>0.5]) * 2
    #
    if(ONEwindow){
      ESTIMATE <- list(TrendCoeff, OutputARorder, OutputWindow, cbind(Window[1], ST, p.value.boot.all, p.value.ass))
      names(ESTIMATE) <- list("common_trend_estimates", "ar.order_used", "Window_used", "all_considered_windows")
      dimnames(ESTIMATE[[4]]) <- list(rep("", NROW(ESTIMATE[[4]])), c("Window", "Statistic", "p-value", "Asympt. p-value"))
    }else{
      ESTIMATE <- list(TrendCoeff, OutputARorder, OutputWindow, cbind(ST, p.value.boot.all, p.value.ass))
      names(ESTIMATE) <- list("common_trend_estimates", "ar.order_used", "Window_used", "all_considered_windows")
      dimnames(ESTIMATE[[4]]) <- list(rep("", NROW(ESTIMATE[[4]])), c("Statistic", "p-value", "Asympt. p-value"))
    }
  }
  
  METHOD <- "Non-parametric test for synchronism of parametric trends"   
  names(STATISTIC) <- "Test statistic"
  ALTERNATIVE <- paste("common trend is not of the form ", frml, ".", sep="")
  structure(list(method = METHOD, data.name = DNAME, statistic = STATISTIC,  p.value = P.VALUE,  alternative = ALTERNATIVE, estimate = ESTIMATE), class = "htest") 	
}
