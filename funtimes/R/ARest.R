ARest <- function(x, ar.order = NULL, ar.method = "HVK", BIC = TRUE)
{
  x <- as.vector(x)
  n <- length(x)
  if(is.null(ar.order)){
    ar.order <- round(10*log10(n))
  }
  bic <- rep(NA, ar.order+1)
  bic[1] <- n*log(var(x)) #no AR-filtering (ar.order==0)
  pheta <- numeric(0) #if no AR-filtering, otherwise will be redefined below  
  if(ar.order > 0){
    if(!BIC){ #BIC==FALSE, use fixed ar.order>0
      if(ar.method == "HVK"){
        pheta <- HVK(x, ar.order=ar.order)
      }else{
        a <- ar(x, aic=FALSE, order.max=ar.order, demean=TRUE, method=ar.method)
        pheta <- a$ar
      }
    }else{ #BIC-based filtering
      for(i in 2:length(bic)){
        if(ar.method == "HVK"){
          pheta0 <- HVK(x, ar.order=i-1)
        }else{
          a <- ar(x, aic=FALSE, order.max=i-1, demean=TRUE, method=ar.method)
          pheta0 <- a$ar
        }
        tmp <- filter(x, pheta0, sides=1)
        et <- x[i:n] - tmp[(i-1):(n-1)]
        bic[i] <- n*log(var(et)) + i*log(n) #here use i = ARorder p + 1 (variance)
      }
      if(which.min(bic)>1){
        if(ar.method == "HVK"){
          pheta <- HVK(x, ar.order=(which.min(bic)-1))
        }else{
          a <- ar(x, aic=FALSE, order.max=(which.min(bic)-1), demean=TRUE, method=ar.method)
          pheta <- a$ar
        }
      }
    }
  }
  return(pheta)
}
