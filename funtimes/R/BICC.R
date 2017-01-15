BICC <- function(X, Alpha=NULL, Beta=NULL, Theta=0.8, p, w, s) 
{
  i <- dim(X)
  N <- i[2] #number of time series
  n <- i[1] #length of time series
  IQRx <- median(apply(X, 2, IQR))
  DELTA <- seq(IQRx/50, IQRx/2, IQRx/10)
  EPS <- seq(1.0/w, 1.0, 1.0/w)
  IC <- array(NA, dim=c(length(DELTA), length(EPS)))
  for (i in 1:(length(DELTA))){ 
    for(j in 1:length(EPS)){
      outputTMP <- CWindowCluster(X, Delta=DELTA[i], Epsilon=EPS[j], Alpha=Alpha, Beta=Beta, Theta=Theta, p=p, w=w, s=s)
      k <- max(outputTMP) 
      Res <- sapply(1:k, function(v) as.matrix(X[,v==outputTMP])-rowMeans(as.matrix(X[,v==outputTMP])))
      VarRes <- var(as.vector(unlist(Res)))*(n*N-1)/(n*N) #biased variance estimate for BIC
      IC[i,j] <- N*log(VarRes) + k*log(N) #BIC
    }
  }
  IC[IC==Inf | IC==-Inf] <- NA
  i <- which(IC==min(IC, na.rm=TRUE), arr.ind=TRUE)
  DELTA_opt <- DELTA[i[1,1]]
  EPSILON_opt <- EPS[i[1,2]]
  output <- CWindowCluster(X, Delta=DELTA_opt, Epsilon=EPSILON_opt, Alpha=Alpha, Beta=Beta, Theta=Theta, p=p, w=w, s=s)
  dimnames(output) <- list(paste("Window", c(1:nrow(output)), sep="_"), colnames(X))
  return(list(Delta.optimal=DELTA_opt, Epsilon.optimal=EPSILON_opt, Clusters=output, IC=IC, Delta.all=DELTA, Epsilon.all=EPS))
}
