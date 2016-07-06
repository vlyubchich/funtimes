BICC <- function(x, w, p, Alpha=NULL, Beta=NULL, Theta=0.8)
  
{
  s <- w
  N <- ncol(x)
  IQRx <- median(apply(x,2,IQR))
  DELTA <- seq(IQRx/100, IQRx, IQRx/N)
  EPS <- seq(0.0, 1.0, 1.0/w)
  bic <- matrix(NA, length(DELTA), length(EPS))
  for (i in 1:length(DELTA)){ 
    for(j in 1:length(EPS)){
      outputTMP <- CWindowCluster(x, Delta=DELTA[i], Theta=Theta, p=p, w=w, s=s, Epsilon=EPS[j], Alpha=Alpha, Beta=Beta)
      k <- max(outputTMP) 
      clusterRes <- sapply(1:k, function(v) as.matrix(x[,v==outputTMP])-rowMeans(as.matrix(x[,v==outputTMP])))
      res <- as.vector(unlist(clusterRes))
      bic[i,j] <- (N-k)*log(var(res))+k*log(N-k)
    }
  }
  i <- which(bic==min(bic, na.rm = T), arr.ind=T)
  DELTA_opt <- DELTA[i[1,1]]
  EPSILON_opt <- EPS[i[1,2]]
  output <- CWindowCluster(x, Delta=DELTA_opt, Theta=Theta, p=p, w=w, s=s, Epsilon=EPSILON_opt, Alpha=Alpha, Beta=Beta)
  print(list(Optimal.Delta = DELTA_opt, Optimal.Epsilon = EPSILON_opt, Clusters = output))
}


