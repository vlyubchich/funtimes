# My comments on GitHub to Ethan:
# - N should be identified from the input x.
# -library() exclude
# -add the DELTA_opt and EPSILON_opt to the output of this function.


BICC <- function(x, N = NULL, WINDOW = NULL, p = NULL, SHIFT = NULL, TREND = NULL)
  
{
  for (window in 1:length(WINDOW)){
    for (shift in 1:length(SHIFT)){
      s <- window  #step to shift a window
      #---------------------------------------------------------------------------------------------
      #Find IQR of Data(x) and Corresponding DELTA and the residual
      
      IQRx <- median(apply(x,2,IQR))
      DELTA <- seq(IQRx/100, IQRx, IQRx/N)
      EPS <- seq(0.0, 1.0, 1.0/WINDOW[window])
      
      #---------------------------------------------------------------------------------------------
      #Apply Clustering to Data(x) using the DELTA values above in order to obtain BIC
      #library(funtimes)
      
      bic <- matrix(NA, length(DELTA), length(EPS))
      for (i in 1:length(DELTA)){ #i=1
        for(j in 1:length(EPS)){
          outputTMP <- CWindowCluster(x, Delta=DELTA[i], Theta=0.8, p=p, w=WINDOW[window], s=s, Epsilon=EPS[j])
          k <- max(outputTMP) #number of clusters for given parameters
          clusterRes <- sapply(1:k, function(v) as.matrix(x[,v==outputTMP])-rowMeans(as.matrix(x[,v==outputTMP])))
          res <- as.vector(unlist(clusterRes))
          bic[i,j] <- (N-k)*log(var(res))+k*log(N-k)
        }
      }
      
      #Determine Optimal Delta and Epsilon for clustering based on BIC
      i <- which(bic==min(bic, na.rm = T), arr.ind=T)
      DELTA_opt <- DELTA[i[1,1]]
      EPSILON_opt <- EPS[i[1,2]]
      
      
      output <- CWindowCluster(x, Delta=DELTA_opt, Theta=0.8, p=p, w=WINDOW[window], s=s, Epsilon=EPSILON_opt)
      print(output)
      
    }
  }
}
