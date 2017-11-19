q.tails <- function(x0, x1, q=0.99){
  n <- length(x0)*(1 - q)
  #Sturges' formula
  k <- ceiling(log2(n) + 1) 
  d <- (1 - q)/k
  m1 <- sapply(1:k, function(i) mean(quantile(x0, probs=c((q+(i-1)*d), (q+i*d)))))
  m2 <- sapply(1:k, function(i) mean(quantile(x1, probs=c((q+(i-1)*d), (q+i*d)))))
  Pk <- m2 - m1	
  return(list(d=d, Pk=Pk))
}
