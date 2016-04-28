q.tails <- function(x, y, q=0.99){
  n <- length(x)*(1-q)
  k <- ceiling(log2(n) + 1)
  d <- (1-q)/k
  m1 <- sapply(1:k, function(i) mean(quantile(x, probs=c((q+(i-1)*d), (q+i*d)))))
  m2 <- sapply(1:k, function(i) mean(quantile(y, probs=c((q+(i-1)*d), (q+i*d)))))
  Pk <- m2 - m1	
  return(list(d=d, Pk=Pk))
}
