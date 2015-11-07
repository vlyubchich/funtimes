q.tails <- function(x, y, q=0.99){
  n <- length(x)*(1-q)
  k <- ceiling(log2(n) + 1)
  d <- (1-q)/k
  m1 <- sapply(1:k, function(i) quantile(x, probs=(q+i*d)) + quantile(x, probs=(q+(i-1)*d)))
  m2 <- sapply(1:k, function(i) quantile(y, probs=(q+i*d)) + quantile(y, probs=(q+(i-1)*d)))
  Pk <- (m2-m1)/2	
  return(list(d=d, Pk=Pk))
}