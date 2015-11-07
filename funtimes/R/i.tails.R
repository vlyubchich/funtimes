i.tails <- function(x, y, d=NULL){
  if(is.null(d)) d <- quantile(x, probs=0.99)
  n <- sum(x >= d)
  #Sturges' formula
  k <- ceiling(log2(n) + 1)
  width <- (max(x) - d) / k
  MAX <- max(c(x,y))
  cutpoints <- seq(from=d, to=MAX+width, by=width)
  nx <- summary(cut(x[x>=d], cutpoints, include.lowest=TRUE))
  ny <- summary(cut(y[y>=d], cutpoints, include.lowest=TRUE))
  Nk <- ny - nx
  Ck <- cutpoints[-1] - width/2
  return(list(Nk=Nk, Ck=Ck))
}