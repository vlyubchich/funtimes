i.tails <- function(x0, x1, d=NULL){
  if(is.null(d)) d <- quantile(x0, probs=0.99)
  n <- sum(x0 >= d)
  #Sturges' formula
  k <- ceiling(log2(n) + 1)
  width <- (max(x0) - d) / k
  MAX <- max(c(x0, x1))
  cutpoints <- seq(from=d, to=MAX+width, by=width)
  nx0 <- summary(cut(x0[x0>=d], cutpoints, include.lowest=TRUE))
  nx1 <- summary(cut(x1[x1>=d], cutpoints, include.lowest=TRUE))
  Nk <- nx1 - nx0
  Ck <- cutpoints[-1] - width/2
  return(list(Nk=Nk, Ck=Ck))
}
