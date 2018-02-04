crad <- function(X, Nbin, StepSize = 1, Rcov = FALSE){

  #require(rPython)
  #require(robustbase)

  if(Rcov){
    cov = robustbase::covMcd(X)$cov
  }else{
    cov = cov(X)
  }

  rPython::python.load("CRAD.py")
  res_cl = rPython::python.call("CRAD", X, StepSize, Nbin, cov)
  return(res_cl)
}

