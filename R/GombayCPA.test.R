#' Change Point Detection in Autoregressive Time Series
#' 
#' The function detects change points in autoregressive (AR) models for time series. Changes 
#' can be detected in any of \code{p+2} (mean, var, phi) autoregressive parameters where \code{p} 
#' is the order of the AR model. The test statistic is based on the efficient score vector \insertCite{Gombay_2008}{funtimes}. 
#' 
#' @details The function allows for
#' testing for a temporary change and for a change in a specific model parameters. 
#' Critical values can be estimated via asymptotic distribution \code{"asymptotic"} (i.e., the
#' default option) or via sieve bootstrap \code{"bootstrap"}. The function employs internal 
#' function \code{change.point} and sieve bootstrap \code{change.point.sieve} function.
#' 
#' @param y A vector that contains univariate time series observations. Missing values are not allowed.
#' @param a.order Order of the autoregressive model which must be a nonnegative integer number. 
#' @param alternatives A string parameter that specifies a type of the test (i.e., "two-sided",
#' "greater", "lesser", and "temporary").  The option "temporary" examines the temporary change
#' in one of the parameters \insertCite{Gombay_2008}{funtimes}.
#' @param crit.type Method of obtaining critical values: "asymptotic" (default) or "bootstrap".
#' @param num.bootstrap Number of bootstrap replications if \code{"crit.type"} ="bootstrap". 
#' Default number is 1,000.

#'
#' 
#' @return A list with the following components:
#' \item{index}{Points of change for each parameter. The value of the \code{"alternatives"}
#'  determines the return: 
#' "temporary - returns max, min and abs.max points;
#' "greater" - returns max points;
#' "lesser" -  returns min points;
#' "two-sided" - returns abs.max. }
#' \item{stats}{Test statistic values for change points in: mean, var, phi.}
#' \item{p.values}{\code{p-value} of the change point test.}
#' 
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @seealso \code{\link{mcusum.test}}  change point test for regression and 
#' \code{\link[Ecdat]{terrorism}} dataset used in the Example 2
#' 
#' @keywords time series, ts
#' 
#' @author Palina Niamkova, Dorcas Ofori-Boateng, Yulia R. Gel
#' 
#' @export
#' @examples 
#' 
#' #Example 1:
#' 
#' #Simulate some time series:
#' 
#' series_1 = arima.sim(n = 100, list(order = c(2,0,0), ar = c(-0.7, -0.1)))
#' series_2 = arima.sim(n = 200, list(order = c(2,0,0), ar = c(0.1, -0.6)))
#' main_series = c(series_1, series_2)
#'
#' result11 = GombayCPA_test(series_1, 2, "two-sided")
#' result11 #== No change point ===#
#'
#' result12 = GombayCPA_test(main_series, 2, "two-sided")
#' result12  #=== One change at phi values ===#
#'
#' result13 = GombayCPA_test(main_series, 2, "two-sided", "bootstrap")
#' result13  #=== One change at phi values ===#
#'
#'
#' 
#' #Example 2:
#' 
#' #From the package 'Ecdat' consider a time series with annual world number of victims of 
#' #terrorism in the US from 1970 till 2016:
#'
#' c.data = (Ecdat::terrorism['nkill.us'])
#' nkill.us.ts<-ts(c.data,start=(1970), end=(2016), frequency(47))
#'
#' #Now let's perform a change point detection with one sided tests:
#'
#' GombayCPA_test(nkill.us.ts, 0, "lesser")
#' GombayCPA_test(nkill.us.ts, 0, "greater")
#' nkill.us.ts[32]
#' year=1970+31
#' print(year)
#' plot(nkill.us.ts)
#'
#' #In both cases we find that the change point is located at the position 31 or 32. We can 
#' # examine it further by checking the value of this position (using: nkill.us.ts[32]) as well as
#' # by plotting the graph (using: plot(nkill.us.ts)). The detected change point corresponds to 
#' #the year of 2001, when the 9/11 attack happened.
#'


GombayCPA_test = function(y, a.order, alternatives = c("two-sided", "greater", "lesser", "temporary"), crit.type = c("asymptotic", "bootstrap"), num.bootstrap=1000)
{
  
  
  change.point = function(y, a.order, alternatives = c("two-sided", "greater", "lesser", "temporary"))
  {
    if(a.order < 0) stop("a.order must be greater than or equal to 0.")
    

    n = length(y)
    
    #=== demeaning ===#
    y.dem = y - mean(na.omit(y))
    ar.p  = arima0(y.dem, order = c(a.order, 0, 0), include.mean = T, method = "ML") 
    # ar.p = ar(y.dem,order=a.order,aic=FALSE) #
    
    mu     = mean(na.omit(y.dem)) #=== this is equal to zero ===#
    sigma2 = ar.p$sigma2
    
    if(a.order > 0)
    {
      phi    = ar.p$coef
      gamma  = matrix(NA, ncol = a.order, nrow = a.order)
      cov    = acf(y, plot = F, type = "covariance")$acf[1:a.order]
      rownum = rep(c(1:a.order), a.order)
      colnum = rep(c(1:a.order), each = a.order)
      abdiff = abs(rownum - colnum) + 1
      gamma  = matrix(cov[abdiff], ncol = a.order, nrow = a.order)
    }
    
    muvec    = double(n)
    muvec[1] = y.dem[1]
    for(i in 2:n)
    {
      if(a.order > 0)
      {
        len      = min(a.order, (i - 1))
        muvec[i] = y.dem[i] - sum(phi[1:len]*y.dem[(i-1):(i-len)])
      }
      if(a.order == 0)
      {
        muvec[i] = y.dem[i]
      }
    }
    
    if(a.order > 0)
    {
      mustat = ((1-sum(phi))/sigma2)*cumsum(muvec)
    }
    if(a.order == 0)
    {
      mustat = 1/sigma2*cumsum(muvec)
    }
    
    sigmastat = -c(1:n)/(2*sigma2) + 1/(2*sigma2^2)*cumsum((muvec^2))
    
    if(a.order > 0)
    {
      phistat = matrix(NA, nrow = a.order, ncol = n)
      
      for(s in 1:a.order)
      {
        phivec.s    = muvec*c(rep(0,s), y.dem[1:(n-s)])
        phistat.s   = 1/sigma2*cumsum(phivec.s)
        phistat[s,] = phistat.s
      }
    }
    
    if(a.order > 0)
    {
      allstat = rbind(mustat, sigmastat, phistat)
    }
    if(a.order == 0)
    {
      allstat = rbind(mustat, sigmastat)
    }
    
    info = matrix(0, nrow = (a.order + 2), ncol = (a.order + 2))
    if(a.order > 0)
    {
      info[c(3:(a.order + 2)),c(3:(a.order + 2))] = 1/sigma2*gamma
      info[1,1]                       = (1/sigma2)*(1 - sum(phi))^2
    }
    if(a.order == 0)
    {
      info[1,1] = 1/sigma2
    }
    info[2,2] = 1/(2*(sigma2^2))
    svd.info  = svd(info)
    info.n    = svd.info$u%*%diag(svd.info$d^(-1/2))%*%t(svd.info$v)
    
    b = n^(-1/2)*info.n%*%allstat
    
    if(alternatives == "two-sided")
    {
      maxstats           = apply(abs(b), 1, max)
      maxstats           = rbind(maxstats)
      rownames(maxstats) = "abs.max"
      if(a.order > 0)
      {
        colnames(maxstats) = c("mean","var",paste("phi",c(1:a.order),sep=""))
      }
      if(a.order == 0)
      {
        colnames(maxstats) = c("mean","var")
      }
      points           = apply(abs(b[,-1]), 1, which.max) + 1 #=== this is for two-sided, with the first one trimmed ===#
      points           = rbind(points)
      rownames(points) = "abs.max"
      if(a.order > 0)
      {
        colnames(points) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(points) = c("mean","var")
      }
    }else if(alternatives == "greater")
    {
      maxstats           = apply(b, 1, max)
      maxstats           = rbind(maxstats)
      rownames(maxstats) = "max"
      if(a.order > 0)
      {
        colnames(maxstats) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(maxstats) = c("mean", "var")
      }
      maxpoints        = apply((b[,-1]), 1, which.max) + 1
      points           = rbind(maxpoints)
      rownames(points) = c("max")
      if(a.order > 0)
      {
        colnames(points) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(points) = c("mean", "var")
      }
    }else if(alternatives == "lesser")
    {
      maxstats           = apply(b, 1, min)
      maxstats           = rbind(maxstats)
      rownames(maxstats) = "min"
      if(a.order > 0)
      {
        colnames(maxstats) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(maxstats) = c("mean", "var")
      }
      maxpoints        = apply((b[,-1]), 1, which.min) + 1
      points           = rbind(maxpoints)
      rownames(points) = c("min")
      if(a.order > 0)
      {
        colnames(points) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(points) = c("mean", "var")
      }
    }else
    {
      alternatives       = "temporary"
      supstats           = apply(b, 1, max)
      infstats           = apply(b, 1, min)
      maxstats           = supstats - infstats
      maxstats           = rbind(maxstats)
      rownames(maxstats) = "diff"
      if(a.order > 0)
      {
        colnames(maxstats) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(maxstats) = c("mean", "var")
      }
      maxpoints        = apply((b[,-1]),1, which.max) + 1
      minpoints        = apply((b[,-1]),1, which.min) + 1
      abspoints        = apply(abs(b[,-1]), 1, which.max) + 1
      points           = rbind(maxpoints, minpoints, abspoints)
      rownames(points) = c("max", "min", "abs.max")
      if(a.order > 0)
      {
        colnames(points) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(points) = c("mean", "var")
      }
    }
    return(list(stats = maxstats, points = points))
  }
  
  ############## change.point function for bootstrapped time series ############## 
  
  change.point.sieve = function(y, a.order, alternatives = c("two-sided", "greater", "lesser", "temporary"), res)
  {
    if(a.order < 0) stop("a.order must be greater than or equal to 0.")
    
    n = length(y)
    
    ### demeaning ###
    y.dem = y - mean(na.omit(y))
    ar.p  = arima0(y.dem, order = c(a.order, 0, 0), include.mean = T, method = "ML") 
    # ar.p = ar(y.dem,order=a.order,aic=FALSE) #
    
    mu     = mean(na.omit(y.dem)) #== this is equal to zero ==#
    # sigma2=ar.p$sigma2 #
    sigma2 = var(res)
    
    if(a.order > 0)
    {
      phi    = ar.p$coef
      gamma  = matrix(NA, ncol = a.order, nrow = a.order)
      cov    = acf(y, plot = F, type = "covariance")$acf[1:a.order]
      rownum = rep(c(1:a.order), a.order)
      colnum = rep(c(1:a.order), each = a.order)
      abdiff = abs(rownum - colnum) + 1
      gamma  = matrix(cov[abdiff], ncol = a.order, nrow = a.order)
    }
    
    muvec = double(n)
    muvec[1] = y.dem[1]
    for(i in 2:n)
    {
      if(a.order > 0)
      {
        len      = min(a.order, (i-1))
        muvec[i] = y.dem[i]-sum(phi[1:len]*y.dem[(i-1):(i-len)])
      }
      if(a.order == 0)
      {
        muvec[i] = y.dem[i]
      }
    }
    
    if(a.order > 0)
    {
      mustat = ((1-sum(phi))/sigma2)*cumsum(muvec)
    }
    if(a.order == 0)
    {
      mustat = 1/sigma2*cumsum(muvec)
    }
    
    sigmastat = -c(1:n)/(2*sigma2) + 1/(2*sigma2^2)*cumsum((muvec^2))
    
    if(a.order > 0)
    {
      phistat = matrix(NA, nrow = a.order, ncol = n)
      
      for(s in 1:a.order)
      {
        phivec.s    = muvec*c(rep(0,s),y.dem[1:(n-s)])
        phistat.s   = 1/sigma2*cumsum(phivec.s)
        phistat[s,] = phistat.s
      }
    }
    
    if(a.order > 0)
    {
      allstat = rbind(mustat, sigmastat, phistat)
    }
    if(a.order == 0)
    {
      allstat = rbind(mustat, sigmastat)
    }
    
    info = matrix(0, nrow = (a.order+2), ncol = (a.order + 2))
    if(a.order > 0)
    {
      info[c(3:(a.order+2)),c(3:(a.order+2))] = 1/sigma2*gamma
      info[1,1]                   = (1/sigma2)*(1-sum(phi))^2
    }
    if(a.order == 0)
    {
      info[1,1] = 1/sigma2
    }
    info[2,2] = 1/(2*(sigma2^2))
    svd.info  = svd(info)
    info.n    = svd.info$u%*%diag(svd.info$d^(-1/2))%*%t(svd.info$v)
    
    b = n^(-1/2)*info.n%*%allstat
    
    if(alternatives == "two-sided")
    {
      maxstats           = apply(abs(b), 1, max)
      maxstats           = rbind(maxstats)
      rownames(maxstats) = "abs.max"
      if(a.order > 0)
      {
        colnames(maxstats) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(maxstats) = c("mean", "var")
      }
      points           = apply(abs(b[,-1]), 1, which.max) + 1 #== this is for two-sided, with the first one trimmed ==#
      points           = rbind(points)
      rownames(points) = "abs.max"
      if(a.order > 0)
      {
        colnames(points) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(points) = c("mean", "var")
      }
    }else if(alternatives == "greater")
    {
      maxstats           = apply(b, 1, max)
      maxstats           = rbind(maxstats)
      rownames(maxstats) = "max"
      if(a.order > 0)
      {
        colnames(maxstats) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(maxstats) = c("mean", "var")
      }
      maxpoints        = apply((b[,-1]), 1, which.max) + 1
      points           = rbind(maxpoints)
      rownames(points) = c("max")
      if(a.order > 0)
      {
        colnames(points) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(points) = c("mean", "var")
      }
    }else if(alternatives == "lesser")
    {
      maxstats           = apply(b, 1, min)
      maxstats           = rbind(maxstats)
      rownames(maxstats) = "min"
      if(a.order > 0)
      {
        colnames(maxstats) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(maxstats) = c("mean", "var")
      }
      maxpoints        = apply((b[,-1]), 1, which.min) + 1
      points           = rbind(maxpoints)
      rownames(points) = c("min")
      if(a.order > 0)
      {
        colnames(points) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(points) = c("mean", "var")
      }
    }else
    {
      alternatives       = "temporary"
      supstats           = apply(b, 1, max)
      infstats           = apply(b, 1, min)
      maxstats           = supstats - infstats
      maxstats           = rbind(maxstats)
      rownames(maxstats) = "diff"
      if(a.order > 0)
      {
        colnames(maxstats) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(maxstats) = c("mean", "var")
      }
      maxpoints        = apply((b[,-1]), 1, which.max) + 1
      minpoints        = apply((b[,-1]), 1, which.min) + 1
      abspoints        = apply(abs(b[,-1]), 1, which.max) + 1
      points           = rbind(maxpoints, minpoints, abspoints)
      rownames(points) = c("max", "min", "abs.max")
      if(a.order > 0)
      {
        colnames(points) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
      }
      if(a.order == 0)
      {
        colnames(points) = c("mean", "var")
      }
    }
    return(list(stats = maxstats, points = points))
  }
  
  
  #=== main part of the change.point.boot ===#
  alternatives = match.arg(alternatives)
  n            = length(y)
  y.dem        = y-mean(na.omit(y))
  ar.p         = arima0(y.dem,order = c(a.order, 0, 0), include.mean = T, method = "ML")
  if(a.order > 0){phi = ar.p$coef[1:a.order]
  }
  if(a.order == 0){phi = 0*ar.p$coef #--- for a.order = 0 ---#  
  }
  
  
  orig       = change.point(y, a.order, alternatives)
  origstats  = orig$stats
  origpoints = orig$points
  
  B = num.bootstrap
  crit.type = match.arg(crit.type)
  
  if(crit.type == "bootstrap"){
    
    residuals = na.omit(ar.p$residuals) - mean(na.omit(ar.p$residuals))
    
    bootstats = matrix(NA, nrow = (a.order + 2), ncol = B)
    
    for(i in 1:B)
    {
      e = sample(residuals, size = n + 100, replace = TRUE)
      
      if(a.order > 0)
      {
        y.sim = arima.sim(list(order = c(a.order, 0, 0), ar = phi), n = n+100, innov=e)[101:(n+100)]
      }
      if(a.order == 0)
      {
        y.sim = e
      }
      bootstats[,i] = change.point.sieve(y.sim, a.order, alternatives, e)$stats
    }
    
    bootpvals = double((a.order + 2))
    for(i in 1:(a.order + 2))
    {
      if(alternatives == "two-sided")
      {
        bootpvals[i]  = length(which(abs(bootstats[i,]) > abs(origstats[i])))/B
      }
      if(alternatives == "greater")
      {
        bootpvals[i]  = length(which(bootstats[i,] > origstats[i]))/B
      }
      if(alternatives == "lesser")
      {
        bootpvals[i]  = length(which(bootstats[i,] < origstats[i]))/B
      }
      if(alternatives == "temporary")
      {
        bootpvals[i] = length(which(bootstats[i,] > origstats[i]))/B
      }
    }
    
    
    bootpvals           = rbind(bootpvals)
    rownames(bootpvals) = "p.value"
    if(a.order > 0)
    {
      colnames(bootpvals) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
    }
    if(a.order == 0)
    {
      colnames(bootpvals) = c("mean", "var")
    }
    return(list(index = origpoints, stats = origstats, p.values = bootpvals))
  }else{
    
    asympvals = double((a.order + 2))
    for(i in 1:(a.order + 2))
    {
      k1   = seq(-1999999, -1, by = 1)
      k2   = seq(1, 1999999, by = 1)
      
      k1sq = k1^2
      k2sq = k2^2
      
      SQtest_stat = origstats[i]^2
      
      if(alternatives == "greater" || alternatives == "lesser"){
        eqna          = exp(-2*SQtest_stat)
        asympvals[i]  = round(eqna,3)
      }
      if(alternatives == "two-sided"){
        eqnb          = sum(((-1)^(k1+1))*exp(-2*SQtest_stat*k1sq))
        eqnc          = sum(((-1)^(k2+1))*exp(-2*SQtest_stat*k2sq))
        asympvals[i]  = round((eqnb + eqnc), 3)
      }
      if(alternatives == "temporary"){
        eqnd          = 1- sum(2*((4*k2sq*SQtest_stat) - 1)*exp(-2*k2sq*SQtest_stat))
        asympvals[i]  = round(eqnd,3)
      }
      
      
    }
    
    asympvals = rbind(asympvals)
    rownames(asympvals) = "p.value"
    if(a.order > 0)
    {
      colnames(asympvals) = c("mean", "var", paste("phi", c(1:a.order), sep = ""))
    }
    if(a.order == 0)
    {
      colnames(asympvals) = c("mean", "var")
    }
    return(list(index = origpoints, stats = origstats, p.values = asympvals))
    
    
    
    
  }
  
}


