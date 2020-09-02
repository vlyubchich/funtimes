#' Change point detection in time series via a linear regression with temporally correlated errors
#' 
#' The function tests for a change point in parameters of a linear regression model with errors exhibiting a general weakly
#' dependent structure. The approach extends earlier methods based on cumulative sums derived under assumption of independent
#' errors. The approach applies smoothing when the time series is dominated
#' by high frequencies.  To detect multiple changes, it is recommended to employ a binary or wild segmentation (see \insertCite{gombay2010change}{funtimes})
#'
#' @param y  A numeric time series vector. Missing values are not allowed.
#' @param a.order Order of the autoregressive model which must be a nonnegative integer number.
#' @param crit.type A string parameter allowing to choose "asymptotic" or "bootstrap" options.
#' @param bootstrap.method A string parameter allowing to choose "nonparametric" or "parametric" method of bootstrapping. 
#' "nonparametric" - resampling of the estimated residuals (with replacement).
#' "parametric"    - sampling innovations from a normal distribution.
#' @param num.bootstrap Number of bootstrap replications if crit.type="bootstrap". Default number is 1,000.
#'
#' @return A list with the following components:
#' \item{index}{Time point where the change has occurred.}
#' \item{stat}{Test statistic.}
#' \item{p.value}{\code{p-value} of the change point test.}
#'  
#' @references 
#' \insertAllCited{}
#' 
#' @seealso \code{\link{mcusum.test}} for change point test for regression
#' 
#' 
#' @author Poli Nemkova, Dorcas Ofori-Boateng, Yulia R. Gel
#' 
#' @export 
#' @examples
#' #Example 1:
#' #Simulate some time series:
#' series_1 = rnorm(157, 2, 1)
#' series_2 = rnorm(43, 7, 10)
#' main_val = c(series_1, series_2)
#' #Now let's perform a change point detection:
#' boot.test.stat(series_1, 1) #- no change -#
#' boot.test.stat(main_val, 1) #- one change, asymptotic critical region -#
#' boot.test.stat(main_val, 1, "bootstrap", "parametric") #- one change, parametric bootstrap -#
#' boot.test.stat(main_val, 1, "bootstrap", "nonparametric") #- one change, nonparametric bootstrap -#
#'
#' 
#' #Example 2:
#' #Consider time series with ratio of real GDP per family to the median income. This is a
#' #skewness and income inequality measure for the US families from 1947 till 2012.       
#' e.data = (Ecdat::incomeInequality['mean.median'])
#' incomeInequality.ts = ts(e.data,start=(1947),end=(2012),frequency = 1)
#' #Now let's perform a change point detection:
#' boot.test.stat(incomeInequality.ts, 0)
#' boot.test.stat(incomeInequality.ts, 0, "bootstrap", "parametric")
#' boot.test.stat(incomeInequality.ts, 0, "bootstrap", "nonparametric")
#' incomeInequality.ts[13] # median income
#' incomeInequality$Year[13] + 1 # year of change point
#' #The first change point occurs at the 13th time point, that is 1960, where the ratio of real GDP per family to the median income is 1.940126.
#' #This ratio shows that in 1960 the national wealth was not distributed equally between all the population
#' #and that most people earn almost twice less than the equal share of the all produced goods and services by the nation
#' #Note: In order to look for the other possible change points run the same function for the segment 
#' #of time series after value #13.



boot.test.stat<-function(y,a.order, crit.type = c("asymptotic", "bootstrap"), bootstrap.method=c("nonparametric","parametric"), num.bootstrap=1000)
{
  test.stat<-function(y)
  {
    n<-length(y)
    x<-c(1:n)/n
    
    inter<-summary(lm(y~x))$coefficients[1,1]
    slope<-summary(lm(y~x))$coefficients[2,1]
    err<-y-inter-slope*x
    #slope<-arima(y, order=c(0,0,2), xreg=x)$coef[(q+2)]
    
    y.bar<-mean(y)
    x.bar<-mean(x)
    
    k.vec<-c(1:n)
    x.bar.k<-cumsum(x)/k.vec
    R.n.sum<-y-y.bar-slope*(x-x.bar)
    R.n<-((n/(k.vec*(n-k.vec)))^(1/2))*cumsum(R.n.sum)
    
    w.n.num<-k.vec*((x.bar-x.bar.k)^2)
    w.n.den<-sum((x-x.bar)^2)*(1-k.vec/n)
    w.n<-(1-w.n.num/w.n.den)^(-1/2)
    
    U.n<-w.n*R.n
    U.n[n]<-0
    
    point<-which.max(abs(U.n))
    tstat<-max(abs(U.n))
    
    n<-length(y)
    T_val = 2*log(n)
    a.T = (2*log(T_val))^0.5
    b.T = 2*log(T_val) + 0.5*log(log(T_val)) - 0.5*log(pi)
    
    p.value<-1-exp(-2*exp(-1*(a.T*tstat - b.T)))
    
    return(list(stat=tstat,point=point, p.value=p.value))
  }
  
  
  orig.test<-test.stat(y)
  orig.stat<-orig.test$stat
  orig.p.value =  orig.test$p.value
  point<-orig.test$point
  
  crit.type = match.arg(crit.type)
  if(crit.type == "bootstrap"){
    n<-length(y)
    x<-c(1:n)/n
    inter<-summary(lm(y~x))$coefficients[1,1]
    slope<-summary(lm(y~x))$coefficients[2,1]
    err<-y-inter-slope*x
    
    a<-arima(err,order=c(a.order,0,0),method="ML")
    if(a.order > 0)
    {
      acoef<-a$coef[1:a.order]
    }
    ainter<-a$coef[(a.order+1)]
    
    residuals<-na.omit(c(a$residuals))-mean(na.omit(c(a$residuals)))
    num.res<-n*num.bootstrap
    
    
    bootstrap.method<-match.arg(bootstrap.method)
    
    ### resample residuals
    if(bootstrap.method=="nonparametric")
    {
      ressample<-sample(residuals,size=num.res,replace=TRUE)
    }
    else
    {
      bootstrap.method<-"parametric"
      res.sd<-sqrt(a$sigma2)
      ressample<-rnorm(num.res,mean=0,sd=res.sd)
    }
    
    resmat<-matrix(ressample,ncol=num.bootstrap,nrow=n)
    
    ### bootstrap samples ###
    boot.stats<-double(num.bootstrap)
    
    for(i in 1:num.bootstrap)
    {
      ### centering the residuals ###
      e<-c(resmat[,i])-mean(c(resmat[,i]))
      
      if(a.order > 0)
      {
        ressim<-arima.sim(list(order = c(a.order,0,0), ar = acoef), n = n, innov=e)
      }
      if(a.order==0)
      {
        ressim<-e
      }
      ### the bootstrap is done under null ###
      ysim<-ressim+inter+slope*x
      ### calculate bootstrap statistic ###
      boot.stats[i]<-test.stat(ysim)$stat
    }
    p.value<-length(which(orig.stat < boot.stats))/num.bootstrap
    return(list(index=point,stat=orig.stat,p.value=p.value))
  }else{
    
    return(list(index=point,stat=orig.stat,p.value = orig.p.value))
  }
  
}




