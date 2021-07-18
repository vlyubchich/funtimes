#' Testing for Change Points in Time Series via Polynomial Regression
#' 
#' The function uses a nonlinear polynomial regression model in which it tests for the null 
#' hypothesis of structural stability in the regression parameters against the alternative of 
#' a break at an unknown time. The method is based on the extreme value distribution of a 
#' maximum-type test statistic which is asymptotically equivalent to the maximally selected 
#' likelihood ratio. The resulting testing approach is easily tractable and delivers accurate 
#' size and power of the test, even in small samples \insertCite{aue2008testing}{funtimes}.
#'
#'
#' @param y a vector that contains univariate time series observations. Missing values are 
#' not allowed.
#' @param a.order order of the autoregressive model which must be a nonnegative integer number.
#' @param alpha significance level for testing hypothesis of no change point. Default value 
#' is 0.05.
#' @param crit.type method of obtaining critical values: "asymptotic" (default) or "bootstrap".
#' @param bootstrap.method type of bootstrap if \code{crit.type = "bootstrap"}: "nonparametric" 
#' (default) or "parametric".
#' @param num.bootstrap number of bootstrap replications if \code{crit.type = "bootstrap"}. 
#' Default number is 1000.
#'
#'
#' @return A list with the following components:
#' \item{index}{time point where the change point has occurred.}
#' \item{stat}{test statistic.} 
#' \item{crit.val}{critical region value (CV(alpha, n)).}
#' \item{p.value}{\code{p-value} of the change point test.}
#'  
#' @references 
#' \insertAllCited{}
#' 
#' @seealso \code{\link{mcusum.test}} change point test for regression
#' 
#' @keywords changepoint ts
#' 
#' @author Palina Niamkova, Dorcas Ofori-Boateng, Yulia R. Gel
#' @export 
#' @examples
#' \dontrun{
#' #Example 1:
#' 
#' #Simulate some time series:
#' set.seed(23450)
#' series_1 = rnorm(137, 3, 5)
#' series_2 = rnorm(213, 0, 1)
#' series_val = c(series_1, series_2)
#' AuePolyReg_test(series_1, 1) # no change (asymptotic)
#' AuePolyReg_test(series_val,1) # one change (asymptotic)
#'
#' #Example 2:
#' 
#' #Consider a time series with annual number of world terrorism incidents from 1970 till 2016:
#' c.data = Ecdat::terrorism["incidents"]
#' incidents.ts <- ts(c.data, start = 1970, end = 2016)
#' 
#' #Run a test for change points:
#' AuePolyReg_test(incidents.ts, 2) # one change (asymptotic)
#' AuePolyReg_test(incidents.ts, 2, 0.05,"bootstrap", "parametric", 200) 
#' # one change (bootstrap)
#' incidents.ts[44] #number of victims at the value of change point
#' year <- 197 + 44 - 1  # year when the change point occurred
#' plot(incidents.ts) # see the visualized data
#' 
#' #The structural change point occurred at the 44th value which corresponds to 2013, 
#' #with 11,990 identified incidents in that year. These findings can be explained with 
#' #a recent rise of nationalism and  extremism due to appearance of the social media, 
#' #Fisher (2019): White Terrorism Shows 'Stunning' Parallels to Islamic State's Rise. 
#' #The New York Times.
#' }
#' 
AuePolyReg_test <- function(y, a.order, alpha = 0.05, 
                            crit.type = c("asymptotic", "bootstrap"), 
                            bootstrap.method = c("nonparametric", "parametric"), 
                            num.bootstrap = 1000)
{
    test.stat <- function(y, alpha)
    {
        n<-length(y)
        x1<-(1:n)/n
        v<-rep(-1, n)
        for(k in 3:(n-3))
        {
            lmk1<-lm(y[1:k]~x1[1:k])
            lmk2<-lm(y[(k+1):n]~x1[(k+1):n])
            sk1<-summary(lmk1)$sigma
            sk2<-summary(lmk2)$sigma
            v[k]<-(k-1)*sk1*sk1+(n-k-1)*sk2*sk2
        }
        v<-v[3:(n-3)]
        lmA<-lm(y~x1)
        sA<-summary(lmA)$sigma
        Tn<--n*(min(log(v))-log(n-1)-2*log(sA))
        
        crit<- -2*log(-0.5*log(1-alpha))+2*log(log(n))+2*log(log(log(n)))
        p.value<-1-exp(-2*exp(-0.5*(Tn-2*log(log(n))-2*log(log(log(n))))))
        index<-which.min(v)+2
        return(list(stat=Tn, crit=crit, p.value=p.value, index=index))
    }
    
    
    orig.test<-test.stat(y, alpha)
    orig.stat<-orig.test$stat
    orig.crit = orig.test$crit
    
    index<-orig.test$index
    
    
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
        
        for (i in 1:num.bootstrap) {
            ### centering the residuals ###
            e <- resmat[,i] - mean(resmat[,i])
            if (a.order > 0) {
                ressim <- arima.sim(list(order = c(a.order,0,0), ar = acoef), n = n, innov = e)
            }
            if (a.order == 0) {
                ressim <- e
            }
            ### the bootstrap is done under null ###
            ysim <- ressim + inter + slope*x
            ### calculate bootstrap statistic ###
            boot.stats[i] <- test.stat(ysim, alpha)$stat
        }
        p.value<-length(which(orig.stat < boot.stats))/num.bootstrap
        return(list(index=index,stat=orig.stat,crit.val = orig.crit, p.value=p.value))
    } else {
        orig.p.value <- orig.test$p.value
        return(list(index = index, stat = orig.stat, crit.val = orig.crit, p.value = orig.p.value))
    }
    
}


