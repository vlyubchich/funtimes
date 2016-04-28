pkgname <- "funtimes"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "funtimes-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('funtimes')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("CExpandSlideCluster")
### * CExpandSlideCluster

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CExpandSlideCluster
### Title: Slide-level time series cluster expansion
### Aliases: CExpandSlideCluster
### Keywords: ts trend

### ** Examples

set.seed(123)
u <- rnorm(10)
Xuncl <- matrix(rt(50, 5), 10, 5)
Alpha <- min(cbind(u,Xuncl))
Beta <- max(cbind(u,Xuncl))
CExpandSlideCluster(u, Xuncl, Alpha, Beta, Delta=0.15, Theta=0.8)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CExpandSlideCluster", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CExpandWindowCluster")
### * CExpandWindowCluster

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CExpandWindowCluster
### Title: Window-level time series cluster expansion
### Aliases: CExpandWindowCluster
### Keywords: ts trend

### ** Examples

set.seed(123)
e <- sample(c(TRUE, FALSE), 5, replace=TRUE)
Euncl <- matrix(sample(c(TRUE, FALSE), 5, replace=TRUE), 5, 5)
CExpandWindowCluster(e, Euncl)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CExpandWindowCluster", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CHomogeneity")
### * CHomogeneity

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CHomogeneity
### Title: Time series cluster homogeneity
### Aliases: CHomogeneity
### Keywords: ts trend

### ** Examples

Bu <- rnorm(10)
Bv <- rnorm(10)
Alpha <- min(c(Bu,Bv))
Beta <- max(c(Bu,Bv))
CHomogeneity(Bu, Bv, Alpha, Beta, Delta=0.5)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CHomogeneity", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CNeighbor")
### * CNeighbor

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CNeighbor
### Title: Neighborhood of time series
### Aliases: CNeighbor
### Keywords: ts trend

### ** Examples

Bu <- rnorm(10)
Bv <- rnorm(10)
Alpha <- min(c(Bu,Bv))
Beta <- max(c(Bu,Bv))
CNeighbor(Bu, Bv, Alpha, Beta, Delta=0.5, Theta=0.8)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CNeighbor", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CSlideCluster")
### * CSlideCluster

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CSlideCluster
### Title: Slide-level time series clustering
### Aliases: CSlideCluster
### Keywords: ts trend

### ** Examples

set.seed(123)
X <- matrix(rnorm(50), 10, 5)
CSlideCluster(X)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CSlideCluster", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CWindowCluster")
### * CWindowCluster

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CWindowCluster
### Title: Window-level time series clustering
### Aliases: CWindowCluster
### Keywords: ts trend

### ** Examples

#For example, weekly data come in slides of 4 weeks
p <- 4 #number of layers in each slide (data come in a slide)

#We want to analyze the trend clusters within a window of 1 year
w <- 13 #number of slides in each window
s <- w  #step to shift a window

#Simulate 26 autoregressive time series with two years of weekly data (52*2 weeks), 
#with a 'burn-in' period of 300.
N <- 26
T <- 2*p*w

set.seed(123) 
phi <- c(0.5) #parameter of autoregression
X <- sapply(1:N, function(x) arima.sim(n=T+300, 
  list(order=c(length(phi),0,0),ar=phi)))[301:(T+300),]
colnames(X) <- paste("TS", c(1:dim(X)[2]), sep="")

tmp <- CWindowCluster(X, Delta=NULL, Theta=0.8, p=p, w=w, s=s, Epsilon=1)

#Time series were simulated with the same parameters, but based on the clustering parameters,
#not all time series join the same cluster. We can plot the main cluster for each window, and 
#time series out of the cluster:
par(mfrow=c(2,2))
ts.plot(X[c(1:(p*w)),tmp[1,]==1], ylim=c(-4,4), 
  main="Time series cluster 1 in window 1")
ts.plot(X[c(1:(p*w)),tmp[1,]!=1], ylim=c(-4,4), 
  main="The rest of the time series in window 1")
ts.plot(X[c(1:(p*w))+s*p,tmp[2,]==1], ylim=c(-4,4), 
  main="Time series cluster 1 in window 2")
ts.plot(X[c(1:(p*w))+s*p,tmp[2,]!=1], ylim=c(-4,4), 
  main="The rest of the time series in window 2")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CWindowCluster", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("HVK")
### * HVK

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: HVK
### Title: HVK estimator
### Aliases: HVK
### Keywords: ts

### ** Examples

X <- arima.sim(n=300, list(order=c(1,0,0), ar=c(0.6)))
HVK(as.vector(X), ar.order=1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("HVK", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("WAVK")
### * WAVK

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: WAVK
### Title: WAVK statistic
### Aliases: WAVK
### Keywords: ts trend

### ** Examples

z <- rnorm(300)
WAVK(z, kn=7)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("WAVK", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sync.test")
### * sync.test

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sync.test
### Title: Time series trend synchronism test
### Aliases: sync.test
### Keywords: ts htest trend

### ** Examples


# Fix seed for reproduceable simulations.
set.seed(123)

# Simulate two autoregressive time series of length n without trend (i.e., with zero trend) 
# and apply the synchronism test.
n <- 200
y1 <- arima.sim(n=n, list(order=c(1,0,0), ar=c(0.6)))
y2 <- arima.sim(n=n, list(order=c(1,0,0), ar=c(-0.2)))
X1 <- cbind(y1, y2)

## Not run: 
##D sync.test(X1, B=1000)
## End(Not run)
# Sample output:
##
##  Non-parametric test for synchronism of parametric linear trends
##
##data:  X1
##Test statistic = -0.0712, p-value = 0.452
##alternative hypothesis: trends are not synchronized.
##sample estimates:
##$common_trend_estimates
##               Estimate Std. Error    t value  Pr(>|t|)
##(Intercept)  0.02944134 0.09871156  0.2982563 0.7658203
##t           -0.05858974 0.17033482 -0.3439681 0.7312353
##
##$ar_order_used
##         y1 y2
##ar_order  1  1
##
##$Window_used
##       y1 y2
##Window  8 15
##
##$all_considered_windows
## Window   Statistic p-value Asympt. p-value
##      8 -0.09419625   0.295       0.3423957
##     11 -0.08168139   0.363       0.4103374
##     15 -0.08831680   0.459       0.3733687
##     20 -0.09337623   0.451       0.3466142


# Add a time series y3 with a different linear trend and apply the synchronism test.
t <- c(1:n)/n
y3 <- 1 + 2*t + arima.sim(n=n, list(order=c(1,0,0), ar=c(-0.2)))
X2 <- cbind(y1, y3)

## Not run: 
##D sync.test(X2, B=1000)
## End(Not run)
# Sample output:
##
##  Non-parametric test for synchronism of parametric linear trends
##
##data:  X2
##Test statistic = 0.3027, p-value < 2.2e-16
##alternative hypothesis: trends are not synchronized.
##sample estimates:
##$common_trend_estimates
##              Estimate Std. Error   t value     Pr(>|t|)
##(Intercept) -0.4047268 0.09862909 -4.103523 5.943524e-05
##t            0.8054264 0.17019251  4.732443 4.215118e-06
##
##$ar_order_used
##         y1 y3
##ar_order  1  1
##
##$Window_used
##       y1 y3
##Window  8  8
##
##$all_considered_windows
## Window Statistic p-value Asympt. p-value
##      8 0.3027026       0    3.464035e-04
##     11 0.3527386       0    3.055570e-05
##     15 0.3608431       0    1.998331e-05
##     20 0.3655885       0    1.552063e-05




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sync.test", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("wavk.test")
### * wavk.test

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: wavk.test
### Title: WAVK trend test
### Aliases: wavk.test
### Keywords: ts htest trend

### ** Examples


# Fix seed for reproduceable simulations.
set.seed(123)

# Simulate autoregressive time series of length n with linear trend 1+2*t, 
# where t is a regular sequence on the interval (0,1].
n <- 100
t <- c(1:n)/n
U <- 1+2*t + arima.sim(n=n, list(order = c(2,0,0), ar = c(-0.7, -0.1)))

# Test for linear trend with output of all results.
## Not run: 
##D wavk.test(U, factor.length = "adaptive.selection", H0="linear", out=TRUE, B=1000)
## End(Not run)
# Sample output:
##
##  Trend test by Wang, Akritas and Van Keilegom
##
##data:  U
##WAVK test statistic = 0.8562, adaptively selected window = 4, p-value = 0.356
##alternative hypothesis: presence of a nonlinear trend
##sample estimates:
##$linear_trend_coefficients
##(Intercept)           t 
##  0.9917251   2.0224272 
##
##$AR_coefficients
##     phi_1      phi_2 
##-0.6814546 -0.2404422 
##
##$all_considered_windows
## Window WAVK-statistic p-value
##      4      0.8561654   0.356
##      5      0.8620023   0.320
##      7      0.8691870   0.288
##     10      0.6837790   0.306


# Test H0 of absence of a trend using asymptotic distribution of statistic.
wavk.test(U, method="asympt")

# Sample output:
##
##        Trend test by Wang, Akritas and Van Keilegom
##
##data:  U
##WAVK test statistic = 18.4712, user-defined window = 10, p-value < 2.2e-16
##alternative hypothesis: presence of a trend
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("wavk.test", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
