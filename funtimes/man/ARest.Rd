\name{ARest}
\alias{ARest}
\title{Estimation of Autoregressive (AR) Parameters}
\description{
Estimate parameters \eqn{\phi} of autoregressive time series model
\deqn{X_t = \sum_{i=1}^p\phi_iX_{t-i} + e_t,} 
by default using robust difference-based estimator and BIC to select the order \eqn{p}. This function is mainly employed for time series filtering in functions \code{\link{sync.test}} and \code{\link{wavk.test}}.
}

\usage{
ARest(x, ar.order = NULL, ar.method = "HVK", BIC = TRUE)
}
\arguments{
  \item{x}{a time series vector.}

  \item{ar.order}{order of autoregressive model when \code{BIC = FALSE}, or the maximal order for BIC-based filtering. Default is \code{round(10*log10(length(x)))}.}
  
    \item{ar.method}{method of estimating autoregression coefficients. Default \code{"HVK"} delivers robust difference-based estimates by Hall and Van Keilegom (2003). Alternatively, options of \command{ar} function can be used, such as \code{"burg"}, \code{"ols"}, \code{"mle"}, and \code{"yw"}.}

  \item{BIC}{logical value indicates whether the order of autoregressive filter should be selected by Bayesian information criterion (BIC). If \code{TRUE} (default), models of orders \eqn{p=} 0,1,...,\code{ar.order} or \eqn{p=} 0,1,...,\code{round(10*log10(length(x)))} are considered, depending on whether \code{ar.order} is defined or not.}
}


\value{A vector of estimated AR coefficients. Returns \code{numeric(0)} if the final \eqn{p=0}.}


\references{
Hall, P. and Van Keilegom, I. (2003). Using difference-based methods for inference in nonparametric regression with time series errors. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 65: 443--456. DOI: \href{http://dx.doi.org/10.1111/1467-9868.00395}{10.1111/1467-9868.00395}
}

\author{Vyacheslav Lyubchich}

\seealso{\code{\link[stats]{ar}}, \code{\link{HVK}}, \code{\link{sync.test}}, \code{\link{wavk.test}}
}


\examples{
# Fix seed for reproduceable simulations:
set.seed(1)

#Simulate some time series, possibly with trend:
n <- 100
Y <- arima.sim(n=n, list(order=c(2,0,0), ar=c(-0.7, -0.1)))
plot.ts(Y)

#Estimate the coefficients:
ARest(Y) #HVK by default
ARest(Y, ar.method="yw") #Yule--Walker
ARest(Y, ar.method="burg") #Burg
}

\keyword{ts}
